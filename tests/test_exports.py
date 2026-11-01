"""Export consistency: no exported value loses its provenance; versions and interpretation levels are always present."""
import json
from pathlib import Path

import pandas as pd

from mxene_xps import __version__
from mxene_xps import sensitivity as sens
from mxene_xps.fitting import fit_core_level
from mxene_xps.identifiability import assess
from mxene_xps.io import parse_two_column
from mxene_xps.models import MeasurementMetadata, SampleSetup
from mxene_xps.multipak import apply_user_mapping, map_to_kb, pair_spin_orbit, provider_component_table, read_multipak_provider
from mxene_xps.provenance import PROVENANCE
from mxene_xps.quantification import QuantOptions, provisional_report, quantify_sample
from mxene_xps.release import OPEN_ISSUES, release_info
from mxene_xps.reporting import calibration_summary, calibration_summary_html, fits_export, project_dict
from mxene_xps.validation import compare_components, compare_grouped, prepare_vendor, summarize

EX = Path(__file__).resolve().parents[1] / "examples" / "synthetic_demo"


def _session(kb):
    comps = kb.components("nbti_m2c_defect_aware", "Nb 3d")
    sp = parse_two_column((EX / "SYNTHETIC_single_region_Nb3d.csv").read_bytes())
    fr = fit_core_level(sp.be, sp.intensity, comps)
    bands = read_multipak_provider((EX / "SYNTHETIC_provider_report.txt").read_bytes(), "demo")
    one = bands[(bands.data_file == "SYNdemo203.spe") & (bands.region == "Nb3d")]
    d, linked = pair_spin_orbit(one)
    auto = map_to_kb(d, comps).table
    did = auto.loc[auto["mapped_component"] == "NbF5_like", "doublet_id"].item()
    mapping, log = apply_user_mapping(auto, {did: None}, comps)          # a user override
    vp = prepare_vendor(provider_component_table(mapping), comps, "total")
    cqc, gqc = compare_components(vp, fr.components, kb.validation_defaults), compare_grouped(vp, fr.components, kb.validation_defaults)
    sres = sens.run(sp.be, sp.intensity, comps, sens.standard_alternatives("Nb 3d", comps, sp.be), spectrum_core_level="Nb 3d")
    fits = {"nbti_m2c_defect_aware | Nb 3d": {
        "fit": fr, "core": "Nb 3d", "ident": assess(fr, alt_spans=sres.component_spans()), "sensitivity": sres,
        "settings": {"profile": "nbti_m2c_defect_aware", "background": "shirley",
                     "model_overrides": [{"component": "Nb2O5", "field": "center", "kb_value": 207.45, "user_value": 207.3}]},
        "provider_raw_bands": linked, "provider_mapping": mapping, "provider_mapping_log": log,
        "component_qc": cqc, "grouped_qc": gqc, "summary": summarize(cqc, gqc).__dict__}}
    return fits


def test_fits_export_preserves_everything(kb):
    ex = fits_export(_session(kb))["nbti_m2c_defect_aware | Nb 3d"]
    assert all(r["provenance"] == "python_fit" for r in ex["components"])
    assert {r["provider_label"] for r in ex["provider_raw_bands"]} >= {"Nb3d band 1", "Nb3d band 2"}
    assert all(r["provenance"] == "provider_raw_fit" for r in ex["provider_raw_bands"])
    pairs = [r for r in ex["provider_mapping"] if r["pairing"] == "paired"]
    assert pairs and all(r["provenance"] == "provider_paired_doublet" and "+" in r["provider_label"] for r in pairs)
    assert any(r["mapping_method"] == "user" for r in ex["provider_mapping"])
    assert ex["provider_mapping_log"] and all(r["provenance"] == "user_override" for r in ex["provider_mapping_log"])
    assert ex["model_overrides"][0]["provenance"] == "user_override"
    assert {r["status"] for r in ex["identifiability"]} and all(r["provenance"] == "python_fit" for r in ex["identifiability"])
    assert all("model_sensitive" in r and r["provenance"] == "derived_quantity" for r in ex["sensitivity"]["spans"])
    assert all(r["provenance_python"] == "python_fit" for r in ex["component_qc"])
    assert all(set(d) == {"severity", "text"} for d in ex["diagnostics"])
    json.dumps(ex, default=str)      # serializable


def test_every_provenance_value_is_in_vocabulary(kb):
    ex = fits_export(_session(kb))
    vals = set()

    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == "provenance" and isinstance(v, str):
                    vals.add(v)
                walk(v)
        elif isinstance(o, list):
            for x in o:
                walk(x)
    walk(ex)
    assert vals and vals <= set(PROVENANCE)


def test_versions_and_levels_in_all_reports(kb):
    fits = _session(kb)
    rel = release_info(kb, commit="abc")
    assert rel["app_version"] == __version__ and rel["kb_version"] == kb.version
    assert set(rel["interpretation_levels"]) == {"fit_succeeded", "statistically_identifiable", "chemically_validated"}
    assert len(rel["scientific_open_issues"]) == len(OPEN_ISSUES) >= 9
    setup = SampleSetup(sample_name="s", structure="(Nb,Ti)2CTx", m_target=2, x_target=1, profile_key="nbti_m2c_defect_aware",
                        m_elements=["Nb", "Ti"], etchants=[], a_element="Al", expected_defects=True)
    proj = project_dict(setup, MeasurementMetadata(), [], kb.version, "nbti_m2c_defect_aware", {}, None, None, None, None, None, None, [])
    assert proj["release"]["app_version"] == __version__ and proj["release"]["kb_version"] == kb.version
    s = calibration_summary("demo", kb, "abc", fits, sensitivity={k: v["sensitivity"] for k, v in fits.items()})
    assert s["release"]["kb_version"] == kb.version and s["levels_not_equivalent"] == [
        "fit succeeded", "statistically identifiable", "chemically validated"]
    assert all({"severity", "text"} == set(w) for f in s["fits"] for w in f["qc_warnings"])
    html = calibration_summary_html(s)
    for txt in (__version__, kb.version, "chemically_validated", "Known scientific open issues"):
        assert txt in html


def test_provisional_alternative_exported():
    t = pd.DataFrame([{"sample": "s", "core_level": c, "component": n, "area": a, "role": r} for c, n, a, r in (
        ("C 1s", "C", 10.0, "framework_C"), ("Nb 3d", "Nb_MXene", 12.0, "framework_metal"),
        ("Nb 3d", "Nb_suboxide_like", 3.0, "metal_suboxide_candidate"))])
    opt = QuantOptions(m_target=2, x_target=1, input_mode="corrected_area")
    rep = provisional_report(quantify_sample(t, {}, opt).values, opt)
    rec = json.loads(rep.to_json(orient="records"))
    m = next(r for r in rec if r["quantity"] == "M_norm")
    assert m["nominal (default)"] == 1.2 and m["alternative (candidate as framework)"] == 1.5 and m["provisional"]
