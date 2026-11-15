"""Identifiability, fit sensitivity, provisional reporting, provenance, metadata, dataset comparison, unsupported regions
and the calibration summary (no change to the frozen scientific state)."""
import copy
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mxene_xps import dataset as ds
from mxene_xps import sample_metadata as sm
from mxene_xps import sensitivity as sens
from mxene_xps.calibration import D4_NB3D_SPLIT_GRID, NO_MODEL_TEXT, UNSUPPORTED_REGIONS
from mxene_xps.fitting import fit_core_level, synthesize_spectrum
from mxene_xps.identifiability import assess, warnings_from
from mxene_xps.knowledge_base import KnowledgeBase
from mxene_xps.models import FitComponent
from mxene_xps.provenance import PROVENANCE, check_kb_provenance, start_value_provenance
from mxene_xps.quantification import QuantOptions, provisional_report, quantify_sample
from mxene_xps.reporting import calibration_summary, calibration_summary_html

ROOT = Path(__file__).resolve().parents[1]
SASAKI = os.environ.get("XPS_SASAKI_DIR")
needs_data = pytest.mark.skipif(not SASAKI, reason="private provider data (set XPS_SASAKI_DIR)")


def _single(name, c, lo, hi, fw=1.0, fwb=(0.5, 2.5)):
    return FitComponent(name=name, core_level="C 1s", kind="single", center=c, bounds=(lo, hi), fwhm=fw, fwhm_bounds=fwb,
                        role="contaminant", source="s", interpretation="i")


# ------------------------------------------------------------------ 1. identifiability
def test_identifiability_statuses():
    x = np.linspace(280, 292, 700)
    sep = [_single("A", 283.0, 282.0, 284.0), _single("B", 287.0, 286.0, 288.0)]
    y = synthesize_spectrum(x, sep, [(283.0, 1.0, 0.2, 5000), (287.0, 1.0, 0.2, 4000)], seed=1)
    t = assess(fit_core_level(x, y, sep, background="linear")).set_index("component")
    assert set(t["status"]) == {"identifiable"}
    ovl = [_single("A", 285.0, 284.5, 285.5, fw=1.4), _single("B", 285.3, 284.8, 285.8, fw=1.4)]
    y2 = synthesize_spectrum(x, ovl, [(285.1, 1.4, 0.2, 5000), (285.2, 1.4, 0.2, 4000)], seed=2)
    t2 = assess(fit_core_level(x, y2, ovl, background="linear")).set_index("component")
    assert set(t2["status"]) <= {"non-identifiable", "weakly identifiable", "boundary-limited"}
    assert t2["strongly_anticorrelated_with"].notna().any()
    pin = [_single("A", 283.9, 283.8, 284.6), _single("B", 287.0, 286.0, 288.0)]   # true peak at 283.0 lies outside A's window
    t3 = assess(fit_core_level(x, y, pin, background="linear")).set_index("component")
    assert t3.at["A", "status"] == "boundary-limited"
    assert all("not a chemical assignment" in n for n in t3["note"])
    assert (t3["provenance"] == "python_fit").all()


def test_identifiability_uses_alternative_spans():
    x = np.linspace(280, 292, 700)
    sep = [_single("A", 283.0, 282.0, 284.0), _single("B", 287.0, 286.0, 288.0)]
    y = synthesize_spectrum(x, sep, [(283.0, 1.0, 0.2, 5000), (287.0, 1.0, 0.2, 4000)], seed=1)
    fr = fit_core_level(x, y, sep, background="linear")
    t = assess(fr, alt_spans={"A": 0.12, "B": 0.06}).set_index("component")
    assert t.at["A", "status"] == "non-identifiable" and t.at["B", "status"] == "weakly identifiable"
    assert any("A: non-identifiable" in w for w in warnings_from(t.reset_index()))


# ------------------------------------------------------------------ 2. sensitivity
def test_standard_alternatives_use_d4_grid_only_for_nb(kb):
    nb = kb.components("nb2c_halim", "Nb 3d")
    alts = sens.standard_alternatives("Nb 3d", nb, np.linspace(199, 219, 161))
    splits = {a.split for a in alts if a.split is not None}
    assert splits == set(D4_NB3D_SPLIT_GRID) - {2.8}
    assert any(a.background == "linear" for a in alts) and any(a.window for a in alts)
    ti = kb.components("nbti_m2c_defect_aware", "Ti 2p")
    assert not [a for a in sens.standard_alternatives("Ti 2p", ti, np.linspace(452, 466, 120)) if a.split is not None]


def test_sensitivity_run_and_model_sensitive_flag(kb):
    comps = kb.components("nb2c_halim", "Nb 3d")
    x = np.linspace(199, 219, 700)
    truth = [(204.5, 1.2, 0.2, 8000), (205.85, 1.5, 0.2, 2500), (207.4, 1.5, 0.2, 4200), (209.4, 1.2, 0.2, 200)]
    y = synthesize_spectrum(x, comps, truth, seed=3)
    alts = sens.standard_alternatives("Nb 3d", comps, x)
    calls = []

    def derived(fr):
        calls.append(fr)
        return {"M_norm": 2 * fr["Nb_MXene"], "apparent_M_deficiency_delta": 2 - 2 * fr["Nb_MXene"]}
    r = sens.run(x, y, comps, alts, derived_fn=derived, spectrum_core_level="Nb 3d")
    assert len(r.table) == len(alts) and len(calls) == len(alts)
    win = r.table[r.table["window"] != "imported"]
    assert win["dAIC_vs_nominal"].isna().all() and not win["AIC_comparable"].any()
    same = r.table[(r.table["window"] == "imported")]
    assert same["dAIC_vs_nominal"].notna().all() and same.iloc[0]["dAIC_vs_nominal"] == 0
    assert set(r.spans["label"]) <= {"model-sensitive", "stable across alternatives"}
    assert (r.spans["model_sensitive"] == (r.spans["span"] >= 0.05)).all()


def test_derived_from_fractions_keeps_total():
    t = pd.DataFrame({"sample": "s", "core_level": ["C 1s", "Nb 3d", "Nb 3d", "Nb 3d"],
                      "component": ["C_MXene", "Nb_MXene", "Nb_suboxide_like", "Nb2O5"], "area": [10.0, 12.0, 3.0, 5.0],
                      "role": ["framework_C", "framework_metal", "metal_suboxide_candidate", "metal_oxide"]})
    opt = QuantOptions(m_target=2, x_target=1, input_mode="corrected_area")
    fn = sens.derived_from_fractions(t, "s", "Nb 3d", lambda tt: quantify_sample(tt, {}, opt).values)
    out = fn({"Nb_MXene": 0.5, "Nb_suboxide_like": 0.25, "Nb2O5": 0.25})
    assert out["M_norm"] == pytest.approx(1.0) and out["M_norm_alt_candidate_as_framework"] == pytest.approx(1.5)


# ------------------------------------------------------------------ 3. provisional reporting
def test_provisional_report_nominal_and_alternative_not_averaged():
    rows = [("C 1s", "C", 10.0, "framework_C"), ("Nb 3d", "Nb_MXene", 9.0, "framework_metal"),
            ("Nb 3d", "Nb_suboxide_like", 3.0, "metal_suboxide_candidate"), ("Ti 2p", "Ti", 3.0, "framework_metal"),
            ("O 1s", "O_low", 5.0, "O_unresolved_low")]
    t = pd.DataFrame([{"sample": "s", "core_level": c, "component": n, "area": a, "role": r} for c, n, a, r in rows])
    opt = QuantOptions(m_target=2, x_target=1, input_mode="corrected_area", nominal_metal_stoich={"Nb": 1.2, "Ti": 0.8},
                       provisional_reasons=("D4 sensitive",))
    v = quantify_sample(t, {}, opt).values
    rep = provisional_report(v, opt).set_index("quantity")
    ti = rep.loc["Ti/(Ti+Nb)"]
    assert ti["nominal (default)"] == pytest.approx(0.25) and ti["alternative (candidate as framework)"] == pytest.approx(0.2)
    m = rep.loc["M_norm"]
    assert m["nominal (default)"] == pytest.approx(1.2) and m["alternative (candidate as framework)"] == pytest.approx(1.5)
    assert m["provisional"] and "D4 sensitive" in m["why provisional"] and "never averaged" in m["why provisional"]
    assert pytest.approx(1.35) not in list(rep["nominal (default)"].dropna()) + list(rep["alternative (candidate as framework)"].dropna())
    assert rep.loc["z_O_total (termination O)", "provisional"] and "lower bound" in rep.loc["z_O_total (termination O)", "why provisional"]
    assert (rep["provenance"] == "derived_quantity").all()


# ------------------------------------------------------------------ 5. provenance
def test_provenance_vocabulary_and_kb_guard(kb):
    assert set(PROVENANCE) >= {"raw_measurement", "provider_raw_fit", "provider_paired_doublet", "python_fit",
                               "literature_reference", "transferred_reference", "calibration_observation",
                               "user_override", "derived_quantity"}
    assert check_kb_provenance(kb) == []
    assert start_value_provenance(next(c for c in kb.components("nb2c_halim", "Nb 3d") if c.name == "Nb_suboxide_like")) == "calibration_observation"
    assert start_value_provenance(next(c for c in kb.components("nbti_211_max", "Ti 2p") if c.name == "Ti_MAX")) == "transferred_reference"
    raw = copy.deepcopy(kb.raw)
    c = next(c for c in raw["profiles"]["nb2alc_max"]["fit_models"]["C 1s"] if c["name"] == "C_MAX")
    c["source"] = "Sasaki2026MultiPak"            # provider observation used as a starting value of a non-candidate
    with pytest.raises(ValueError, match="provenance violation"):
        KnowledgeBase(raw)


# ------------------------------------------------------------------ 6. sample metadata
def test_metadata_schema_and_example():
    assert list(sm.template().columns) == list(sm.FIELDS)
    d = sm.read((ROOT / "templates" / "next_dataset_example.csv").read_bytes())
    assert list(d["phase"]) == ["MAX", "MAX", "MAX", "MXene", "MXene"]
    assert [sm.infer_structure(r) for r in d.to_dict("records")] == ["Nb2AlC", "(Nb,Ti)2AlC", "(Nb,Ti)2AlC", "Nb2CTx", "(Nb,Ti)2CTx"]
    assert d.loc[3, "treatment"] == "H2SO4" and d.loc[3, "treatment_temperature_C"] == "140"
    assert d.loc[4, "treatment_temperature_C"] == "100"
    with pytest.raises(ValueError, match="phase"):
        sm.normalize(pd.DataFrame({"sample_id": ["x"], "phase": ["film"]}))
    with pytest.raises(ValueError, match="numeric"):
        sm.normalize(pd.DataFrame({"sample_id": ["x"], "phase": ["MAX"], "treatment_temperature_C": ["hot"]}))


def test_metadata_not_in_kb(kb):
    txt = json.dumps(kb.raw)
    assert "H2SO4" not in txt and "140" not in json.dumps(kb.structures)


# ------------------------------------------------------------------ 4. dataset / cross-sample
ACT_TXT = """--------------------------
File Information
--------------------------
File Name:  demo108.spe
File Type:  XPS SPECTRUM
Comment:  4_Nb2CTx
--------------------------
Atomic Concentration Table
--------------------------
Abscissa        C1s        O1s        F1s        Cl2p       Nb3d
               0.314      0.733      1.000      0.954      3.127    RSF
              16.528     39.914     55.540     43.107    140.359    CorrectedRSF
1.00           32.17      39.37       2.80       3.56      21.37
"""


def test_parse_atomic_concentration():
    a = ds.parse_atomic_concentration(ACT_TXT, "demo")
    assert dict(zip(a["core_level"], a["at_percent"])) == {"C 1s": 32.17, "O 1s": 39.37, "F 1s": 2.80, "Cl 2p": 3.56, "Nb 3d": 21.37}
    assert (a["data_file"] == "demo108.spe").all() and (a["provenance"] == "provider_raw_fit").all()


def test_unresolved_tags_are_neutral(kb):
    comps = kb.components("nb2alc_max", "Nb 3d")
    assert ds.unresolved_tag("Nb 3d", 201.8, comps) == "unresolved_low_BE_Nb"
    o = ds.known_observation("Nb 3d", 201.8, "MAX")
    assert o["identifiability"] == "non-identifiable" and "Not assigned" in o["note"]
    for word in ("metal", "NbC", "Nb-Al", "aluminide"):
        assert word not in o["label"]


# ------------------------------------------------------------------ 7. unsupported regions
def test_unsupported_regions_have_no_model(kb):
    assert set(UNSUPPORTED_REGIONS) == {"Cl 2p", "S 2p", "Li 1s"}
    for p in kb.profile_keys:
        assert not set(kb.fit_models(p)) & set(UNSUPPORTED_REGIONS)


def test_fit_refuses_unrelated_core_level(kb):
    x = np.linspace(196, 204, 200)
    y = 100 + np.exp(-0.5 * ((x - 199.3) / 0.6) ** 2)
    with pytest.raises(ValueError, match="unrelated model"):
        fit_core_level(x, y, kb.components("nb2c_halim", "Nb 3d"), spectrum_core_level="Cl 2p")


def test_app_unsupported_region_and_new_tab():
    st_testing = pytest.importorskip("streamlit.testing.v1")
    at = st_testing.AppTest.from_file(str(ROOT / "app.py"), default_timeout=120)
    at.session_state["lang"] = "en"
    at.run()
    assert not at.exception and len(at.tabs) == 6
    at.radio(key="v_src").set_value("upload").run()
    at.selectbox(key="v_core").select("Cl 2p").run()
    assert not at.exception
    assert any(NO_MODEL_TEXT in w.value and "CAUTION" in w.value for w in at.warning)   # not an error: nothing is wrong
    assert at.button(key="btn_fit").disabled


def test_request_s2p_for_h2so4(kb):
    from mxene_xps.models import SampleSetup
    from mxene_xps.request_generator import required_regions
    s = SampleSetup(sample_name="x", structure="Nb2CTx", m_target=2, x_target=1, profile_key="nb2c_halim", m_elements=["Nb"],
                    etchants=["LiF/HCl", "H2SO4 (post-treatment)"], a_element="Al", expected_defects=False)
    regs = {r["region"]: r for r in required_regions(s)}
    assert "S 2p" in regs and "no validated KB model" in regs["S 2p"]["reason"] and "Li 1s" in regs


# ------------------------------------------------------------------ 8. calibration summary
def test_calibration_summary_distinguishes_fit_from_interpretation(kb):
    comps = kb.components("nbti_m2c_defect_aware", "Nb 3d")
    x = np.linspace(202, 214, 600)
    y = synthesize_spectrum(x, comps, [(204.9, 1.3, 0.2, 2800), (205.8, 1.5, 0.2, 600), (207.4, 1.5, 0.2, 800),
                                       (209.4, 1.2, 0.2, 100)], seed=5)
    fr = fit_core_level(x, y, comps)
    fits = {"nbti_m2c_defect_aware | Nb 3d": {"fit": fr, "core": "Nb 3d", "ident": assess(fr),
                                              "settings": {"profile": "nbti_m2c_defect_aware", "background": "shirley",
                                                           "model_overrides": [{"component": "Nb2O5", "field": "center"}]}}}
    s = calibration_summary("demo", kb, "abc123", fits)
    f = s["fits"][0]
    assert f["fit_succeeded"] is True and f["chemical_interpretation_validated"] is False
    assert s["kb_version"] == kb.version and s["git_commit"] == "abc123"
    assert s["frozen_scientific_state"]["nb3d_suboxide_like_validated"] is False
    assert any(u["component"] == "Nb_suboxide_like" for u in s["unresolved_components"])
    assert "Halim2018Vacancy" in s["references_used"] and f["model_overrides"]
    html = calibration_summary_html(s)
    assert "fit succeeded = True" in html and "chemical interpretation validated = False" in html


# ------------------------------------------------------------------ real data (private)
BASE = {  # KB 1.5.0 outputs before this feature work (captured 2026-11-02); must be unchanged
    ("sasaki108", "Nb 3d", "shirley"): (91.479226, {"Nb_MXene": 0.44574, "Nb_suboxide_like": 0.291902, "Nb2O5": 0.257522}),
    ("sasaki110", "Nb 3d", "shirley"): (58.087969, {"Nb_MXene": 0.438476, "Nb_suboxide_like": 0.308838, "Nb2O5": 0.245367}),
    ("sasaki102", "Nb 3d", "shirley"): (72.643666, {"Nb_MAX": 0.183926, "Nb_MAX_2nd": 0.210947, "Nb_suboxide_like": 0.365841}),
    ("sasaki106", "Ti 2p", "shirley"): (125.87808, {"Ti_MAX": 0.625253, "Ti_oxide": 0.374747}),
    ("sasaki110", "O 1s", "linear"): (69.717262, {"O_low": 0.152382, "O_mid": 0.230623, "OH_like": 0.548417}),
}
PROF = {"sasaki102": "nb2alc_max", "sasaki104": "nbti_211_max", "sasaki106": "nbti_211_max", "sasaki108": "nb2c_halim",
        "sasaki110": "nbti_m2c_defect_aware"}


@needs_data
@pytest.mark.parametrize("key", list(BASE))
def test_real_fit_outputs_unchanged(kb, key):
    from mxene_xps.io import parse_two_column
    stem, core, bg = key
    f = next(Path(SASAKI).glob(f"{stem}_*{core.replace(' ', '')}*curve.csv"))
    sp = parse_two_column(f.read_bytes())
    fr = fit_core_level(sp.be, sp.intensity, kb.components(PROF[stem], core), background=bg)
    rmse, fr_ref = BASE[key]
    assert fr.metrics["RMSE"] == pytest.approx(rmse, abs=1e-5)
    got = dict(zip(fr.components["component"], fr.components["area_fraction"]))
    for k, v in fr_ref.items():
        assert got[k] == pytest.approx(v, abs=1e-6)


@needs_data
def test_real_dataset_comparison(kb):
    root = Path(SASAKI)
    bands, act = ds.read_provider_report(next(root.glob("*.pptx")).read_bytes(), "p")
    assert act["data_file"].nunique() == 5
    S = [ds.build_sample(f, bands, act, kb) for f in dict.fromkeys(bands["data_file"])]
    assert [s.structure for s in S] == ["Nb2AlC", "(Nb,Ti)2AlC", "(Nb,Ti)2AlC", "Nb2CTx", "(Nb,Ti)2CTx"]
    q = {s.data_file: s.quant_values for s in S}
    assert q["sasaki108.spe"]["M_norm"] == pytest.approx(1.543774, abs=1e-6)
    assert q["sasaki108.spe"]["M_norm_alt_candidate_as_framework"] == pytest.approx(1.990561, abs=1e-6)
    assert q["sasaki110.spe"]["apparent_M_deficiency_delta"] == pytest.approx(0.83408, abs=1e-6)
    for s in S:
        assert s.quant_values["apparent_M_deficiency_status"] == "provisional"
    s1 = S[0].mappings["Nb 3d"]
    row = s1[s1["comparison_label"] == "unresolved_low_BE_Nb"]
    assert len(row) == 1 and row["identifiability"].item() == "non-identifiable" and row["mapped_component"].isna().all()
    for core in ("Cl 2p", "S 2p"):
        m = S[3].mappings[core]
        assert m["comparison_label"].str.contains(NO_MODEL_TEXT).all() and (m["pairing"] == "paired").all()
    long = ds.long_table(S)
    assert set(long["provenance"]) <= set(PROVENANCE)
    assert not (long["provenance"] == "literature_reference").any()


@needs_data
def test_real_identifiability_known_cases(kb):
    from mxene_xps.io import parse_two_column
    root = Path(SASAKI)
    sp = parse_two_column(next(root.glob("sasaki108_*Nb3dcurve.csv")).read_bytes())
    t = assess(fit_core_level(sp.be, sp.intensity, kb.components("nb2c_halim", "Nb 3d"))).set_index("component")
    assert t.at["Nb_MXene", "strongly_anticorrelated_with"] == "Nb_suboxide_like" or \
        t.at["Nb_suboxide_like", "strongly_anticorrelated_with"] == "Nb_MXene"
    # Sample-1 low-BE Nb candidate (script-local definition from the calibration study) stays non-identifiable
    spec = importlib.util.spec_from_file_location("s1", ROOT / "tools" / "calibration" / "s1_lowbe_nb_sensitivity.py")
    s1 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s1)
    sp1 = parse_two_column(next(root.glob("sasaki102_*Nb3dcurve.csv")).read_bytes())
    for bg in ("shirley", "linear"):
        fr = fit_core_level(sp1.be, sp1.intensity, kb.components("nb2alc_max", "Nb 3d") + s1.candidates("doublet"), background=bg)
        assert assess(fr).set_index("component").at["Nb_lowBE_cand", "status"] == "non-identifiable"
    w = kb.profile("nbti_m2c_defect_aware")["core_level_warnings"]["Ti 2p"][0]
    assert "Nb 3s" in w
