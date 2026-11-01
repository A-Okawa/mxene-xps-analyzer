"""End-to-end software validation on the SYNTHETIC example dataset (examples/synthetic_demo; no private data)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mxene_xps import dataset as ds
from mxene_xps import sample_metadata as sm
from mxene_xps.calibration import NO_MODEL_TEXT
from mxene_xps.fitting import fit_core_level
from mxene_xps.identifiability import assess
from mxene_xps.io import parse_phi_multiregion, parse_two_column
from mxene_xps.multipak import map_to_kb, pair_spin_orbit, provider_component_table, read_multipak_provider
from mxene_xps.reporting import calibration_summary, calibration_summary_html, fits_export
from mxene_xps.validation import compare_components, compare_grouped, prepare_vendor, summarize

EX = Path(__file__).resolve().parents[1] / "examples" / "synthetic_demo"


def test_example_is_labelled_synthetic():
    assert (EX / "README.md").read_text().lower().count("synthetic") >= 3
    assert "SYNTHETIC" in (EX / "SYNTHETIC_provider_report.txt").read_text().splitlines()[0]
    assert json.loads((EX / "SYNTHETIC_truth.json").read_text())["note"].startswith("synthetic")
    assert all("SYNTHETIC" in p.name or p.name == "README.md" for p in EX.iterdir())


def test_single_region_fit_identifiability_and_provider_compare(kb):
    sp = parse_two_column((EX / "SYNTHETIC_single_region_Nb3d.csv").read_bytes())
    comps = kb.components("nbti_m2c_defect_aware", "Nb 3d")
    fr = fit_core_level(sp.be, sp.intensity, comps)
    assert fr.metrics["success"] and fr.metrics["R2"] > 0.98
    ident = assess(fr)
    assert set(ident["status"]) <= {"identifiable", "weakly identifiable", "boundary-limited", "non-identifiable"}
    vend = pd.read_csv(EX / "SYNTHETIC_provider_components_Nb3d.csv")
    vp = prepare_vendor(vend, comps, "total")
    summ = summarize(compare_components(vp, fr.components, kb.validation_defaults), compare_grouped(vp, fr.components, kb.validation_defaults))
    assert summ.status in ("PASS", "REVIEW", "CHECK", "CHECK-MISSING", "FAIL")


def test_phi_multiregion_and_provider_pairing(kb):
    regs = parse_phi_multiregion((EX / "SYNdemo203_3_SYNTHETIC_multi.csv").read_bytes())
    assert [r.core_level for r in regs] == ["C 1s", "O 1s", "Nb 3d", "Ti 2p", "F 1s"]
    with pytest.raises(ValueError, match="multi-region"):
        parse_two_column((EX / "SYNdemo203_3_SYNTHETIC_multi.csv").read_bytes())
    bands = read_multipak_provider((EX / "SYNTHETIC_provider_report.txt").read_bytes(), "demo")
    nb = bands[(bands.data_file == "SYNdemo203.spe") & (bands.region == "Nb3d")]
    d, _ = pair_spin_orbit(nb)
    assert (d["pairing"] == "paired").all() and len(d) == 4
    m = map_to_kb(d, kb.components("nbti_m2c_defect_aware", "Nb 3d")).table
    assert m["mapped_component"].notna().all()
    assert set(provider_component_table(m)["provider_label"]) == set(d["provider_label"])


def test_cross_sample_quant_and_report(kb):
    bands, act = ds.read_provider_report((EX / "SYNTHETIC_provider_report.txt").read_bytes(), "demo")
    assert act["data_file"].nunique() == 3
    meta = sm.read((EX / "SYNTHETIC_sample_metadata.csv").read_bytes()).set_index("provider_file")
    samples = []
    for f in dict.fromkeys(bands["data_file"]):
        stem = f.rsplit(".", 1)[0]
        raw = parse_phi_multiregion(next(EX.glob(f"{stem}_*multi.csv")).read_bytes())
        samples.append(ds.build_sample(f, bands, act, kb, meta=meta.loc[stem].to_dict() | {"provider_file": stem}, raw_regions=raw))
    assert [s.structure for s in samples] == ["Nb2AlC", "Nb2CTx", "(Nb,Ti)2CTx"]
    for s in samples:
        assert s.quant_values["apparent_M_deficiency_status"] == "provisional"
        assert len(s.provisional) and (s.provisional["provenance"] == "derived_quantity").all()
        assert s.python and all("ident" in v for v in s.python.values())
    cl = samples[1].mappings["Cl 2p"]
    assert cl["comparison_label"].str.contains(NO_MODEL_TEXT).all() and (cl["pairing"] == "paired").all()
    long = ds.long_table(samples)
    assert {"provider_paired_doublet", "provider_raw_fit", "python_fit", "derived_quantity"} <= set(long["provenance"])
    s = calibration_summary("synthetic demo", kb, "test", {}, samples=samples)
    assert s["release"]["app_version"] and s["release"]["kb_version"] == kb.version
    html = calibration_summary_html(s)
    assert "Interpretation levels (not equivalent)" in html and "Known scientific open issues" in html
