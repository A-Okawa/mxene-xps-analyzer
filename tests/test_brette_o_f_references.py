"""Brette SI O 1s / F 1s experimental values (HF10, HF48): reference metadata only (2026-10-16 correction).

Values transcribed from the Brette2024 SI (downloadSupplement.pdf), Tables S5-A (HF48), S5-B (HF10),
S8 (HF10 before sputtering) and S4-B (HF10 before sputtering).
"""
import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import ROOT
from mxene_xps.fitting import fit_core_level, synthesize_spectrum
from mxene_xps.knowledge_base import KnowledgeBase, with_overrides
from mxene_xps.models import FitComponent, SampleSetup
from mxene_xps.quantification import QuantOptions, quantify_sample
from mxene_xps.request_generator import generate_request
from mxene_xps.validation import compare_components, compare_grouped, prepare_vendor

P = "ti3c2_brette"
SEED = json.loads((ROOT / "tools" / "kb_seed_from_handoff_s27.json").read_text())

# (component, table, dataset) -> (center, FWHM, SI label)
EXPECTED = {
    ("O_ideal", "S5-B", "HF10"): (529.7, 1.0, "MXene (P1)"),
    ("O_disordered", "S5-B", "HF10"): (531.0, 2.0, "MXene (P2)"),
    ("OH_like", "S5-B", "HF10"): (532.5, 1.6, "MXene (P3)"),
    ("O_adsorbate", "S5-B", "HF10"): (534.0, 1.6, "H2O + OR (P4)"),
    ("O_ideal", "S8", "HF10"): (529.6, 0.9, "MXene (P1)"),
    ("O_disordered", "S8", "HF10"): (531.1, 2.0, "MXene (P2)"),
    ("OH_like", "S8", "HF10"): (532.5, 1.5, "MXene (P3)"),
    ("O_adsorbate", "S8", "HF10"): (534.0, 2.1, "H2O + OR (P4)"),
    ("O_ideal", "S5-A", "HF48"): (529.8, 1.0, "MXene (P1)"),
    ("O_disordered", "S5-A", "HF48"): (531.0, 2.0, "MXene (P2)"),
    ("OH_like", "S5-A", "HF48"): (532.0, 1.7, "MXene (P3)"),
    ("F_MXene", "S5-B", "HF10"): (684.9, 1.3, "MXene"),
    ("F_MXene", "S4-B", "HF10"): (684.8, 1.3, "MXene"),
    ("F_MXene", "S5-A", "HF48"): (685.1, 1.6, "MXene"),
}
HF48_ADSORBATE = {(533.8, 1.7, "OR (P4)"), (535.8, 2.0, "H2O (P4)")}


def comp(kb, core, name):
    return next(c for c in kb.components(P, core) if c.name == name)


def core_of(name):
    return "F 1s" if name.startswith("F_") else "O 1s"


@pytest.mark.parametrize("key", list(EXPECTED))
def test_exact_si_values(kb, key):
    name, table, dataset = key
    refs = [r for r in comp(kb, core_of(name), name).literature_reference
            if r["table"] == table and r["dataset"] == dataset]
    assert len(refs) == 1
    r = refs[0]
    assert (r["center"], r["fwhm"], r["si_label"]) == EXPECTED[key]
    assert r["source"] == "Brette2024" and "not a starting value" in r["note"]


def test_hf48_two_adsorbate_contributions_kept(kb):
    refs = [r for r in comp(kb, "O 1s", "O_adsorbate").literature_reference if r["dataset"] == "HF48"]
    assert {(r["center"], r["fwhm"], r["si_label"]) for r in refs} == HF48_ADSORBATE


def test_no_sputtered_values_and_no_f_residue_value(kb):
    for core in ("O 1s", "F 1s"):
        for c in kb.components(P, core):
            for r in c.literature_reference or []:
                assert "after" not in r["sample_state"].lower()
    assert comp(kb, "F 1s", "F_residue").literature_reference is None


def test_s27_starting_values_and_windows_not_overwritten(kb):
    seed = {(core, c["name"]): c for core in ("O 1s", "F 1s") for c in SEED["profiles"][P]["fit_models"][core]}
    for (core, name), sc in seed.items():
        c = comp(kb, core, name)
        assert (c.center, c.fwhm) == (sc["center"], sc["fwhm"])
        assert c.bounds == tuple(sc["bounds"]) and c.fwhm_bounds == tuple(sc["fwhm_bounds"])
        assert c.role == sc["role"]
    # e.g. HF10 OH-like is 532.5 eV in the SI but the §27 starting value 532.2 eV is kept
    assert comp(kb, "O 1s", "OH_like").center == 532.2


def test_reference_outside_window_does_not_invalidate_component():
    d = {"name": "X", "kind": "single", "center": 530.0, "bounds": [529.5, 530.5], "fwhm": 1.0, "fwhm_bounds": [0.5, 1.5],
         "role": "termination_O", "source": "Brette2024", "interpretation": "test",
         "literature_reference": [{"source": "Brette2024", "line": "O 1s", "center": 540.0, "fwhm": 9.0}]}
    c = FitComponent.from_dict(d, "O 1s")
    assert c.bounds == (529.5, 530.5)


def test_fit_is_independent_of_literature_reference(kb):
    comps = kb.components(P, "O 1s")
    x = np.linspace(526, 538, 1201)
    y = synthesize_spectrum(x, comps, [(529.75, 1.0, 0.3, 3000), (531.0, 1.9, 0.3, 1800), (532.3, 1.6, 0.3, 900),
                                       (534.2, 1.7, 0.3, 300)], seed=7)
    stripped = [dataclasses.replace(c, literature_reference=None, literature_interpretation=None) for c in comps]
    bogus = [dataclasses.replace(c, literature_reference=[{"source": "Brette2024", "line": "O 1s", "center": 1e3, "fwhm": 50.0}])
             for c in comps]
    cols = ["BE", "FWHM", "eta", "area"]
    ref = fit_core_level(x, y, comps, background="shirley").components[cols].to_numpy()
    for variant in (stripped, bogus):
        np.testing.assert_array_equal(fit_core_level(x, y, variant, background="shirley").components[cols].to_numpy(), ref)


def test_validation_is_independent_of_literature_reference(kb):
    comps = kb.components(P, "F 1s")
    fit = pd.DataFrame({"component": ["F_MXene", "F_residue"], "BE": [684.95, 686.4], "FWHM": [1.4, 1.8],
                        "area": [9.0, 1.0], "role": ["termination_F", "residue_F"], "area_fraction": [0.9, 0.1]})
    vendor = pd.DataFrame({"component": ["F_MXene", "F_residue"], "BE": [685.0, 686.5], "FWHM": [1.35, 1.9], "area": [8.8, 1.2]})
    th = kb.validation_defaults
    stripped = [dataclasses.replace(c, literature_reference=None) for c in comps]
    for a, b in ((compare_components, compare_components), (compare_grouped, compare_grouped)):
        pd.testing.assert_frame_equal(a(prepare_vendor(vendor, comps), fit, th), b(prepare_vendor(vendor, stripped), fit, th))


def test_engine_modules_never_read_literature_reference():
    src = ROOT / "src" / "mxene_xps"
    for mod in ("fitting.py", "validation.py", "quantification.py", "io.py"):
        assert "literature_reference" not in (src / mod).read_text(), mod


def test_literature_reference_not_editable(kb):
    comps = kb.components(P, "O 1s")
    with pytest.raises(ValueError, match="not editable"):
        with_overrides(comps, {"O_disordered": {"literature_reference": []}})
    new, _ = with_overrides(comps, {"O_disordered": {"bounds": (530.0, 532.0)}})
    assert next(c for c in new if c.name == "O_disordered").literature_reference == comp(kb, "O 1s", "O_disordered").literature_reference


def test_kb_json_values_change_nothing_but_metadata(kb):
    """Removing every literature_reference from the KB JSON leaves all fit models identical otherwise."""
    raw = json.loads(kb.to_json())
    for prof in raw["profiles"].values():
        for comps in prof["fit_models"].values():
            for c in comps:
                c.pop("literature_reference", None)
                c.pop("literature_interpretation", None)
    kb2 = KnowledgeBase(raw)
    for p in kb.profile_keys:
        for core, comps in kb.fit_models(p).items():
            for a, b in zip(comps, kb2.fit_models(p)[core]):
                assert dataclasses.replace(a, literature_reference=None, literature_interpretation=None) == b


def test_request_labels_values_as_reference_only(kb):
    setup = SampleSetup("s", "Ti3C2Tx", 3, 2, P, ["Ti"], ["HF"])
    txt, _ = generate_request(setup, kb)
    line = next(l for l in txt.splitlines() if "Table S5-A" in l and "535.80" in l)
    assert "reference value, not a fixed position" in line and "HF48" in line
    assert "start BE 531.00 eV" in txt and "start BE 532.20 eV" in txt   # §27 starting values still used


# ------------------------------------------------ ~531 eV interpretation (Brette)
def test_531_contribution_is_disorder_vacancy_not_unique_oh(kb):
    od = comp(kb, "O 1s", "O_disordered")
    assert od.role == "termination_O_disordered" and od.role != "termination_OH_like"
    assert "vacanc" in od.interpretation and "not automatically OH" in od.interpretation
    li = od.literature_interpretation
    assert "disordered" in li and "vacancies on T-groups" in li and "NOT a uniquely assigned OH" in li
    assert all("P2" in r["si_label"] for r in od.literature_reference)
    oh = comp(kb, "O 1s", "OH_like")
    assert all("P3" in r["si_label"] for r in oh.literature_reference)
    assert not any("P2" in r["si_label"] for r in oh.literature_reference)


def test_531_component_quantified_separately_from_oh(kb):
    rows = [("C 1s", "C_MXene", 1.0, "framework_C"), ("O 1s", "O_disordered", 0.5, "termination_O_disordered"),
            ("O 1s", "OH_like", 0.2, "termination_OH_like")]
    df = pd.DataFrame([{"sample": "s", "core_level": c, "component": n, "area": a, "role": r} for c, n, a, r in rows])
    v = quantify_sample(df, {}, QuantOptions(m_target=3, x_target=2, input_mode="corrected_area")).values
    assert v["z_O_disordered"] == pytest.approx(1.0) and v["z_OH_like"] == pytest.approx(0.4)
