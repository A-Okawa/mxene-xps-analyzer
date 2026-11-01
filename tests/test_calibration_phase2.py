"""Calibration phase 2 (KB 1.3.0): candidates are present, separate, flagged, and not used as final values."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mxene_xps.fitting import fit_core_level, synthesize_spectrum
from mxene_xps.multipak import map_to_kb
from mxene_xps.quantification import QuantOptions, quantify_sample

MX = ("nb2c_halim", "nbti_m2c_defect_aware")


def _c(kb, p, core, name):
    return next(c for c in kb.components(p, core) if c.name == name)


def test_kb_version_and_status(kb):
    assert kb.version in ("1.3.1", "1.4.0", "1.5.0")
    cs = kb.raw["calibration_status"]
    assert cs["nb3d_suboxide_like_validated"] is False
    assert cs["o1s_crosscheck_nb_suboxide_like"]["result"] == "inconclusive"
    assert cs["d4_nb3d_dSO_calibration_sensitive"] is True and len(cs["delta_finalization_blocked_by"]) == 2
    d4 = next(d for d in kb.raw["spec_discrepancies"] if d["id"] == "D4")
    assert d4["status"].startswith("calibration_sensitive")


def test_d4_kb_splittings_not_changed(kb):
    """D4 was resolved by sensitivity analysis only: KB splittings stay at the literature defaults."""
    for p in ("nb2c_halim", "nbti_m2c_defect_aware"):
        assert {c.split for c in kb.components(p, "Nb 3d")} == {2.8}
    assert {c.split for c in kb.components("nbti_m2c_defect_aware", "Ti 2p")} == {5.8, 5.7, 5.6}
    assert 6.17 not in {c.split for p in kb.profile_keys for cs in kb.fit_models(p).values() for c in cs}


def test_provisional_reasons_reported():
    r = quantify_sample(_table(), {}, QuantOptions(m_target=2, x_target=1, input_mode="corrected_area",
                                                   provisional_reasons=("D4 sensitive",)))
    assert any(w == "delta not finalized: D4 sensitive" for w in r.warnings)


@pytest.mark.parametrize("p", MX)
def test_nb_suboxide_candidate(kb, p):
    c = _c(kb, p, "Nb 3d", "Nb_suboxide_like")
    assert c.kind == "doublet" and c.role == "metal_suboxide_candidate"
    assert c.bounds[0] <= 205.6 and c.bounds[1] >= 205.9
    assert 208.3 <= c.center + c.split <= 208.6
    assert c.calibration_candidate and not c.validated and c.assignment_status == "calibration_candidate"
    assert c.confidence == "medium" and c.cross_check == ["O 1s"]
    assert c.alternative_role == "framework_metal" and "termination-modified framework nb" in c.alternative_interpretation.lower()
    names = [x.name for x in kb.components(p, "Nb 3d")]
    assert {"Nb_MXene", "Nb2O5", "NbF5_like", "Nb_suboxide_like"} == set(names)   # separate, not merged


def test_ti_oxide_bound_and_nb3s_warning(kb):
    t = _c(kb, "nbti_m2c_defect_aware", "Ti 2p", "Ti_oxide")
    assert t.bounds[0] <= 458.17 < t.center and t.bounds[0] >= 457.9      # widened slightly only
    assert t.role == "metal_oxide"
    for p in ("nbti_m2c_defect_aware", "nbti_211_max"):
        w = kb.profile(p)["core_level_warnings"]["Ti 2p"][0]
        assert "Nb 3s" in w and "Ti 2p1/2" in w
    # Ti3C2Tx (no Nb) is untouched
    assert _c(kb, "ti3c2_brette", "Ti 2p", "Ti_oxide_oxyfluoride").bounds == (458.3, 460.1)


@pytest.mark.parametrize("p", MX)
def test_low_be_f_unresolved(kb, p):
    f = _c(kb, p, "F 1s", "F_lowBE_unresolved")
    assert f.role == "unresolved_low_BE_F" and f.bounds[0] <= 682.66 <= f.bounds[1]
    assert f.assignment_status == "unresolved"


def test_max_profiles_are_separate(kb):
    assert kb.structures["Nb2AlC"]["default_profile"] == "nb2alc_max"
    assert kb.structures["(Nb,Ti)2AlC"]["default_profile"] == "nbti_211_max"
    for p, cores in (("nb2alc_max", {"Nb 3d", "C 1s", "O 1s", "Al 2p"}),
                     ("nbti_211_max", {"Ti 2p", "Nb 3d", "C 1s", "O 1s", "Al 2p"})):
        prof = kb.profile(p)
        assert prof["phase"] == "MAX" and set(kb.fit_models(p)) == cores
        nb = {c.name for c in kb.components(p, "Nb 3d")}
        assert "Nb_MAX" in nb and "Nb_MXene" not in nb
        assert _c(kb, p, "Nb 3d", "Nb_MAX").bounds[1] < _c(kb, "nb2c_halim", "Nb 3d", "Nb_MXene").bounds[0]
        # KB 1.4.0: MAX components are literature/transferred; only the unchanged Nb candidate stays a calibration candidate
        assert {c.name for cs in kb.fit_models(p).values() for c in cs if c.calibration_candidate} == {"Nb_suboxide_like"}
        assert {c.role for c in kb.components(p, "Al 2p")} >= {"A_site_metal", "A_site_oxide"}


def _table(sample="s", extra=()):
    rows = [("C 1s", "C_MXene", 10.0, "framework_C"), ("Nb 3d", "Nb_MXene", 12.0, "framework_metal"),
            ("Nb 3d", "Nb_suboxide_like", 3.5, "metal_suboxide_candidate"), ("Nb 3d", "Nb2O5", 6.0, "metal_oxide"),
            ("F 1s", "F_lowBE_unresolved", 0.3, "unresolved_low_BE_F"), ("F 1s", "Nb_F_term", 2.0, "termination_F"),
            *extra]
    return pd.DataFrame([{"sample": sample, "core_level": c, "component": n, "area": a, "role": r} for c, n, a, r in rows])


def test_quant_candidate_not_folded_and_provisional():
    r = quantify_sample(_table(), {}, QuantOptions(m_target=2, x_target=1, input_mode="corrected_area"))
    v = r.values
    assert v["M_norm"] == pytest.approx(1.2)                        # framework Nb only (12/10)
    assert v["M_suboxide_candidate_norm"] == pytest.approx(0.35)
    assert v["M_norm_alt_candidate_as_framework"] == pytest.approx(1.55)
    assert v["apparent_M_deficiency_delta"] == pytest.approx(0.8)
    assert v["apparent_M_deficiency_delta_alt"] == pytest.approx(0.45)
    assert v["apparent_M_deficiency_status"] == "provisional"
    assert v["z_F_raw"] == pytest.approx(0.2) and v["F_unresolved_lowBE_norm"] == pytest.approx(0.03)
    assert any("PROVISIONAL" in w for w in r.warnings) and any("Calibration-candidate" in w for w in r.warnings)


def test_quant_status_after_validation_and_without_nb():
    t = _table()
    t = t[t["role"] != "metal_suboxide_candidate"]
    r = quantify_sample(t, {}, QuantOptions(m_target=2, x_target=1, input_mode="corrected_area", nb_suboxide_validated=True))
    assert r.values["apparent_M_deficiency_status"] == "standard"
    ti = pd.DataFrame([{"sample": "t", "core_level": "C 1s", "component": "C", "area": 2.0, "role": "framework_C"},
                       {"sample": "t", "core_level": "Ti 2p", "component": "Ti", "area": 3.0, "role": "framework_metal"}])
    r = quantify_sample(ti, {}, QuantOptions(m_target=3, x_target=2, input_mode="corrected_area"))
    assert r.values["apparent_M_deficiency_status"] == "standard"


def test_quant_max_phase():
    t = _table(extra=[("Al 2p", "Al_MAX", 5.0, "A_site_metal"), ("Al 2p", "Al_oxide", 2.0, "A_site_oxide")])
    r = quantify_sample(t, {}, QuantOptions(m_target=2, x_target=1, input_mode="corrected_area", phase="MAX", a_target=1))
    assert r.values["A_site_metal_norm"] == pytest.approx(0.5) and r.values["A_target"] == 1
    assert any("MAX-phase" in w for w in r.warnings)


def test_provider_components_map_to_candidates(kb):
    d = pd.DataFrame({"doublet_id": ["a", "b", "c"], "provider_label": ["Nb3d band 1 + Nb3d band 3",
                      "Nb3d band 2 + Nb3d band 5", "Nb3d band 4 + Nb3d band 6"], "BE": [204.50, 205.89, 207.43],
                      "pairing": "paired"})
    m = map_to_kb(d, kb.components("nb2c_halim", "Nb 3d")).table.set_index("doublet_id")
    assert list(m["mapped_component"]) == ["Nb_MXene", "Nb_suboxide_like", "Nb2O5"]
    assert m.at["b", "provider_label"] == "Nb3d band 2 + Nb3d band 5"
    ti = pd.DataFrame({"doublet_id": ["x"], "provider_label": ["Ti2p band 2 + Ti2p band 4"], "BE": [458.17], "pairing": "paired"})
    assert map_to_kb(ti, kb.components("nbti_m2c_defect_aware", "Ti 2p")).table["mapped_component"].item() == "Ti_oxide"
    f = pd.DataFrame({"doublet_id": ["f"], "provider_label": ["F1s band 1"], "BE": [682.66], "pairing": "singlet"})
    assert map_to_kb(f, kb.components("nb2c_halim", "F 1s")).table["mapped_component"].item() == "F_lowBE_unresolved"


def test_fit_keeps_candidate_separate(kb):
    """Synthetic Nb 3d with a 205.8 eV doublet: recovered as its own component, not absorbed by Nb2O5."""
    comps = kb.components("nb2c_halim", "Nb 3d")
    names = [c.name for c in comps]
    truth = {"Nb_MXene": (204.5, 1.2, 0.2, 8000), "Nb_suboxide_like": (205.85, 1.5, 0.2, 2500),
             "Nb2O5": (207.4, 1.5, 0.2, 4200), "NbF5_like": (209.4, 1.2, 0.2, 200)}
    x = np.linspace(202, 214, 900)
    y = synthesize_spectrum(x, comps, [truth[n] for n in names], seed=11)
    r = fit_core_level(x, y, comps, background="linear").components.set_index("component")
    assert abs(r.at["Nb_suboxide_like", "BE"] - 205.85) < 0.15
    frac_true = 2500 / sum(t[3] for t in truth.values())
    assert r.at["Nb_suboxide_like", "area_fraction"] == pytest.approx(frac_true, abs=0.06)
    assert r.at["Nb2O5", "FWHM"] < 2.0


def test_max_structure_renders():
    st_testing = pytest.importorskip("streamlit.testing.v1")
    at = st_testing.AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=120)
    at.run()
    for s in ("Nb2AlC", "(Nb,Ti)2AlC"):
        at.sidebar.selectbox[0].select(s).run()
        assert not at.exception
