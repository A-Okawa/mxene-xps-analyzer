"""Regression tests for the user's decisions on spec discrepancies D1–D3 (2026-10-16)."""
import json

import numpy as np
import pytest

from conftest import FIX, ROOT
from mxene_xps.fitting import fit_core_level
from mxene_xps.io import parse_two_column
from mxene_xps.knowledge_base import with_overrides
from mxene_xps.models import SampleSetup
from mxene_xps.request_generator import generate_request

SEED = json.loads((ROOT / "tools" / "kb_seed_from_handoff_s27.json").read_text())


def _comp(kb, profile, core, name):
    return next(c for c in kb.components(profile, core) if c.name == name)


def _ref(c, table=None):
    refs = c.literature_reference or []
    return refs[0] if table is None else next(r for r in refs if r.get("table") == table)


def _decision(kb, did):
    return next(d for d in kb.raw["spec_discrepancies"] if d["id"] == did)


# ---------------------------------------------------------------- D1
D1_S27_TI2P = {  # §27 ti3c2_brette: (center, bounds, fwhm, fwhm_bounds)
    "Ti_lowBE": (454.9, (454.3, 455.4), 1.0, (0.6, 1.8)),
    "Ti_midBE": (455.6, (455.0, 456.2), 1.6, (0.8, 2.4)),
    "Ti_highBE": (456.8, (456.1, 457.5), 1.7, (0.8, 2.5)),
    "Ti_oxide_oxyfluoride": (459.3, (458.3, 460.1), 1.8, (1.0, 3.2)),
}
D1_BRETTE_SI = {  # Brette SI HF10 % before sputtering (handoff §12): (2p3/2 center, FWHM, 2p1/2 center, FWHM)
    "Ti_lowBE": (454.9, 1.0, 460.8, 1.1),
    "Ti_midBE": (455.6, 1.6, 461.5, 1.6),
    "Ti_highBE": (456.8, 1.7, 462.4, 1.6),
    "Ti_oxide_oxyfluoride": (459.3, 1.8, 464.8, 3.2),
}
D1_BRETTE_C1S = {"C_MXene": (281.8, 0.6), "C_contam_1": (284.1, 1.9), "C_contam_2": (286.4, 1.8), "C_contam_3": (288.9, 2.0)}


def test_d1_ti3c2_uses_s27_json_not_s12_table(kb):
    for name, (center, bounds, fwhm, fb) in D1_S27_TI2P.items():
        c = _comp(kb, "ti3c2_brette", "Ti 2p", name)
        assert (c.center, c.bounds, c.fwhm, c.fwhm_bounds) == (center, bounds, fwhm, fb)
    # the §12-table values (= (Nb,Ti) profile) must not leak into Ti3C2Tx
    low = _comp(kb, "ti3c2_brette", "Ti 2p", "Ti_lowBE")
    assert low.bounds != (454.2, 455.5) and low.fwhm != 1.1


def test_d1_brette_si_stored_as_literature_reference(kb):
    for name, (c32, f32, c12, f12) in D1_BRETTE_SI.items():
        lr = _ref(_comp(kb, "ti3c2_brette", "Ti 2p", name), "S3-A")
        assert lr["source"] == "Brette2024" and lr["line"] == "Ti 2p3/2" and lr["dataset"] == "HF10"
        assert (lr["center"], lr["fwhm"]) == (c32, f32)
        assert (lr["secondary"]["line"], lr["secondary"]["center"], lr["secondary"]["fwhm"]) == ("Ti 2p1/2", c12, f12)
    for name, (cen, fw) in D1_BRETTE_C1S.items():
        lr = _ref(_comp(kb, "ti3c2_brette", "C 1s", name), "S4-A")
        assert (lr["center"], lr["fwhm"]) == (cen, fw)
    # reference values are separate data: they differ from the starting value where Brette SI differs
    c = _comp(kb, "ti3c2_brette", "C 1s", "C_MXene")
    assert c.center == 281.85 and _ref(c, "S4-A")["center"] == 281.8


def test_d1_windows_are_editable_implementation_bounds(kb):
    comps = kb.components("ti3c2_brette", "Ti 2p")
    new, log = with_overrides(comps, {"Ti_lowBE": {"bounds": (454.0, 455.8), "fwhm_bounds": (0.5, 2.2)}})
    low = next(c for c in new if c.name == "Ti_lowBE")
    assert low.bounds == (454.0, 455.8) and low.fwhm_bounds == (0.5, 2.2)
    assert low.literature_reference == comps[0].literature_reference          # reference untouched
    assert {r["field"] for r in log} == {"bounds", "fwhm_bounds"}             # edits are logged
    assert _comp(kb, "ti3c2_brette", "Ti 2p", "Ti_lowBE").bounds == (454.3, 455.4)   # KB unchanged
    with pytest.raises(ValueError):
        with_overrides(comps, {"Ti_lowBE": {"bounds": (455.0, 455.8)}})     # start value outside window
    with pytest.raises(ValueError, match="not editable"):
        with_overrides(comps, {"Ti_lowBE": {"role": "metal_oxide"}})        # roles cannot be changed here
    with pytest.raises(ValueError, match="not editable"):
        with_overrides(comps, {"Ti_lowBE": {"split": 6.0}})


def test_d1_overridden_window_is_used_by_fit(kb):
    comps = kb.components("nbti_m2c_defect_aware", "Nb 3d")
    sp = parse_two_column((FIX / "synthetic_Nb3d_raw.csv").read_bytes())   # true Nb_MXene at 204.95 eV
    narrow, _ = with_overrides(comps, {"Nb_MXene": {"bounds": (204.5, 204.7), "center": 204.6}})
    r = fit_core_level(sp.be, sp.intensity, narrow, background="linear")
    be = r.components.set_index("component").loc["Nb_MXene", "BE"]
    assert 204.5 <= be <= 204.7
    assert any("Nb_MXene" in d and "edge" in d for d in r.diagnostics)


def test_d1_resolved_in_kb(kb):
    d = _decision(kb, "D1")
    assert d["status"] == "resolved" and d["decided_by"] == "user" and "§27" in d["decision"]


# ---------------------------------------------------------------- D2
D2_HALIM = {"Nb_MXene": (204.8, 1.1), "Nb2O5": (207.3, 1.7), "NbF5_like": (209.2, 1.0)}


def test_d2_nb2c_reference_values(kb):
    seed = {c["name"]: c for c in SEED["profiles"]["nb2c_halim"]["fit_models"]["Nb 3d"]}
    for name, (cen, fw) in D2_HALIM.items():
        c = _comp(kb, "nb2c_halim", "Nb 3d", name)
        assert (c.center, c.fwhm) == (cen, fw)
        lr = _ref(c)
        assert lr["center"] == cen and lr["fwhm"] == fw
        assert lr["source"] == "Halim2018Vacancy" and lr["line"] == "Nb 3d5/2"
        # windows and doublet constraints are unchanged from §27
        assert c.bounds == tuple(seed[name]["bounds"]) and c.fwhm_bounds == tuple(seed[name]["fwhm_bounds"])
        assert (c.split, c.ratio_primary_secondary, c.role) == (seed[name]["split"], seed[name]["ratio_primary_secondary"], seed[name]["role"])


def test_d2_does_not_change_mixed_metal_profile(kb):
    seed = {c["name"]: c for c in SEED["profiles"]["nbti_m2c_defect_aware"]["fit_models"]["Nb 3d"]}
    for c in kb.components("nbti_m2c_defect_aware", "Nb 3d"):
        assert (c.center, c.fwhm) == (seed[c.name]["center"], seed[c.name]["fwhm"])
        assert c.literature_reference is None


def test_d2_values_are_starting_values_not_fixed(kb):
    """The fit must be free to move away from 204.8 eV when the data say otherwise."""
    comps = kb.components("nb2c_halim", "Nb 3d")
    sp = parse_two_column((FIX / "synthetic_Nb3d_raw.csv").read_bytes())   # true positions 204.95 / 207.40 / 209.45
    r = fit_core_level(sp.be, sp.intensity, comps, background="linear").components.set_index("component")
    assert abs(r.loc["Nb_MXene", "BE"] - 204.95) < 0.05 and abs(r.loc["Nb_MXene", "BE"] - 204.8) > 0.1
    assert abs(r.loc["NbF5_like", "BE"] - 209.45) < 0.1


def test_d2_request_shows_reference_values(kb):
    setup = SampleSetup("s", "Nb2CTx", 2, 1, "nb2c_halim", ["Nb"], ["HF"])
    txt, _ = generate_request(setup, kb)
    assert "start BE 204.80 eV" in txt and "start FWHM 1.10 eV" in txt
    assert "literature reference (Halim2018Vacancy, Nb2CTx, Nb 3d5/2): 207.30 eV, FWHM 1.70 eV" in txt
    assert "not a fixed position" in txt


def test_d2_resolved_in_kb(kb):
    d = _decision(kb, "D2")
    assert d["status"] == "resolved" and "204.8/207.3/209.2" in d["decision"]


# ---------------------------------------------------------------- D3
def test_d3_reference_ids(kb):
    refs = kb.references
    assert "Halim2016" in refs and "Halim2016XPS" not in refs
    assert "Chen2024" in refs and "Chen2023" not in refs
    assert refs["Chen2024"]["year"] == 2024 and refs["Chen2024"]["copyright_or_online_year"] == 2023
    assert refs["Chen2024"]["doi"] == "10.1002/smll.202307408"
    assert _comp(kb, "ti3c2_brette", "F 1s", "F_residue").source == "Halim2016"


def test_d3_no_stale_ids_anywhere(kb):
    stale = {"Chen2023", "Halim2016XPS"}
    assert not stale & set(kb.references)
    sources = {c.source for p in kb.profile_keys for cs in kb.fit_models(p).values() for c in cs}
    assert not stale & sources and sources <= set(kb.references)
    for path in [*(ROOT / "src").rglob("*.py"), ROOT / "app.py"]:
        text = path.read_text(encoding="utf-8")
        assert not any(k in text for k in stale), path


def test_d3_resolved_in_kb(kb):
    d = _decision(kb, "D3")
    assert d["status"] == "resolved" and "Chen2024" in d["decision"] and "Halim2016" in d["decision"]


def test_all_decisions_recorded_and_version_bumped(kb):
    assert kb.version >= "1.1.0"
    assert {d["id"] for d in kb.raw["spec_discrepancies"] if d["status"] == "resolved"} == {"D1", "D2", "D3"}
    assert any(e["version"] == "1.1.0" for e in kb.raw["changelog"])
