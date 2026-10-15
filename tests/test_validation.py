import pytest

from conftest import FIX
from mxene_xps.fitting import fit_core_level
from mxene_xps.io import parse_two_column
from mxene_xps.validation import compare_components, compare_grouped, prepare_vendor, role_disagreements, summarize


@pytest.fixture(scope="module")
def ti(kb):
    comps = kb.components("nbti_m2c_defect_aware", "Ti 2p")
    sp = parse_two_column((FIX / "synthetic_Ti2p_raw.csv").read_bytes())
    return comps, fit_core_level(sp.be, sp.intensity, comps, background="linear").components


def vendor_from(fit, dbe=(0.05, -0.08, 0.10, -0.05), ffw=(1.05, 0.95, 1.08, 1.02), far=(1.03, 0.98, 1.04, 0.96)):
    v = fit[["component", "BE", "FWHM", "area"]].copy()
    v["BE"] += list(dbe); v["FWHM"] *= list(ffw); v["area"] *= list(far)
    return v


def run(kb, comps, fit, vendor, **kw):
    th = kb.validation_defaults
    vp = prepare_vendor(vendor, comps, **kw)
    c = compare_components(vp, fit, th)
    g = compare_grouped(vp, fit, th)
    return c, g, summarize(c, g)


def test_small_perturbation_pass(kb, ti):
    comps, fit = ti
    c, g, s = run(kb, comps, fit, vendor_from(fit))
    assert set(c["QC"]) == {"PASS"} and set(g["QC"]) == {"PASS"} and s.status == "PASS"
    d = c.set_index("component")["delta_BE_eV"]   # signed: Python - Vendor
    assert d.to_dict() == pytest.approx({"Ti_lowBE": -0.05, "Ti_midBE": 0.08, "Ti_highBE": -0.10, "Ti_oxide": 0.05})


def test_forced_1eV_shift_fail(kb, ti):
    comps, fit = ti
    c, _, _ = run(kb, comps, fit, vendor_from(fit, dbe=(1.0, 0, 0, 0)))
    assert c.loc[c["component"] == "Ti_lowBE", "QC"].item() == "FAIL"


def test_missing_vendor_component(kb, ti):
    comps, fit = ti
    c, g, s = run(kb, comps, fit, vendor_from(fit).iloc[:-1])
    assert c.loc[c["component"] == "Ti_oxide", "QC"].item() == "CHECK-MISSING"
    assert "CHECK-MISSING" in set(g["QC"])   # metal_oxide role absent from vendor
    assert s.level == "warning"


def test_grouped_pass_but_internal_ti_differs(kb, ti):
    comps, fit = ti
    v = vendor_from(fit, dbe=(0, 0, 0, 0), ffw=(1, 1, 1, 1), far=(1, 1, 1, 1))
    shift = 0.6 * v.loc[1, "area"]
    v.loc[1, "area"] -= shift; v.loc[2, "area"] += shift   # move area midBE -> highBE (both framework_metal)
    c, g, s = run(kb, comps, fit, v)
    assert set(g["QC"]) == {"PASS"}
    assert set(c["QC"]) - {"PASS"}
    assert s.level == "warning" and "not identifiable" in s.message


def test_grouped_fail_is_error(kb, ti):
    comps, fit = ti
    v = vendor_from(fit, dbe=(0, 0, 0, 0), ffw=(1, 1, 1, 1), far=(1, 1, 1, 4))
    _, g, s = run(kb, comps, fit, v)
    assert "FAIL" in set(g["QC"]) and s.level == "error"


def test_missing_fwhm_is_check_not_pass(kb, ti):
    comps, fit = ti
    v = vendor_from(fit).drop(columns="FWHM")
    c, _, _ = run(kb, comps, fit, v)
    assert set(c["QC_FWHM"]) == {"CHECK"} and "PASS" not in set(c["QC"])


def test_primary_only_vendor_areas_converted(kb, ti):
    comps, fit = ti
    v = vendor_from(fit, dbe=(0, 0, 0, 0), ffw=(1, 1, 1, 1), far=(1, 1, 1, 1))
    v["area"] = v["area"] / 1.5
    vp = prepare_vendor(v, comps, doublet_area="primary")
    assert vp["area"].to_numpy() == pytest.approx(fit["area"].to_numpy())


def test_vendor_role_disagreement_reported(kb, ti):
    comps, fit = ti
    v = vendor_from(fit); v["role"] = ["framework_metal", "framework_metal", "metal_oxide", "metal_oxide"]
    msgs = role_disagreements(prepare_vendor(v, comps))
    assert len(msgs) == 1 and "Ti_highBE" in msgs[0]
