import numpy as np
import pandas as pd
import pytest

from conftest import FIX
from mxene_xps.fitting import fit_core_level, shirley_background
from mxene_xps.io import parse_two_column

PROFILE = "nbti_m2c_defect_aware"


def _fit(kb, truth, tag, fname, bg="linear"):
    # Synthetic truth was generated with the pre-calibration model (no calibration candidates).
    comps = [c for c in kb.components(PROFILE, truth[tag]["core"]) if not c.calibration_candidate]
    sp = parse_two_column((FIX / fname).read_bytes())
    r = fit_core_level(sp.be, sp.intensity, comps, background=bg)
    t = pd.DataFrame(truth[tag]["truth"])
    t["area_true"] = t["area_primary"] * [1 + 1 / c.ratio_primary_secondary for c in comps]
    m = t.merge(r.components, on="component")
    return r, m


@pytest.mark.parametrize("bg", ["linear", "shirley"])
def test_ti2p_refit_29_4(kb, truth, bg):
    r, m = _fit(kb, truth, "Ti2p", "synthetic_Ti2p_raw.csv", bg)
    assert r.metrics["success"] and r.metrics["R2"] > 0.999
    assert (m["center"] - m["BE"]).abs().max() < 0.2
    assert (m["fwhm"] - m["FWHM"]).abs().max() < 0.3
    assert (r.components["area"] >= 0).all()
    fw = m["role"] == "framework_metal"
    assert abs(m.loc[fw, "area"].sum() / m.loc[fw, "area_true"].sum() - 1) < 0.10   # role-grouped partition stable
    # identifiability: internal Ti sub-components are allowed to redistribute, but must be flagged
    assert any("strongly correlated" in d and "Ti_" in d for d in r.diagnostics)


def test_ti2p_matches_prototype_numbers(kb, truth):
    """Prototype (v0.1) reference: R2~0.99984, max|dBE|~0.110, max|dFWHM|~0.217, framework error ~7.3 %."""
    r, m = _fit(kb, truth, "Ti2p", "synthetic_Ti2p_raw.csv", "linear")
    assert r.metrics["R2"] == pytest.approx(0.99984, abs=5e-5)
    assert (m["center"] - m["BE"]).abs().max() == pytest.approx(0.110, abs=0.03)
    assert (m["fwhm"] - m["FWHM"]).abs().max() == pytest.approx(0.217, abs=0.05)


@pytest.mark.parametrize("bg", ["linear", "shirley"])
def test_nb3d_refit_29_5(kb, truth, bg):
    r, m = _fit(kb, truth, "Nb3d", "synthetic_Nb3d_raw.csv", bg)
    assert r.metrics["success"] and r.metrics["R2"] > 0.999
    assert (m["center"] - m["BE"]).abs().max() < 0.05
    assert (m["fwhm"] - m["FWHM"]).abs().max() < 0.1
    assert (m["area"] / m["area_true"] - 1).abs().max() < 0.25


def test_descending_headerless_input_gives_identical_fit(kb, truth):
    r1, _ = _fit(kb, truth, "Ti2p", "synthetic_Ti2p_raw.csv")
    r2, _ = _fit(kb, truth, "Ti2p", "synthetic_Ti2p_raw_desc_noheader.txt")
    np.testing.assert_allclose(r1.components[["BE", "FWHM", "area"]].to_numpy(),
                               r2.components[["BE", "FWHM", "area"]].to_numpy(), rtol=1e-6, atol=1e-6)


def test_metrics_and_curves(kb, truth):
    r, _ = _fit(kb, truth, "Nb3d", "synthetic_Nb3d_raw.csv")
    for k in ("RMSE", "R2", "AIC", "BIC", "success", "message"):
        assert k in r.metrics
    for col in ("BE", "observed", "background", "fit_total", "residual", "Nb_MXene"):
        assert col in r.curves
    assert np.all(np.diff(r.curves["BE"]) > 0)


def test_too_few_points(kb):
    with pytest.raises(ValueError, match="20"):
        fit_core_level(np.arange(10.0), np.ones(10), kb.components(PROFILE, "Nb 3d"))


def test_component_outside_range_flagged(kb):
    x = np.linspace(202, 206, 200)
    y = 100 + 1000 * np.exp(-0.5 * ((x - 204.9) / 0.6) ** 2)
    r = fit_core_level(x, y, kb.components(PROFILE, "Nb 3d"), background="linear")
    assert any("outside the measured BE range" in d for d in r.diagnostics)


def test_shirley_endpoints_and_monotonic():
    x = np.linspace(450, 470, 400)
    y = 100 + 50 * (x > 458) + 800 * np.exp(-0.5 * ((x - 458) / 0.8) ** 2)
    b = shirley_background(x, y)
    assert b[0] == pytest.approx(y[0]) and b[-1] == pytest.approx(y[-1])
    assert np.all(np.diff(b) >= -1e-9)
