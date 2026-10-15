"""Headless UI smoke tests via streamlit.testing (no browser)."""
from pathlib import Path

import pytest

st_testing = pytest.importorskip("streamlit.testing.v1")
APP = str(Path(__file__).resolve().parents[1] / "app.py")


def _app():
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.run()
    assert not at.exception, at.exception
    return at


def test_app_renders_all_tabs():
    at = _app()
    assert len(at.tabs) == 5
    req = at.text_area[0].value
    assert "Ti 2p: 4 candidate doublet(s)" in req and "do NOT force a vacancy-specific" in req


def test_raw_area_without_rsf_is_blocked():
    at = _app()
    next(b for b in at.button if b.label == "Calculate").click().run()   # empty template, raw_area mode
    assert not at.exception
    assert any("no component areas entered" in e.value for e in at.error)


@pytest.mark.parametrize("structure", ["Nb2CTx", "Ti3C2Tx", "Custom"])
def test_structures_render(structure):
    at = _app()
    at.sidebar.selectbox[0].select(structure).run()
    assert not at.exception


def test_demo_fit_runs():
    at = _app()
    demo = next(c for c in at.checkbox if "synthetic demo" in c.label)
    demo.check().run()
    next(b for b in at.button if b.label == "Run Python fit").click().run()
    assert not at.exception
    assert not at.error
    assert any(m.label == "R²" for m in at.metric)


def test_switching_profile_after_fit_does_not_crash():
    """Regression: a Ti 2p fit from the (Nb,Ti) profile must not be re-plotted with the Ti3C2Tx model."""
    at = _app()
    next(c for c in at.checkbox if "synthetic demo" in c.label).check().run()
    next(b for b in at.button if b.label == "Run Python fit").click().run()
    assert not at.exception
    at.sidebar.selectbox[0].select("Ti3C2Tx").run()
    assert not at.exception
