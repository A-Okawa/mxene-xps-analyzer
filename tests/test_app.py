"""Headless UI smoke tests via streamlit.testing (no browser)."""
from pathlib import Path

import pytest

st_testing = pytest.importorskip("streamlit.testing.v1")
APP = str(Path(__file__).resolve().parents[1] / "app.py")


def _app(lang="en"):
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.session_state["lang"] = lang
    at.run()
    assert not at.exception, at.exception
    return at


def test_app_renders_all_tabs():
    at = _app()
    assert len(at.tabs) == 6
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
    at.selectbox(key="sb_structure").select(structure).run()
    assert not at.exception


def test_demo_fit_runs():
    at = _app()
    at.checkbox(key="v_demo").check().run()
    at.button(key="btn_fit").click().run()
    assert not at.exception
    assert not at.error
    assert any(m.label == "R²" for m in at.metric)


def test_switching_profile_after_fit_does_not_crash():
    """Regression: a Ti 2p fit from the (Nb,Ti) profile must not be re-plotted with the Ti3C2Tx model."""
    at = _app()
    at.checkbox(key="v_demo").check().run()
    at.button(key="btn_fit").click().run()
    assert not at.exception
    at.selectbox(key="sb_structure").select("Ti3C2Tx").run()
    assert not at.exception


def _gated(password="lab-secret"):
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.secrets["password"] = password
    at.session_state["lang"] = "en"
    return at


def test_password_gate_blocks_until_correct_password():
    at = _gated().run()
    assert not at.exception
    assert len(at.tabs) == 0 and len(at.file_uploader if hasattr(at, "file_uploader") else []) == 0
    at.text_input(key="_auth_pw").set_value("wrong").run()
    assert len(at.tabs) == 0 and any("Incorrect password" in e.value for e in at.error)
    at.text_input(key="_auth_pw").set_value("lab-secret").run()
    assert not at.exception and len(at.tabs) == 6


def test_no_password_secret_means_no_gate():
    assert len(_app().tabs) == 6


# ---------------------------------------------------------------- screen language (日本語 / English)
def test_japanese_is_the_default_language():
    at = st_testing.AppTest.from_file(APP, default_timeout=120).run()
    assert not at.exception
    assert at.session_state["lang"] == "ja"
    assert [t.label for t in at.tabs] == ["1. 測定依頼", "2. 一括読込・試料比較", "3. フィット検証", "4. 定量",
                                          "5. レポート", "6. 知識ベース"]
    assert at.sidebar.header[0].value == "① 試料の設定"


def test_english_tabs():
    assert [t.label for t in _app().tabs] == ["1. Request", "2. Import & compare", "3. Fit check", "4. Quantification",
                                              "5. Report", "6. Knowledge base"]


def test_query_parameter_selects_english():
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.query_params["lang"] = "en"
    at.run()
    assert at.session_state["lang"] == "en" and at.tabs[0].label == "1. Request"


def test_switching_language_keeps_inputs():
    at = _app("ja")
    at.text_input(key="sb_name_(Nb,Ti)2CTx").set_value("S-42").run()
    at.checkbox(key="v_demo").check().run()
    at.radio(key="lang").set_value("en").run()
    assert not at.exception
    assert at.text_input(key="sb_name_(Nb,Ti)2CTx").value == "S-42"
    assert at.checkbox(key="v_demo").value is True
    assert at.tabs[0].label == "1. Request"


def test_japanese_demo_fit_and_messages():
    at = _app("ja")
    at.checkbox(key="v_demo").check().run()
    at.button(key="btn_fit").click().run()
    assert not at.exception and not at.error
    assert any(m.label == "収束" for m in at.metric)


def test_japanese_password_gate():
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.secrets["password"] = "pw"
    at.run()
    at.text_input(key="_auth_pw").set_value("x").run()
    assert any("パスワードが違います" in e.value for e in at.error)


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_synthetic_demo_dataset_in_compare_tab(lang):
    at = _app(lang)
    at.checkbox(key="xs_demo").check().run()
    assert not at.exception
    assert not at.error, [e.value for e in at.error]
    samples = at.session_state["xs_samples"]
    assert len(samples) == 3
    assert at.multiselect(key="xs_pick").value == [s.sample for s in samples]
