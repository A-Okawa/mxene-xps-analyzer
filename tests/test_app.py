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
    at.selectbox(key="sb_structure").select("(Nb,Ti)2CTx").run()     # request content of the mixed-metal profile
    req = next(t for t in at.text_area if t.label == "Request text").value
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
    at.radio(key="v_src").set_value("upload").run()
    at.checkbox(key="v_demo").check().run()
    at.button(key="btn_fit").click().run()
    assert not at.exception
    assert not at.error
    assert any(m.label == "R²" for m in at.metric)


def test_switching_profile_after_fit_does_not_crash():
    """Regression: a Ti 2p fit from the (Nb,Ti) profile must not be re-plotted with the Ti3C2Tx model."""
    at = _app()
    at.radio(key="v_src").set_value("upload").run()
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
    assert [t.label for t in at.tabs] == ["1. 一括読込・試料比較", "2. フィット検証", "3. 定量", "4. レポート",
                                          "5. 知識ベース", "6. 測定依頼"]
    assert at.sidebar.header[0].value == "① 試料の設定"


def test_english_tabs():
    assert [t.label for t in _app().tabs] == ["1. Import & compare", "2. Fit check", "3. Quantification", "4. Report",
                                              "5. Knowledge base", "6. Request"]


def test_query_parameter_selects_english():
    at = st_testing.AppTest.from_file(APP, default_timeout=120)
    at.query_params["lang"] = "en"
    at.run()
    assert at.session_state["lang"] == "en" and at.tabs[0].label == "1. Import & compare"


def test_switching_language_keeps_inputs():
    at = _app("ja")
    at.text_input(key="sb_name_Ti3C2Tx").set_value("S-42").run()
    at.radio(key="v_src").set_value("upload").run()
    at.checkbox(key="v_demo").check().run()
    at.radio(key="lang").set_value("en").run()
    assert not at.exception
    assert at.text_input(key="sb_name_Ti3C2Tx").value == "S-42"
    assert at.checkbox(key="v_demo").value is True
    assert at.tabs[0].label == "1. Import & compare"


def test_japanese_demo_fit_and_messages():
    at = _app("ja")
    at.radio(key="v_src").set_value("upload").run()
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
    at.radio(key="xs_src").set_value("demo").run()
    assert not at.exception
    assert not at.error, [e.value for e in at.error]
    samples = at.session_state["xs_samples"]
    assert len(samples) == 3
    assert at.multiselect(key="xs_pick").value == [s.sample for s in samples]


def test_hf_hcl_requests_f1s_and_cl2p():
    at = _app()
    at.selectbox(key="sb_structure").select("Ti3C2Tx").run()
    at.multiselect(key="sb_etch_Ti3C2Tx").set_value(["HF/HCl"]).run()
    assert not at.exception
    regions = set(at.dataframe[0].value["region"]) if False else None
    req = next(t for t in at.text_area if t.label == "Request text").value
    assert "HF/HCl" in req and "F 1s" in req and "Cl 2p" in req


def test_fit_check_uses_the_data_loaded_in_tab_1():
    """Tab 1 practice data -> tab 2 'use tab-1 data': pick sample and core, fit, compare with the report automatically."""
    at = _app("ja")
    at.radio(key="xs_src").set_value("demo").run()
    assert at.radio(key="v_src").value == "tab1"
    smp = at.selectbox(key="v_smp")
    i = next(k for k, lbl in enumerate(smp.options) if "(Nb0.6Ti0.4)2CTx" in lbl)
    smp.set_value(i).run()
    core = at.selectbox(key="v_core_SYNdemo203")
    assert core.value == "Ti 2p" and "Nb 3d" in core.options
    core.set_value("Nb 3d").run()
    at.button(key="btn_fit").click().run()
    assert not at.exception and not at.error, [e.value for e in at.error]
    assert any(m.label == "収束" for m in at.metric)
    fits = at.session_state["fits"]
    key = next(k for k in fits if k.startswith("3_(Nb0.6Ti0.4)2CTx | nbti_m2c_defect_aware | Nb 3d"))
    assert "grouped_qc" in fits[key]           # provider report from tab 1 was compared automatically
    assert at.session_state["lang"] == "ja"


def test_fit_check_tab1_without_data_explains_what_to_do():
    at = _app("ja")
    assert any("タブ 1" in i.value for i in at.info)
    assert not at.exception


def test_default_structure_is_ti3c2tx():
    assert _app().selectbox(key="sb_structure").value == "Ti3C2Tx"
