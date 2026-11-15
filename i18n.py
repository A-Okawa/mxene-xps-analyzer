# -*- coding: utf-8 -*-
"""Screen languages: 日本語, English (same scheme as figapp / XRD Maker).

app.py writes each screen text as ui('日本語', 'English'). The language of the current screen is
st.session_state['lang'] (a code, 'ja' or 'en'); ?lang=en in the URL chooses the first language.
Texts with {placeholders} are filled after the language is chosen: ui('{n} 件', '{n} items', n=5).

Only the screen changes. Scientific messages written by src/mxene_xps (fit diagnostics, quantification
warnings, knowledge-base notes) and the provider request text stay in English, as do exported files.
"""

LANGUAGES = {"日本語": "ja", "English": "en"}
CODES = list(LANGUAGES.values())
NAMES = {code: name for name, code in LANGUAGES.items()}


def _japanese():
    return "ja"


lang_getter = _japanese


def use_streamlit_session(key="lang"):
    """Take the screen language from st.session_state[key]."""
    import streamlit as st
    global lang_getter

    def from_session():
        code = st.session_state.get(key)
        return code if code in CODES else "ja"
    lang_getter = from_session


def lang():
    try:
        code = lang_getter()
    except Exception:  # noqa: BLE001  (outside a Streamlit run)
        return "ja"
    return code if code in CODES else "ja"


def ui(ja, en, **kw):
    """The text for the current screen language; {placeholders} are filled with kw."""
    text = ja if lang() == "ja" else en
    return text.format(**kw) if kw else text
