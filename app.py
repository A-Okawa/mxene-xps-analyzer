"""MAX/MXene XPS Analyzer – Streamlit UI (handoff §4, §5).

All scientific logic lives in src/mxene_xps; this file only collects input and displays results.
No data are sent to any external service.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from mxene_xps import __version__  # noqa: E402
from mxene_xps.fitting import fit_core_level  # noqa: E402
from mxene_xps.io import MULTIREGION_IMPORTERS, SPECTRUM_IMPORTERS, crop, read_csv_bytes  # noqa: E402
from mxene_xps.multipak import (apply_user_mapping, map_to_kb, pair_spin_orbit, provider_component_table,  # noqa: E402
                                read_multipak_provider)
from mxene_xps.knowledge_base import default_kb, with_overrides  # noqa: E402
from mxene_xps.models import (BE_REFERENCE_METHODS, ORTHOGONAL_METHODS, ROLE_LABELS, ROLES, SURFACE_STATES,  # noqa: E402
                              MeasurementMetadata, OrthogonalEvidence, SampleSetup)
from mxene_xps.quantification import (LEGACY_CONTROLS_DEFAULT, LEGACY_LABEL, QuantOptions,  # noqa: E402
                                      legacy_quantify_table, legacy_template, quantify_table, role_reclassifications)
from mxene_xps.reporting import evidence_summary, html_report, metadata_warnings, project_dict  # noqa: E402
from mxene_xps.request_generator import generate_request, make_peak_template, metal_core  # noqa: E402
from mxene_xps import dataset as ds  # noqa: E402
from mxene_xps import sample_metadata as smeta  # noqa: E402
from mxene_xps.calibration import NO_MODEL_TEXT, UNSUPPORTED_REGIONS  # noqa: E402
from mxene_xps.identifiability import assess as assess_identifiability, warnings_from as ident_warnings  # noqa: E402
from mxene_xps.provenance import PROVENANCE, start_value_provenance  # noqa: E402
from mxene_xps.quantification import provisional_report  # noqa: E402
from mxene_xps.reporting import calibration_summary, calibration_summary_html, fits_export, git_commit  # noqa: E402
from mxene_xps.release import release_info  # noqa: E402
from mxene_xps import sensitivity as sens  # noqa: E402
from mxene_xps.validation import (IDENTIFIABILITY_WARNING, compare_components, compare_grouped,  # noqa: E402
                                  prepare_vendor, role_disagreements, summarize)

import i18n  # noqa: E402
from i18n import ui  # noqa: E402
from upload_match import match_uploads  # noqa: E402

KB = default_kb()
FIX = ROOT / "tests" / "fixtures"
DEMO = ROOT / "examples" / "synthetic_demo"
STATUS_COLOR = {"PASS": "#d8f0d8", "CHECK": "#fff0c2", "REVIEW": "#fff0c2", "CHECK-MISSING": "#ececec", "FAIL": "#f7d0d0"}
APP_TITLE = "🧪 MXene XPS Analyzer"

st.set_page_config(page_title="MXene XPS Analyzer", page_icon="🧪", layout="wide")

# 画面の言語（コード ja・en）は st.session_state['lang']。URL の ?lang=en で最初の言語を選べる
i18n.use_streamlit_session("lang")
if "lang" not in st.session_state:
    _ql = str(st.query_params.get("lang", "")).lower()
    st.session_state["lang"] = _ql if _ql in i18n.CODES else "ja"
st.sidebar.radio("🌐 Language / 言語", i18n.CODES, format_func=i18n.NAMES.get, key="lang", horizontal=True)


def require_password() -> None:
    """Optional access gate for shared deployments.

    Active only when a ``password`` entry exists in Streamlit secrets (e.g. on Streamlit Community Cloud).
    Without a secrets file (local use, tests) the app opens directly.
    """
    import hmac

    try:
        expected = st.secrets.get("password")
    except Exception:  # no secrets file configured
        expected = None
    if not expected or st.session_state.get("_auth_ok"):
        return
    st.title(APP_TITLE)
    entered = st.text_input(ui("パスワード", "Password"), type="password", key="_auth_pw")
    if entered:
        if hmac.compare_digest(entered.encode(), str(expected).encode()):
            st.session_state["_auth_ok"] = True
            st.rerun()
        st.error(ui("パスワードが違います。", "Incorrect password."))
    st.stop()


require_password()


def style_status(df: pd.DataFrame):
    cols = [c for c in df.columns if c == "QC" or c.startswith("QC_") or c == "metal_QC"]
    fmt = {c: "{:.4g}" for c in df.columns if pd.api.types.is_float_dtype(df[c])}
    return df.style.format(fmt, na_rep="").map(lambda v: f"background-color:{STATUS_COLOR.get(v, '')}", subset=cols)


from mxene_xps.severity import QC_STATUS_SEVERITY, classify as sev_classify  # noqa: E402

_SEV_UI = {"INFO": st.info, "CAUTION": st.warning, "WARNING": st.warning, "INVALID": st.error}
_SEV_ICON = {"INFO": "ℹ️", "CAUTION": "🟡", "WARNING": "⚠️", "INVALID": "⛔"}
_SEV_NAME = {"INFO": lambda: ui("情報", "INFO"), "CAUTION": lambda: ui("注意", "CAUTION"),
             "WARNING": lambda: ui("警告", "WARNING"), "INVALID": lambda: ui("無効", "INVALID")}


def say(a: str, b: str | None = None) -> None:
    """Show a message with a consistent severity: say(text[, LEVEL]) or say(LEVEL, text)."""
    if a in _SEV_UI and b is not None:
        level, text = a, b
    else:
        text, level = a, (b or sev_classify(a))
    _SEV_UI[level](f"**{_SEV_NAME[level]()}** — {text}", icon=_SEV_ICON[level])


def csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8-sig")


def _cal_text(x) -> str:
    if not x.calibration_candidate:
        return ""
    t = f"{'validated' if x.validated else 'CANDIDATE'} (confidence {x.confidence or 'n/a'})"
    if x.alternative_role:
        t += f"; alt role {x.alternative_role}"
    if x.cross_check:
        t += f"; cross-check {', '.join(x.cross_check)}"
    return t


def split_list(text: str) -> list[str]:
    return [t.strip() for t in text.split(",") if t.strip()]


# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.header(ui("① 試料の設定", "① Sample setup"))
    st.caption(ui("最初に試料の構造を選んでください。フィットモデルと定量の基準がここで決まります。",
                  "Choose the sample structure first; it sets the fit model and the quantification basis."))
    _order = ["(Nb,Ti)2CTx", "Nb2CTx", "Ti3C2Tx"]
    _custom = {"Custom": ui("その他（手動設定）", "Custom")}
    structure = st.selectbox(ui("試料の構造", "Nominal structure"),
                             _order + [k for k in KB.structures if k not in _order] + ["Custom"], key="sb_structure",
                             format_func=lambda k: _custom.get(k, k))
    if structure == "Custom":
        m_target = st.number_input(ui("理想の M サイト数", "Ideal M-site count"), min_value=0.1, value=2.0, step=0.1, key="sb_m")
        x_target = st.number_input(ui("理想の X サイト数", "Ideal X-site count"), min_value=0.1, value=1.0, step=0.1, key="sb_x")
        profile_key = st.selectbox(ui("参照プロファイル", "Reference profile"), KB.profile_keys, key="sb_profile")
    else:
        spec = KB.structures[structure]
        m_target, x_target, profile_key = float(spec["m_target"]), float(spec["x_target"]), spec["default_profile"]
        st.caption(f"M_target = {m_target:g}, X_target = {x_target:g}, profile = `{profile_key}`")
    phase = (KB.structures.get(structure) or KB.profile(profile_key)).get("phase", "MXene")
    a_target = (KB.structures.get(structure) or {}).get("a_target")
    if phase == "MAX":
        st.caption(ui("MAX 相のプロファイル：初期値は文献ではなく測定データ由来の校正候補です。",
                      "MAX-phase profile: data-derived calibration candidates (not literature)."))
    default_m = {"(Nb,Ti)2CTx": "Nb,Ti", "Nb2CTx": "Nb", "Ti3C2Tx": "Ti", "Nb2AlC": "Nb", "(Nb,Ti)2AlC": "Nb,Ti"}.get(structure, "")
    m_elements = split_list(st.text_input(ui("M サイトの元素（コンマ区切り）", "M-site elements (comma-separated)"), default_m,
                                          key=f"sb_mel_{structure}"))
    sample_name = st.text_input(ui("試料名", "Sample name"), "NbTi_MXene_01" if structure == "(Nb,Ti)2CTx" else
                                ("MAX_01" if phase == "MAX" else "MXene_01"), key=f"sb_name_{structure}",
                                help=ui("ダウンロードするファイル名にも使います。", "Also used in the names of downloaded files."))
    etchants = st.multiselect(ui("エッチング・処理の薬品", "Etching / processing chemistry"),
                              ["HF", "LiF/HCl", "NaF/HCl", "KF/HCl", "NH4F/HCl", "TBAOH", "molten chloride", "H2SO4 (post-treatment)", "other"],
                              default=["LiF/HCl"] if structure == "(Nb,Ti)2CTx" else [], key=f"sb_etch_{structure}",
                              help=ui("測定を依頼する領域（F 1s、Cl 2p、Li 1s など）の判断に使います。",
                                      "Used to decide which regions to request (F 1s, Cl 2p, Li 1s, ...)."))
    a_element = st.text_input(ui("確認する A サイト元素", "A-site element to check"), "Al", key="sb_a")
    expected_defects = st.checkbox(ui("M サイト欠損・選択溶出が予想される", "M-site vacancies / selective dissolution expected"),
                                   value=structure == "(Nb,Ti)2CTx", key=f"sb_def_{structure}")
    cl_interest = st.checkbox(ui("Cl 終端を見たい（塩化物を使っていなくても Cl 2p を依頼）",
                                 "Cl termination of interest (request Cl 2p even without chloride chemistry)"),
                              value=False, key="sb_cl")
    dopants = split_list(st.text_input(ui("合成後に加えたドーパント（コンマ区切り）", "Post-synthesis dopant(s), comma-separated"), "",
                                       key="sb_dop"))
    nominal: dict[str, float] = {}
    if len(m_elements) > 1:
        st.subheader(ui("M サイトの仕込み組成", "Nominal M-site stoichiometry"))
        for el in m_elements:
            nominal[el] = st.number_input(ui("{el} の仕込み量", "{el} nominal amount", el=el), min_value=0.0,
                                          value=float(m_target / len(m_elements)), step=0.05, key=f"nom_{el}")
    with st.expander(ui("試料の詳細（任意）", "More sample details (optional)")):
        parent = st.text_input(ui("前駆体（MAX 相）の組成", "Parent/MAX composition"), "", key="sb_parent")
        post = st.text_input(ui("合成後の処理", "Post-synthesis treatment"), "", key="sb_post")
        term_family = st.text_input(ui("予想される終端基", "Expected termination family"), "", key="sb_term")
        notes = st.text_area(ui("メモ", "Notes"), "", key="sb_notes")

    with st.expander(ui("測定条件（任意・レポートに記録）", "Measurement conditions (optional; recorded in reports)"), expanded=False):
        md = MeasurementMetadata(
            provider=st.text_input(ui("測定の依頼先", "Provider"), "", key="md_provider") or None,
            instrument=st.text_input(ui("装置", "Instrument"), "", key="md_instr") or None,
            xray_source=st.text_input(ui("X 線源", "X-ray source"), "", key="md_xray") or None,
            pass_energy_eV=st.number_input(ui("パスエネルギー（eV、0 = 不明）", "Pass energy (eV, 0 = unknown)"), min_value=0.0, value=0.0,
                                           key="md_pe") or None,
            step_eV=st.number_input(ui("ステップ幅（eV、0 = 不明）", "Step size (eV, 0 = unknown)"), min_value=0.0, value=0.0,
                                    format="%.3f", key="md_step") or None,
            spot_size=st.text_input(ui("スポット径", "Spot size"), "", key="md_spot") or None,
            charge_neutralization=st.text_input(ui("帯電中和", "Charge neutralization"), "", key="md_cn") or None,
            be_reference_method=st.selectbox(ui("結合エネルギーの基準", "BE reference method"), BE_REFERENCE_METHODS, key="md_beref"),
            be_reference_detail=st.text_input(ui("基準の詳細（どの金属・何 eV など）", "BE reference detail (e.g. which metal / value)"), "",
                                              key="md_bedet") or None,
            surface_state=st.selectbox(ui("表面の状態", "Surface state"), SURFACE_STATES, key="md_surf",
                                       help=ui("スパッタ後のデータは、未処理の表面と同じ扱いにはなりません。",
                                               "Data measured after sputtering are not treated as equivalent to the as-prepared surface.")),
        )
        if md.sputtered:
            md.sputter_species = st.text_input(ui("スパッタのイオン種", "Sputter ion"), "Ar+", key="md_sp_ion") or None
            md.sputter_energy_eV = st.number_input(ui("スパッタのエネルギー（eV）", "Sputter energy (eV)"), min_value=0.0, value=0.0,
                                                   key="md_sp_e") or None
            md.sputter_time_s = st.number_input(ui("スパッタ時間（s）", "Sputter time (s)"), min_value=0.0, value=0.0, key="md_sp_t") or None
    st.divider()
    st.caption(ui("App v{v} · 知識ベース v{k}。アップロードしたファイルはこのアプリを動かすサーバーのメモリ上で処理し、"
                  "ほかのサービスには送りません。",
                  "App v{v} · knowledge base v{k}. Uploaded files are processed in memory by the server running this app "
                  "and are not sent to any other service.", v=__version__, k=KB.version))

setup = SampleSetup(sample_name=sample_name, structure=structure, m_target=m_target, x_target=x_target,
                    profile_key=profile_key, m_elements=m_elements, etchants=etchants, a_element=a_element,
                    expected_defects=expected_defects, dopants=dopants, nominal_metal_stoich=nominal,
                    parent_composition=parent or None, post_treatment=post or None,
                    termination_family=term_family or None, cl_interest=cl_interest, notes=notes or None)
profile = KB.profile(profile_key)
models = KB.fit_models(profile_key)
kb_roles = {(core, c.name): c.role for core, comps in models.items() for c in comps}
for d in dopants:
    core_d, _ = metal_core(d)
    kb_roles.update({(core_d, f"{d}_lattice_or_vacancy"): "dopant_lattice_candidate",
                     (core_d, f"{d}_adsorbed"): "dopant_adsorbed", (core_d, f"{d}_oxide"): "dopant_oxide"})

st.title(APP_TITLE)
st.caption(ui("MAX 相・MXene の XPS：測定依頼の作成、業者データの読込と比較、独立フィットによる検証、欠損・終端の定量。"
              "結合エネルギーから化学状態を一対一で決めつけない設計です。",
              "XPS for MAX phases and MXenes: measurement requests, provider-data import and comparison, independent "
              "refits for QC, and defect/termination-aware quantification. Binding energy is not treated as a one-to-one "
              "chemical-state lookup."))
st.caption(ui("検証の状況：リリース候補版。校正に使ったのは 1 データセットのみで、独立データでの検証はまだです。"
              "フィットに基づく値は暫定値として扱ってください。",
              "Validation status: release candidate. Calibrated on one development dataset; prospective hold-out "
              "validation not yet performed. Treat fit-based results as provisional."))

_started = any(st.session_state.get(k) for k in ("quant", "fits", "xs_samples"))
with st.expander(ui("📖 使い方", "📖 How to use"), expanded=not _started):
    st.markdown(ui(
        "まず左の **① 試料の設定** で試料の構造（例：(Nb,Ti)2CTx）を選び、目的に合うタブへ進んでください。",
        "First choose the sample structure (e.g. (Nb,Ti)2CTx) under **① Sample setup** on the left, then go to the tab "
        "for what you want to do."))
    st.markdown(ui(
        "| やりたいこと | タブ | 用意するもの |\n|---|---|---|\n"
        "| 業者から届いたデータを試料ごとに並べて見る | **1. 一括読込・試料比較** | 業者の報告書（MultiPak の .pptx か Curve Fit Summary）、同じフォルダの CSV（任意） |\n"
        "| 1 本のスペクトルを自分でフィットし、業者の結果と照合する | **2. フィット検証** | 生データ（`*_multi.csv` か 2 列の CSV）、業者のフィット結果（任意） |\n"
        "| ピーク面積から M サイト欠損・終端基を計算する | **3. 定量** | 成分表の CSV（テンプレートあり） |\n"
        "| 結果をまとめて保存する | **4. レポート** | 1〜3 で出した結果 |\n"
        "| 初期値・帰属の根拠と文献を確かめる | **5. 知識ベース** | なし |\n"
        "| 次の測定・フィットを業者に依頼する | **6. 測定依頼** | なし（依頼文をダウンロード） |",
        "| I want to… | Tab | What you need |\n|---|---|---|\n"
        "| see a provider's data side by side, sample by sample | **1. Import & compare** | provider report (MultiPak .pptx or Curve Fit Summary), the CSV files from the same folder (optional) |\n"
        "| fit one spectrum myself and check it against the provider | **2. Fit check** | raw data (`*_multi.csv` or a two-column CSV), provider fit (optional) |\n"
        "| compute M-site deficiency and terminations from peak areas | **3. Quantification** | component CSV (template provided) |\n"
        "| save everything in one report | **4. Report** | results from tabs 1–3 |\n"
        "| check starting values, assignments and references | **5. Knowledge base** | nothing |\n"
        "| ask a provider for the next measurement and fit | **6. Request** | nothing (download the request text) |"))
    st.markdown(ui(
        "データが手元になくても、タブ 1 と 2 の **合成デモデータで試す** で一通りの動きを確かめられます。"
        "ファイルは業者から届いた名前のまま入れてください（報告書と生データをファイル名で対応付けます）。",
        "Without data of your own, **Try with synthetic demo data** in tabs 1 and 2 shows the whole flow. Keep the "
        "provider's file names unchanged (the report and the raw files are matched by name)."))

ctx_warn = metadata_warnings(md)
if ctx_warn:
    with st.expander(ui("⚠ 測定条件についての注意（{n} 件）", "⚠ Measurement-context warnings ({n})", n=len(ctx_warn)),
                     expanded=md.sputtered):
        for w in ctx_warn:
            say(w)

tab_x, tab_v, tab_q, tab_rep, tab_kb, tab_req = st.tabs(
    [ui("1. 一括読込・試料比較", "1. Import & compare"), ui("2. フィット検証", "2. Fit check"),
     ui("3. 定量", "3. Quantification"), ui("4. レポート", "4. Report"),
     ui("5. 知識ベース", "5. Knowledge base"), ui("6. 測定依頼", "6. Request")])


def steps(ja: str, en: str) -> None:
    """The short 'what this tab does' box at the top of each tab."""
    st.info(ui(ja, en), icon="🧭")


# ------------------------------------------------------------------ 1 request
with tab_req:
    steps("① サイドバーで試料の構造と薬品を選ぶ → ② 測定すべき領域を確認 → ③ 依頼文を必要なら直してダウンロードし、業者に送る。",
          "① Set the structure and chemistry in the sidebar → ② check the regions to measure → ③ edit the request if "
          "needed, download it and send it to the provider.")
    req_text, regions = generate_request(setup, KB)
    st.subheader(ui("測定すべき領域", "Recommended regions"))
    st.dataframe(regions, hide_index=True)
    st.subheader(ui("測定・フィットの依頼文", "Measurement / fitting request"))
    st.caption(ui("そのまま編集できます。依頼文は英語です（業者に渡す文書のため）。各成分の初期値は知識ベースから入ります。",
                  "Editable. The request is written in English. Candidate fit numbers come from the active knowledge-base profile."))
    edited_req = st.text_area(ui("依頼文", "Request text"), req_text, height=600, key=f"req_{hash(req_text)}")
    st.download_button(ui("依頼文をダウンロード（.txt）", "Download request (.txt)"), edited_req,
                       file_name=f"{sample_name}_XPS_request.txt", mime="text/plain", type="primary")

# ------------------------------------------------------------------ 2 quant
with tab_q:
    steps("① テンプレート（成分表 CSV）をダウンロード → ② 業者のフィット結果の面積を `area` 列に入れてアップロード"
          "（画面の表に直接入力しても可） → ③ 面積の種類と感度係数（RSF）を設定 → ④「計算」。",
          "① Download the component CSV template → ② fill the `area` column from the provider fit and upload it (or "
          "type into the table) → ③ set the area type and RSFs → ④ press Calculate.")
    use_legacy = False
    if structure == "Ti3C2Tx":
        _qm = {"role": ui("役割ベース（Brette2024 プロファイル）", "Role-based (Brette2024 profile)"),
               "legacy": ui("研究室の旧 Excel モデル（Ti3C2Tx）", "Legacy Ti3C2Tx workbook model")}
        use_legacy = st.radio(ui("定量モデル", "Quantification model"), list(_qm), horizontal=True, key="q_model",
                              format_func=_qm.get) == "legacy"

    if use_legacy:
        say(LEGACY_LABEL, "INFO")
        tmpl = legacy_template(sample_name)
        st.download_button(ui("旧モデル用のテンプレートをダウンロード", "Download legacy component template"), csv_bytes(tmpl),
                           file_name=f"{sample_name}_legacy_template.csv")
        up = st.file_uploader(ui("旧モデル用の成分表 CSV（sample, core_level, component, area）",
                                 "Upload legacy component CSV (sample, core_level, component, area)"), type=["csv"], key="leg_up")
        base = read_csv_bytes(up.getvalue()) if up else tmpl
        leg_df = st.data_editor(base, num_rows="dynamic", key=f"leg_ed_{up.name if up else sample_name}")
        c1, c2 = st.columns(2)
        with c1:
            ctrl = {
                "F_oxide_per_Ooxide": st.number_input(ui("F_oxide_per_Ooxide（旧モデル・モデル依存）", "F_oxide_per_Ooxide (legacy, model-dependent)"), value=LEGACY_CONTROLS_DEFAULT["F_oxide_per_Ooxide"], format="%.4f"),
                "Include_O_highBE": int(st.checkbox("Include_O_highBE", True)),
                "Include_Cl": int(st.checkbox("Include_Cl", True)),
                "Include_F1": int(st.checkbox("Include_F1", False)),
                "Include_F2": int(st.checkbox("Include_F2", False)),
            }
        with c2:
            st.caption(ui("RSF は旧 Excel の既定値です。新しいデータでは業者に確認してください。",
                          "Historical RSFs are legacy defaults only – confirm with the provider for new data."))
            rsf_df = st.data_editor(pd.DataFrame({"core_level": list(KB.legacy["default_rsf"]),
                                                  "RSF": list(KB.legacy["default_rsf"].values())}),
                                    hide_index=True, key="leg_rsf")
        rsf = {r.core_level: float(r.RSF) for r in rsf_df.itertuples() if pd.notna(r.RSF)}
        if st.button(ui("計算（旧モデル）", "Calculate (legacy model)"), type="primary", key="btn_calc_legacy"):
            try:
                res, warns = legacy_quantify_table(leg_df, rsf, ctrl)
                st.session_state["quant"] = {"results": res, "warnings": warns, "model": "legacy", "profile": profile_key, "options": ctrl,
                                             "rsf": rsf, "table": leg_df, "reclass": None}
            except Exception as e:
                say(str(e), "INVALID")
    else:
        st.subheader(ui("成分表（業者のフィット結果）", "Provider component table"))
        st.markdown(ui("必須の列：`sample, core_level, component, area, role`。あるとよい列：`element, BE, FWHM, notes`。"
                       "`area` が空欄の成分は「フィットしていない」とみなします。テンプレートには試料の構造に合った成分が並んでいます。",
                       "Required columns: `sample, core_level, component, area, role`. Recommended: "
                       "`element, BE, FWHM, notes`. Blank areas = component not fitted. The template lists the components "
                       "of the chosen structure."))
        tmpl = make_peak_template(setup, KB)
        c_t, c_u = st.columns([1, 2])
        c_t.download_button(ui("テンプレートをダウンロード", "Download component CSV template"), csv_bytes(tmpl),
                            file_name=f"{sample_name}_peak_template.csv")
        up = c_u.file_uploader(ui("成分表 CSV をアップロード", "Upload fitted component CSV"), type=["csv"], key="peak_up")
        base = read_csv_bytes(up.getvalue()) if up else tmpl
        ed_key = f"peak_ed_{up.name if up else ''}_{profile_key}_{sample_name}_{','.join(dopants)}"
        peaks = st.data_editor(
            base, num_rows="dynamic", key=ed_key,
            column_config={"role": st.column_config.SelectboxColumn("role", options=list(ROLES), required=True,
                                                                    help=ui("定量での役割。変更は記録されます。", "Quantification role (§6). Changing a role is logged.")),
                           "area": st.column_config.NumberColumn("area", min_value=0.0)})

        reclass = role_reclassifications(peaks, kb_roles) if {"core_level", "component", "role", "sample"} <= set(peaks.columns) else pd.DataFrame()
        if len(reclass):
            say("CAUTION", ui("知識ベースの既定から役割を変えた成分があります。利用者による割り当てとして出力に記録されます。"
                              "それぞれに根拠（金属内殻との照合や他の測定）が必要です。",
                              "Roles changed from the knowledge-base default. These are user assignments and are recorded in "
                              "the export; each must be justified (e.g. by metal-core cross-check or orthogonal evidence)."))
            st.dataframe(reclass, hide_index=True)

        st.subheader(ui("定量の設定", "Quantification controls"))
        c1, c2, c3 = st.columns(3)
        with c1:
            _im = {"raw_area": ui("生の面積（RSF で補正する）", "raw_area (corrected with RSFs)"),
                   "corrected_area": ui("補正済みの面積", "corrected_area"), "atomic_percent": ui("原子 %", "atomic_percent")}
            input_mode = st.selectbox(ui("面積の種類", "Area input mode"), list(_im), key="q_mode", format_func=_im.get)
            inc_dis = st.checkbox(ui("乱れた O 終端（termination_O_disordered）を Tx に数える", "Count disordered O (termination_O_disordered) as Tx"),
                                  value=True, key="q_dis")
            inc_oh = st.checkbox(ui("OH 様の O を Tx に数える", "Count OH-like O as Tx"), value=True, key="q_oh")
            st.caption(ui("O_unresolved_low / O_unresolved_mid は Tx に数えません。照合の根拠があるときだけ、成分表で役割を明示的に変えてください。",
                          "O_unresolved_low / O_unresolved_mid are never counted as Tx; reclassify a component's role "
                          "explicitly if a cross-check justifies it."))
        with c2:
            inc_cl = st.checkbox(ui("Cl を Tx に数える", "Count Cl as Tx"), value=True, key="q_cl")
            inc_fres = st.checkbox(ui("高 BE 側の F 残渣も Tx に数える（例外設定）", "Count high-BE F residue as Tx (override)"), value=False, key="q_fres")
            fcorr = st.number_input(ui("酸化物由来の F 補正係数（0 = 補正なし）", "Oxide-derived F correction factor (0 = off)"), min_value=0.0,
                                    value=0.6667 if structure == "Ti3C2Tx" else 0.0, step=0.05, format="%.4f", key=f"q_fcorr_{structure}",
                                    help=ui("0.6667 は旧 Ti3C2Tx モデルの値です。Nb/Ti 混合系には既定では使いません。",
                                            "0.6667 is the legacy Ti3C2Tx model; not transferred to mixed Nb/Ti by default."))
            tol = st.number_input(ui("金属 QC の許容幅 |M_target − M_norm|（ソフトの既定値）", "Metal-QC tolerance |M_target − M_norm| (software default)"),
                                  min_value=0.0, value=0.35, step=0.05, key="q_tol")
        with c3:
            cores = sorted(set(peaks["core_level"].dropna().astype(str))) if "core_level" in peaks else []
            prefill = st.checkbox(ui("研究室の旧 RSF を入れる（要確認）", "Prefill legacy laboratory RSFs where available (must be confirmed)"),
                                  value=False, disabled=input_mode != "raw_area", key="q_prefill")
            leg = KB.legacy["default_rsf"]
            rsf_in = pd.DataFrame({"core_level": cores,
                                   "RSF": [leg.get(c, np.nan) if prefill else np.nan for c in cores]})
            rsf_df = st.data_editor(rsf_in, hide_index=True, key=f"rsf_{prefill}_{'|'.join(cores)}",
                                    disabled=input_mode != "raw_area")
            rsf = {r.core_level: float(r.RSF) for r in rsf_df.itertuples() if pd.notna(r.RSF)}
            if input_mode == "raw_area":
                miss = [c for c in cores if c not in rsf]
                if miss:
                    say("CAUTION", ui("計算の前に RSF を入れてください（1.0 で代用はしません）：{m}",
                                      "RSF required before calculation (no 1.0 default): {m}", m=", ".join(miss)))

        with st.expander(ui("XPS 以外の裏付け（ICP・STEM など。表示のみで点数化はしない）", "Orthogonal evidence (§2.2, §24.4) – displayed, not scored")):
            if "evidence" not in st.session_state:
                st.session_state["evidence"] = pd.DataFrame({"method": list(ORTHOGONAL_METHODS), "available": False,
                                                             "result_summary": "", "supports": ""})
            ev_df = st.data_editor(st.session_state["evidence"], hide_index=True, num_rows="dynamic", key="ev_ed",
                                   column_config={"supports": st.column_config.TextColumn(
                                       "supports", help=ui("コンマ区切り（例：Ti depletion, M vacancy, dopant incorporation）",
                                                           "comma-separated, e.g. Ti depletion, M vacancy, dopant incorporation"))})
            st.session_state["evidence_current"] = ev_df

        if st.button(ui("計算", "Calculate"), type="primary", key="btn_calc"):
            opt = QuantOptions(m_target=m_target, x_target=x_target, input_mode=input_mode, include_o_disordered=inc_dis,
                               include_oh_like=inc_oh, include_cl=inc_cl, include_f_residue=inc_fres,
                               oxide_f_correction_factor=fcorr, metal_qc_tolerance=tol, nominal_metal_stoich=nominal,
                               nb_suboxide_validated=bool(KB.raw.get("calibration_status", {}).get("nb3d_suboxide_like_validated", False)),
                               phase=phase, a_target=a_target,
                               provisional_reasons=tuple(KB.raw.get("calibration_status", {}).get("delta_finalization_blocked_by", [])))
            try:
                res, warns = quantify_table(peaks, rsf, opt)
                if len(reclass):
                    for s in warns:
                        warns[s].append(f"{int((reclass['sample'] == s).sum())} component role(s) were reclassified by the user (see reclassification log).")
                st.session_state["quant"] = {"results": res, "warnings": warns, "model": "role-based", "profile": profile_key,
                                             "options": opt.__dict__, "rsf": rsf, "table": peaks, "reclass": reclass}
            except Exception as e:
                say(str(e), "INVALID")

    q = st.session_state.get("quant")
    if q and q["model"] == ("legacy" if use_legacy else "role-based") and q.get("profile") == profile_key:
        st.subheader(ui("結果", "Results"))
        res = q["results"]
        show = (res.set_index("sample").T
                .map(lambda x: f"{x:.6g}" if isinstance(x, (float, np.floating)) else x)
                .reset_index().rename(columns={"index": "quantity"}))
        st.dataframe(show, hide_index=True, height=min(900, 36 * (len(show) + 1)))
        st.download_button(ui("結果をダウンロード（.csv）", "Download quantification results (.csv)"), csv_bytes(res.assign(provenance="derived_quantity",
                                                                                         app_version=__version__, kb_version=KB.version)), file_name=f"{sample_name}_quantification.csv")
        if q["model"] == "role-based":
            st.markdown(ui("**基準値と代替値（暫定値）** — 基準値は校正候補の成分を骨格に数えず、代替値は数えます。両者は平均しません。",
                           "**Nominal vs alternative (provisional quantities)** — default excludes the calibration candidate "
                           "from the framework; the alternative counts it. The two are never averaged."))
            qopt = QuantOptions(**q["options"])
            prov_tabs = []
            for rec in res.to_dict("records"):
                pr = provisional_report(rec, qopt)
                pr.insert(0, "sample", rec["sample"])
                prov_tabs.append(pr)
            prov_df = pd.concat(prov_tabs, ignore_index=True) if prov_tabs else pd.DataFrame()
            st.dataframe(prov_df, hide_index=True)
            q["provisional"] = prov_df
        for s, ws in q["warnings"].items():
            for w in ws:
                say(f"{s}: {w}")
        if q["model"] == "role-based" and "apparent_M_deficiency_status" in res and (res["apparent_M_deficiency_status"] == "provisional").any():
            say("CAUTION", ui("apparent_M_deficiency_delta は暫定値です（Nb 3d ~205.6–205.9 eV の候補成分が未検証）。"
                              "`apparent_M_deficiency_delta_alt`（候補を骨格に数えた値）と比べ、どちらも確定値として報告しないでください。",
                              "apparent_M_deficiency_delta is PROVISIONAL (Nb 3d ~205.6–205.9 eV candidate not validated). "
                              "Compare with `apparent_M_deficiency_delta_alt` (candidate counted as framework); do not report "
                              "either as final."))
        if q["model"] == "role-based":
            say("INFO", ui("`apparent_M_deficiency_delta` は炭化物 C で規格化した表面敏感な XPS の指標で、空孔の絶対濃度ではありません。"
                           "選択溶出は `<M>_fraction_shift` / `<M>_retention_ratio_fractional` で別に示します。",
                           "`apparent_M_deficiency_delta` is a surface-sensitive XPS index after carbide-C normalization, not "
                           "an absolute vacancy concentration. Selective depletion is reported separately via "
                           "`<M>_fraction_shift` / `<M>_retention_ratio_fractional`."))
            ev_lines, ev_status = evidence_summary([
                OrthogonalEvidence(r.method, bool(r.available), r.result_summary or None, split_list(str(r.supports or "")))
                for r in st.session_state.get("evidence_current", pd.DataFrame(columns=["method", "available", "result_summary", "supports"])).itertuples()])
            st.caption(ev_status)
            for line in ev_lines:
                st.write("• " + line)

# ------------------------------------------------------------------ 3 validation
with tab_v:
    steps("① コア準位を選ぶ → ② 生データを入れる（`*_multi.csv` ならファイル内の領域も選ぶ） → ③ 業者のフィット結果を入れる（任意） "
          "→ ④「Python でフィット」。業者の結果を入れると、成分ごとと役割ごとに一致度（PASS / REVIEW / FAIL）を判定します。"
          "Python のフィットは独立した検算で、正解ではありません。",
          "① Choose the core level → ② load the raw data (for `*_multi.csv`, also pick the region) → ③ load the provider "
          "fit (optional) → ④ Run Python fit. With a provider fit, agreement is judged per component and per role "
          "(PASS / REVIEW / FAIL). The Python fit is an independent check, not ground truth.")
    cores = list(models) + [r for r in UNSUPPORTED_REGIONS if r not in models]
    _core_lbl = {c: c if c in models else ui("{c}（検証済みモデルなし）", "{c} ({t})", c=c, t=NO_MODEL_TEXT) for c in cores}
    core = st.selectbox(ui("① コア準位", "① Core level"), cores, key="v_core", format_func=_core_lbl.get)
    kb_comps = models.get(core, [])
    unsupported = not kb_comps
    if unsupported:
        say("CAUTION", ui("{c}：検証済みのモデルがありません。読み込み・表示・出力と業者ピークの取り込みはできますが、Python のフィットは行いません。",
                          "{c}: {t}. The region can be imported, displayed and exported, and provider peaks can be "
                          "imported, but no Python fit is run.", c=core, t=NO_MODEL_TEXT) + " " + UNSUPPORTED_REGIONS.get(core, ""))
    for w in profile.get("core_level_warnings", {}).get(core, []):
        say(w, "CAUTION")
    importer_names = list(MULTIREGION_IMPORTERS) + list(SPECTRUM_IMPORTERS)
    importer_label = {list(MULTIREGION_IMPORTERS)[0]: ui("PHI 複数領域・サーベイ（*_multi.csv, *_su.csv）", "PHI multi-region / survey (*_multi.csv, *_su.csv)"),
                      list(SPECTRUM_IMPORTERS)[0]: ui("2 列（結合エネルギー, 強度）", "generic two-column (BE, intensity)")}
    c1, c2 = st.columns(2)
    with c1:
        st.markdown(ui("**② 生データ**", "**② Raw spectrum**"))
        importer = st.selectbox(ui("形式", "Format"), importer_names, key="v_importer", format_func=importer_label.get)
        raw_file = st.file_uploader(ui("生データのファイル", "Raw spectrum file"), type=["csv", "txt", "dat"], key="raw",
                                    help=ui("2 列の場合：見出しの有無・並び順・区切り文字は問いません。",
                                            "Two-column files: header optional, any order, any delimiter."))
        region_sp = None
        if raw_file is not None and importer in MULTIREGION_IMPORTERS:
            try:
                regs = MULTIREGION_IMPORTERS[importer](raw_file.getvalue())
                labels = [f"{r.region} → {r.core_level or ui('内殻ではない', 'not a core level')} ({r.be.min():.1f}–{r.be.max():.1f} eV, "
                          f"{len(r.be)} pts)" for r in regs]
                default = next((i for i, r in enumerate(regs) if r.core_level == core), 0)
                pick = st.selectbox(ui("ファイル内の領域", "Region in file"), list(range(len(regs))), index=default,
                                    format_func=lambda i: labels[i], key=f"region_{raw_file.name}_{core}")
                region_sp = regs[pick]
                if region_sp.core_level != core:
                    say("INVALID", ui("選んだ領域 {r} は {c} ではないため、{c} のモデルではフィットしません。",
                                      "Selected region {r} is not {c}; fitting it with the {c} model is refused.", r=region_sp.region, c=core))
            except Exception as e:
                say("INVALID", ui("複数領域ファイルを読めませんでした：{e}", "Multi-region import failed: {e}", e=e))
        demo_tag = {"Ti 2p": "Ti2p", "Nb 3d": "Nb3d"}.get(core)
        use_demo = st.checkbox(ui("合成デモデータで試す（Ti 2p・Nb 3d のみ。生データ未指定のとき）",
                                  "Try with synthetic demo spectrum (Ti 2p, Nb 3d; used when no raw file is loaded)"),
                               value=False, disabled=demo_tag is None, key="v_demo")
    with c2:
        st.markdown(ui("**③ 業者のフィット結果（任意）**", "**③ Provider fit (optional)**"))
        _pf = {"PHI MultiPak Curve Fit Summary": ui("PHI MultiPak（.pptx / Curve Fit Summary）", "PHI MultiPak (.pptx / Curve Fit Summary)"),
               "Generic component CSV": ui("成分表 CSV", "Component CSV")}
        prov_fmt = st.radio(ui("形式", "Format"), list(_pf), horizontal=True, key="v_prov_fmt", format_func=_pf.get)
        if prov_fmt == "Generic component CSV":
            vendor_file = st.file_uploader(ui("成分表 CSV（component, BE, FWHM, area[, role, core_level]）",
                                              "Provider component CSV (component, BE, FWHM, area[, role, core_level])"), type=["csv"], key="vend")
            _dc = {"total": ui("2 本の合計", "total (both spin-orbit lines)"), "primary": ui("主線 1 本のみ", "primary line only")}
            doublet_conv = st.radio(ui("ダブレットの面積は…", "Provider doublet area is…"), list(_dc), horizontal=True, key="v_dconv",
                                    format_func=_dc.get)
        else:
            vendor_file = st.file_uploader(ui("MultiPak の報告書（.pptx）または Curve Fit Summary（.txt / .csv）",
                                              "MultiPak report (.pptx) or Curve Fit Summary text (.txt/.csv)"),
                                           type=["pptx", "txt", "csv"], key="vend_mp")
            doublet_conv = "total"
            st.caption(ui("別々に報告されたスピン軌道の 2 本は、領域・分裂幅・面積比・近さで対にします。業者のラベルはそのまま残し、"
                          "知識ベースへの対応付けは仮のもので、下の表で直せます。",
                          "Spin-orbit lines reported as separate bands are paired by region, splitting, area ratio and "
                          "proximity. Provider labels are kept; the KB mapping below is tentative and editable."))
    with st.expander(ui("フィットの条件（バックグラウンド・範囲）", "Fit settings (background, window)"), expanded=False):
        f1, f2 = st.columns(2)
        bg = f1.selectbox(ui("バックグラウンド", "Background"), ["shirley", "linear"], key="v_bg",
                          help=ui("Shirley は未校正の近似です。", "The Shirley background is an uncalibrated approximation."))
        endavg = f2.number_input(ui("端点の平均点数", "Background endpoint averaging (points)"), min_value=1, value=1, step=1, key="v_endavg")
        restrict = st.checkbox(ui("フィットする結合エネルギーの範囲を絞る", "Restrict fit BE window"), value=False, key="v_restrict",
                               help=ui("隣の内殻を除くときなど。設定はフィット条件に記録されます。",
                                       "E.g. exclude a neighbouring core level. Recorded in the fit settings."))
        win = None
        if restrict and kb_comps:
            w1, w2 = st.columns(2)
            win = (w1.number_input(ui("低 BE 側（eV）", "low BE (eV)"), value=float(min(x.bounds[0] for x in models[core]) - 3.0),
                                   format="%.2f", key=f"v_lo_{core}"),
                   w2.number_input(ui("高 BE 側（eV）", "high BE (eV)"), value=float(max(x.bounds[1] + (x.split or 0) for x in models[core]) + 3.0),
                                   format="%.2f", key=f"v_hi_{core}"))
    with st.expander(ui("QC の判定基準（ソフトの既定値で、物理定数ではない）", "QC thresholds (software defaults, not physical constants)")):
        vd = KB.validation_defaults
        a, b, c = st.columns(3)
        th = {
            "delta_be_pass_eV": a.number_input("|ΔBE| PASS ≤ (eV)", 0.0, value=float(vd["delta_be_pass_eV"]), step=0.05, key="th1"),
            "delta_be_review_eV": a.number_input("|ΔBE| REVIEW ≤ (eV)", 0.0, value=float(vd["delta_be_review_eV"]), step=0.05, key="th2"),
            "fwhm_rel_pass": b.number_input("rel. ΔFWHM PASS ≤", 0.0, value=float(vd["fwhm_rel_pass"]), step=0.05, key="th3"),
            "fwhm_rel_review": b.number_input("rel. ΔFWHM REVIEW ≤", 0.0, value=float(vd["fwhm_rel_review"]), step=0.05, key="th4"),
            "area_fraction_abs_pass": c.number_input("Δarea-fraction PASS ≤", 0.0, value=float(vd["area_fraction_abs_pass"]), step=0.02, key="th5"),
            "area_fraction_abs_review": c.number_input("Δarea-fraction REVIEW ≤", 0.0, value=float(vd["area_fraction_abs_review"]), step=0.02, key="th6"),
        }
        corr_th = st.number_input(ui("面積の相関で警告する |r| ≥（同定可能性の診断）", "Area-correlation flag |r| ≥ (identifiability diagnostic)"),
                                  0.0, 1.0, 0.8, 0.05, key="th_corr")
        st.caption(vd.get("note", ""))
    if kb_comps:
        with st.expander(ui("開始モデル（初期値と範囲）を確認・編集", "Starting model (start values and windows): view / edit"), expanded=False):
            st.markdown(ui("初期値とフィット範囲は編集できます（変更は記録されます）。文献値は参考表示だけで、拘束には使いません。",
                           "Starting values and fitting windows are editable implementation bounds (edits are logged); "
                           "literature reference values are shown separately and are never used as constraints."))

            model_df = pd.DataFrame([{"name": x.name, "kind": x.kind, "role": x.role, "status": x.assignment_status,
                                      "center": x.center, "BE_low": x.bounds[0], "BE_high": x.bounds[1],
                                      "FWHM": x.fwhm, "FWHM_low": x.fwhm_bounds[0], "FWHM_high": x.fwhm_bounds[1],
                                      "ΔSO": x.split, "ratio": x.ratio_primary_secondary, "ΔSO basis": x.split_basis or "",
                                      "calibration": _cal_text(x),
                                      "literature reference (display only)": x.literature_summary(),
                                      "transferred reference (display only)": x.transferred_summary(),
                                      "provider observation (display only)": x.observation_summary(),
                                      "basis": x.source_type or "", "start provenance": start_value_provenance(x), "source": x.source,
                                      "DOI": KB.references[x.source].get("doi")} for x in kb_comps])
            ver = st.session_state.setdefault("model_reset", 0)
            edited_model = st.data_editor(
                model_df, hide_index=True, key=f"model_{profile_key}_{core}_{ver}",
                disabled=["name", "kind", "role", "status", "ΔSO", "ratio", "ΔSO basis", "calibration", "literature reference (display only)",
                          "transferred reference (display only)", "provider observation (display only)", "basis", "start provenance", "source", "DOI"])
            if st.button(ui("知識ベースの値に戻す", "Reset model to knowledge-base values"), key="btn_reset_model"):
                st.session_state["model_reset"] = ver + 1
                st.rerun()
        overrides = {r["name"]: {"center": r["center"], "bounds": (r["BE_low"], r["BE_high"]), "fwhm": r["FWHM"],
                                 "fwhm_bounds": (r["FWHM_low"], r["FWHM_high"])} for r in edited_model.to_dict("records")}
        try:
            comps, model_log = with_overrides(kb_comps, overrides)
        except ValueError as e:
            say("INVALID", ui("開始モデルの編集が無効です：{e}", "Invalid model edit: {e}", e=e))
            comps, model_log = kb_comps, []
        if model_log:
            say("CAUTION", ui("知識ベース v{k} から初期値・範囲を {n} か所変えています。フィット条件と出力に記録されます。",
                              "{n} starting value/window edit(s) relative to KB v{k}; recorded in the fit settings and the project export.",
                              n=len(model_log), k=KB.version))
            st.dataframe(pd.DataFrame(model_log).astype(str), hide_index=True)

    else:
        comps, model_log = [], []

    mp_state = None
    if vendor_file is not None and prov_fmt != "Generic component CSV":
        try:
            mp_bands = read_multipak_provider(vendor_file.getvalue(), vendor_file.name)
            keys = (mp_bands.drop_duplicates("summary_index")
                    .assign(lbl=lambda d: d["comment"].astype(str) + " | " + d["region"].astype(str) + " | "
                            + d["data_file"].astype(str) + " | " + d["origin"].astype(str)))
            cands = keys[keys["core_level"] == core]
            st.markdown(ui("**業者のフィット（MultiPak）** — ファイル内 {n} 件、うち {c} は {m} 件",
                           "**Provider fit (MultiPak)** — {n} summaries in file, {m} for {c}", n=len(keys), m=len(cands), c=core))
            if len(cands):
                opts = list(cands["summary_index"])
                raw_name = raw_file.name if raw_file is not None else ""
                stems = [str(f).rsplit(".", 1)[0] for f in cands["data_file"]]
                default = next((i for i, (st_, cm) in enumerate(zip(stems, cands["comment"].astype(str)))
                                if (st_ and st_ in raw_name) or (cm and cm == sample_name)), 0)
                sidx = st.selectbox(ui("照合する業者のフィット", "Provider summary"), opts, index=default,
                                    format_func=lambda i: cands.set_index("summary_index").at[i, "lbl"],
                                    key=f"mp_sum_{core}_{raw_name}")
                if raw_name and not any(st_ and st_ in raw_name for st_ in stems):
                    say("WARNING", ui("生データのファイル名に合う業者のフィットがありません。選んだフィットが同じ試料・領域のものか確認してください。",
                                      "No provider summary matches the raw-file name; confirm the selected summary belongs to "
                                      "the same sample/region."))
                one = mp_bands[mp_bands["summary_index"] == sidx]
                doublets, linked = pair_spin_orbit(one)
                auto = map_to_kb(doublets, comps)
                with st.expander(ui("業者のフィットピーク（報告どおり・未加工）", "Provider raw fitted peaks (as reported; never modified)")):
                    st.dataframe(linked.drop(columns=["origin"]), hide_index=True)
                st.markdown(ui("スピン軌道の対と、知識ベースの成分への仮の対応（`mapped_component` を直せます。変更は記録されます）",
                               "Spin-orbit pairing and tentative KB mapping (edit `mapped_component`; edits are logged)"))
                ed = st.data_editor(
                    auto.table[["doublet_id", "provider_label", "pairing", "area_convention", "BE", "split_obs", "ratio_obs",
                                "FWHM", "area_reported", "area_total_estimated", "mapped_component", "mapped_role",
                                "mapping_method"]],
                    hide_index=True, key=f"mp_map_{sidx}_{profile_key}",
                    disabled=["doublet_id", "provider_label", "pairing", "area_convention", "BE", "split_obs", "ratio_obs",
                              "FWHM", "area_reported", "area_total_estimated", "mapped_role", "mapping_method"],
                    column_config={"mapped_component": st.column_config.SelectboxColumn(
                        "mapped_component", options=[x.name for x in comps], required=False)})
                edits = {r["doublet_id"]: (r["mapped_component"] or None) for r in ed.to_dict("records")}
                mapping, map_log = apply_user_mapping(auto.table, edits, comps)
                for n in auto.notes:
                    say(n, "CAUTION")
                est = mapping[mapping["area_convention"] == "primary_line_only"]
                if len(est):
                    say("CAUTION", ui("スピン軌道の 1 本だけが報告された成分：{l}。合計面積は縮退比から推定しています。",
                                      "Only one spin-orbit line was reported for: {l}. Total areas are ESTIMATED with the degeneracy ratio.",
                                      l=", ".join(est["provider_label"])))
                chk = mapping[mapping["area_convention"] == "unpaired_check"]
                if len(chk):
                    say("WARNING", ui("対になっていないピークがダブレットの隣にあります：{l}。確認してください。",
                                      "Unpaired band(s) next to paired doublets: {l} — check.", l=", ".join(chk["provider_label"])))
                if map_log:
                    say("CAUTION", ui("自動の仮対応から {n} か所変えています（記録済み）。",
                                      "{n} user mapping change(s) relative to the automatic tentative mapping (logged).", n=len(map_log)))
                    st.dataframe(pd.DataFrame(map_log).astype(str), hide_index=True)
                mp_state = {"bands": linked, "mapping": mapping, "mapping_log": map_log,
                            "vendor": provider_component_table(mapping)}
        except Exception as e:
            say("INVALID", ui("MultiPak のファイルを読めませんでした：{e}", "MultiPak import failed: {e}", e=e))

    if unsupported:
        sp_show = region_sp if region_sp is not None else None
        if sp_show is None and raw_file is not None and importer in SPECTRUM_IMPORTERS:
            try:
                sp_show = SPECTRUM_IMPORTERS[importer](raw_file.getvalue())
            except Exception as e:
                say("INVALID", ui("読み込めませんでした：{e}", "Import failed: {e}", e=e))
        if sp_show is not None:
            fig, ax = plt.subplots(figsize=(8, 3.5))
            ax.plot(sp_show.be, sp_show.intensity, ".-", ms=2, lw=0.6, color="0.3")
            ax.invert_xaxis(); ax.set_xlabel("Binding energy (eV)"); ax.set_ylabel("Intensity")
            ax.set_title(f"{core} — imported spectrum ({NO_MODEL_TEXT})")
            st.pyplot(fig); plt.close(fig)
            st.download_button(ui("読み込んだスペクトル（.csv）", "Imported spectrum .csv"), csv_bytes(pd.DataFrame({"BE": sp_show.be, "intensity": sp_show.intensity,
                                                                                 "provenance": "raw_measurement"})),
                               file_name=f"{sample_name}_{core.replace(' ', '')}_raw.csv")
        if mp_state is not None:
            st.download_button(ui("業者のピーク（対にしたもの、.csv）", "Provider peaks (paired) .csv"), csv_bytes(mp_state["mapping"].assign(provenance="provider_paired_doublet")),
                               file_name=f"{sample_name}_{core.replace(' ', '')}_provider_peaks.csv")

    if st.button(ui("④ Python でフィット", "④ Run Python fit"), type="primary", disabled=unsupported, key="btn_fit"):
        try:
            if region_sp is not None:
                sp = region_sp
            elif raw_file is not None and importer in MULTIREGION_IMPORTERS:
                raise ValueError(ui("複数領域ファイルの中から有効な領域を選んでください。", "Select a valid region of the multi-region file."))
            elif raw_file is not None:
                sp = SPECTRUM_IMPORTERS[importer](raw_file.getvalue())
            elif use_demo and demo_tag:
                sp = list(SPECTRUM_IMPORTERS.values())[0]((FIX / f"synthetic_{demo_tag}_raw.csv").read_bytes())
            else:
                raise ValueError(ui("生データを入れるか、合成デモデータにチェックを入れてください。", "Upload a raw spectrum or tick the demo option."))
            if win is not None:
                sp = crop(sp, *win)
            fr = fit_core_level(sp.be, sp.intensity, comps, background=bg, endpoint_avg=int(endavg), corr_threshold=corr_th,
                                spectrum_core_level=sp.core_level)
            st.session_state.setdefault("fits", {})[f"{profile_key} | {core}"] = {"fit": fr, "import_notes": sp.notes,
                                                             "core": core, "spectrum": (sp.be, sp.intensity, sp.core_level),
                                                             "ident": assess_identifiability(fr),
                                                             "settings": {"background": bg, "endpoint_avg": int(endavg),
                                                                          "profile": profile_key, "kb_version": KB.version,
                                                                          "model_overrides": model_log,
                                                                          "source_region": sp.region,
                                                                          "fit_window_eV": list(win) if win else None}}
        except Exception as e:
            say("INVALID", ui("フィットできませんでした：{e}", "Fit failed: {e}", e=e))

    fitstate = st.session_state.get("fits", {}).get(f"{profile_key} | {core}")
    if fitstate:
        fr = fitstate["fit"]
        for n in fitstate["import_notes"]:
            st.caption(ui("読み込み：", "Import: ") + n)
        m = fr.metrics
        st.subheader(ui("フィットの結果", "Fit result"))
        k1, k2, k3, k4, k5 = st.columns(5)
        k1.metric(ui("収束", "converged"), ui("はい", "yes") if m["success"] else ui("いいえ", "no"))
        k2.metric("R²", f"{m['R2']:.5f}")
        k3.metric("RMSE", f"{m['RMSE']:.4g}")
        k4.metric("AIC", f"{m['AIC']:.1f}")
        k5.metric("BIC", f"{m['BIC']:.1f}")
        st.dataframe(style_status(fr.components.assign(provenance="python_fit")), hide_index=True)
        ident = fitstate.get("ident")
        if ident is not None:
            st.markdown(ui("**同定可能性**（統計的な状態のみ。化学的な帰属ではありません）",
                           "**Identifiability / fit confidence** (statistical status only — not a chemical assignment)"))
            st.dataframe(ident.style.format({c: "{:.4g}" for c in ident.columns if pd.api.types.is_float_dtype(ident[c])}, na_rep="")
                         .map(lambda v: {"identifiable": "background-color:#d8f0d8", "weakly identifiable": "background-color:#fff0c2",
                                         "boundary-limited": "background-color:#fde2c2", "non-identifiable": "background-color:#f7d0d0"}.get(v, ""),
                              subset=["status"]), hide_index=True)
            for w in ident_warnings(ident):
                st.caption(f"• [{sev_classify(w)}] {w}")
        if fr.diagnostics:
            for lv in ("INVALID", "WARNING", "CAUTION", "INFO"):
                ds_ = [d for d in fr.diagnostics if sev_classify(d) == lv]
                if ds_:
                    say(lv, ui("フィットの診断：", "Fit diagnostics:") + "\n\n" + "\n".join(f"- {d}" for d in ds_))
        cand = [x for x in comps if x.calibration_candidate and not x.validated]
        if cand:
            say("CAUTION", ui("このモデルには校正候補の成分が含まれます（フィットはしたが未検証）：", "Calibration candidate(s) in this model (fitted, not validated): ") + "; ".join(
                f"{x.name} [{x.role}, confidence {x.confidence}" + (f", cross-check {', '.join(x.cross_check)}" if x.cross_check else "")
                + "]" for x in cand))
        cv = fr.curves
        fig, (ax, axr) = plt.subplots(2, 1, figsize=(9, 6), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
        ax.plot(cv["BE"], cv["observed"], ".", ms=2, color="0.3", label="observed")
        ax.plot(cv["BE"], cv["background"], "--", lw=0.9, color="0.5", label=f"background ({bg})")
        for x in comps:
            ax.plot(cv["BE"], cv[x.name] + cv["background"], lw=0.9, label=f"{x.name} [{x.role}]")
        ax.plot(cv["BE"], cv["fit_total"], lw=1.4, color="crimson", label="Python fit")
        ax.set_ylabel("Intensity"); ax.legend(fontsize=7, ncol=2); ax.set_title(f"{sample_name} — {core}")
        axr.plot(cv["BE"], cv["residual"], lw=0.7, color="0.2"); axr.axhline(0, color="0.6", lw=0.6)
        axr.set_xlabel("Binding energy (eV)"); axr.set_ylabel("resid.")
        ax.invert_xaxis()
        fig.tight_layout()
        st.pyplot(fig)
        dl = st.columns(5)
        for i, fmt in enumerate(("png", "svg", "pdf")):
            buf = io.BytesIO(); fig.savefig(buf, format=fmt, dpi=200)
            dl[i].download_button(ui("図 .{f}", "Plot .{f}", f=fmt), buf.getvalue(), file_name=f"{sample_name}_{core.replace(' ', '')}_fit.{fmt}")
        plt.close(fig)
        dl[3].download_button(ui("成分 .csv", "Fit components .csv"), csv_bytes(fr.components.assign(provenance="python_fit")), file_name=f"{sample_name}_{core.replace(' ', '')}_python_fit.csv")
        dl[4].download_button(ui("曲線 .csv", "Curves .csv"), csv_bytes(cv.assign(provenance="observed=raw_measurement; other columns=python_fit")), file_name=f"{sample_name}_{core.replace(' ', '')}_curves.csv")
        if fr.area_correlation is not None:
            with st.expander(ui("面積パラメータの相関行列", "Area-parameter correlation matrix")):
                st.dataframe(fr.area_correlation.style.format("{:+.2f}"))

        with st.expander(ui("条件を変えたときの感度（あらかじめ決めた代替条件。「最良の化学」は自動で選ばない）",
                            "Fit-sensitivity comparison (predefined alternatives; no automatic 'best chemistry')")):
            be_s, y_s, core_s = fitstate["spectrum"]
            alts = sens.standard_alternatives(core, comps, be_s)
            st.caption(ui("代替条件：", "Alternatives: ") + "; ".join(f"{a.label} [{a.basis}]" for a in alts)
                       + ui("。ΔAIC/ΔBIC はデータ点が基準のフィットと同じ場合だけ表示します。",
                            ". ΔAIC/ΔBIC are shown only where the data points are identical to the nominal fit."))
            qd = st.session_state.get("quant")
            derived_fn, dsample = None, None
            if qd and qd.get("model") == "role-based" and "core_level" in qd["table"]:
                tab = qd["table"]
                samples_q = sorted(set(tab.loc[tab["core_level"] == core, "sample"].astype(str)))
                if samples_q:
                    dsample = st.selectbox(ui("定量タブの成分表から派生量を出す試料", "Derived quantities from the quantification table of sample"),
                                           samples_q, key=f"sens_s_{core}")
                    qopt = QuantOptions(**qd["options"])
                    derived_fn = sens.derived_from_fractions(tab, dsample, core,
                                                             lambda t: quantify_table(t, qd["rsf"], qopt)[0].iloc[0].to_dict())
            if st.button(ui("代替条件でフィット", "Run fitting alternatives"), key=f"sens_btn_{core}"):
                with st.spinner(ui("代替条件でフィット中…", "Fitting alternatives…")):
                    sres = sens.run(be_s, y_s, comps, alts, endpoint_avg=int(endavg), derived_fn=derived_fn, spectrum_core_level=core_s)
                fitstate["sensitivity"] = sres
                fitstate["ident"] = assess_identifiability(fr, alt_spans=sres.component_spans())
                st.rerun()
            sres = fitstate.get("sensitivity")
            if sres is not None:
                st.dataframe(sres.table, hide_index=True)
                st.markdown(ui("代替条件ごとの面積分率", "Area fraction by alternative"))
                st.dataframe(sres.fractions.style.format("{:.3f}"))
                if len(sres.derived):
                    st.markdown(ui("派生量（試料 {s}）", "Derived quantities (sample {s})", s=dsample))
                    st.dataframe(sres.derived)
                st.markdown(ui("代替条件間の幅（0.05 以上は **model-sensitive**）", "Spans across alternatives (≥ 0.05 → **model-sensitive**, D4 materiality criteria)"))
                st.dataframe(sres.spans, hide_index=True)
                ms = sres.spans[sres.spans["model_sensitive"]]["quantity"].tolist() if len(sres.spans) else []
                if ms:
                    say("WARNING", "model-sensitive: " + ", ".join(ms))

        if vendor_file is not None and (prov_fmt == "Generic component CSV" or mp_state is not None):
            try:
                if mp_state is not None:
                    vend = mp_state["vendor"]
                else:
                    vend = read_csv_bytes(vendor_file.getvalue())
                    if "core_level" in vend.columns:
                        vend = vend[vend["core_level"].astype(str).str.strip() == core]
                vp = prepare_vendor(vend, comps, "primary" if doublet_conv.startswith("primary") else "total")
                cqc = compare_components(vp, fr.components, th)
                if "provider_label" in vp.columns:
                    cqc.insert(1, "provider_label", cqc["component"].map(dict(zip(vp["component"], vp["provider_label"]))))
                gqc = compare_grouped(vp, fr.components, th)
                summ = summarize(cqc, gqc)
                say(QC_STATUS_SEVERITY.get(summ.status, "CAUTION"), f"provider-vs-Python QC {summ.status}: {summ.message}")
                for msg in role_disagreements(vp):
                    say(msg, "CAUTION")
                st.markdown(ui("**役割ごとの照合**（骨格金属・酸化物などの合計で比べる。結論はこちらを優先）", "**Role-grouped QC**"))
                st.dataframe(style_status(gqc), hide_index=True)
                st.markdown(ui("**成分ごとの照合**（ΔBE = Python − 業者）", "**Component-level QC** (ΔBE = Python − Provider)"))
                st.dataframe(style_status(cqc), hide_index=True)
                st.caption(IDENTIFIABILITY_WARNING)
                d1, d2 = st.columns(2)
                d1.download_button(ui("成分ごとの照合 .csv", "Component QC .csv"), csv_bytes(cqc.assign(provenance_vendor="provider", provenance_python="python_fit")), file_name=f"{sample_name}_{core.replace(' ', '')}_validation.csv")
                d2.download_button(ui("役割ごとの照合 .csv", "Grouped QC .csv"), csv_bytes(gqc.assign(provenance="derived_quantity")), file_name=f"{sample_name}_{core.replace(' ', '')}_validation_grouped.csv")
                fitstate.update({"component_qc": cqc, "grouped_qc": gqc, "summary": summ.__dict__, "thresholds": th})
                if mp_state is not None:
                    fitstate.update({"provider_raw_bands": mp_state["bands"], "provider_mapping": mp_state["mapping"],
                                     "provider_mapping_log": mp_state["mapping_log"]})
            except Exception as e:
                say("INVALID", ui("業者の結果と照合できませんでした：{e}", "Provider comparison failed: {e}", e=e))

# ------------------------------------------------------------------ 4 KB
with tab_kb:
    steps("このアプリが使う初期値・範囲・帰属の根拠と文献です。内容は英語のままです。数値はここではなく、タブ 2 の「開始モデル」で編集できます。",
          "Starting values, windows, assignment rationale and references used by the app. Values are edited in tab 2 "
          "(Starting model), not here.")
    st.subheader(ui("解釈の原則", "Interpretation principles"))
    for p in KB.raw["principles"]:
        st.markdown(f"- {p}")
    st.subheader(ui("使用中のプロファイル：`{k}`（{m}）", "Active profile: `{k}` ({m})", k=profile_key, m=profile["material"]))
    for n in profile.get("notes", []):
        say(n, "INFO")
    for c_core, cs in models.items():
        st.markdown(f"#### {c_core}")
        st.dataframe(pd.DataFrame([{"name": x.name, "center": x.center, "bounds": f"{x.bounds[0]}–{x.bounds[1]}",
                                    "FWHM": x.fwhm, "FWHM bounds": f"{x.fwhm_bounds[0]}–{x.fwhm_bounds[1]}",
                                    "literature reference (display only)": x.literature_summary(),
                                    "transferred reference (display only)": x.transferred_summary(),
                                    "provider observation (display only)": x.observation_summary(),
                                    "start-value basis": x.source_type or "",
                                    "literature interpretation": x.literature_interpretation or "",
                                    "role": x.role, "status": x.assignment_status, "calibration": _cal_text(x),
                                    "alternative interpretation": x.alternative_interpretation or "",
                                    "interpretation": x.interpretation, "source": x.source,
                                    "reference": KB.reference_text(x.source)} for x in cs]), hide_index=True)
    for c_core, ws in profile.get("core_level_warnings", {}).items():
        for w in ws:
            say("CAUTION", f"{c_core}: {w}")
    if "calibration_status" in KB.raw:
        st.subheader(ui("校正の状況", "Calibration status"))
        st.json(KB.raw["calibration_status"])
    st.subheader(ui("役割の一覧", "Role ontology"))
    st.dataframe(pd.DataFrame({"role": list(ROLE_LABELS), "meaning": list(ROLE_LABELS.values())}), hide_index=True)
    st.subheader(ui("文献", "References"))
    st.dataframe(pd.DataFrame([{"key": k, **v} for k, v in KB.references.items()]), hide_index=True)
    st.subheader("Ti₂₋ᵧNbᵧCTx XPS atomic % (Han2020SolidSolution)")
    ssr = KB.raw["solid_solution_reference"]["Ti2-yNbyCTx"]
    st.caption(ssr["note"] + " XPS in that work was measured after Ar sputtering (2 kV, 2 µA, 2 min).")
    st.dataframe(pd.DataFrame(ssr["xps_atomic_percent"]), hide_index=True)
    st.subheader(ui("仕様の食い違いと決定事項", "Specification discrepancies and decisions"))
    st.dataframe(pd.DataFrame(KB.raw.get("spec_discrepancies", [])), hide_index=True)
    st.caption(f"KB v{KB.version} · created {KB.raw.get('created_at')} · reviewed_by: {KB.raw.get('reviewed_by') or 'not yet reviewed'}")
    st.download_button(ui("knowledge_base.json をダウンロード", "Download knowledge_base.json"), KB.to_json().encode("utf-8"),
                       file_name="knowledge_base.json", mime="application/json")

# ------------------------------------------------------------------ 5 report
with tab_rep:
    steps("タブ 1〜3 で出した結果（比較・フィット・照合・定量）を、試料の設定・測定条件・文献・未解決の課題と一緒に 1 つにまとめて保存します。",
          "Saves the results of tabs 1–3 (comparison, fits, QC, quantification) together with the sample setup, measurement "
          "conditions, references and open issues.")
    st.subheader(ui("レポートと出力", "Report & export"))
    q = st.session_state.get("quant")
    fits = st.session_state.get("fits", {})
    ev_df = st.session_state.get("evidence_current")
    evidence = [] if ev_df is None else [OrthogonalEvidence(r.method, bool(r.available), r.result_summary or None,
                                                            split_list(str(r.supports or ""))) for r in ev_df.itertuples()]
    ev_lines, ev_status = evidence_summary(evidence)
    sources = sorted({x.source for cs in models.values() for x in cs})
    st.write(ui("プロファイル `{p}`、知識ベース v{k}。定量：{q}。フィット済みの内殻：{f}。",
                "Profile `{p}`, KB v{k}; quantification: {q}; fitted core levels: {f}.",
                p=profile_key, k=KB.version, q=q["model"] if q else ui("未実行", "not run"), f=", ".join(fits) or ui("なし", "none")))
    proj = project_dict(setup, md, evidence, KB.version, profile_key, q["rsf"] if q else {},
                        q["options"] if q else None, q["table"] if q else None, q["results"] if q else None,
                        q["warnings"] if q else None,
                        {k: v["settings"] for k, v in fits.items()} or None,
                        next((v.get("thresholds") for v in fits.values() if v.get("thresholds")), None), sources,
                        q.get("reclass") if q else None)
    def _rec(x):
        return None if x is None else json.loads(x.to_json(orient="records"))
    proj["provenance_vocabulary"] = list(PROVENANCE)
    proj["fits"] = fits_export(fits)
    if q and q.get("warnings"):
        proj["warnings_with_severity"] = {s_: [{"severity": sev_classify(w), "text": w} for w in ws] for s_, ws in q["warnings"].items()}
    if q is not None and q.get("results") is not None:
        proj["quantification_results"] = _rec(q["results"].assign(provenance="derived_quantity"))
    if q and q.get("provisional") is not None:
        proj["provisional_quantities"] = _rec(q["provisional"])
    if q and q.get("reclass") is not None and len(q["reclass"]):
        proj["role_reclassifications"] = _rec(q["reclass"].assign(provenance="user_override"))
    rel = release_info(KB)
    sections: list = [("Release / versions", pd.DataFrame([{"app_version": rel["app_version"], "kb_version": rel["kb_version"],
                                                            "git_commit": rel["git_commit"] or "unavailable",
                                                            "python": rel["python"], "generated_at": rel["generated_at"]}])),
                      ("Interpretation levels (not equivalent)", [f"{k}: {v}" for k, v in rel["interpretation_levels"].items()]),
                      ("Sample setup", pd.DataFrame([{k: str(v) for k, v in proj["sample_setup"].items()}]).T.reset_index()
                       .rename(columns={"index": "field", 0: "value"})),
                      ("Measurement metadata", pd.DataFrame([{k: str(v) for k, v in md.to_dict().items()}]).T.reset_index()
                       .rename(columns={"index": "field", 0: "value"})),
                      ("Measurement-context warnings", ctx_warn),
                      ("Orthogonal evidence", [ev_status] + ev_lines)]
    if q:
        if q.get("provisional") is not None:
            sections.append(("Provisional quantities (nominal vs alternative; never averaged)", q["provisional"]))
        sections += [("Quantification results", q["results"].set_index("sample").T.reset_index()),
                     ("Quantification warnings", [f"{s}: {w}" for s, ws in q["warnings"].items() for w in ws])]
        if q.get("reclass") is not None and len(q["reclass"]):
            sections.append(("User role reclassifications", q["reclass"]))
    for k, v in fits.items():
        sections += [(f"Python fit – {k}", v["fit"].components), (f"Fit diagnostics – {k}", v["fit"].diagnostics)]
        if v.get("ident") is not None:
            sections.append((f"Identifiability – {k} (statistical status only)", v["ident"]))
        if v.get("provider_raw_bands") is not None:
            sections += [(f"Provider raw fitted peaks – {k}", v["provider_raw_bands"]),
                         (f"Provider pairing / tentative KB mapping – {k}", v["provider_mapping"])]
        if "grouped_qc" in v:
            sections += [(f"Role-grouped QC – {k}", v["grouped_qc"]), (f"Component QC – {k}", v["component_qc"]),
                         (f"QC summary – {k}", [f"{v['summary']['status']}: {v['summary']['message']}"])]
    sections.append(("Literature sources (active profile)", [f"{s}: {KB.reference_text(s)}" for s in sources]))
    sections.append(("Known scientific open issues", rel["scientific_open_issues"]))
    c1, c2, c3 = st.columns(3)
    c1.download_button(ui("作業の記録（.json）", "Project/session (.json)"), json.dumps(proj, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
                       file_name=f"{sample_name}_xps_project.json", mime="application/json")
    c2.download_button(ui("まとめのレポート（.html）", "Consolidated report (.html)"),
                       html_report(f"XPS analysis – {sample_name}", sections).encode("utf-8"),
                       file_name=f"{sample_name}_xps_report.html", mime="text/html", type="primary")
    c3.download_button("knowledge_base.json", KB.to_json().encode("utf-8"), file_name="knowledge_base.json")
    st.caption(ui("レポートの中身は英語です。出力はこの画面の作業内容から作り、リポジトリには書き込みません。",
                  "Exports are generated in the browser session only. Unpublished data are never written to the repository."))

    st.subheader(ui("校正のまとめ", "Calibration summary report"))
    st.caption(ui("「フィットが収束した」と「化学的な解釈が検証された」を区別して記録します（後者をこのアプリが主張することはありません）。",
                  "Distinguishes 'fit succeeded' from 'chemical interpretation validated' (the latter is never claimed by the app)."))
    ds_name = st.text_input(ui("データセット名", "Dataset name"), value=st.session_state.get("xs_name", sample_name), key="cal_ds_name")
    overrides_all = {"role_reclassifications": json.loads(q["reclass"].to_json(orient="records")) if q and q.get("reclass") is not None and len(q["reclass"]) else []}
    cal = calibration_summary(ds_name, KB, git_commit(ROOT), fits, samples=st.session_state.get("xs_samples"),
                              sensitivity={k: v.get("sensitivity") for k, v in fits.items()}, overrides=overrides_all)
    cal["provenance_vocabulary"] = list(PROVENANCE)
    r1, r2 = st.columns(2)
    r1.download_button(ui("校正のまとめ（.json）", "Calibration summary (.json)"), json.dumps(cal, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
                       file_name=f"{ds_name}_calibration_summary.json", mime="application/json")
    r2.download_button(ui("校正のまとめ（.html）", "Calibration summary (.html)"), calibration_summary_html(cal).encode("utf-8"),
                       file_name=f"{ds_name}_calibration_summary.html", mime="text/html")
    st.caption(f"KB v{KB.version} · git {cal['git_commit'] or 'unavailable'}")

# ------------------------------------------------------------------ 6 cross-sample
class _DemoFile:
    """A file of examples/synthetic_demo shaped like an uploaded file (name + getvalue)."""

    def __init__(self, path: Path):
        self.name, self._data = path.name, path.read_bytes()

    def getvalue(self) -> bytes:
        return self._data


with tab_x:
    steps("① 業者の報告書（MultiPak の .pptx か Curve Fit Summary）を入れる → ② 同じフォルダの CSV をまとめて入れる（任意。"
          "入れると Python の再フィットも並べます） → ③ 試料メタデータ（任意） → 下に「ファイルと試料の対応」と試料ごとの比較表が出ます。"
          "ファイル名は業者から届いたままにしてください。",
          "① Load the provider report (MultiPak .pptx or Curve Fit Summary) → ② add the CSV files from the same folder "
          "(optional; adds Python refits) → ③ sample metadata (optional) → the file-to-sample matching and the comparison "
          "tables appear below. Keep the provider's file names unchanged.")
    st.caption(ui("業者のフィット（スピン軌道の対にして、各試料の知識ベースのモデルに仮対応）、サーベイの at.%、暫定の定量値を並べます。"
                  "ここで化学的な解釈を加えることはなく、すべての値に出どころを付けています。",
                  "Side-by-side view of provider fits (paired, tentatively mapped to each sample's KB model), survey at.% and "
                  "provisional quantities. Nothing here changes or adds a chemical interpretation; every value carries its "
                  "provenance."))
    use_demo_x = st.checkbox(ui("合成デモデータで試す（3 試料分の合成データ。実データではありません）",
                                "Try with synthetic demo data (three synthetic samples; not real data)"), value=False, key="xs_demo")
    c1, c2, c3 = st.columns(3)
    with c1:
        rep_file = st.file_uploader(ui("① 業者の報告書（.pptx / Curve Fit Summary）", "① Provider report (MultiPak .pptx or Curve Fit Summary text)"),
                                    type=["pptx", "txt", "csv"], key="xs_rep", disabled=use_demo_x)
    with c2:
        raw_files = st.file_uploader(ui("② 生データの CSV（任意・まとめて選べます）", "② Raw-data CSV files (optional; select them all at once)"),
                                     type=["csv", "txt"], accept_multiple_files=True, key="xs_raw", disabled=use_demo_x,
                                     help=ui("業者から届いたフォルダの CSV をまとめて選んでください。*_multi.csv は再フィットに、"
                                             "*curve.csv は複数領域ファイルがない内殻の再フィットに使い、*_su.csv（サーベイ）はフィットしません。",
                                             "Select all CSV files from the provider's folder. *_multi.csv files are refitted, *curve.csv "
                                             "files are used for core levels without a multi-region file, and *_su.csv (survey) is not fitted."))
    with c3:
        meta_file = st.file_uploader(ui("③ 試料メタデータ（任意、.csv）", "③ Sample metadata (optional, .csv)"), type=["csv"], key="xs_meta",
                                     disabled=use_demo_x)
        st.download_button(ui("メタデータのテンプレート（.csv）", "Metadata template (.csv)"), csv_bytes(smeta.template()),
                           file_name="sample_metadata_template.csv")
        ex = ROOT / "templates" / "next_dataset_example.csv"
        if ex.exists():
            st.download_button(ui("記入例：次の 5 試料（.csv）", "Example: next five-sample dataset (.csv)"), ex.read_bytes(),
                               file_name="next_dataset_example.csv")
    if use_demo_x:
        rep_file = _DemoFile(DEMO / "SYNTHETIC_provider_report.txt")
        raw_files = [_DemoFile(f) for f in sorted(DEMO.glob("SYNdemo*_SYNTHETIC_multi.csv"))]
        meta_file = _DemoFile(DEMO / "SYNTHETIC_sample_metadata.csv")
        say("INFO", ui("合成デモデータを使っています（examples/synthetic_demo）。", "Using the synthetic demo data (examples/synthetic_demo)."))
    st.session_state["xs_name"] = st.text_input(ui("データセット名", "Dataset name"), value=st.session_state.get("xs_name", "dataset"),
                                                key="xs_name_in")
    if rep_file is not None:
        try:
            meta_df = smeta.read(meta_file.getvalue()) if meta_file is not None else None
            bands_x, act_x = ds.read_provider_report(rep_file.getvalue(), rep_file.name)
            files = list(dict.fromkeys(bands_x["data_file"]))
            matched = match_uploads(raw_files, files)
            raw_by_stem = matched.raw_regions
            st.session_state["xs_match"] = matched.rows
            samples_x = []
            for f in files:
                stem = str(f).rsplit(".", 1)[0]
                meta_row = None
                if meta_df is not None:
                    hit = meta_df[(meta_df["provider_file"] == stem) | (meta_df["sample_id"] == str(
                        bands_x.loc[bands_x["data_file"] == f, "comment"].iloc[0]))]
                    meta_row = hit.iloc[0].to_dict() if len(hit) else None
                samples_x.append(ds.build_sample(f, bands_x, act_x, KB, meta=meta_row, raw_regions=raw_by_stem.get(stem)))
            st.session_state["xs_samples"] = samples_x
        except Exception as e:
            say("INVALID", ui("試料を組み立てられませんでした：{e}", "Dataset assembly failed: {e}", e=e))
            samples_x = []
    else:
        samples_x = st.session_state.get("xs_samples") or []
    _rows = st.session_state.get("xs_match") if samples_x else None
    if _rows:
        _kind = {"multi": ui("複数領域（再フィットに使用）", "multi-region (refit)"),
                 "curve": ui("業者のフィット範囲", "provider curve window"),
                 "survey": ui("サーベイ", "survey"), "other": ui("その他", "other")}
        _status = {"used": ui("✅ 使用", "✅ used"),
                   "superseded": ui("— 同じ内殻の複数領域ファイルを優先", "— multi-region file used instead"),
                   "survey": ui("— フィットしない", "— not fitted"),
                   "no_sample": ui("⚠ 報告書に対応する試料がない", "⚠ no matching sample in the report"),
                   "unknown_region": ui("⚠ 領域名を読めない", "⚠ region not recognised"),
                   "no_core_level": ui("⚠ 内殻の領域がない", "⚠ no core-level region"),
                   "failed": ui("⛔ 読めない", "⛔ import failed")}
        _labels = {str(s.data_file).rsplit(".", 1)[0]: s.sample for s in samples_x}
        _n_used = sum(r["status"] == "used" for r in _rows)
        _n_bad = sum(r["status"] not in ("used", "superseded", "survey") for r in _rows)
        _n_smp = len({r["stem"] for r in _rows if r["status"] == "used"})
        (st.warning if _n_bad else st.success)(ui(
            "CSV {n} ファイルのうち {u} ファイルを {k} 試料の再フィットに使いました。対応できなかったファイル：{b}。詳しくは下の表を開いてください。",
            "{u} of {n} CSV files are used for the refits of {k} sample(s); files not matched or not read: {b}. Open the table below for details.",
            n=len(_rows), u=_n_used, k=_n_smp, b=_n_bad))
        with st.expander(ui("読み込んだファイルと試料の対応（{n} ファイル）", "Uploaded files and the samples they were matched to ({n} files)",
                            n=len(_rows)), expanded=any(r["status"] not in ("used", "superseded", "survey") for r in _rows)):
            st.caption(ui("報告書のデータファイル名（例：okawa121.spe → okawa121_…）で照合します。同じ内殻に複数領域ファイルと curve ファイルの"
                          "両方があるときは、複数領域ファイルで再フィットします。",
                          "Files are matched by the provider data-file name (e.g. okawa121.spe → okawa121_…). When both a multi-region "
                          "file and a curve file supply a core level, the multi-region file is refitted."))
            st.dataframe(pd.DataFrame([{ui("ファイル", "file"): r["file"], ui("種類", "kind"): _kind[r["kind"]],
                                        ui("試料", "sample"): _labels.get(r["stem"], "—"), ui("内殻", "core levels"): r["regions"],
                                        ui("状態", "status"): _status[r["status"]] + (f"（{r['detail']}）" if r["detail"] else "")}
                                       for r in _rows]), hide_index=True)
        if not any(r["kind"] == "multi" and r["status"] == "used" for r in _rows):
            say("CAUTION", ui("複数領域ファイル（*_multi.csv）がありません。curve ファイルは業者のフィット範囲だけを含むため、"
                              "その範囲で再フィットしています。",
                              "No multi-region file (*_multi.csv). Curve files cover only the provider's fit window, so the "
                              "refit uses that window."))
    if samples_x:
        info = pd.DataFrame([{"sample": s.sample, "data_file": s.data_file, "phase": s.meta.get("phase"),
                              "composition": s.meta.get("nominal_composition"), "treatment": s.meta.get("treatment"),
                              "treatment_temperature_C": s.meta.get("treatment_temperature_C"), "structure": s.structure,
                              "profile": s.profile, "notes": "; ".join(s.notes)} for s in samples_x])
        st.dataframe(info, hide_index=True)
        g1, g2 = st.columns(2)
        _gb = {"none": ui("まとめない", "none"), "phase": ui("相（MAX / MXene）", "phase"),
               "composition": ui("組成", "composition"), "treatment": ui("処理", "treatment")}
        group_by = g1.selectbox(ui("まとめ方", "Group by"), list(_gb), key="xs_group", format_func=_gb.get)
        chosen = g2.multiselect(ui("表示する試料", "Samples"), list(info["sample"]), default=list(info["sample"]), key="xs_pick")
        sel = [s for s in samples_x if s.sample in chosen]
        long = ds.long_table(sel)

        def _w(quantity, core=None):
            w = ds.wide(long, quantity, core)
            if w.empty or group_by == "none":
                return w
            gmap = long.drop_duplicates("sample").set_index("sample")[group_by].astype(str).to_dict()
            w.columns = pd.MultiIndex.from_tuples([(gmap.get(c, ""), c) for c in w.columns], names=[group_by, "sample"])
            return w.sort_index(axis=1)

        st.markdown(ui("#### 暫定の定量値（基準値＝校正候補を骨格に数えない、代替値＝数える。平均はしない）",
                       "#### Provisional quantities (nominal = candidate excluded; alternative = candidate as framework; never averaged)"))
        qs = ds.quant_summary(sel)
        st.dataframe(qs, hide_index=True)
        cores_x = sorted({c for s in sel for c in s.mappings})
        core_pick = st.selectbox(ui("下の表で見るコア準位", "Core level for BE / area-fraction tables"), cores_x, key="xs_core")
        if core_pick in UNSUPPORTED_REGIONS:
            say("CAUTION", ui("{c}：検証済みモデルなし。業者のピークを読み込んだまま表示しています（対にできるものは対にしています）。",
                              "{c}: {t}; provider peaks shown as imported (paired where applicable).", c=core_pick, t=NO_MODEL_TEXT))
        st.markdown(ui("#### 成分ごとの結合エネルギー — {c}（業者のフィット）", "#### Fitted BE by component — {c} (provider; provenance provider_paired_doublet / provider_raw_fit)", c=core_pick))
        st.dataframe(_w("BE (eV)", core_pick))
        st.markdown(ui("#### 成分ごとの面積分率 — {c}", "#### Area fraction by component — {c}", c=core_pick))
        st.dataframe(_w("area fraction", core_pick))
        if any(s.python for s in sel):
            st.markdown(ui("#### Python の再フィット — {c}（面積分率。同定可能性は縦長の表に）", "#### Python refit — {c} (area fraction; identifiability in the long table)", c=core_pick))
            st.dataframe(_w("area fraction [python]", core_pick))
        for oc in ("O 1s", "F 1s"):
            if oc in cores_x and oc != core_pick:
                st.markdown(ui("#### {o} の分率", "#### {o} fractions", o=oc))
                st.dataframe(_w("area fraction", oc))
        st.markdown(ui("#### サーベイの at.%（業者の原子濃度表。測っていれば残留 Al・Cl・S・Li も）",
                       "#### Survey at.% (provider atomic concentration table) — incl. Al residual, Cl, S, Li where measured"))
        st.dataframe(_w("at.%"))
        unres = long[(long["quantity"] == "BE (eV)") & long["identifiability"].astype(str).ne("")]
        if len(unres):
            st.markdown(ui("#### 未解決・モデル外の業者成分（ラベルは中立。帰属はしない）",
                           "#### Unresolved / unmodelled provider components (labels are neutral; no assignment)"))
            st.dataframe(unres[["sample", "core_level", "component", "provider_label", "value", "identifiability"]]
                         .rename(columns={"value": "BE (eV)"}), hide_index=True)
        st.download_button(ui("比較表（縦長・出どころ付き、.csv）", "Comparison long table with provenance (.csv)"), csv_bytes(long.assign(app_version=__version__, kb_version=KB.version)),
                           file_name=f"{st.session_state['xs_name']}_comparison_long.csv")
        st.download_button(ui("暫定の定量値（.csv）", "Provisional quantities (.csv)"), csv_bytes(qs.assign(app_version=__version__, kb_version=KB.version)), file_name=f"{st.session_state['xs_name']}_provisional.csv")
