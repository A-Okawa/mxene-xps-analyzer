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

KB = default_kb()
FIX = ROOT / "tests" / "fixtures"
STATUS_COLOR = {"PASS": "#d8f0d8", "CHECK": "#fff0c2", "REVIEW": "#fff0c2", "CHECK-MISSING": "#ececec", "FAIL": "#f7d0d0"}

st.set_page_config(page_title="MAX/MXene XPS Analyzer", layout="wide")


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
    st.title("MAX/MXene XPS Analyzer")
    entered = st.text_input("Password", type="password", key="_auth_pw")
    if entered:
        if hmac.compare_digest(entered.encode(), str(expected).encode()):
            st.session_state["_auth_ok"] = True
            st.rerun()
        st.error("Incorrect password.")
    st.stop()


require_password()


def style_status(df: pd.DataFrame):
    cols = [c for c in df.columns if c == "QC" or c.startswith("QC_") or c == "metal_QC"]
    fmt = {c: "{:.4g}" for c in df.columns if pd.api.types.is_float_dtype(df[c])}
    return df.style.format(fmt, na_rep="").map(lambda v: f"background-color:{STATUS_COLOR.get(v, '')}", subset=cols)


from mxene_xps.severity import QC_STATUS_SEVERITY, classify as sev_classify  # noqa: E402

_SEV_UI = {"INFO": st.info, "CAUTION": st.warning, "WARNING": st.warning, "INVALID": st.error}
_SEV_ICON = {"INFO": "ℹ️", "CAUTION": "🟡", "WARNING": "⚠️", "INVALID": "⛔"}


def say(a: str, b: str | None = None) -> None:
    """Show a message with a consistent severity: say(text[, LEVEL]) or say(LEVEL, text)."""
    if a in _SEV_UI and b is not None:
        level, text = a, b
    else:
        text, level = a, (b or sev_classify(a))
    _SEV_UI[level](f"**{level}** — {text}", icon=_SEV_ICON[level])


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
    st.header("Sample setup")
    _order = ["(Nb,Ti)2CTx", "Nb2CTx", "Ti3C2Tx"]
    structure = st.selectbox("Nominal structure", _order + [k for k in KB.structures if k not in _order] + ["Custom"])
    if structure == "Custom":
        m_target = st.number_input("Ideal M-site count", min_value=0.1, value=2.0, step=0.1)
        x_target = st.number_input("Ideal X-site count", min_value=0.1, value=1.0, step=0.1)
        profile_key = st.selectbox("Reference profile", KB.profile_keys)
    else:
        spec = KB.structures[structure]
        m_target, x_target, profile_key = float(spec["m_target"]), float(spec["x_target"]), spec["default_profile"]
        st.caption(f"M_target = {m_target:g}, X_target = {x_target:g}, profile = `{profile_key}`")
    phase = (KB.structures.get(structure) or KB.profile(profile_key)).get("phase", "MXene")
    a_target = (KB.structures.get(structure) or {}).get("a_target")
    if phase == "MAX":
        st.caption("MAX-phase profile: data-derived calibration candidates (not literature).")
    default_m = {"(Nb,Ti)2CTx": "Nb,Ti", "Nb2CTx": "Nb", "Ti3C2Tx": "Ti", "Nb2AlC": "Nb", "(Nb,Ti)2AlC": "Nb,Ti"}.get(structure, "")
    m_elements = split_list(st.text_input("M-site elements (comma-separated)", default_m))
    sample_name = st.text_input("Sample name", "NbTi_MXene_01" if structure == "(Nb,Ti)2CTx" else
                                ("MAX_01" if phase == "MAX" else "MXene_01"))
    etchants = st.multiselect("Etching / processing chemistry",
                              ["HF", "LiF/HCl", "NaF/HCl", "KF/HCl", "NH4F/HCl", "TBAOH", "molten chloride", "H2SO4 (post-treatment)", "other"],
                              default=["LiF/HCl"] if structure == "(Nb,Ti)2CTx" else [])
    a_element = st.text_input("A-site element to check", "Al")
    expected_defects = st.checkbox("M-site vacancies / selective dissolution expected", value=structure == "(Nb,Ti)2CTx")
    cl_interest = st.checkbox("Cl termination of interest (request Cl 2p even without chloride chemistry)", value=False)
    dopants = split_list(st.text_input("Post-synthesis dopant(s), comma-separated", ""))
    nominal: dict[str, float] = {}
    if len(m_elements) > 1:
        st.subheader("Nominal M-site stoichiometry")
        for el in m_elements:
            nominal[el] = st.number_input(f"{el} nominal amount", min_value=0.0, value=float(m_target / len(m_elements)),
                                          step=0.05, key=f"nom_{el}")
    with st.expander("More sample details"):
        parent = st.text_input("Parent/MAX composition", "")
        post = st.text_input("Post-synthesis treatment", "")
        term_family = st.text_input("Expected termination family", "")
        notes = st.text_area("Notes", "")

    with st.expander("Measurement metadata (§2.8, §21, §22)", expanded=False):
        md = MeasurementMetadata(
            provider=st.text_input("Provider", "") or None,
            instrument=st.text_input("Instrument", "") or None,
            xray_source=st.text_input("X-ray source", "") or None,
            pass_energy_eV=st.number_input("Pass energy (eV, 0 = unknown)", min_value=0.0, value=0.0) or None,
            step_eV=st.number_input("Step size (eV, 0 = unknown)", min_value=0.0, value=0.0, format="%.3f") or None,
            spot_size=st.text_input("Spot size", "") or None,
            charge_neutralization=st.text_input("Charge neutralization", "") or None,
            be_reference_method=st.selectbox("BE reference method", BE_REFERENCE_METHODS),
            be_reference_detail=st.text_input("BE reference detail (e.g. which metal / value)", "") or None,
            surface_state=st.selectbox("Surface state", SURFACE_STATES),
        )
        if md.sputtered:
            md.sputter_species = st.text_input("Sputter ion", "Ar+") or None
            md.sputter_energy_eV = st.number_input("Sputter energy (eV)", min_value=0.0, value=0.0) or None
            md.sputter_time_s = st.number_input("Sputter time (s)", min_value=0.0, value=0.0) or None
    st.divider()
    st.caption(f"App v{__version__} · knowledge base v{KB.version}. Uploaded files are processed in memory by the server "
               "running this app and are not sent to any other service.")

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

st.title("MAX/MXene XPS Analyzer")
st.caption("Defect/dopant-aware XPS request generation, quantification, and provider-vs-Python validation. "
           "Binding energy is not treated as a one-to-one chemical-state lookup.")
st.caption("Validation status: release candidate. Calibrated on one development dataset; prospective hold-out "
           "validation not yet performed. Treat fit-based results as provisional (see tab 4 and the README).")
ctx_warn = metadata_warnings(md)
if ctx_warn:
    with st.expander(f"⚠ Measurement-context warnings ({len(ctx_warn)})", expanded=md.sputtered):
        for w in ctx_warn:
            say(w)

tab_req, tab_q, tab_v, tab_kb, tab_rep, tab_x = st.tabs(
    ["1. Request generator", "2. Peak table & quantification", "3. Provider ↔ Python validation",
     "4. Knowledge base", "5. Report & export", "6. Cross-sample comparison"])

# ------------------------------------------------------------------ 1 request
with tab_req:
    req_text, regions = generate_request(setup, KB)
    st.subheader("Recommended regions")
    st.dataframe(regions, hide_index=True)
    st.subheader("Measurement / fitting request")
    st.caption("Editable. English only in v1. Candidate fit numbers come from the active knowledge-base profile.")
    edited_req = st.text_area("Request text", req_text, height=600, key=f"req_{hash(req_text)}")
    st.download_button("Download request (.txt)", edited_req, file_name=f"{sample_name}_XPS_request.txt", mime="text/plain")

# ------------------------------------------------------------------ 2 quant
with tab_q:
    use_legacy = False
    if structure == "Ti3C2Tx":
        use_legacy = st.radio("Quantification model", ["Role-based (Brette2024 profile)", "Legacy Ti3C2Tx workbook model"],
                              horizontal=True) == "Legacy Ti3C2Tx workbook model"

    if use_legacy:
        say(LEGACY_LABEL, "INFO")
        tmpl = legacy_template(sample_name)
        st.download_button("Download legacy component template", csv_bytes(tmpl), file_name=f"{sample_name}_legacy_template.csv")
        up = st.file_uploader("Upload legacy component CSV (sample, core_level, component, area)", type=["csv"], key="leg_up")
        base = read_csv_bytes(up.getvalue()) if up else tmpl
        leg_df = st.data_editor(base, num_rows="dynamic", key=f"leg_ed_{up.name if up else sample_name}")
        c1, c2 = st.columns(2)
        with c1:
            ctrl = {
                "F_oxide_per_Ooxide": st.number_input("F_oxide_per_Ooxide (legacy, model-dependent)", value=LEGACY_CONTROLS_DEFAULT["F_oxide_per_Ooxide"], format="%.4f"),
                "Include_O_highBE": int(st.checkbox("Include_O_highBE", True)),
                "Include_Cl": int(st.checkbox("Include_Cl", True)),
                "Include_F1": int(st.checkbox("Include_F1", False)),
                "Include_F2": int(st.checkbox("Include_F2", False)),
            }
        with c2:
            st.caption("Historical RSFs are legacy defaults only – confirm with the provider for new data.")
            rsf_df = st.data_editor(pd.DataFrame({"core_level": list(KB.legacy["default_rsf"]),
                                                  "RSF": list(KB.legacy["default_rsf"].values())}),
                                    hide_index=True, key="leg_rsf")
        rsf = {r.core_level: float(r.RSF) for r in rsf_df.itertuples() if pd.notna(r.RSF)}
        if st.button("Calculate (legacy model)", type="primary"):
            try:
                res, warns = legacy_quantify_table(leg_df, rsf, ctrl)
                st.session_state["quant"] = {"results": res, "warnings": warns, "model": "legacy", "profile": profile_key, "options": ctrl,
                                             "rsf": rsf, "table": leg_df, "reclass": None}
            except Exception as e:
                say(str(e), "INVALID")
    else:
        st.subheader("Provider component table")
        st.markdown("Required columns: `sample, core_level, component, area, role`. Recommended: "
                    "`element, BE, FWHM, notes`. Blank areas = component not fitted.")
        tmpl = make_peak_template(setup, KB)
        st.download_button("Download component CSV template", csv_bytes(tmpl), file_name=f"{sample_name}_peak_template.csv")
        up = st.file_uploader("Upload fitted component CSV", type=["csv"], key="peak_up")
        base = read_csv_bytes(up.getvalue()) if up else tmpl
        ed_key = f"peak_ed_{up.name if up else ''}_{profile_key}_{sample_name}_{','.join(dopants)}"
        peaks = st.data_editor(
            base, num_rows="dynamic", key=ed_key,
            column_config={"role": st.column_config.SelectboxColumn("role", options=list(ROLES), required=True,
                                                                    help="Quantification role (§6). Changing a role is logged."),
                           "area": st.column_config.NumberColumn("area", min_value=0.0)})

        reclass = role_reclassifications(peaks, kb_roles) if {"core_level", "component", "role", "sample"} <= set(peaks.columns) else pd.DataFrame()
        if len(reclass):
            say("CAUTION", "Roles changed from the knowledge-base default. These are user assignments and are recorded in "
                       "the export; each must be justified (e.g. by metal-core cross-check or orthogonal evidence).")
            st.dataframe(reclass, hide_index=True)

        st.subheader("Quantification controls")
        c1, c2, c3 = st.columns(3)
        with c1:
            input_mode = st.selectbox("Area input mode", ["raw_area", "corrected_area", "atomic_percent"])
            inc_dis = st.checkbox("Count disordered O (termination_O_disordered) as Tx", value=True)
            inc_oh = st.checkbox("Count OH-like O as Tx", value=True)
            st.caption("O_unresolved_low / O_unresolved_mid are never counted as Tx; reclassify a component's role "
                       "explicitly if a cross-check justifies it.")
        with c2:
            inc_cl = st.checkbox("Count Cl as Tx", value=True)
            inc_fres = st.checkbox("Count high-BE F residue as Tx (override)", value=False)
            fcorr = st.number_input("Oxide-derived F correction factor (0 = off)", min_value=0.0,
                                    value=0.6667 if structure == "Ti3C2Tx" else 0.0, step=0.05, format="%.4f",
                                    help="0.6667 is the legacy Ti3C2Tx model; not transferred to mixed Nb/Ti by default.")
            tol = st.number_input("Metal-QC tolerance |M_target − M_norm| (software default)", min_value=0.0, value=0.35, step=0.05)
        with c3:
            cores = sorted(set(peaks["core_level"].dropna().astype(str))) if "core_level" in peaks else []
            prefill = st.checkbox("Prefill legacy laboratory RSFs where available (must be confirmed)", value=False,
                                  disabled=input_mode != "raw_area")
            leg = KB.legacy["default_rsf"]
            rsf_in = pd.DataFrame({"core_level": cores,
                                   "RSF": [leg.get(c, np.nan) if prefill else np.nan for c in cores]})
            rsf_df = st.data_editor(rsf_in, hide_index=True, key=f"rsf_{prefill}_{'|'.join(cores)}",
                                    disabled=input_mode != "raw_area")
            rsf = {r.core_level: float(r.RSF) for r in rsf_df.itertuples() if pd.notna(r.RSF)}
            if input_mode == "raw_area":
                miss = [c for c in cores if c not in rsf]
                if miss:
                    say("CAUTION", f"RSF required before calculation (no 1.0 default): {', '.join(miss)}")

        with st.expander("Orthogonal evidence (§2.2, §24.4) – displayed, not scored"):
            if "evidence" not in st.session_state:
                st.session_state["evidence"] = pd.DataFrame({"method": list(ORTHOGONAL_METHODS), "available": False,
                                                             "result_summary": "", "supports": ""})
            ev_df = st.data_editor(st.session_state["evidence"], hide_index=True, num_rows="dynamic", key="ev_ed",
                                   column_config={"supports": st.column_config.TextColumn(
                                       "supports", help="comma-separated, e.g. Ti depletion, M vacancy, dopant incorporation")})
            st.session_state["evidence_current"] = ev_df

        if st.button("Calculate", type="primary"):
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
        st.subheader("Results")
        res = q["results"]
        show = (res.set_index("sample").T
                .map(lambda x: f"{x:.6g}" if isinstance(x, (float, np.floating)) else x)
                .reset_index().rename(columns={"index": "quantity"}))
        st.dataframe(show, hide_index=True, height=min(900, 36 * (len(show) + 1)))
        st.download_button("Download quantification results (.csv)", csv_bytes(res.assign(provenance="derived_quantity",
                                                                                         app_version=__version__, kb_version=KB.version)), file_name=f"{sample_name}_quantification.csv")
        if q["model"] == "role-based":
            st.markdown("**Nominal vs alternative (provisional quantities)** — default excludes the calibration candidate "
                        "from the framework; the alternative counts it. The two are never averaged.")
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
            say("CAUTION", "apparent_M_deficiency_delta is PROVISIONAL (Nb 3d ~205.6–205.9 eV candidate not validated). "
                     "Compare with `apparent_M_deficiency_delta_alt` (candidate counted as framework); do not report "
                     "either as final.")
        if q["model"] == "role-based":
            say("INFO", "`apparent_M_deficiency_delta` is a surface-sensitive XPS index after carbide-C normalization, not "
                    "an absolute vacancy concentration. Selective depletion is reported separately via "
                    "`<M>_fraction_shift` / `<M>_retention_ratio_fractional`.")
            ev_lines, ev_status = evidence_summary([
                OrthogonalEvidence(r.method, bool(r.available), r.result_summary or None, split_list(str(r.supports or "")))
                for r in st.session_state.get("evidence_current", pd.DataFrame(columns=["method", "available", "result_summary", "supports"])).itertuples()])
            st.caption(ev_status)
            for line in ev_lines:
                st.write("• " + line)

# ------------------------------------------------------------------ 3 validation
with tab_v:
    st.subheader("Independent Python refit")
    st.caption("The Python fit is an independent QC fit, not ground truth. Pseudo-Voigt components with fixed spin-orbit "
               "split and area ratio; Shirley background is an uncalibrated approximation.")
    cores = list(models) + [r for r in UNSUPPORTED_REGIONS if r not in models]
    core = st.selectbox("Core level", cores, format_func=lambda c: c if c in models else f"{c} ({NO_MODEL_TEXT})")
    kb_comps = models.get(core, [])
    unsupported = not kb_comps
    if unsupported:
        say("CAUTION", f"{core}: {NO_MODEL_TEXT}. The region can be imported, displayed and exported, and provider peaks can be "
                 f"imported, but no Python fit is run. {UNSUPPORTED_REGIONS.get(core, '')}")
    for w in profile.get("core_level_warnings", {}).get(core, []):
        say(w, "CAUTION")
    c1, c2 = st.columns(2)
    with c1:
        importer = st.selectbox("Raw-spectrum format", list(SPECTRUM_IMPORTERS) + list(MULTIREGION_IMPORTERS))
        raw_file = st.file_uploader("Raw spectrum (BE, intensity; header optional; any order)", type=["csv", "txt", "dat"], key="raw")
        region_sp = None
        if raw_file is not None and importer in MULTIREGION_IMPORTERS:
            try:
                regs = MULTIREGION_IMPORTERS[importer](raw_file.getvalue())
                labels = [f"{r.region} → {r.core_level or 'not a core level'} ({r.be.min():.1f}–{r.be.max():.1f} eV, "
                          f"{len(r.be)} pts)" for r in regs]
                default = next((i for i, r in enumerate(regs) if r.core_level == core), 0)
                pick = st.selectbox("Region in file", list(range(len(regs))), index=default, format_func=lambda i: labels[i],
                                    key=f"region_{raw_file.name}_{core}")
                region_sp = regs[pick]
                if region_sp.core_level != core:
                    say("INVALID", f"Selected region {region_sp.region} is not {core}; fitting it with the {core} model is refused.")
            except Exception as e:
                say("INVALID", f"Multi-region import failed: {e}")
        restrict = st.checkbox("Restrict fit BE window", value=False,
                               help="E.g. exclude a neighbouring core level. Recorded in the fit settings.")
        win = None
        if restrict and kb_comps:
            w1, w2 = st.columns(2)
            win = (w1.number_input("low BE (eV)", value=float(min(x.bounds[0] for x in models[core]) - 3.0), format="%.2f"),
                   w2.number_input("high BE (eV)", value=float(max(x.bounds[1] + (x.split or 0) for x in models[core]) + 3.0), format="%.2f"))
        demo_tag = {"Ti 2p": "Ti2p", "Nb 3d": "Nb3d"}.get(core)
        use_demo = st.checkbox("Use synthetic demo spectrum (public test fixture)", value=False, disabled=demo_tag is None)
    with c2:
        prov_fmt = st.radio("Provider fit format", ["Generic component CSV", "PHI MultiPak Curve Fit Summary"], horizontal=True)
        if prov_fmt == "Generic component CSV":
            vendor_file = st.file_uploader("Provider component CSV (component, BE, FWHM, area[, role, core_level])", type=["csv"], key="vend")
            doublet_conv = st.radio("Provider doublet area is…", ["total (both spin-orbit lines)", "primary line only"], horizontal=True)
        else:
            vendor_file = st.file_uploader("MultiPak report (.pptx) or Curve Fit Summary text (.txt/.csv)",
                                           type=["pptx", "txt", "csv"], key="vend_mp")
            doublet_conv = "total (both spin-orbit lines)"
            st.caption("Spin-orbit lines reported as separate bands are paired by region, splitting, area ratio and "
                       "proximity. Provider labels are kept; the KB mapping below is tentative and editable.")
        bg = st.selectbox("Background", ["shirley", "linear"])
        endavg = st.number_input("Background endpoint averaging (points)", min_value=1, value=1, step=1)
    with st.expander("QC thresholds (software defaults, not physical constants)"):
        vd = KB.validation_defaults
        a, b, c = st.columns(3)
        th = {
            "delta_be_pass_eV": a.number_input("|ΔBE| PASS ≤ (eV)", 0.0, value=float(vd["delta_be_pass_eV"]), step=0.05),
            "delta_be_review_eV": a.number_input("|ΔBE| REVIEW ≤ (eV)", 0.0, value=float(vd["delta_be_review_eV"]), step=0.05),
            "fwhm_rel_pass": b.number_input("rel. ΔFWHM PASS ≤", 0.0, value=float(vd["fwhm_rel_pass"]), step=0.05),
            "fwhm_rel_review": b.number_input("rel. ΔFWHM REVIEW ≤", 0.0, value=float(vd["fwhm_rel_review"]), step=0.05),
            "area_fraction_abs_pass": c.number_input("Δarea-fraction PASS ≤", 0.0, value=float(vd["area_fraction_abs_pass"]), step=0.02),
            "area_fraction_abs_review": c.number_input("Δarea-fraction REVIEW ≤", 0.0, value=float(vd["area_fraction_abs_review"]), step=0.02),
        }
        corr_th = st.number_input("Area-correlation flag |r| ≥ (identifiability diagnostic)", 0.0, 1.0, 0.8, 0.05)
        st.caption(vd.get("note", ""))
    if kb_comps:
        st.markdown("**Active starting model** — starting values and fitting windows are editable implementation "
                    "bounds; literature reference values are shown separately and are never used as constraints.")

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
        if st.button("Reset model to knowledge-base values"):
            st.session_state["model_reset"] = ver + 1
            st.rerun()
        overrides = {r["name"]: {"center": r["center"], "bounds": (r["BE_low"], r["BE_high"]), "fwhm": r["FWHM"],
                                 "fwhm_bounds": (r["FWHM_low"], r["FWHM_high"])} for r in edited_model.to_dict("records")}
        try:
            comps, model_log = with_overrides(kb_comps, overrides)
        except ValueError as e:
            say("INVALID", f"Invalid model edit: {e}")
            comps, model_log = kb_comps, []
        if model_log:
            say("CAUTION", f"{len(model_log)} starting value/window edit(s) relative to KB v{KB.version}; recorded in the "
                       "fit settings and the project export.")
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
            st.markdown(f"**Provider fit (MultiPak)** — {len(keys)} summaries in file, {len(cands)} for {core}")
            if len(cands):
                opts = list(cands["summary_index"])
                raw_name = raw_file.name if raw_file is not None else ""
                stems = [str(f).rsplit(".", 1)[0] for f in cands["data_file"]]
                default = next((i for i, (st_, cm) in enumerate(zip(stems, cands["comment"].astype(str)))
                                if (st_ and st_ in raw_name) or (cm and cm == sample_name)), 0)
                sidx = st.selectbox("Provider summary", opts, index=default,
                                    format_func=lambda i: cands.set_index("summary_index").at[i, "lbl"],
                                    key=f"mp_sum_{core}_{raw_name}")
                if raw_name and not any(st_ and st_ in raw_name for st_ in stems):
                    say("WARNING", "No provider summary matches the raw-file name; confirm the selected summary belongs to "
                               "the same sample/region.")
                one = mp_bands[mp_bands["summary_index"] == sidx]
                doublets, linked = pair_spin_orbit(one)
                auto = map_to_kb(doublets, comps)
                with st.expander("Provider raw fitted peaks (as reported; never modified)"):
                    st.dataframe(linked.drop(columns=["origin"]), hide_index=True)
                st.markdown("Spin-orbit pairing and tentative KB mapping (edit `mapped_component`; edits are logged)")
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
                    say("CAUTION", "Only one spin-orbit line was reported for: " + ", ".join(est["provider_label"])
                               + ". Total areas are ESTIMATED with the degeneracy ratio.")
                chk = mapping[mapping["area_convention"] == "unpaired_check"]
                if len(chk):
                    say("WARNING", "Unpaired band(s) next to paired doublets: " + ", ".join(chk["provider_label"]) + " — check.")
                if map_log:
                    say("CAUTION", f"{len(map_log)} user mapping change(s) relative to the automatic tentative mapping (logged).")
                    st.dataframe(pd.DataFrame(map_log).astype(str), hide_index=True)
                mp_state = {"bands": linked, "mapping": mapping, "mapping_log": map_log,
                            "vendor": provider_component_table(mapping)}
        except Exception as e:
            say("INVALID", f"MultiPak import failed: {e}")

    if unsupported:
        sp_show = region_sp if region_sp is not None else None
        if sp_show is None and raw_file is not None and importer in SPECTRUM_IMPORTERS:
            try:
                sp_show = SPECTRUM_IMPORTERS[importer](raw_file.getvalue())
            except Exception as e:
                say("INVALID", f"Import failed: {e}")
        if sp_show is not None:
            fig, ax = plt.subplots(figsize=(8, 3.5))
            ax.plot(sp_show.be, sp_show.intensity, ".-", ms=2, lw=0.6, color="0.3")
            ax.invert_xaxis(); ax.set_xlabel("Binding energy (eV)"); ax.set_ylabel("Intensity")
            ax.set_title(f"{core} — imported spectrum ({NO_MODEL_TEXT})")
            st.pyplot(fig); plt.close(fig)
            st.download_button("Imported spectrum .csv", csv_bytes(pd.DataFrame({"BE": sp_show.be, "intensity": sp_show.intensity,
                                                                                 "provenance": "raw_measurement"})),
                               file_name=f"{sample_name}_{core.replace(' ', '')}_raw.csv")
        if mp_state is not None:
            st.download_button("Provider peaks (paired) .csv", csv_bytes(mp_state["mapping"].assign(provenance="provider_paired_doublet")),
                               file_name=f"{sample_name}_{core.replace(' ', '')}_provider_peaks.csv")

    if st.button("Run Python fit", type="primary", disabled=unsupported):
        try:
            if region_sp is not None:
                sp = region_sp
            elif raw_file is not None and importer in MULTIREGION_IMPORTERS:
                raise ValueError("Select a valid region of the multi-region file.")
            elif raw_file is not None:
                sp = SPECTRUM_IMPORTERS[importer](raw_file.getvalue())
            elif use_demo and demo_tag:
                sp = SPECTRUM_IMPORTERS[importer]((FIX / f"synthetic_{demo_tag}_raw.csv").read_bytes())
            else:
                raise ValueError("Upload a raw spectrum or tick the demo option.")
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
            say("INVALID", f"Fit failed: {e}")

    fitstate = st.session_state.get("fits", {}).get(f"{profile_key} | {core}")
    if fitstate:
        fr = fitstate["fit"]
        for n in fitstate["import_notes"]:
            st.caption("Import: " + n)
        m = fr.metrics
        k1, k2, k3, k4, k5 = st.columns(5)
        k1.metric("converged", "yes" if m["success"] else "no")
        k2.metric("R²", f"{m['R2']:.5f}")
        k3.metric("RMSE", f"{m['RMSE']:.4g}")
        k4.metric("AIC", f"{m['AIC']:.1f}")
        k5.metric("BIC", f"{m['BIC']:.1f}")
        st.dataframe(style_status(fr.components.assign(provenance="python_fit")), hide_index=True)
        ident = fitstate.get("ident")
        if ident is not None:
            st.markdown("**Identifiability / fit confidence** (statistical status only — not a chemical assignment)")
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
                    say(lv, "Fit diagnostics:\n\n" + "\n".join(f"- {d}" for d in ds_))
        cand = [x for x in comps if x.calibration_candidate and not x.validated]
        if cand:
            say("CAUTION", "Calibration candidate(s) in this model (fitted, not validated): " + "; ".join(
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
            dl[i].download_button(f"Plot .{fmt}", buf.getvalue(), file_name=f"{sample_name}_{core.replace(' ', '')}_fit.{fmt}")
        plt.close(fig)
        dl[3].download_button("Fit components .csv", csv_bytes(fr.components.assign(provenance="python_fit")), file_name=f"{sample_name}_{core.replace(' ', '')}_python_fit.csv")
        dl[4].download_button("Curves .csv", csv_bytes(cv.assign(provenance="observed=raw_measurement; other columns=python_fit")), file_name=f"{sample_name}_{core.replace(' ', '')}_curves.csv")
        if fr.area_correlation is not None:
            with st.expander("Area-parameter correlation matrix"):
                st.dataframe(fr.area_correlation.style.format("{:+.2f}"))

        with st.expander("Fit-sensitivity comparison (predefined alternatives; no automatic 'best chemistry')"):
            be_s, y_s, core_s = fitstate["spectrum"]
            alts = sens.standard_alternatives(core, comps, be_s)
            st.caption("Alternatives: " + "; ".join(f"{a.label} [{a.basis}]" for a in alts)
                       + ". ΔAIC/ΔBIC are shown only where the data points are identical to the nominal fit.")
            qd = st.session_state.get("quant")
            derived_fn, dsample = None, None
            if qd and qd.get("model") == "role-based" and "core_level" in qd["table"]:
                tab = qd["table"]
                samples_q = sorted(set(tab.loc[tab["core_level"] == core, "sample"].astype(str)))
                if samples_q:
                    dsample = st.selectbox("Derived quantities from the quantification table of sample", samples_q, key=f"sens_s_{core}")
                    qopt = QuantOptions(**qd["options"])
                    derived_fn = sens.derived_from_fractions(tab, dsample, core,
                                                             lambda t: quantify_table(t, qd["rsf"], qopt)[0].iloc[0].to_dict())
            if st.button("Run fitting alternatives", key=f"sens_btn_{core}"):
                with st.spinner("Fitting alternatives…"):
                    sres = sens.run(be_s, y_s, comps, alts, endpoint_avg=int(endavg), derived_fn=derived_fn, spectrum_core_level=core_s)
                fitstate["sensitivity"] = sres
                fitstate["ident"] = assess_identifiability(fr, alt_spans=sres.component_spans())
                st.rerun()
            sres = fitstate.get("sensitivity")
            if sres is not None:
                st.dataframe(sres.table, hide_index=True)
                st.markdown("Area fraction by alternative")
                st.dataframe(sres.fractions.style.format("{:.3f}"))
                if len(sres.derived):
                    st.markdown(f"Derived quantities (sample {dsample})")
                    st.dataframe(sres.derived)
                st.markdown("Spans across alternatives (≥ 0.05 → **model-sensitive**, D4 materiality criteria)")
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
                st.markdown("**Role-grouped QC**")
                st.dataframe(style_status(gqc), hide_index=True)
                st.markdown("**Component-level QC** (ΔBE = Python − Provider)")
                st.dataframe(style_status(cqc), hide_index=True)
                st.caption(IDENTIFIABILITY_WARNING)
                d1, d2 = st.columns(2)
                d1.download_button("Component QC .csv", csv_bytes(cqc.assign(provenance_vendor="provider", provenance_python="python_fit")), file_name=f"{sample_name}_{core.replace(' ', '')}_validation.csv")
                d2.download_button("Grouped QC .csv", csv_bytes(gqc.assign(provenance="derived_quantity")), file_name=f"{sample_name}_{core.replace(' ', '')}_validation_grouped.csv")
                fitstate.update({"component_qc": cqc, "grouped_qc": gqc, "summary": summ.__dict__, "thresholds": th})
                if mp_state is not None:
                    fitstate.update({"provider_raw_bands": mp_state["bands"], "provider_mapping": mp_state["mapping"],
                                     "provider_mapping_log": mp_state["mapping_log"]})
            except Exception as e:
                say("INVALID", f"Provider comparison failed: {e}")

# ------------------------------------------------------------------ 4 KB
with tab_kb:
    st.subheader("Interpretation principles")
    for p in KB.raw["principles"]:
        st.markdown(f"- {p}")
    st.subheader(f"Active profile: `{profile_key}` ({profile['material']})")
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
        st.subheader("Calibration status")
        st.json(KB.raw["calibration_status"])
    st.subheader("Role ontology")
    st.dataframe(pd.DataFrame({"role": list(ROLE_LABELS), "meaning": list(ROLE_LABELS.values())}), hide_index=True)
    st.subheader("References")
    st.dataframe(pd.DataFrame([{"key": k, **v} for k, v in KB.references.items()]), hide_index=True)
    st.subheader("Ti₂₋ᵧNbᵧCTx XPS atomic % (Han2020SolidSolution)")
    ssr = KB.raw["solid_solution_reference"]["Ti2-yNbyCTx"]
    st.caption(ssr["note"] + " XPS in that work was measured after Ar sputtering (2 kV, 2 µA, 2 min).")
    st.dataframe(pd.DataFrame(ssr["xps_atomic_percent"]), hide_index=True)
    st.subheader("Specification discrepancies and decisions")
    st.dataframe(pd.DataFrame(KB.raw.get("spec_discrepancies", [])), hide_index=True)
    st.caption(f"KB v{KB.version} · created {KB.raw.get('created_at')} · reviewed_by: {KB.raw.get('reviewed_by') or 'not yet reviewed'}")
    st.download_button("Download knowledge_base.json", KB.to_json().encode("utf-8"), file_name="knowledge_base.json", mime="application/json")

# ------------------------------------------------------------------ 5 report
with tab_rep:
    st.subheader("Report & export")
    q = st.session_state.get("quant")
    fits = st.session_state.get("fits", {})
    ev_df = st.session_state.get("evidence_current")
    evidence = [] if ev_df is None else [OrthogonalEvidence(r.method, bool(r.available), r.result_summary or None,
                                                            split_list(str(r.supports or ""))) for r in ev_df.itertuples()]
    ev_lines, ev_status = evidence_summary(evidence)
    sources = sorted({x.source for cs in models.values() for x in cs})
    st.write(f"Profile `{profile_key}`, KB v{KB.version}; quantification: {q['model'] if q else 'not run'}; "
             f"fitted core levels: {', '.join(fits) or 'none'}.")
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
    c1.download_button("Project/session (.json)", json.dumps(proj, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
                       file_name=f"{sample_name}_xps_project.json", mime="application/json")
    c2.download_button("Consolidated report (.html)", html_report(f"XPS analysis – {sample_name}", sections).encode("utf-8"),
                       file_name=f"{sample_name}_xps_report.html", mime="text/html")
    c3.download_button("knowledge_base.json", KB.to_json().encode("utf-8"), file_name="knowledge_base.json")
    st.caption("Exports are generated in the browser session only. Unpublished data are never written to the repository.")

    st.subheader("Calibration summary report")
    st.caption("Distinguishes 'fit succeeded' from 'chemical interpretation validated' (the latter is never claimed by the app).")
    ds_name = st.text_input("Dataset name", value=st.session_state.get("xs_name", sample_name), key="cal_ds_name")
    overrides_all = {"role_reclassifications": json.loads(q["reclass"].to_json(orient="records")) if q and q.get("reclass") is not None and len(q["reclass"]) else []}
    cal = calibration_summary(ds_name, KB, git_commit(ROOT), fits, samples=st.session_state.get("xs_samples"),
                              sensitivity={k: v.get("sensitivity") for k, v in fits.items()}, overrides=overrides_all)
    cal["provenance_vocabulary"] = list(PROVENANCE)
    r1, r2 = st.columns(2)
    r1.download_button("Calibration summary (.json)", json.dumps(cal, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
                       file_name=f"{ds_name}_calibration_summary.json", mime="application/json")
    r2.download_button("Calibration summary (.html)", calibration_summary_html(cal).encode("utf-8"),
                       file_name=f"{ds_name}_calibration_summary.html", mime="text/html")
    st.caption(f"KB v{KB.version} · git {cal['git_commit'] or 'unavailable'}")

# ------------------------------------------------------------------ 6 cross-sample
with tab_x:
    st.subheader("Cross-sample comparison")
    st.caption("Side-by-side view of provider fits (paired, tentatively mapped to each sample's KB model), survey at.% and "
               "provisional quantities. Nothing here changes or adds a chemical interpretation; every value carries its "
               "provenance.")
    c1, c2, c3 = st.columns(3)
    with c1:
        rep_file = st.file_uploader("Provider report (MultiPak .pptx or Curve Fit Summary text)", type=["pptx", "txt", "csv"], key="xs_rep")
    with c2:
        raw_files = st.file_uploader("Optional raw multi-region files (*_multi.csv) for Python refits", type=["csv", "txt"],
                                     accept_multiple_files=True, key="xs_raw")
    with c3:
        meta_file = st.file_uploader("Optional sample metadata (.csv)", type=["csv"], key="xs_meta")
        st.download_button("Metadata template (.csv)", csv_bytes(smeta.template()), file_name="sample_metadata_template.csv")
        ex = ROOT / "templates" / "next_dataset_example.csv"
        if ex.exists():
            st.download_button("Example: next five-sample dataset (.csv)", ex.read_bytes(), file_name="next_dataset_example.csv")
    st.session_state["xs_name"] = st.text_input("Dataset name", value=st.session_state.get("xs_name", "dataset"), key="xs_name_in")
    if rep_file is not None:
        try:
            meta_df = smeta.read(meta_file.getvalue()) if meta_file is not None else None
            bands_x, act_x = ds.read_provider_report(rep_file.getvalue(), rep_file.name)
            files = list(dict.fromkeys(bands_x["data_file"]))
            raw_by_stem = {}
            for f in raw_files or []:
                stem = f.name.split("_", 1)[0]
                try:
                    raw_by_stem[stem] = MULTIREGION_IMPORTERS[list(MULTIREGION_IMPORTERS)[0]](f.getvalue())
                except Exception as e:
                    say("INVALID", f"{f.name}: {e}")
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
            say("INVALID", f"Dataset assembly failed: {e}")
            samples_x = []
    else:
        samples_x = st.session_state.get("xs_samples") or []
    if samples_x:
        info = pd.DataFrame([{"sample": s.sample, "data_file": s.data_file, "phase": s.meta.get("phase"),
                              "composition": s.meta.get("nominal_composition"), "treatment": s.meta.get("treatment"),
                              "treatment_temperature_C": s.meta.get("treatment_temperature_C"), "structure": s.structure,
                              "profile": s.profile, "notes": "; ".join(s.notes)} for s in samples_x])
        st.dataframe(info, hide_index=True)
        g1, g2 = st.columns(2)
        group_by = g1.selectbox("Group by", ["none", "phase", "composition", "treatment"], key="xs_group")
        chosen = g2.multiselect("Samples", list(info["sample"]), default=list(info["sample"]), key="xs_pick")
        sel = [s for s in samples_x if s.sample in chosen]
        long = ds.long_table(sel)

        def _w(quantity, core=None):
            w = ds.wide(long, quantity, core)
            if w.empty or group_by == "none":
                return w
            gmap = long.drop_duplicates("sample").set_index("sample")[group_by].astype(str).to_dict()
            w.columns = pd.MultiIndex.from_tuples([(gmap.get(c, ""), c) for c in w.columns], names=[group_by, "sample"])
            return w.sort_index(axis=1)

        st.markdown("#### Provisional quantities (nominal = candidate excluded; alternative = candidate as framework; never averaged)")
        qs = ds.quant_summary(sel)
        st.dataframe(qs, hide_index=True)
        cores_x = sorted({c for s in sel for c in s.mappings})
        core_pick = st.selectbox("Core level for BE / area-fraction tables", cores_x, key="xs_core")
        if core_pick in UNSUPPORTED_REGIONS:
            say("CAUTION", f"{core_pick}: {NO_MODEL_TEXT}; provider peaks shown as imported (paired where applicable).")
        st.markdown(f"#### Fitted BE by component — {core_pick} (provider; provenance provider_paired_doublet / provider_raw_fit)")
        st.dataframe(_w("BE (eV)", core_pick))
        st.markdown(f"#### Area fraction by component — {core_pick}")
        st.dataframe(_w("area fraction", core_pick))
        if any(s.python for s in sel):
            st.markdown(f"#### Python refit — {core_pick} (area fraction; identifiability in the long table)")
            st.dataframe(_w("area fraction [python]", core_pick))
        for oc in ("O 1s", "F 1s"):
            if oc in cores_x and oc != core_pick:
                st.markdown(f"#### {oc} fractions")
                st.dataframe(_w("area fraction", oc))
        st.markdown("#### Survey at.% (provider atomic concentration table) — incl. Al residual, Cl, S, Li where measured")
        st.dataframe(_w("at.%"))
        unres = long[(long["quantity"] == "BE (eV)") & long["identifiability"].astype(str).ne("")]
        if len(unres):
            st.markdown("#### Unresolved / unmodelled provider components (labels are neutral; no assignment)")
            st.dataframe(unres[["sample", "core_level", "component", "provider_label", "value", "identifiability"]]
                         .rename(columns={"value": "BE (eV)"}), hide_index=True)
        st.download_button("Comparison long table with provenance (.csv)", csv_bytes(long.assign(app_version=__version__, kb_version=KB.version)),
                           file_name=f"{st.session_state['xs_name']}_comparison_long.csv")
        st.download_button("Provisional quantities (.csv)", csv_bytes(qs.assign(app_version=__version__, kb_version=KB.version)), file_name=f"{st.session_state['xs_name']}_provisional.csv")
