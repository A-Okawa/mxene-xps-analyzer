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
from mxene_xps.io import SPECTRUM_IMPORTERS, read_csv_bytes  # noqa: E402
from mxene_xps.knowledge_base import default_kb, with_overrides  # noqa: E402
from mxene_xps.models import (BE_REFERENCE_METHODS, ORTHOGONAL_METHODS, ROLE_LABELS, ROLES, SURFACE_STATES,  # noqa: E402
                              MeasurementMetadata, OrthogonalEvidence, SampleSetup)
from mxene_xps.quantification import (LEGACY_CONTROLS_DEFAULT, LEGACY_LABEL, QuantOptions,  # noqa: E402
                                      legacy_quantify_table, legacy_template, quantify_table, role_reclassifications)
from mxene_xps.reporting import evidence_summary, html_report, metadata_warnings, project_dict  # noqa: E402
from mxene_xps.request_generator import generate_request, make_peak_template, metal_core  # noqa: E402
from mxene_xps.validation import (IDENTIFIABILITY_WARNING, compare_components, compare_grouped,  # noqa: E402
                                  prepare_vendor, role_disagreements, summarize)

KB = default_kb()
FIX = ROOT / "tests" / "fixtures"
STATUS_COLOR = {"PASS": "#d8f0d8", "CHECK": "#fff0c2", "REVIEW": "#fff0c2", "CHECK-MISSING": "#ececec", "FAIL": "#f7d0d0"}

st.set_page_config(page_title="MAX/MXene XPS Analyzer", layout="wide")


def style_status(df: pd.DataFrame):
    cols = [c for c in df.columns if c == "QC" or c.startswith("QC_") or c == "metal_QC"]
    fmt = {c: "{:.4g}" for c in df.columns if pd.api.types.is_float_dtype(df[c])}
    return df.style.format(fmt, na_rep="").map(lambda v: f"background-color:{STATUS_COLOR.get(v, '')}", subset=cols)


def csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8-sig")


def split_list(text: str) -> list[str]:
    return [t.strip() for t in text.split(",") if t.strip()]


# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.header("Sample setup")
    structure = st.selectbox("Nominal structure", ["(Nb,Ti)2CTx", "Nb2CTx", "Ti3C2Tx", "Custom"])
    if structure == "Custom":
        m_target = st.number_input("Ideal M-site count", min_value=0.1, value=2.0, step=0.1)
        x_target = st.number_input("Ideal X-site count", min_value=0.1, value=1.0, step=0.1)
        profile_key = st.selectbox("Reference profile", KB.profile_keys)
    else:
        spec = KB.structures[structure]
        m_target, x_target, profile_key = float(spec["m_target"]), float(spec["x_target"]), spec["default_profile"]
        st.caption(f"M_target = {m_target:g}, X_target = {x_target:g}, profile = `{profile_key}`")
    default_m = {"(Nb,Ti)2CTx": "Nb,Ti", "Nb2CTx": "Nb", "Ti3C2Tx": "Ti"}.get(structure, "")
    m_elements = split_list(st.text_input("M-site elements (comma-separated)", default_m))
    sample_name = st.text_input("Sample name", "NbTi_MXene_01" if structure == "(Nb,Ti)2CTx" else "MXene_01")
    etchants = st.multiselect("Etching / processing chemistry",
                              ["HF", "LiF/HCl", "NaF/HCl", "KF/HCl", "NH4F/HCl", "TBAOH", "molten chloride", "other"],
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
    st.caption(f"App v{__version__} · knowledge base v{KB.version}. Runs locally; no data leave this server.")

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
ctx_warn = metadata_warnings(md)
if ctx_warn:
    with st.expander(f"⚠ Measurement-context warnings ({len(ctx_warn)})", expanded=md.sputtered):
        for w in ctx_warn:
            st.warning(w)

tab_req, tab_q, tab_v, tab_kb, tab_rep = st.tabs(
    ["1. Request generator", "2. Peak table & quantification", "3. Provider ↔ Python validation",
     "4. Knowledge base", "5. Report & export"])

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
        st.info(LEGACY_LABEL)
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
                st.error(str(e))
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
            st.warning("Roles changed from the knowledge-base default. These are user assignments and are recorded in "
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
                    st.warning(f"RSF required (no 1.0 default): {', '.join(miss)}")

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
                               oxide_f_correction_factor=fcorr, metal_qc_tolerance=tol, nominal_metal_stoich=nominal)
            try:
                res, warns = quantify_table(peaks, rsf, opt)
                if len(reclass):
                    for s in warns:
                        warns[s].append(f"{int((reclass['sample'] == s).sum())} component role(s) were reclassified by the user (see reclassification log).")
                st.session_state["quant"] = {"results": res, "warnings": warns, "model": "role-based", "profile": profile_key,
                                             "options": opt.__dict__, "rsf": rsf, "table": peaks, "reclass": reclass}
            except Exception as e:
                st.error(str(e))

    q = st.session_state.get("quant")
    if q and q["model"] == ("legacy" if use_legacy else "role-based") and q.get("profile") == profile_key:
        st.subheader("Results")
        res = q["results"]
        show = (res.set_index("sample").T
                .map(lambda x: f"{x:.6g}" if isinstance(x, (float, np.floating)) else x)
                .reset_index().rename(columns={"index": "quantity"}))
        st.dataframe(show, hide_index=True, height=min(900, 36 * (len(show) + 1)))
        st.download_button("Download quantification results (.csv)", csv_bytes(res), file_name=f"{sample_name}_quantification.csv")
        for s, ws in q["warnings"].items():
            for w in ws:
                st.warning(f"{s}: {w}")
        if q["model"] == "role-based":
            st.info("`apparent_M_deficiency_delta` is a surface-sensitive XPS index after carbide-C normalization, not "
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
    cores = list(models)
    core = st.selectbox("Core level", cores)
    kb_comps = models[core]
    c1, c2 = st.columns(2)
    with c1:
        importer = st.selectbox("Raw-spectrum format", list(SPECTRUM_IMPORTERS))
        raw_file = st.file_uploader("Raw spectrum (BE, intensity; header optional; any order)", type=["csv", "txt", "dat"], key="raw")
        demo_tag = {"Ti 2p": "Ti2p", "Nb 3d": "Nb3d"}.get(core)
        use_demo = st.checkbox("Use synthetic demo spectrum (public test fixture)", value=False, disabled=demo_tag is None)
    with c2:
        vendor_file = st.file_uploader("Provider component CSV (component, BE, FWHM, area[, role, core_level])", type=["csv"], key="vend")
        doublet_conv = st.radio("Provider doublet area is…", ["total (both spin-orbit lines)", "primary line only"], horizontal=True)
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
    st.markdown("**Active starting model** — starting values and fitting windows are editable implementation "
                "bounds; literature reference values are shown separately and are never used as constraints.")

    model_df = pd.DataFrame([{"name": x.name, "kind": x.kind, "role": x.role, "status": x.assignment_status,
                              "center": x.center, "BE_low": x.bounds[0], "BE_high": x.bounds[1],
                              "FWHM": x.fwhm, "FWHM_low": x.fwhm_bounds[0], "FWHM_high": x.fwhm_bounds[1],
                              "ΔSO": x.split, "ratio": x.ratio_primary_secondary,
                              "literature reference (display only)": x.literature_summary(),
                              "source": x.source, "DOI": KB.references[x.source].get("doi")} for x in kb_comps])
    ver = st.session_state.setdefault("model_reset", 0)
    edited_model = st.data_editor(
        model_df, hide_index=True, key=f"model_{profile_key}_{core}_{ver}",
        disabled=["name", "kind", "role", "status", "ΔSO", "ratio", "literature reference (display only)", "source", "DOI"])
    if st.button("Reset model to knowledge-base values"):
        st.session_state["model_reset"] = ver + 1
        st.rerun()
    overrides = {r["name"]: {"center": r["center"], "bounds": (r["BE_low"], r["BE_high"]), "fwhm": r["FWHM"],
                             "fwhm_bounds": (r["FWHM_low"], r["FWHM_high"])} for r in edited_model.to_dict("records")}
    try:
        comps, model_log = with_overrides(kb_comps, overrides)
    except ValueError as e:
        st.error(f"Invalid model edit: {e}")
        comps, model_log = kb_comps, []
    if model_log:
        st.warning(f"{len(model_log)} starting value/window edit(s) relative to KB v{KB.version}; recorded in the "
                   "fit settings and the project export.")
        st.dataframe(pd.DataFrame(model_log).astype(str), hide_index=True)

    if st.button("Run Python fit", type="primary"):
        try:
            if raw_file is not None:
                sp = SPECTRUM_IMPORTERS[importer](raw_file.getvalue())
            elif use_demo and demo_tag:
                sp = SPECTRUM_IMPORTERS[importer]((FIX / f"synthetic_{demo_tag}_raw.csv").read_bytes())
            else:
                raise ValueError("Upload a raw spectrum or tick the demo option.")
            fr = fit_core_level(sp.be, sp.intensity, comps, background=bg, endpoint_avg=int(endavg), corr_threshold=corr_th)
            st.session_state.setdefault("fits", {})[f"{profile_key} | {core}"] = {"fit": fr, "import_notes": sp.notes,
                                                             "settings": {"background": bg, "endpoint_avg": int(endavg),
                                                                          "profile": profile_key, "kb_version": KB.version,
                                                                          "model_overrides": model_log}}
        except Exception as e:
            st.error(f"Fit failed: {e}")

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
        st.dataframe(style_status(fr.components), hide_index=True)
        if fr.diagnostics:
            st.warning("Fit diagnostics:\n\n" + "\n".join(f"- {d}" for d in fr.diagnostics))
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
        dl[3].download_button("Fit components .csv", csv_bytes(fr.components), file_name=f"{sample_name}_{core.replace(' ', '')}_python_fit.csv")
        dl[4].download_button("Curves .csv", csv_bytes(cv), file_name=f"{sample_name}_{core.replace(' ', '')}_curves.csv")
        if fr.area_correlation is not None:
            with st.expander("Area-parameter correlation matrix"):
                st.dataframe(fr.area_correlation.style.format("{:+.2f}"))

        if vendor_file is not None:
            try:
                vend = read_csv_bytes(vendor_file.getvalue())
                if "core_level" in vend.columns:
                    vend = vend[vend["core_level"].astype(str).str.strip() == core]
                vp = prepare_vendor(vend, comps, "primary" if doublet_conv.startswith("primary") else "total")
                cqc = compare_components(vp, fr.components, th)
                gqc = compare_grouped(vp, fr.components, th)
                summ = summarize(cqc, gqc)
                getattr(st, summ.level)(f"**{summ.status}** — {summ.message}")
                for msg in role_disagreements(vp):
                    st.warning(msg)
                st.markdown("**Role-grouped QC**")
                st.dataframe(style_status(gqc), hide_index=True)
                st.markdown("**Component-level QC** (ΔBE = Python − Provider)")
                st.dataframe(style_status(cqc), hide_index=True)
                st.caption(IDENTIFIABILITY_WARNING)
                d1, d2 = st.columns(2)
                d1.download_button("Component QC .csv", csv_bytes(cqc), file_name=f"{sample_name}_{core.replace(' ', '')}_validation.csv")
                d2.download_button("Grouped QC .csv", csv_bytes(gqc), file_name=f"{sample_name}_{core.replace(' ', '')}_validation_grouped.csv")
                fitstate.update({"component_qc": cqc, "grouped_qc": gqc, "summary": summ.__dict__, "thresholds": th})
            except Exception as e:
                st.error(f"Provider comparison failed: {e}")

# ------------------------------------------------------------------ 4 KB
with tab_kb:
    st.subheader("Interpretation principles")
    for p in KB.raw["principles"]:
        st.markdown(f"- {p}")
    st.subheader(f"Active profile: `{profile_key}` ({profile['material']})")
    for n in profile.get("notes", []):
        st.info(n)
    for c_core, cs in models.items():
        st.markdown(f"#### {c_core}")
        st.dataframe(pd.DataFrame([{"name": x.name, "center": x.center, "bounds": f"{x.bounds[0]}–{x.bounds[1]}",
                                    "FWHM": x.fwhm, "FWHM bounds": f"{x.fwhm_bounds[0]}–{x.fwhm_bounds[1]}",
                                    "literature reference (display only)": x.literature_summary(),
                                    "literature interpretation": x.literature_interpretation or "",
                                    "role": x.role, "status": x.assignment_status,
                                    "interpretation": x.interpretation, "source": x.source,
                                    "reference": KB.reference_text(x.source)} for x in cs]), hide_index=True)
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
    proj["fits"] = {k: {"metrics": v["fit"].metrics, "components": json.loads(v["fit"].components.to_json(orient="records")),
                        "diagnostics": v["fit"].diagnostics, "summary": v.get("summary")} for k, v in fits.items()}
    sections: list = [("Sample setup", pd.DataFrame([{k: str(v) for k, v in proj["sample_setup"].items()}]).T.reset_index()
                       .rename(columns={"index": "field", 0: "value"})),
                      ("Measurement metadata", pd.DataFrame([{k: str(v) for k, v in md.to_dict().items()}]).T.reset_index()
                       .rename(columns={"index": "field", 0: "value"})),
                      ("Measurement-context warnings", ctx_warn),
                      ("Orthogonal evidence", [ev_status] + ev_lines)]
    if q:
        sections += [("Quantification results", q["results"].set_index("sample").T.reset_index()),
                     ("Quantification warnings", [f"{s}: {w}" for s, ws in q["warnings"].items() for w in ws])]
        if q.get("reclass") is not None and len(q["reclass"]):
            sections.append(("User role reclassifications", q["reclass"]))
    for k, v in fits.items():
        sections += [(f"Python fit – {k}", v["fit"].components), (f"Fit diagnostics – {k}", v["fit"].diagnostics)]
        if "grouped_qc" in v:
            sections += [(f"Role-grouped QC – {k}", v["grouped_qc"]), (f"Component QC – {k}", v["component_qc"]),
                         (f"QC summary – {k}", [f"{v['summary']['status']}: {v['summary']['message']}"])]
    sections.append(("Literature sources (active profile)", [f"{s}: {KB.reference_text(s)}" for s in sources]))
    c1, c2, c3 = st.columns(3)
    c1.download_button("Project/session (.json)", json.dumps(proj, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
                       file_name=f"{sample_name}_xps_project.json", mime="application/json")
    c2.download_button("Consolidated report (.html)", html_report(f"XPS analysis – {sample_name}", sections).encode("utf-8"),
                       file_name=f"{sample_name}_xps_report.html", mime="text/html")
    c3.download_button("knowledge_base.json", KB.to_json().encode("utf-8"), file_name="knowledge_base.json")
    st.caption("Exports are generated in the browser session only. Unpublished data are never written to the repository.")
