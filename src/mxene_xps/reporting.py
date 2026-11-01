"""Context warnings, project/session export, and consolidated HTML report (handoff §21, §22, §24.4, §32.1)."""
from __future__ import annotations

import html
import json
from dataclasses import asdict
from datetime import datetime, timezone

import pandas as pd

from .models import MeasurementMetadata, OrthogonalEvidence, SampleSetup

LITERATURE_REFERENCING = "Fermi edge"  # Brette2024 and Halim2018Vacancy references were Fermi-edge calibrated

W_SPUTTER = ("Spectra were acquired after Ar+ sputtering / depth profiling. Sputtering can remove F groups, modify the "
             "O/OH balance, reduce or alter oxide layers, damage C/M sites, and change apparent chemical states. Do not "
             "pool these data with unsputtered surface-termination measurements.")


def metadata_warnings(md: MeasurementMetadata) -> list[str]:
    w: list[str] = []
    ref = md.be_reference_method or "unknown/other"
    if ref == "unknown/other":
        w.append("Binding-energy referencing method is unknown. Absolute BE comparison with literature windows or "
                 "other samples is unreliable until the method is recorded.")
    elif ref != LITERATURE_REFERENCING:
        w.append(f"BE reference is '{ref}', whereas the seeded literature models were Fermi-edge referenced. Absolute "
                 "BE offsets relative to the starting windows may reflect referencing, not chemistry.")
    if ref == "adventitious C 1s":
        w.append("Adventitious-C referencing is not universally reliable for MXenes; do not force all spectra to 284.8 eV.")
    if md.sputtered:
        w.append(W_SPUTTER)
    missing = [k for k in ("xray_source", "pass_energy_eV", "step_eV", "charge_neutralization") if getattr(md, k) in (None, "")]
    if missing:
        w.append(f"Measurement metadata not yet recorded: {missing}.")
    return w


def evidence_summary(evidence: list[OrthogonalEvidence]) -> tuple[list[str], str]:
    avail = [e for e in evidence if e.available]
    lines = [f"{e.method}: {e.result_summary or '(no summary)'}" + (f" — supports: {', '.join(e.supports)}" if e.supports else "")
             for e in avail]
    status = ("No orthogonal evidence recorded: vacancy / depletion / dopant-incorporation statements remain XPS-only "
              "hypotheses." if not avail else
              f"{len(avail)} orthogonal method(s) recorded (displayed only; no automated confidence scoring in v1).")
    return lines, status


def project_dict(setup: SampleSetup, md: MeasurementMetadata, evidence: list[OrthogonalEvidence], kb_version: str,
                 profile_key: str, rsf: dict, quant_options: dict | None, peak_table: pd.DataFrame | None,
                 quant_results: pd.DataFrame | None, warnings: dict | None, fit_settings: dict | None,
                 qc_thresholds: dict | None, sources_used: list[str], reclassifications: pd.DataFrame | None = None) -> dict:
    def df(x):
        return None if x is None else json.loads(x.to_json(orient="records"))
    return {
        "app": "mxene-xps-analyzer", "project_schema": 2, "release": _release(kb_version),
        "exported_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "kb_version": kb_version, "profile": profile_key,
        "sample_setup": asdict(setup), "measurement_metadata": md.to_dict(),
        "orthogonal_evidence": [asdict(e) for e in evidence],
        "rsf": rsf, "quantification_options": quant_options, "peak_table": df(peak_table),
        "role_reclassifications": df(reclassifications),
        "quantification_results": df(quant_results), "warnings": warnings,
        "fit_settings": fit_settings, "qc_thresholds": qc_thresholds, "literature_sources": sources_used,
    }


def html_report(title: str, sections: list[tuple[str, object]]) -> str:
    """Minimal self-contained HTML report. sections: (heading, str | list[str] | DataFrame)."""
    parts = [f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(title)}</title>",
             "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 16px;color:#222}"
             "table{border-collapse:collapse;font-size:12px;margin:.5em 0}td,th{border:1px solid #ccc;padding:3px 6px}"
             "th{background:#f3f3f3}.w{background:#fff6e0;border-left:4px solid #e0a000;padding:6px 10px;margin:4px 0}</style>",
             f"</head><body><h1>{html.escape(title)}</h1>"]
    for head, body in sections:
        parts.append(f"<h2>{html.escape(head)}</h2>")
        if isinstance(body, pd.DataFrame):
            parts.append(body.to_html(index=False, float_format=lambda v: f"{v:.4g}", na_rep=""))
        elif isinstance(body, list):
            parts += [f"<div class='w'>{html.escape(str(x))}</div>" for x in body] or ["<p>(none)</p>"]
        else:
            parts.append(f"<pre>{html.escape(str(body))}</pre>")
    parts.append("</body></html>")
    return "".join(parts)


# ---------------------------------------------------------------------------------------------------------------
# Calibration summary report
def git_commit(repo_dir) -> str | None:
    import subprocess
    try:
        out = subprocess.run(["git", "-C", str(repo_dir), "rev-parse", "--short", "HEAD"], capture_output=True, text=True, timeout=5)
        dirty = subprocess.run(["git", "-C", str(repo_dir), "status", "--porcelain"], capture_output=True, text=True, timeout=5)
        return (out.stdout.strip() + ("+dirty" if dirty.stdout.strip() else "")) if out.returncode == 0 else None
    except Exception:
        return None


INTERPRETATION_STATUS = ("chemical interpretation NOT validated: component roles are literature/transferred candidates or "
                         "calibration candidates; no orthogonal evidence has validated them.")


def fits_export(fits: dict) -> dict:
    """JSON-ready export of session fits; every value keeps its provenance."""
    def rec(x):
        return None if x is None else json.loads(x.to_json(orient="records"))
    out = {}
    for k, v in (fits or {}).items():
        fr = v["fit"]
        sres = v.get("sensitivity")
        out[k] = {
            "settings": {kk: vv for kk, vv in v.get("settings", {}).items() if kk != "model_overrides"},
            "metrics": fr.metrics,
            "components": rec(fr.components.assign(provenance="python_fit")),
            "diagnostics": [{"severity": _sev(d), "text": d} for d in fr.diagnostics],
            "identifiability": rec(v.get("ident")),
            "sensitivity": None if sres is None else {
                "alternatives": rec(sres.table.assign(provenance="python_fit")),
                "fractions": rec(sres.fractions.reset_index()),
                "derived": rec(sres.derived.reset_index()) if len(sres.derived) else [],
                "spans": rec(sres.spans.assign(provenance="derived_quantity")) if len(sres.spans) else []},
            "model_overrides": [{**o, "provenance": "user_override"} for o in v.get("settings", {}).get("model_overrides") or []],
            "provider_raw_bands": rec(v.get("provider_raw_bands")),
            "provider_mapping": rec(v.get("provider_mapping")),
            "provider_mapping_log": [{**o, "provenance": "user_override"} for o in v.get("provider_mapping_log") or []],
            "component_qc": rec(v["component_qc"].assign(provenance_vendor="provider", provenance_python="python_fit"))
            if v.get("component_qc") is not None else None,
            "grouped_qc": rec(v["grouped_qc"].assign(provenance="derived_quantity")) if v.get("grouped_qc") is not None else None,
            "qc_summary": v.get("summary"),
        }
    return out


def _sev(text):
    from .severity import classify
    return classify(text)


def calibration_summary(dataset: str, kb, commit: str | None, fits: dict, samples: list | None = None,
                        sensitivity: dict | None = None, overrides: dict | None = None) -> dict:
    """Machine-readable calibration summary.

    fits: {key: {"fit": FitResult, "settings": {...}, "summary": {...}|None, "ident": DataFrame|None, ...}}
    samples: dataset.SampleResult list (cross-sample tab) or None.
    sensitivity: {key: SensitivityResult}.
    """
    from .provenance import start_value_provenance
    from .release import release_info
    cal = kb.raw.get("calibration_status", {})
    out = {
        "dataset": dataset, "kb_version": kb.version, "git_commit": commit,
        "release": release_info(kb, commit=commit),
        "frozen_scientific_state": {
            "nb3d_suboxide_like_validated": cal.get("nb3d_suboxide_like_validated"),
            "d4_nb3d_dSO_calibration_sensitive": cal.get("d4_nb3d_dSO_calibration_sensitive"),
            "delta_finalization_blocked_by": cal.get("delta_finalization_blocked_by", []),
        },
        "levels_not_equivalent": ["fit succeeded", "statistically identifiable", "chemically validated"],
        "fits": [], "model_sensitive": [], "unresolved_components": [], "user_overrides": overrides or {},
        "references_used": {}, "samples": [],
    }
    refs = {}
    for key, f in (fits or {}).items():
        fr = f["fit"]
        prof = f["settings"].get("profile")
        comps = kb.fit_models(prof).get(f.get("core") or key.split(" | ")[-1], []) if prof else []
        for c in comps:
            for r in (c.literature_reference or []):
                refs.setdefault(r["source"], set()).add("literature_reference")
            for r in (c.transferred_reference or []):
                refs.setdefault(r["source"], set()).add("transferred_reference")
            refs.setdefault(c.source, set()).add(f"start value ({start_value_provenance(c)})")
            if c.assignment_status in ("unresolved", "calibration_candidate") or (c.calibration_candidate and not c.validated):
                out["unresolved_components"].append({"fit": key, "component": c.name, "status": c.assignment_status,
                                                     "role": c.role})
        ident = f.get("ident")
        out["fits"].append({
            "key": key, "profile": prof, "settings": {k: v for k, v in f["settings"].items() if k != "model_overrides"},
            "fit_succeeded": bool(fr.metrics.get("success")), "chemical_interpretation_validated": False,
            "interpretation_status": INTERPRETATION_STATUS,
            "metrics": {k: fr.metrics[k] for k in ("RMSE", "R2", "AIC", "BIC", "n_points")},
            "qc_summary": f.get("summary"), "qc_warnings": [{"severity": _sev(d), "text": d} for d in fr.diagnostics],
            "identifiability": None if ident is None else ident[["component", "status", "status_reasons",
                                                                  "strongly_anticorrelated_with"]].to_dict("records"),
            "model_overrides": f["settings"].get("model_overrides") or [],
            "provider_mapping_log": f.get("provider_mapping_log") or [],
        })
    for key, sres in (sensitivity or {}).items():
        if sres is None or not len(sres.spans):
            continue
        for r in sres.spans[sres.spans["model_sensitive"]].to_dict("records"):
            out["model_sensitive"].append({"fit": key, **r})
    for s in samples or []:
        out["samples"].append({"sample": s.sample, "structure": s.structure, "profile": s.profile, "meta": s.meta,
                               "notes": s.notes, "quant_warnings": [{"severity": _sev(w), "text": w} for w in s.quant_warnings],
                               "provisional": json.loads(s.provisional.to_json(orient="records")) if len(s.provisional) else [],
                               "unresolved_provider_components": [
                                   {"core": core, "label": r.comparison_label, "provider_label": r.provider_label,
                                    "BE": r.BE, "identifiability": r.identifiability}
                                   for core, m in s.mappings.items() for r in m.itertuples()
                                   if not isinstance(r.mapped_component, str)]})
    out["references_used"] = {k: {"uses": sorted(v), "citation": kb.reference_text(k),
                                  "doi": kb.references.get(k, {}).get("doi"),
                                  "citation_status": kb.references.get(k, {}).get("citation_status")} for k, v in refs.items()}
    return out


def calibration_summary_html(summary: dict) -> str:
    rel = summary.get("release", {})
    secs: list = [("Dataset / versions", pd.DataFrame([{"dataset": summary["dataset"], "app_version": rel.get("app_version"),
                                                        "kb_version": summary["kb_version"],
                                                        "git_commit": summary["git_commit"] or "unavailable",
                                                        "python": rel.get("python"), "generated_at": rel.get("generated_at")}])),
                  ("Interpretation levels (not equivalent)", [f"{k}: {v}" for k, v in rel.get("interpretation_levels", {}).items()]),
                  ("Frozen scientific state", [f"{k}: {v}" for k, v in summary["frozen_scientific_state"].items()]),
                  ("Fit succeeded vs chemical interpretation validated",
                   [f"{f['key']}: fit succeeded = {f['fit_succeeded']}; chemical interpretation validated = "
                    f"{f['chemical_interpretation_validated']}" for f in summary["fits"]] or ["no fits in this session"])]
    for f in summary["fits"]:
        secs.append((f"Fit {f['key']} — metrics / settings", pd.DataFrame([{**f["metrics"], **{k: str(v) for k, v in f["settings"].items()}}])))
        if f["qc_summary"]:
            secs.append((f"Provider-vs-Python QC — {f['key']}", [f"{f['qc_summary'].get('status')}: {f['qc_summary'].get('message')}"]))
        secs.append((f"QC warnings — {f['key']}", [f"[{w['severity']}] {w['text']}" for w in f["qc_warnings"]]))
        if f["identifiability"]:
            secs.append((f"Identifiability — {f['key']}", pd.DataFrame(f["identifiability"])))
        if f["model_overrides"] or f["provider_mapping_log"]:
            secs.append((f"User overrides — {f['key']}", [str(x) for x in f["model_overrides"] + f["provider_mapping_log"]]))
    secs.append(("Model-sensitive quantities", pd.DataFrame(summary["model_sensitive"]) if summary["model_sensitive"] else ["none flagged"]))
    secs.append(("Unresolved / calibration-candidate components",
                 pd.DataFrame(summary["unresolved_components"]) if summary["unresolved_components"] else ["none"]))
    for s in summary["samples"]:
        secs.append((f"Sample {s['sample']} — provisional quantities", pd.DataFrame(s["provisional"]) if s["provisional"] else ["none"]))
        if s["unresolved_provider_components"]:
            secs.append((f"Sample {s['sample']} — unresolved provider components", pd.DataFrame(s["unresolved_provider_components"])))
    if summary["user_overrides"]:
        secs.append(("User overrides (session)", [f"{k}: {v}" for k, v in summary["user_overrides"].items()]))
    secs.append(("Literature / transferred references used",
                 pd.DataFrame([{"key": k, **{kk: (", ".join(vv) if isinstance(vv, list) else vv) for kk, vv in v.items()}}
                               for k, v in summary["references_used"].items()]) if summary["references_used"] else ["none"]))
    secs.append(("Known scientific open issues", rel.get("scientific_open_issues", [])))
    return html_report(f"Calibration summary – {summary['dataset']}", secs)


def _release(kb_version):
    from .release import release_info
    r = release_info(None)
    r["kb_version"] = kb_version
    return r
