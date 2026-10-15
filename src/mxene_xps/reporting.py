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
        "app": "mxene-xps-analyzer", "project_schema": 1,
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
