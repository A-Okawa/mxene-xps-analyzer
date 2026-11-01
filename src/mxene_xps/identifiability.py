"""Identifiability / fit-confidence layer (no chemical meaning).

For every fitted component: BE, FWHM, area, area fraction, bound proximity, parameter uncertainty, strongest
correlations, strong anticorrelation partner, stability across fitting alternatives, and a status:

  non-identifiable   area rel. std. error >= rel_stderr_non, |r_area| >= corr_non, alternative span >= alt_span_non,
                     or area fraction < min_area_fraction
  boundary-limited   BE or FWHM on (or near) a fit-window edge
  weakly identifiable  rel. std. error >= rel_stderr_weak, |r_area| >= corr_weak, or alternative span >= alt_span_weak
  identifiable       none of the above

Precedence: non-identifiable > boundary-limited > weakly identifiable > identifiable.  All individual flags are kept
as separate columns.  The status is a statistical statement about the fit only; it is never a chemical assignment.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .calibration import IDENT_THRESHOLDS

STATUSES = ("identifiable", "weakly identifiable", "boundary-limited", "non-identifiable")


def _near(value, lo, hi, frac):
    span = hi - lo
    if not np.isfinite(value) or span <= 0:
        return False
    return min(value - lo, hi - value) <= frac * span


def assess(fit, alt_spans: dict | None = None, thresholds: dict | None = None) -> pd.DataFrame:
    """Return one row per component. ``alt_spans`` = {component: area-fraction span over fitting alternatives}."""
    th = {**IDENT_THRESHOLDS, **(thresholds or {})}
    comp = fit.components.set_index("component")
    ac = fit.area_correlation
    pc = fit.param_correlation
    rows = []
    for name, r in comp.iterrows():
        b = (fit.bounds or {}).get(name, {})
        at_bound = str(r.get("at_bound") or "")
        near = [p for p, col in (("BE", "BE"), ("FWHM", "FWHM")) if p in b and _near(r[col], *b[p], th["near_bound_frac"])]
        partner, r_part = None, np.nan
        if ac is not None and name in ac.index:
            others = ac.loc[name].drop(name, errors="ignore").dropna()
            if len(others):
                partner = others.abs().idxmax()
                r_part = float(others[partner])
        p_lbl, p_r = None, np.nan
        if pc is not None:
            own = [i for i in pc.index if i.startswith(name + ":")]
            other = [i for i in pc.index if not i.startswith(name + ":")]
            if own and other:
                sub = pc.loc[own, other].abs().stack().dropna()
                if len(sub):
                    (a, bb) = sub.idxmax()
                    p_lbl, p_r = f"{a} ~ {bb}", float(pc.loc[a, bb])
        rel = float(r.get("area_rel_stderr", np.nan))
        frac = float(r["area_fraction"])
        span = (alt_spans or {}).get(name, np.nan)
        reasons_non, reasons_weak = [], []
        if np.isfinite(rel) and rel >= th["rel_stderr_non"]:
            reasons_non.append(f"area rel. SE {rel:.2f}")
        elif np.isfinite(rel) and rel >= th["rel_stderr_weak"]:
            reasons_weak.append(f"area rel. SE {rel:.2f}")
        if np.isfinite(r_part) and abs(r_part) >= th["corr_non"]:
            reasons_non.append(f"|r| {abs(r_part):.2f} with {partner}")
        elif np.isfinite(r_part) and abs(r_part) >= th["corr_weak"]:
            reasons_weak.append(f"|r| {abs(r_part):.2f} with {partner}")
        if np.isfinite(span) and span >= th["alt_span_non"]:
            reasons_non.append(f"area-fraction span {span:.3f} across alternatives")
        elif np.isfinite(span) and span >= th["alt_span_weak"]:
            reasons_weak.append(f"area-fraction span {span:.3f} across alternatives")
        if frac < th["min_area_fraction"]:
            reasons_non.append(f"area fraction {frac:.3f} < {th['min_area_fraction']:g}")
        if reasons_non:
            status, why = "non-identifiable", reasons_non
        elif at_bound or near:
            status, why = "boundary-limited", [f"at bound: {at_bound}" if at_bound else f"near bound: {','.join(near)}"]
        elif reasons_weak:
            status, why = "weakly identifiable", reasons_weak
        else:
            status, why = "identifiable", []
        rows.append({
            "component": name, "BE": r["BE"], "BE_stderr": r.get("BE_stderr", np.nan), "FWHM": r["FWHM"],
            "FWHM_stderr": r.get("FWHM_stderr", np.nan), "area": r["area"], "area_fraction": frac,
            "area_rel_stderr": rel, "at_bound": at_bound, "near_bound": ",".join(near),
            "strongest_area_corr_with": partner, "strongest_area_r": r_part,
            "strongly_anticorrelated_with": partner if np.isfinite(r_part) and r_part <= -th["corr_weak"] else None,
            "strongest_param_corr": p_lbl, "strongest_param_r": p_r,
            "alt_area_fraction_span": span, "status": status, "status_reasons": "; ".join(why),
            "provenance": "python_fit",
            "note": "Statistical fit-confidence status only; not a chemical assignment.",
        })
    return pd.DataFrame(rows)


def warnings_from(table: pd.DataFrame) -> list[str]:
    out = []
    for r in table.itertuples():
        if r.status != "identifiable":
            out.append(f"{r.component}: {r.status}" + (f" ({r.status_reasons})" if r.status_reasons else ""))
        if r.strongly_anticorrelated_with:
            out.append(f"{r.component} is strongly anticorrelated with {r.strongly_anticorrelated_with} "
                       f"(r = {r.strongest_area_r:+.2f}); report their sum rather than individual areas.")
    return out
