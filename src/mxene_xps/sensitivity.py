"""Fit-sensitivity comparison over a small predefined set of fitting alternatives.

The comparison never selects "the best chemistry".  It reports dAIC/dBIC (flagging when they are not comparable),
component area-fraction changes and, if a quantification callback is supplied, derived-quantity changes.  Spans at or
above the D4 materiality criteria mark a quantity ``model-sensitive``.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import pandas as pd

from .calibration import D4_NB3D_SPLIT_GRID, MATERIAL_AREA_FRACTION_SPAN, MATERIAL_DELTA_SPAN
from .fitting import fit_core_level
from .models import FitComponent

DERIVED_KEYS = ("M_norm", "apparent_M_deficiency_delta", "apparent_M_deficiency_delta_alt", "M_norm_alt_candidate_as_framework")


@dataclass(frozen=True)
class FitAlternative:
    label: str
    background: str = "shirley"
    window: Optional[tuple] = None          # (lo, hi) eV; None = as imported
    split: Optional[float] = None           # fixed dSO applied to every doublet of the model (calibration grid)
    basis: str = ""


def model_window(components: list[FitComponent], pad: float = 2.0) -> tuple[float, float]:
    lo = min(c.bounds[0] for c in components) - pad
    hi = max(c.bounds[1] + (c.split or 0.0) for c in components) + pad
    return lo, hi


def standard_alternatives(core: str, components: list[FitComponent], be: np.ndarray) -> list[FitAlternative]:
    """Nominal (Shirley, imported range), linear background, model-window crop (both backgrounds), and for Nb 3d the
    fixed D4 splitting grid (calibration metadata).  No other splittings are invented."""
    alts = [FitAlternative("nominal (Shirley, imported range)", basis="reference for deltas"),
            FitAlternative("linear background", background="linear", basis="background alternative")]
    lo, hi = model_window(components)
    lo, hi = max(lo, float(np.min(be))), min(hi, float(np.max(be)))
    if (lo > float(np.min(be)) + 0.2 or hi < float(np.max(be)) - 0.2) and hi - lo > 3:
        alts += [FitAlternative(f"model window {lo:.1f}-{hi:.1f} eV (Shirley)", window=(lo, hi), basis="window alternative"),
                 FitAlternative(f"model window {lo:.1f}-{hi:.1f} eV (linear)", background="linear", window=(lo, hi),
                                basis="window alternative")]
    if core == "Nb 3d" and any(c.kind == "doublet" for c in components):
        cur = {c.split for c in components if c.kind == "doublet"}
        for s in D4_NB3D_SPLIT_GRID:
            if s not in cur:
                alts.append(FitAlternative(f"fixed dSO {s:.2f} eV (D4 grid)", split=s, basis="D4 calibration grid (KB unchanged)"))
    return alts


@dataclass
class SensitivityResult:
    table: pd.DataFrame                  # one row per alternative
    fractions: pd.DataFrame              # alternatives x components (area fraction)
    derived: pd.DataFrame                # alternatives x derived quantities (may be empty)
    spans: pd.DataFrame                  # quantity, span, model_sensitive
    fits: dict = field(default_factory=dict)

    def component_spans(self) -> dict:
        return {k: v for k, v in zip(self.spans["quantity"], self.spans["span"]) if k in self.fractions.columns}


def run(be, y, components: list[FitComponent], alternatives: list[FitAlternative], endpoint_avg: int = 1,
        derived_fn: Optional[Callable[[dict], dict]] = None, spectrum_core_level: Optional[str] = None) -> SensitivityResult:
    be, y = np.asarray(be, float), np.asarray(y, float)
    rows, fr_rows, der_rows, fits = [], [], [], {}
    for a in alternatives:
        m = np.ones_like(be, bool) if a.window is None else (be >= a.window[0]) & (be <= a.window[1])
        comps = [dataclasses.replace(c, split=a.split) if (a.split is not None and c.kind == "doublet") else c for c in components]
        try:
            fr = fit_core_level(be[m], y[m], comps, background=a.background, endpoint_avg=endpoint_avg,
                                spectrum_core_level=spectrum_core_level)
        except ValueError as e:
            rows.append({"alternative": a.label, "basis": a.basis, "error": str(e)})
            continue
        fits[a.label] = fr
        fracs = dict(zip(fr.components["component"], fr.components["area_fraction"]))
        rows.append({"alternative": a.label, "basis": a.basis, "background": a.background,
                     "window": "imported" if a.window is None else f"{a.window[0]:.2f}-{a.window[1]:.2f}",
                     "dSO": a.split, "n_points": fr.metrics["n_points"], "RMSE": fr.metrics["RMSE"],
                     "AIC": fr.metrics["AIC"], "BIC": fr.metrics["BIC"],
                     "at_bound": "; ".join(f"{c}:{b}" for c, b in zip(fr.components["component"], fr.components["at_bound"]) if b),
                     "error": ""})
        fr_rows.append({"alternative": a.label, **fracs})
        if derived_fn is not None:
            try:
                der_rows.append({"alternative": a.label, **derived_fn(fracs)})
            except Exception as e:  # derived quantities are optional
                der_rows.append({"alternative": a.label, "derived_error": str(e)})
    t = pd.DataFrame(rows)
    ok = t[t["error"] == ""] if "error" in t else t
    if len(ok):
        ref = ok.iloc[0]
        t["dAIC_vs_nominal"] = t["AIC"] - ref["AIC"]
        t["dBIC_vs_nominal"] = t["BIC"] - ref["BIC"]
        t["AIC_comparable"] = t["n_points"] == ref["n_points"]
        t.loc[~t["AIC_comparable"], ["dAIC_vs_nominal", "dBIC_vs_nominal"]] = np.nan
    fr_df = pd.DataFrame(fr_rows).set_index("alternative") if fr_rows else pd.DataFrame()
    der_df = pd.DataFrame(der_rows).set_index("alternative") if der_rows else pd.DataFrame()
    spans = []
    for col in fr_df.columns:
        s = float(np.ptp(fr_df[col].dropna())) if fr_df[col].notna().any() else np.nan
        spans.append({"quantity": col, "kind": "area fraction", "span": s,
                      "model_sensitive": bool(np.isfinite(s) and s >= MATERIAL_AREA_FRACTION_SPAN)})
    for col in der_df.columns:
        vals = pd.to_numeric(der_df[col], errors="coerce").dropna()
        if not len(vals):
            continue
        s = float(np.ptp(vals))
        spans.append({"quantity": col, "kind": "derived", "span": s,
                      "model_sensitive": bool(s >= MATERIAL_DELTA_SPAN)})
    sp = pd.DataFrame(spans)
    if len(sp):
        sp["label"] = np.where(sp["model_sensitive"], "model-sensitive", "stable across alternatives")
    return SensitivityResult(table=t, fractions=fr_df, derived=der_df, spans=sp, fits=fits)


def derived_from_fractions(peak_table: pd.DataFrame, sample: str, core: str, quantify: Callable[[pd.DataFrame], dict]):
    """Return a callback: redistribute the sample's ``core`` area over components by new fractions, re-quantify.

    The total area of the core level in the peak table is kept; only components present in the fit are rescaled.
    """
    base = peak_table[peak_table["sample"].astype(str) == str(sample)].copy()
    base["area"] = pd.to_numeric(base["area"], errors="coerce")

    def fn(fracs: dict) -> dict:
        t = base.copy()
        m = (t["core_level"] == core) & t["component"].isin(list(fracs))
        tot = float(t.loc[m, "area"].sum())
        if not tot > 0:
            raise ValueError(f"no {core} areas for {sample} in the peak table")
        sub = {k: v for k, v in fracs.items() if k in set(t.loc[m, "component"])}
        norm = sum(sub.values()) or 1.0
        t.loc[m, "area"] = t.loc[m, "component"].map(lambda c: tot * sub.get(c, 0.0) / norm)
        v = quantify(t)
        return {k: v[k] for k in v if k in DERIVED_KEYS or k.endswith("_observed_fraction") or k.endswith("_observed_fraction_alt")}
    return fn
