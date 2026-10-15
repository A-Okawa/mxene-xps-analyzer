"""Independent Python refit (handoff §17, §19).

The Python fit is an independent QC fit, not ground truth.  The Shirley
background is an iterative approximation and has not been calibrated against
CasaXPS/provider software.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, nnls

from .models import FitComponent

ETA0 = 0.35
_trapz = getattr(np, "trapezoid", None) or np.trapz


def pseudo_voigt_unit_area(x: np.ndarray, center: float, fwhm: float, eta: float) -> np.ndarray:
    fwhm = max(float(fwhm), 1e-6)
    eta = float(np.clip(eta, 0.0, 1.0))
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    g = np.exp(-0.5 * ((x - center) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))
    gamma = fwhm / 2.0
    lor = (gamma / np.pi) / ((x - center) ** 2 + gamma ** 2)
    return eta * lor + (1.0 - eta) * g


def component_curve(x: np.ndarray, c: FitComponent, center: float, fwhm: float, eta: float, area_primary: float) -> np.ndarray:
    cur = area_primary * pseudo_voigt_unit_area(x, center, fwhm, eta)
    if c.kind == "doublet":
        cur = cur + (area_primary / c.ratio_primary_secondary) * pseudo_voigt_unit_area(x, center + c.split, fwhm, eta)
    return cur


def total_area_factor(c: FitComponent) -> float:
    return 1.0 + 1.0 / c.ratio_primary_secondary if c.kind == "doublet" else 1.0


def _end_mean(y: np.ndarray, n: int, side: str) -> float:
    n = max(1, min(int(n), len(y) // 4 or 1))
    return float(np.mean(y[:n] if side == "low" else y[-n:]))


def linear_background(x: np.ndarray, y: np.ndarray, endpoint_avg: int = 1) -> np.ndarray:
    lo, hi = _end_mean(y, endpoint_avg, "low"), _end_mean(y, endpoint_avg, "high")
    return lo + (hi - lo) * (x - x[0]) / (x[-1] - x[0])


def shirley_background(x: np.ndarray, y: np.ndarray, endpoint_avg: int = 1, max_iter: int = 100, tol: float = 1e-6) -> np.ndarray:
    """Iterative Shirley approximation on ascending BE.

    B(E) = I_low + (I_high - I_low) * A(E_low..E) / A_total, with A computed from
    (y - B) clipped at zero.  Uses a trapezoidal cumulative integral (non-uniform grids OK).
    """
    if len(y) < 3:
        return linear_background(x, y, endpoint_avg)
    lo, hi = _end_mean(y, endpoint_avg, "low"), _end_mean(y, endpoint_avg, "high")
    b = linear_background(x, y, endpoint_avg)
    for _ in range(max_iter):
        s = np.maximum(y - b, 0.0)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(x))])
        if cum[-1] <= 0:
            return linear_background(x, y, endpoint_avg)
        new_b = lo + (hi - lo) * cum / cum[-1]
        if np.max(np.abs(new_b - b)) < tol * max(1.0, float(np.max(np.abs(y)))):
            return new_b
        b = 0.5 * b + 0.5 * new_b
    return b


BACKGROUNDS = {"shirley": shirley_background, "linear": linear_background}


@dataclass
class FitResult:
    components: pd.DataFrame
    curves: pd.DataFrame
    metrics: dict
    diagnostics: list[str] = field(default_factory=list)
    area_correlation: pd.DataFrame | None = None


def fit_core_level(x, y, components: list[FitComponent], background: str = "shirley",
                   endpoint_avg: int = 1, corr_threshold: float = 0.8, max_nfev: int = 5000) -> FitResult:
    if not components:
        raise ValueError("No fit components defined for this core level.")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 20:
        raise ValueError("At least 20 finite data points are required for fitting.")
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if background not in BACKGROUNDS:
        raise ValueError(f"background must be one of {list(BACKGROUNDS)}")

    diag: list[str] = []
    lo_x, hi_x = x[0], x[-1]
    for c in components:
        hi_line = c.center + (c.split or 0.0) if c.kind == "doublet" else c.center
        if c.center < lo_x or hi_line > hi_x:
            diag.append(f"{c.name}: starting position(s) lie outside the measured BE range "
                        f"({lo_x:.2f}–{hi_x:.2f} eV); the component may be unconstrained by data.")

    bg = BACKGROUNDS[background](x, y, endpoint_avg=endpoint_avg)
    ys = np.maximum(y - bg, 0.0)
    total_area = max(float(_trapz(ys, x)), 1.0)

    # NNLS initialization of areas at literature starting centers/FWHMs (§17.4).
    basis = np.column_stack([component_curve(x, c, c.center, c.fwhm, ETA0, 1.0) for c in components])
    try:
        a0, _ = nnls(basis, ys)
    except Exception:
        a0 = np.full(len(components), total_area / len(components))
    if not np.all(np.isfinite(a0)) or a0.sum() <= 0:
        a0 = np.full(len(components), total_area / len(components))

    p0, lb, ub = [], [], []
    for i, c in enumerate(components):
        p0 += [c.center, c.fwhm, ETA0, max(float(a0[i]), 1e-6)]
        lb += [c.bounds[0], c.fwhm_bounds[0], 0.0, 0.0]
        ub += [c.bounds[1], c.fwhm_bounds[1], 1.0, total_area * 3.0]
    p0 = np.clip(np.array(p0), np.array(lb), np.array(ub))

    def model(p):
        curves = [component_curve(x, c, *p[4 * i: 4 * i + 4]) for i, c in enumerate(components)]
        return np.sum(curves, axis=0), curves

    yscale = max(float(np.std(ys)), float(np.max(ys)) * 0.05, 1.0)
    res = least_squares(lambda p: (model(p)[0] - ys) / yscale, p0, bounds=(np.array(lb), np.array(ub)), max_nfev=max_nfev)
    sig, curves = model(res.x)
    fit_total = sig + bg
    resid = y - fit_total
    n, k = len(y), len(res.x)
    sse = float(np.sum(resid ** 2))
    sst = float(np.sum((y - y.mean()) ** 2))
    metrics = {
        "success": bool(res.success), "message": str(res.message), "n_points": n, "n_params": k,
        "RMSE": float(np.sqrt(sse / n)), "R2": 1.0 - sse / sst if sst > 0 else np.nan,
        "AIC": float(n * np.log(max(sse / n, 1e-300)) + 2 * k),
        "BIC": float(n * np.log(max(sse / n, 1e-300)) + k * np.log(n)),
        "background": background, "endpoint_avg": endpoint_avg,
        "note": "Independent QC fit; Shirley background is an uncalibrated approximation." if background == "shirley"
        else "Independent QC fit.",
    }

    rows, tot = [], 0.0
    span = np.array(ub) - np.array(lb)
    for i, c in enumerate(components):
        ctr, fw, eta, ap = res.x[4 * i: 4 * i + 4]
        area = ap * total_area_factor(c)
        tot += area
        at_bound = []
        for j, nm in enumerate(("BE", "FWHM")):
            q = 4 * i + j
            if min(res.x[q] - lb[q], ub[q] - res.x[q]) <= 1e-3 * span[q]:
                at_bound.append(nm)
        rows.append({"component": c.name, "BE": ctr, "FWHM": fw, "eta": eta, "area": area, "area_primary": ap,
                     "role": c.role, "assignment_status": c.assignment_status, "source": c.source,
                     "at_bound": ",".join(at_bound)})
    comp = pd.DataFrame(rows)
    comp["area_fraction"] = comp["area"] / tot if tot > 0 else np.nan

    for _, r in comp.iterrows():
        if r["at_bound"]:
            diag.append(f"{r['component']}: {r['at_bound']} at the edge of its allowed window – the model window "
                        "or component choice may be inappropriate for this spectrum.")
        if tot > 0 and r["area"] / tot < 0.01:
            diag.append(f"{r['component']}: area < 1 % of the fitted envelope – the data may not require this "
                        "component (do not report it as detected without justification).")

    corr_df = None
    try:
        J = res.jac
        dof = max(n - k, 1)
        cov = np.linalg.pinv(J.T @ J) * (2 * res.cost / dof)
        ai = [4 * i + 3 for i in range(len(components))]
        sd = np.sqrt(np.clip(np.diag(cov)[ai], 0, None))
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = cov[np.ix_(ai, ai)] / np.outer(sd, sd)
        names = [c.name for c in components]
        corr_df = pd.DataFrame(corr, index=names, columns=names)
        rel_sd = sd * total_area_factor_vec(components) / np.where(comp["area"].to_numpy() > 0, comp["area"].to_numpy(), np.nan)
        comp["area_rel_stderr"] = rel_sd
        for a in range(len(names)):
            for b in range(a + 1, len(names)):
                r_ab = corr[a, b]
                if np.isfinite(r_ab) and abs(r_ab) >= corr_threshold:
                    diag.append(f"Areas of {names[a]} and {names[b]} are strongly correlated (r = {r_ab:+.2f}); "
                                "their individual fractions are not uniquely identifiable – prefer role-grouped totals.")
    except Exception as e:  # diagnostics must never break the fit
        diag.append(f"Parameter covariance could not be estimated ({e}).")

    curves_df = pd.DataFrame({"BE": x, "observed": y, "background": bg, "fit_total": fit_total, "residual": resid})
    for c, cur in zip(components, curves):
        curves_df[c.name] = cur
    return FitResult(components=comp, curves=curves_df, metrics=metrics, diagnostics=diag, area_correlation=corr_df)


def total_area_factor_vec(components: list[FitComponent]) -> np.ndarray:
    return np.array([total_area_factor(c) for c in components])


def synthesize_spectrum(x: np.ndarray, components: list[FitComponent], truth: list[tuple[float, float, float, float]],
                        baseline=(120.0, 180.0), noise_frac: float = 0.0015, seed: int = 0) -> np.ndarray:
    """Synthetic spectrum: truth = [(center, fwhm, eta, area_primary), ...]; linear baseline + Gaussian noise."""
    rng = np.random.default_rng(seed)
    y = np.linspace(baseline[0], baseline[1], len(x))
    for c, t in zip(components, truth):
        y = y + component_curve(x, c, *t)
    return y + rng.normal(0, float(np.max(y)) * noise_frac, len(x))
