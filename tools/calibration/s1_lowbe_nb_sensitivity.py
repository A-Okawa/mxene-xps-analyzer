"""Sensitivity of the Sample-1 low-BE Nb 3d5/2 feature (~201.8 eV) — calibration candidate only (KB unchanged).

Usage:
  XPS_SASAKI_DIR=~/Downloads/261030_sasakiXPS PYTHONPATH=src python -I tools/calibration/s1_lowbe_nb_sensitivity.py OUT_DIR

Candidate (script-local, never written to the KB): Nb_metal_like / unresolved_low_BE_Nb, confidence low,
source_type = calibration_observation (provider band at 201.80 eV, 3.6 % of Nb 3d in sample 1).

Grid
  sample     : 1_Nb2AlC (nb2alc_max) + controls 2/3 (nbti_211_max; provider did not fit such a component)
  window     : provider curve window; raw region full; raw 199.5-212.5; raw 200.5-213.5 (tight low side)
  background : shirley, linear
  line shape : symmetric pseudo-Voigt only (the fitter has no asymmetric/metallic profile). Proxy for the metallic
               asymmetry of the MAX line: with vs without the Nb_MAX_2nd tail surrogate.
  candidate  : none | doublet (dSO 2.8, ratio 1.5, KB convention) | 3d5/2 singlet only |
               two free singlets 201.8 + 204.6 eV (tests whether the 3d3/2 partner is independently required)
Reported: RMSE, AIC, BIC, dAIC/dBIC vs the same setting without candidate, candidate BE/FWHM/area fraction/rel. SE,
window-edge hits, area correlations, residual structure in 200.5-203.0 eV (mean, RMS, lag-1 autocorrelation, runs z).

Decision rule (declared before running): the candidate is called "stable" only if, for sample 1, dAIC <= -10 AND
dBIC <= -10 in every window/background (doublet, with tail surrogate), its area fraction spans < 0.02, it is not on a
window edge, |r| with Nb_MAX < 0.8, AND the controls (samples 2/3) do NOT show a comparable improvement. Even if stable,
no metallic assignment without orthogonal support (e.g. XRD bcc Nb in the same specimen).
"""
from __future__ import annotations

import dataclasses
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from mxene_xps.fitting import fit_core_level  # noqa: E402
from mxene_xps.io import crop, parse_phi_multiregion, parse_two_column  # noqa: E402
from mxene_xps.knowledge_base import default_kb  # noqa: E402
from mxene_xps.models import FitComponent  # noqa: E402

SAMPLES = {"1_Nb2AlC": ("sasaki102", "nb2alc_max"), "2_(Nb0.9Ti0.1)2AlC": ("sasaki104", "nbti_211_max"),
           "3_(Nb0.6Ti0.4)2AlC": ("sasaki106", "nbti_211_max")}
CAND_KW = dict(core_level="Nb 3d", bounds=(201.3, 202.3), fwhm=1.2, fwhm_bounds=(0.6, 2.0), role="other",
               source="Sasaki2026MultiPak", assignment_status="unresolved", confidence="low",
               calibration_candidate=True, validated=False, source_type="provider_derived",
               interpretation="Nb_metal_like / unresolved_low_BE_Nb (script-local calibration candidate; not in KB).")
RES_LO, RES_HI = 200.5, 203.0


def candidates(kind: str) -> list[FitComponent]:
    if kind == "none":
        return []
    if kind == "doublet":
        return [FitComponent(name="Nb_lowBE_cand", kind="doublet", center=201.8, split=2.8, ratio_primary_secondary=1.5, **CAND_KW)]
    if kind == "singlet_5/2":
        return [FitComponent(name="Nb_lowBE_cand", kind="single", center=201.8, **CAND_KW)]
    if kind == "free_pair":
        return [FitComponent(name="Nb_lowBE_cand", kind="single", center=201.8, **CAND_KW),
                FitComponent(name="Nb_lowBE_partner32", kind="single", center=204.6,
                             **{**CAND_KW, "bounds": (204.1, 205.1)})]
    raise ValueError(kind)


def residual_stats(cv: pd.DataFrame) -> dict:
    m = (cv["BE"] >= RES_LO) & (cv["BE"] <= RES_HI)
    r = cv.loc[m, "residual"].to_numpy()
    full = cv["residual"].to_numpy()
    if len(r) < 5:
        return {}
    ac1 = float(np.corrcoef(r[:-1], r[1:])[0, 1])
    s = np.sign(r - 0)
    s = s[s != 0]
    runs = 1 + int(np.sum(s[1:] != s[:-1]))
    n1, n2 = int(np.sum(s > 0)), int(np.sum(s < 0))
    mu = 2 * n1 * n2 / (n1 + n2) + 1
    var = 2 * n1 * n2 * (2 * n1 * n2 - n1 - n2) / ((n1 + n2) ** 2 * (n1 + n2 - 1)) if n1 + n2 > 1 else np.nan
    return {"res_mean_lowBE": float(r.mean()), "res_rms_lowBE": float(np.sqrt(np.mean(r ** 2))),
            "res_rms_ratio_lowBE_vs_all": float(np.sqrt(np.mean(r ** 2)) / np.sqrt(np.mean(full ** 2))),
            "res_ac1_lowBE": ac1, "res_runs_z_lowBE": float((runs - mu) / np.sqrt(var)) if var and var > 0 else np.nan}


def windows(root: Path, stem: str) -> dict:
    curve = parse_two_column(next(root.glob(f"{stem}_*Nb3dcurve.csv")).read_bytes())
    raw = next(r for r in parse_phi_multiregion(next(root.glob(f"{stem}_*multi.csv")).read_bytes()) if r.core_level == "Nb 3d")
    return {"provider_window": curve, "raw_full": raw, "raw_199.5-212.5": crop(raw, 199.5, 212.5),
            "raw_200.5-213.5": crop(raw, 200.5, 213.5)}


def run(root: Path, out: Path):
    kb = default_kb()
    rows, curves = [], {}
    for sample, (stem, prof) in SAMPLES.items():
        base_all = kb.components(prof, "Nb 3d")
        for wname, sp in windows(root, stem).items():
            for bg in ("shirley", "linear"):
                for tail in (True, False):
                    base = [c for c in base_all if tail or c.name != "Nb_MAX_2nd"]
                    for kind in ("none", "doublet", "singlet_5/2", "free_pair"):
                        comps = base + candidates(kind)
                        fr = fit_core_level(sp.be, sp.intensity, comps, background=bg)
                        t = fr.components.set_index("component")
                        corr = fr.area_correlation
                        row = {"sample": sample, "window": wname, "background": bg, "tail_surrogate": tail, "candidate": kind,
                               "BE_lo": float(sp.be.min()), "BE_hi": float(sp.be.max()), "n_params": fr.metrics["n_params"],
                               "RMSE": fr.metrics["RMSE"], "AIC": fr.metrics["AIC"], "BIC": fr.metrics["BIC"],
                               "f_Nb_MAX": float(t.at["Nb_MAX", "area_fraction"]),
                               "f_Nb_suboxide_like": float(t.at["Nb_suboxide_like", "area_fraction"]),
                               "at_bound": "; ".join(f"{i}:{v}" for i, v in t["at_bound"].items() if v),
                               **residual_stats(fr.curves)}
                        if kind != "none":
                            row.update({"cand_BE": float(t.at["Nb_lowBE_cand", "BE"]), "cand_FWHM": float(t.at["Nb_lowBE_cand", "FWHM"]),
                                        "cand_frac": float(t.at["Nb_lowBE_cand", "area_fraction"]),
                                        "cand_relSE": float(t.at["Nb_lowBE_cand", "area_rel_stderr"]),
                                        "r_cand_NbMAX": float(corr.at["Nb_lowBE_cand", "Nb_MAX"]) if corr is not None else np.nan,
                                        "r_cand_tail": float(corr.at["Nb_lowBE_cand", "Nb_MAX_2nd"]) if corr is not None and tail else np.nan})
                            if kind == "free_pair":
                                a5 = float(t.at["Nb_lowBE_cand", "area"]); a3 = float(t.at["Nb_lowBE_partner32", "area"])
                                row.update({"partner_BE": float(t.at["Nb_lowBE_partner32", "BE"]),
                                            "partner_frac": float(t.at["Nb_lowBE_partner32", "area_fraction"]),
                                            "partner_relSE": float(t.at["Nb_lowBE_partner32", "area_rel_stderr"]),
                                            "ratio_5/2_to_3/2": a5 / a3 if a3 > 0 else np.inf})
                        rows.append(row)
                        if sample == "1_Nb2AlC" and wname == "provider_window" and bg == "shirley" and tail and kind in ("none", "doublet"):
                            curves[kind] = fr.curves
    df = pd.DataFrame(rows)
    key = ["sample", "window", "background", "tail_surrogate"]
    ref = df[df.candidate == "none"].set_index(key)[["AIC", "BIC", "RMSE", "res_rms_lowBE", "res_ac1_lowBE"]]
    df = df.join(ref, on=key, rsuffix="_nocand")
    df["dAIC"] = df["AIC"] - df["AIC_nocand"]
    df["dBIC"] = df["BIC"] - df["BIC_nocand"]
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "s1_lowbe_grid.csv", index=False)
    for k, cv in curves.items():
        cv.to_csv(out / f"s1_curves_{k.replace('/', '')}.csv", index=False)
    return df


def verdict(df: pd.DataFrame) -> dict:
    s1 = df[(df["sample"] == "1_Nb2AlC") & (df.candidate == "doublet") & df.tail_surrogate]
    ctrl = df[(df["sample"] != "1_Nb2AlC") & (df.candidate == "doublet") & df.tail_surrogate]
    v = {
        "s1_all_dAIC<=-10": bool((s1.dAIC <= -10).all()), "s1_all_dBIC<=-10": bool((s1.dBIC <= -10).all()),
        "s1_frac_span": float(np.ptp(s1.cand_frac)), "s1_any_edge": bool(s1.at_bound.str.contains("Nb_lowBE_cand").any()),
        "s1_max_abs_r_NbMAX": float(s1.r_cand_NbMAX.abs().max()),
        "ctrl_median_dAIC": float(ctrl.dAIC.median()), "s1_median_dAIC": float(s1.dAIC.median()),
    }
    v["stable"] = (v["s1_all_dAIC<=-10"] and v["s1_all_dBIC<=-10"] and v["s1_frac_span"] < 0.02 and not v["s1_any_edge"]
                   and v["s1_max_abs_r_NbMAX"] < 0.8 and v["ctrl_median_dAIC"] > -10)
    return v


if __name__ == "__main__":
    root = Path(os.environ["XPS_SASAKI_DIR"]).expanduser()
    out = Path(sys.argv[1]).expanduser()
    if ROOT in out.resolve().parents or out.resolve() == ROOT:
        raise SystemExit("Write outputs outside the repository (private data).")
    df = run(root, out)
    pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 400)
    cols = ["sample", "window", "background", "tail_surrogate", "candidate", "RMSE", "dAIC", "dBIC", "cand_BE", "cand_FWHM",
            "cand_frac", "cand_relSE", "r_cand_NbMAX", "r_cand_tail", "res_rms_lowBE", "res_ac1_lowBE", "at_bound"]
    print(df[cols].round(3).to_string(index=False))
    fp = df[df.candidate == "free_pair"][["sample", "window", "background", "tail_surrogate", "dAIC", "cand_frac", "partner_BE",
                                          "partner_frac", "partner_relSE", "ratio_5/2_to_3/2"]]
    print(fp.round(3).to_string(index=False))
    print(verdict(df))
