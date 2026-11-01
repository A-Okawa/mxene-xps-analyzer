"""D4 sensitivity: refit real Nb 3d spectra with a fixed grid of spin-orbit splittings (KB unchanged).

Usage:
  XPS_SASAKI_DIR=~/Downloads/261030_sasakiXPS PYTHONPATH=src python -I tools/calibration/d4_nb3d_split_sensitivity.py OUT_DIR

Private data are read from XPS_SASAKI_DIR; outputs go to OUT_DIR (keep it outside the repository).

Design (user instruction 2026-11-02): the KB splitting is NOT changed. Every Nb 3d doublet of the active profile is
refitted with dSO fixed at 2.70, 2.72, 2.75, 2.80 eV (no continuous dSO fit). For each case: RMSE, AIC/BIC, Nb_MXene
and Nb_suboxide_like areas, their area correlation and relative standard error, and the resulting M_norm / delta.

delta is computed in the same way as the app's atomic_percent route:
  M_norm = (Nb_at% x f_Nb_MXene [+ Ti_at% x f_Ti_framework]) / (C_at% x f_C_framework) x X_target
Nb fractions come from the Python refit (they change across the grid); Nb/Ti/C at.% and the C and Ti framework
fractions come from the provider report (held fixed: they are not affected by the Nb 3d splitting).

Materiality criteria (declared before running): D4 is calibration-sensitive if, across the grid (primary window,
Shirley), any of f_Nb_MXene / f_Nb_suboxide_like spans >= 0.05 absolute, or delta spans >= 0.05.
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
from mxene_xps.io import parse_phi_multiregion, parse_two_column  # noqa: E402
from mxene_xps.knowledge_base import default_kb  # noqa: E402
from mxene_xps.multipak import pair_spin_orbit, read_multipak_provider  # noqa: E402

GRID = (2.70, 2.72, 2.75, 2.80)
AREA_TOL, DELTA_TOL = 0.05, 0.05
SAMPLES = {
    # sample id: (profile, data_file stem, nominal M, provider at.% Nb, C, Ti)
    "4": dict(profile="nb2c_halim", stem="sasaki108", label="4_Nb2CTx", at={"Nb": 21.37, "C": 32.17, "Ti": 0.0}),
    "5": dict(profile="nbti_m2c_defect_aware", stem="sasaki110", label="5_(Nb0.6Ti0.4)2CTx",
              at={"Nb": 16.58, "C": 35.56, "Ti": 5.31}),
}
M_TARGET, X_TARGET = 2.0, 1.0


def provider_fixed_fractions(bands: pd.DataFrame, stem: str) -> dict:
    """C framework fraction (lowest-BE C 1s band) and Ti framework fraction (doublet below 457 eV) from the provider."""
    out = {}
    c = bands[(bands.data_file == f"{stem}.spe") & (bands.region == "C1s")]
    out["f_C_framework"] = float(c.loc[c.BE.idxmin(), "area"] / c.area.sum())
    ti = bands[(bands.data_file == f"{stem}.spe") & (bands.region == "Ti2p")]
    if len(ti):
        d, _ = pair_spin_orbit(ti)
        out["f_Ti_framework"] = float(d.loc[d.BE < 457.0, "area_reported"].sum() / d.area_reported.sum())
    else:
        out["f_Ti_framework"] = 0.0
    return out


def run(root: Path, out: Path) -> pd.DataFrame:
    kb = default_kb()
    bands = read_multipak_provider(next(root.glob("*.pptx")).read_bytes(), "provider.pptx")
    rows = []
    for sid, cfg in SAMPLES.items():
        fixed = provider_fixed_fractions(bands, cfg["stem"])
        comps0 = kb.components(cfg["profile"], "Nb 3d")
        curve = next(root.glob(f"{cfg['stem']}_*Nb3dcurve.csv"))
        multi = next(root.glob(f"{cfg['stem']}_*multi.csv"))
        windows = {
            "provider_window": parse_two_column(curve.read_bytes()),
            "raw_full_region": next(r for r in parse_phi_multiregion(multi.read_bytes()) if r.core_level == "Nb 3d"),
        }
        for wname, sp in windows.items():
            for bg in ("shirley", "linear"):
                for so in GRID:
                    comps = [dataclasses.replace(c, split=so) if c.kind == "doublet" else c for c in comps0]
                    fr = fit_core_level(sp.be, sp.intensity, comps, background=bg)
                    t = fr.components.set_index("component")
                    corr = fr.area_correlation
                    f_mx, f_sub = float(t.at["Nb_MXene", "area_fraction"]), float(t.at["Nb_suboxide_like", "area_fraction"])
                    c_fw = cfg["at"]["C"] * fixed["f_C_framework"]
                    ti_fw = cfg["at"]["Ti"] * fixed["f_Ti_framework"]
                    m = (cfg["at"]["Nb"] * f_mx + ti_fw) / c_fw * X_TARGET
                    m_alt = (cfg["at"]["Nb"] * (f_mx + f_sub) + ti_fw) / c_fw * X_TARGET
                    rows.append({
                        "sample": cfg["label"], "window": wname, "background": bg, "dSO_eV": so,
                        "BE_lo": float(sp.be.min()), "BE_hi": float(sp.be.max()),
                        "RMSE": fr.metrics["RMSE"], "AIC": fr.metrics["AIC"], "BIC": fr.metrics["BIC"],
                        "BE_Nb_MXene": float(t.at["Nb_MXene", "BE"]), "BE_Nb_suboxide_like": float(t.at["Nb_suboxide_like", "BE"]),
                        "BE_Nb2O5": float(t.at["Nb2O5", "BE"]),
                        "area_Nb_MXene": float(t.at["Nb_MXene", "area"]), "area_Nb_suboxide_like": float(t.at["Nb_suboxide_like", "area"]),
                        "f_Nb_MXene": f_mx, "f_Nb_suboxide_like": f_sub, "f_Nb2O5": float(t.at["Nb2O5", "area_fraction"]),
                        "f_NbF5_like": float(t.at["NbF5_like", "area_fraction"]),
                        "r_MXene_suboxide": float(corr.at["Nb_MXene", "Nb_suboxide_like"]) if corr is not None else np.nan,
                        "r_suboxide_Nb2O5": float(corr.at["Nb_suboxide_like", "Nb2O5"]) if corr is not None else np.nan,
                        "relSE_Nb_MXene": float(t.at["Nb_MXene", "area_rel_stderr"]),
                        "relSE_Nb_suboxide_like": float(t.at["Nb_suboxide_like", "area_rel_stderr"]),
                        "at_bound": "; ".join(f"{i}:{v}" for i, v in t["at_bound"].items() if v),
                        "M_norm": m, "delta": M_TARGET - m, "M_norm_alt": m_alt, "delta_alt": M_TARGET - m_alt,
                        **fixed,
                    })
    df = pd.DataFrame(rows)
    df["dAIC_vs_best"] = df["AIC"] - df.groupby(["sample", "window", "background"])["AIC"].transform("min")
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "d4_nb3d_split_grid.csv", index=False)

    spans = (df.groupby(["sample", "window", "background"])
             .agg(span_f_MXene=("f_Nb_MXene", np.ptp), span_f_suboxide=("f_Nb_suboxide_like", np.ptp),
                  span_delta=("delta", np.ptp), span_delta_alt=("delta_alt", np.ptp),
                  best_dSO_AIC=("dAIC_vs_best", lambda s: df.loc[s.idxmin(), "dSO_eV"]),
                  max_dAIC=("dAIC_vs_best", "max"))
             .reset_index())
    spans["material"] = (spans[["span_f_MXene", "span_f_suboxide"]].max(axis=1) >= AREA_TOL) | (spans["span_delta"] >= DELTA_TOL)
    spans.to_csv(out / "d4_nb3d_split_spans.csv", index=False)
    return df, spans


if __name__ == "__main__":
    root = Path(os.environ["XPS_SASAKI_DIR"]).expanduser()
    out = Path(sys.argv[1]).expanduser()
    if ROOT in out.resolve().parents or out.resolve() == ROOT:
        raise SystemExit("Write outputs outside the repository (private data).")
    df, spans = run(root, out)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    cols = ["sample", "window", "background", "dSO_eV", "RMSE", "AIC", "BIC", "dAIC_vs_best", "f_Nb_MXene",
            "f_Nb_suboxide_like", "f_Nb2O5", "r_MXene_suboxide", "relSE_Nb_suboxide_like", "M_norm", "delta", "delta_alt", "at_bound"]
    print(df[cols].round(4).to_string(index=False))
    print(spans.round(4).to_string(index=False))
