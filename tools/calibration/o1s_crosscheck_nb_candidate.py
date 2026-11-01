"""O 1s cross-check of the Nb 3d ~205.6-205.9 eV calibration candidate (Nb_suboxide_like). KB is not changed.

Usage:
  XPS_SASAKI_DIR=/path/to/private/provider_folder PYTHONPATH=src python -I tools/calibration/o1s_crosscheck_nb_candidate.py OUT_DIR [D4_GRID_CSV]

Question: does the O 1s spectrum contain the oxide-type O that a sub-oxide (NbOx) assignment of the candidate
requires, or not (which would favour the alternative: termination-modified framework Nb)?

Stoichiometric O budget in the low-BE (oxide-type) O 1s window, in at.% (provider survey quantification):
  O_low_obs   = O_at% x (fraction of O 1s area in components with BE <= O_LOW_MAX_EV)
  O_known     = 2.5 x Nb(Nb2O5) + 2.0 x Ti(TiO2)              (fully oxidized species)
  R           = O_low_obs - O_known                          (O left for a sub-oxide)
  rho_x       = R / (x x Nb(candidate)),  x = 1 (NbO), 1.5 (Nb2O3 / Nb(III)), 2 (NbO2)

Decision rule (declared before running): over ALL decomposition combinations (Nb 3d: provider + Python D4 grid on
the provider window; O 1s: provider + Python refit, Shirley and linear),
  min rho_1 >= 1.0  -> O 1s SUPPORTS a sub-oxide assignment
  max rho_1 <  0.5  -> O 1s ARGUES AGAINST a sub-oxide assignment (favours termination-modified framework Nb)
  otherwise         -> INCONCLUSIVE.
Known biases (reported, not corrected): RSF/transmission set of the provider; attenuation (an oxide overlayer
enhances its own O and Nb signals differently, O 1s KE ~950 eV vs Nb 3d ~1280 eV); low-BE O may also contain ordered
O terminations (~530 eV), which would make R an OVER-estimate of oxide O; Nb2O5 surface O may be O-deficient.
"""
from __future__ import annotations

import os
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from mxene_xps.fitting import fit_core_level  # noqa: E402
from mxene_xps.io import parse_two_column  # noqa: E402
from mxene_xps.knowledge_base import default_kb  # noqa: E402
from mxene_xps.multipak import pair_spin_orbit, read_multipak_provider  # noqa: E402

O_LOW_MAX_EV = 530.8
X_STOICH = {"NbO (x=1)": 1.0, "Nb2O3 / Nb(III) (x=1.5)": 1.5, "NbO2 (x=2)": 2.0}
SAMPLES = {
    "4": dict(profile="nb2c_halim", stem="sasaki108", label="4_Nb2CTx", at={"Nb": 21.37, "O": 39.37, "Ti": 0.0}),
    "5": dict(profile="nbti_m2c_defect_aware", stem="sasaki110", label="5_(Nb0.6Ti0.4)2CTx",
              at={"Nb": 16.58, "O": 35.02, "Ti": 5.31}),
}


def nb_provider(bands, stem):
    d, _ = pair_spin_orbit(bands[(bands.data_file == f"{stem}.spe") & (bands.region == "Nb3d")])
    f = d.set_index("BE")["area_reported"] / d["area_reported"].sum()
    get = lambda lo, hi: float(f[(f.index >= lo) & (f.index < hi)].sum())   # noqa: E731
    return {"f_MXene": get(203.9, 205.2), "f_sub": get(205.2, 206.4), "f_ox": get(206.4, 208.5), "src": "provider"}


def ti_ox_fraction(bands, stem):
    t = bands[(bands.data_file == f"{stem}.spe") & (bands.region == "Ti2p")]
    if not len(t):
        return 0.0
    d, _ = pair_spin_orbit(t)
    return float(d.loc[d.BE >= 457.0, "area_reported"].sum() / d["area_reported"].sum())


def o_low_fractions(bands, root, stem, kb, profile):
    out = []
    p = bands[(bands.data_file == f"{stem}.spe") & (bands.region == "O1s")]
    out.append({"o_src": "provider", "f_O_low": float(p.loc[p.BE <= O_LOW_MAX_EV, "area"].sum() / p.area.sum()),
                "O_low_components": ", ".join(f"{b:.2f}" for b in p.loc[p.BE <= O_LOW_MAX_EV, "BE"])})
    sp = parse_two_column(next(root.glob(f"{stem}_*O1scurve.csv")).read_bytes())
    for bg in ("shirley", "linear"):
        c = fit_core_level(sp.be, sp.intensity, kb.components(profile, "O 1s"), background=bg).components
        low = c[c.BE <= O_LOW_MAX_EV]
        out.append({"o_src": f"python_{bg}", "f_O_low": float(low.area.sum() / c.area.sum()),
                    "O_low_components": ", ".join(f"{r.component}@{r.BE:.2f}" for r in low.itertuples())})
    return out


def run(root: Path, out: Path, d4_csv: Path | None):
    kb = default_kb()
    bands = read_multipak_provider(next(root.glob("*.pptx")).read_bytes(), "provider.pptx")
    grid = pd.read_csv(d4_csv) if d4_csv and d4_csv.exists() else None
    rows = []
    for sid, cfg in SAMPLES.items():
        nb_sets = [nb_provider(bands, cfg["stem"])]
        if grid is not None:
            g = grid[(grid["sample"] == cfg["label"]) & (grid["window"] == "provider_window")]
            for r in g.itertuples():
                nb_sets.append({"f_MXene": r.f_Nb_MXene, "f_sub": r.f_Nb_suboxide_like, "f_ox": r.f_Nb2O5,
                                "src": f"python_{r.background}_dSO{r.dSO_eV:.2f}"})
        o_sets = o_low_fractions(bands, root, cfg["stem"], kb, cfg["profile"])
        f_ti_ox = ti_ox_fraction(bands, cfg["stem"])
        for nb, o in product(nb_sets, o_sets):
            nb_ox, nb_sub = cfg["at"]["Nb"] * nb["f_ox"], cfg["at"]["Nb"] * nb["f_sub"]
            o_low = cfg["at"]["O"] * o["f_O_low"]
            o_known = 2.5 * nb_ox + 2.0 * cfg["at"]["Ti"] * f_ti_ox
            r = o_low - o_known
            row = {"sample": cfg["label"], "nb_src": nb["src"], "o_src": o["o_src"], "O_low_components": o["O_low_components"],
                   "f_Nb_sub": nb["f_sub"], "f_Nb_ox": nb["f_ox"], "f_O_low": o["f_O_low"], "f_Ti_ox": f_ti_ox,
                   "Nb_sub_at": nb_sub, "Nb_ox_at": nb_ox, "O_low_obs_at": o_low, "O_known_at": o_known, "R_at": r}
            for k, x in X_STOICH.items():
                row[f"rho {k}"] = r / (x * nb_sub) if nb_sub > 0 else np.nan
            rows.append(row)
    df = pd.DataFrame(rows)
    verdict = []
    for s, g in df.groupby("sample", sort=False):
        r1 = g["rho NbO (x=1)"]
        v = "SUPPORTS sub-oxide" if r1.min() >= 1.0 else ("ARGUES AGAINST sub-oxide" if r1.max() < 0.5 else "INCONCLUSIVE")
        verdict.append({"sample": s, "n_combinations": len(g), "rho1_min": r1.min(), "rho1_median": r1.median(),
                        "rho1_max": r1.max(), "rho2_median": g["rho NbO2 (x=2)"].median(),
                        "R_at_min": g.R_at.min(), "R_at_max": g.R_at.max(), "verdict": verdict_text(v)})
    vd = pd.DataFrame(verdict)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "o1s_crosscheck_combinations.csv", index=False)
    vd.to_csv(out / "o1s_crosscheck_verdict.csv", index=False)
    return df, vd


def verdict_text(v: str) -> str:
    return v


if __name__ == "__main__":
    root = Path(os.environ["XPS_SASAKI_DIR"]).expanduser()
    out = Path(sys.argv[1]).expanduser()
    if ROOT in out.resolve().parents or out.resolve() == ROOT:
        raise SystemExit("Write outputs outside the repository (private data).")
    d4 = Path(sys.argv[2]).expanduser() if len(sys.argv) > 2 else None
    df, vd = run(root, out, d4)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 60)
    show = df.groupby(["sample", "o_src"]).agg(O_low=("O_low_components", "first"), f_O_low=("f_O_low", "first"),
                                               O_low_at=("O_low_obs_at", "first"), R_min=("R_at", "min"), R_max=("R_at", "max"),
                                               rho1_min=("rho NbO (x=1)", "min"), rho1_max=("rho NbO (x=1)", "max")).reset_index()
    print(show.round(3).to_string(index=False))
    print(df[df.nb_src == "provider"].round(3).drop(columns=["O_low_components"]).to_string(index=False))
    print(vd.round(3).to_string(index=False))
