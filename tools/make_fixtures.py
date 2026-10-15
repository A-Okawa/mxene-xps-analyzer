"""Generate synthetic (public, non-experimental) test spectra in tests/fixtures/.

Usage:  PYTHONPATH=src python -I tools/make_fixtures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from mxene_xps.fitting import synthesize_spectrum
from mxene_xps.knowledge_base import default_kb

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures"

CASES = {
    # truth per component: (center, fwhm, eta, primary area)
    "Ti2p": {"core": "Ti 2p", "x": (452.0, 467.0, 1601), "seed": 20261006,
             "truth": [(454.82, 1.05, 0.30, 2500), (455.75, 1.55, 0.38, 2100), (456.92, 1.72, 0.42, 1600), (459.15, 1.95, 0.35, 900)]},
    "Nb3d": {"core": "Nb 3d", "x": (202.0, 214.0, 1601), "seed": 20261007,
             "truth": [(204.95, 1.35, 0.25, 2850), (207.40, 1.55, 0.20, 860), (209.45, 1.25, 0.15, 210)]},
}


def main() -> None:
    FIX.mkdir(parents=True, exist_ok=True)
    kb = default_kb()
    meta = {}
    for tag, case in CASES.items():
        comps = kb.components("nbti_m2c_defect_aware", case["core"])
        x = np.linspace(*case["x"])
        y = synthesize_spectrum(x, comps, case["truth"], seed=case["seed"])
        pd.DataFrame({"BE": x, "intensity": y}).to_csv(FIX / f"synthetic_{tag}_raw.csv", index=False)
        # descending, header-less, whitespace-delimited variant with a junk row
        lines = [f"{a:.6f}\t{b:.6f}" for a, b in zip(x[::-1], y[::-1])]
        lines.insert(10, "nan\tnan")
        lines.insert(500, "not a number row")
        (FIX / f"synthetic_{tag}_raw_desc_noheader.txt").write_text("\n".join(lines) + "\n")
        meta[tag] = {"core": case["core"], "profile": "nbti_m2c_defect_aware",
                     "truth": [dict(zip(("component", "center", "fwhm", "eta", "area_primary"), (c.name, *t)))
                               for c, t in zip(comps, case["truth"])]}
    (FIX / "synthetic_truth.json").write_text(json.dumps(meta, indent=2))
    print("fixtures written to", FIX)


if __name__ == "__main__":
    main()
