"""Generate the SYNTHETIC example dataset in examples/synthetic_demo/ (software validation / onboarding only).

Usage:  python -I tools/make_example_dataset.py

Everything is generated from KB starting values plus arbitrary areas and Gaussian noise. No real or private
measurement is used or derived. The chemistry is not meant to be realistic.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mxene_xps.fitting import synthesize_spectrum  # noqa: E402
from mxene_xps.knowledge_base import default_kb  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "synthetic_demo"
KB = default_kb()

REGION = {"C 1s": ("C1s", 280.0, 292.0), "O 1s": ("O1s", 526.0, 538.0), "Nb 3d": ("Nb3d", 200.0, 214.0),
          "Ti 2p": ("Ti2p", 452.0, 468.0), "F 1s": ("F1s", 681.0, 691.0), "Al 2p": ("Al2p", 69.0, 79.0),
          "Cl 2p": ("Cl2p", 195.0, 205.0)}

SAMPLES = [
    {"stem": "SYNdemo201", "comment": "1_Nb2AlC", "profile": "nb2alc_max", "phase": "MAX", "composition": "Nb2AlC",
     "cores": ["C 1s", "O 1s", "Nb 3d", "Al 2p"], "seed": 201,
     "at": {"C 1s": 30.0, "O 1s": 45.0, "Nb 3d": 17.0, "Al 2p": 8.0}},
    {"stem": "SYNdemo202", "comment": "2_Nb2CTx", "profile": "nb2c_halim", "phase": "MXene", "composition": "Nb2CTx",
     "cores": ["C 1s", "O 1s", "Nb 3d", "F 1s", "Cl 2p"], "seed": 202,
     "at": {"C 1s": 33.0, "O 1s": 38.0, "Nb 3d": 20.0, "F 1s": 5.0, "Cl 2p": 4.0}},
    {"stem": "SYNdemo203", "comment": "3_(Nb0.6Ti0.4)2CTx", "profile": "nbti_m2c_defect_aware", "phase": "MXene",
     "composition": "(Nb0.6Ti0.4)2CTx", "cores": ["C 1s", "O 1s", "Nb 3d", "Ti 2p", "F 1s"], "seed": 203,
     "at": {"C 1s": 34.0, "O 1s": 36.0, "Nb 3d": 15.0, "Ti 2p": 8.0, "F 1s": 7.0}},
]


def truth_for(comps, rng):
    t = []
    for c in comps:
        area = float(rng.uniform(800, 4000)) if c.role not in ("metal_fluoride",) else 150.0
        t.append((c.center + float(rng.normal(0, 0.08)), c.fwhm, 0.25, area))
    return t


def phi_block(region, x, y):
    s = "Full\r\n \r\n" + f"{region}\r\n1\r\n"
    return s + "".join(f"{a:.4f},{b:.4f}\r\n" for a, b in zip(x[::-1], y[::-1]))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summaries, acts, meta, truth_all = [], [], [], {}
    for smp in SAMPLES:
        rng = np.random.default_rng(smp["seed"])
        models = KB.fit_models(smp["profile"])
        multi = ""
        truth_all[smp["comment"]] = {}
        for core in smp["cores"]:
            reg, lo, hi = REGION[core]
            x = np.linspace(lo, hi, 161)
            comps = models.get(core)
            if comps:
                tr = truth_for(comps, rng)
                y = synthesize_spectrum(x, comps, tr, baseline=(800.0, 1100.0), noise_frac=0.006, seed=smp["seed"])
                truth_all[smp["comment"]][core] = [{"component": c.name, "center": t[0], "fwhm": t[1], "area_primary": t[3]}
                                                   for c, t in zip(comps, tr)]
                bands = []
                for c, t in zip(comps, tr):
                    bands.append((t[0], t[1], t[3]))
                    if c.kind == "doublet":
                        bands.append((t[0] + c.split, t[1], t[3] / c.ratio_primary_secondary))
            else:   # unsupported region (Cl 2p): one synthetic doublet, no KB model
                g = lambda c, w: np.exp(-0.5 * ((x - c) / w) ** 2)   # noqa: E731
                y = 900 + 600 * g(199.3, 0.55) + 300 * g(200.9, 0.55) + rng.normal(0, 8, len(x))
                bands = [(199.3, 1.3, 600 * 0.55 * 2.5066), (200.9, 1.3, 300 * 0.55 * 2.5066)]
            multi += phi_block(reg, x, y)
            bands.sort()
            tot = sum(b[2] for b in bands)
            lines = [" Curve Fit Summary", "=" * 78, "Summary File: CfSpecSummary.csv 01-Nov-2026", f"Setup File: {reg}.fit",
                     "=" * 78, f"Data File: {smp['stem']}.spe (XPS SPECTRUM)   Region: {reg}",
                     f"Data File Comment: {smp['comment']}", "=" * 78,
                     "Spec Band    Pos PosSep B_FWHM  FWHM Height %Gauss  Area %Area ChiSquared"]
            for i, (be, fw, ar) in enumerate(bands, 1):
                pre = "   1 " if i == 1 else "     "
                row = (f"{pre}{i:4d} {be:7.2f} {be - bands[0][0]:6.2f} {fw:6.2f} {fw:5.2f} {ar / fw:7.0f} {80:6d} "
                       f"{ar:7.0f} {100 * ar / tot:5.2f}")
                lines.append(row + ("       1.00" if i == 1 else ""))
            summaries.append("\n".join(lines))
            if core == "Nb 3d" and smp["stem"] == "SYNdemo203":
                pd.DataFrame({"BE": x, "intensity": y}).to_csv(OUT / "SYNTHETIC_single_region_Nb3d.csv", index=False)
                prov = pd.DataFrame([{"core_level": core, "component": c.name, "BE": t[0], "FWHM": t[1],
                                      "area": t[3] * (1 + 1 / c.ratio_primary_secondary) if c.kind == "doublet" else t[3]}
                                     for c, t in zip(comps, tr)])
                prov.to_csv(OUT / "SYNTHETIC_provider_components_Nb3d.csv", index=False)
        (OUT / f"{smp['stem']}_{smp['comment'].split('_', 1)[0]}_SYNTHETIC_multi.csv").write_bytes(multi.encode())
        regs = [REGION[c][0] for c in smp["cores"]]
        acts.append("\n".join(["--------------------------", "     File Information", "--------------------------",
                               f"File Name:  {smp['stem']}.spe", "File Type:  XPS SPECTRUM", f"Comment:  {smp['comment']}",
                               "--------------------------", "Atomic Concentration Table", "--------------------------",
                               "Abscissa " + " ".join(f"{r:>10s}" for r in regs),
                               "         " + " ".join(f"{1.0:10.3f}" for _ in regs) + "    RSF",
                               "         " + " ".join(f"{10.0:10.3f}" for _ in regs) + "    CorrectedRSF",
                               "1.00     " + " ".join(f"{smp['at'][c]:10.2f}" for c in smp["cores"])]))
        meta.append({"sample_id": smp["comment"], "provider_file": smp["stem"], "nominal_composition": smp["composition"],
                     "phase": smp["phase"], "synthesis_route": "synthetic", "etchant": "", "treatment": "",
                     "treatment_temperature_C": "", "measurement_date": "", "be_reference_method": "synthetic (none)",
                     "instrument_operator_notes": "", "sputtered": "no", "comments": "SYNTHETIC demo data - not a measurement"})
    header = ("SYNTHETIC EXAMPLE DATA - generated by tools/make_example_dataset.py for software validation only.\n"
              "Not a measurement; chemistry is not realistic.\n\n")
    (OUT / "SYNTHETIC_provider_report.txt").write_text(header + "\n\n".join(summaries) + "\n\n" + "\n\n".join(acts) + "\n")
    pd.DataFrame(meta).to_csv(OUT / "SYNTHETIC_sample_metadata.csv", index=False)
    (OUT / "SYNTHETIC_truth.json").write_text(json.dumps({"note": "synthetic ground truth (KB start values + noise)",
                                                         "samples": truth_all}, indent=1))
    print("written", sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()
