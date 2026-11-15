# -*- coding: utf-8 -*-
"""Match uploaded raw exports to the samples of a provider report (Import & compare tab).

Same rules as tools/validation/holdout_first_pass.py: a file belongs to the sample whose provider data file
(e.g. ``okawa121.spe``) it starts with, as ``<stem>_...``. File kinds are told apart by name:

* ``*multi.csv``  PHI multi-region export: every core-level region is used for the Python refit
* ``*curve.csv``  a single region exported over the provider's curve-fit window (``<stem>_..._<Region>curve.csv``);
                  used for a core level only when no multi-region file supplies it
* ``*_su.csv``    survey: kept out of peak fitting
* anything else   read as a multi-region file if possible

UI-only glue: no change to the scientific model. Kinds and statuses are codes; app.py shows them in the screen
language.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from mxene_xps.io import parse_phi_multiregion, parse_two_column, phi_region_to_core  # noqa: E402

KINDS = ("multi", "curve", "survey", "other")


@dataclass
class Matched:
    raw_regions: dict = field(default_factory=dict)   # stem -> list[Spectrum] for dataset.build_sample
    rows: list = field(default_factory=list)          # one row per uploaded file (for the matching table)


def file_kind(name: str) -> str:
    low = name.lower()
    if low.endswith(("curve.csv", "curve.txt")):
        return "curve"
    if "multi" in low:
        return "multi"
    if re.search(r"_su\.(csv|txt)$", low) or re.search(r"_survey\.(csv|txt)$", low):
        return "survey"
    return "other"


def match_stem(name: str, stems: list[str]) -> str | None:
    """The provider data-file stem this file belongs to: the longest stem with name.startswith(stem + '_')."""
    hits = [s for s in stems if s and name.startswith(s + "_")]
    if hits:
        return max(hits, key=len)
    head = name.split("_", 1)[0]          # earlier rule (prefix before the first underscore)
    return head if head in stems else None


def match_uploads(files, data_files) -> Matched:
    """files: objects with .name and .getvalue(); data_files: provider data-file names (e.g. 'x.spe')."""
    stems = [str(f).rsplit(".", 1)[0] for f in dict.fromkeys(data_files)]
    out = Matched()
    multi: dict[str, dict] = {s: {} for s in stems}
    curve: dict[str, dict] = {s: {} for s in stems}
    for f in files or []:
        name, kind = f.name, file_kind(f.name)
        stem = match_stem(name, stems)
        row = {"file": name, "kind": kind, "stem": stem, "regions": "", "status": "", "detail": ""}
        if stem is None:
            row["status"] = "no_sample"
            out.rows.append(row)
            continue
        try:
            if kind == "survey":
                row["status"] = "survey"
            elif kind == "curve":
                sp = parse_two_column(f.getvalue())
                reg = name[len(stem) + 1:].rsplit("curve", 1)[0].split("_")[-1]
                core = phi_region_to_core(reg)
                if core is None:
                    row["status"], row["detail"] = "unknown_region", reg
                else:
                    sp.region, sp.core_level = reg, core
                    curve[stem][core] = sp
                    row["regions"], row["status"] = core, "used"
            else:
                regs = [r for r in parse_phi_multiregion(f.getvalue()) if r.core_level]
                if not regs:
                    row["status"] = "no_core_level"
                for r in regs:
                    multi[stem][r.core_level] = r
                row["regions"] = ", ".join(r.core_level for r in regs)
                row["status"] = row["status"] or "used"
        except Exception as e:  # noqa: BLE001
            row["status"], row["detail"] = "failed", str(e)
        out.rows.append(row)
    for s in stems:
        regs = dict(multi[s])
        for core, sp in curve[s].items():
            regs.setdefault(core, sp)
        out.raw_regions[s] = list(regs.values())
    # a curve file whose core level is already supplied by a multi-region file is not used for the refit
    for row in out.rows:
        if row["kind"] == "curve" and row["status"] == "used" and row["regions"] in multi.get(row["stem"], {}):
            row["status"] = "superseded"
    return out
