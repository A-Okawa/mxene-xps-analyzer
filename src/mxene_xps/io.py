"""Raw-spectrum and peak-table import (handoff §3, §17.1, §29.7, §33).

Importers are kept in a registry so provider-specific export formats can be
added later without touching the fitting code.
"""
from __future__ import annotations

import io as _io
import re
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pandas as pd

from .models import ROLE_LABELS, infer_element

MIN_POINTS = 20


@dataclass
class Spectrum:
    be: np.ndarray          # ascending
    intensity: np.ndarray
    n_dropped: int = 0
    was_descending: bool = False
    had_header: bool = False
    notes: list[str] = field(default_factory=list)


_NUM = re.compile(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?")


def _split(line: str) -> list[str]:
    return [t for t in re.split(r"[,\t;\s]+", line.strip()) if t]


def _is_num(tok: str) -> bool:
    return bool(_NUM.fullmatch(tok))


def parse_two_column(data: bytes | str) -> Spectrum:
    """Generic two-column BE/intensity reader.

    Accepts comma/tab/semicolon/whitespace delimiters, optional header lines and
    non-numeric/NaN rows (dropped and counted). Sorted to ascending BE.
    """
    text = data.decode("utf-8-sig", errors="replace") if isinstance(data, bytes) else data
    xs, ys = [], []
    dropped, header = 0, False
    seen_numeric = False
    for line in text.splitlines():
        if not line.strip():
            continue
        toks = _split(line)
        if len(toks) >= 2 and _is_num(toks[0]) and _is_num(toks[1]):
            x, y = float(toks[0]), float(toks[1])
            if np.isfinite(x) and np.isfinite(y):
                xs.append(x); ys.append(y); seen_numeric = True
            else:
                dropped += 1
        else:
            if seen_numeric:
                dropped += 1
            else:
                header = True
    x, y = np.asarray(xs, float), np.asarray(ys, float)
    if len(x) < MIN_POINTS:
        raise ValueError(f"Only {len(x)} finite BE/intensity points found; at least {MIN_POINTS} are required.")
    descending = bool(x[0] > x[-1])
    order = np.argsort(x, kind="stable")
    sp = Spectrum(be=x[order], intensity=y[order], n_dropped=dropped, was_descending=descending, had_header=header)
    if dropped:
        sp.notes.append(f"{dropped} non-numeric/NaN row(s) were ignored.")
    if len(np.unique(sp.be)) != len(sp.be):
        sp.notes.append("Duplicate BE values present; check the export (e.g. concatenated sweeps).")
    return sp


SPECTRUM_IMPORTERS: dict[str, Callable[[bytes | str], Spectrum]] = {
    "generic two-column (BE, intensity)": parse_two_column,
}


# --- Peak tables ----------------------------------------------------------------
REQUIRED_PEAK_COLUMNS = ("sample", "core_level", "component", "area", "role")
RECOMMENDED_PEAK_COLUMNS = ("sample", "core_level", "element", "component", "BE", "FWHM", "area", "role", "notes")


@dataclass
class PeakTable:
    df: pd.DataFrame                 # only rows with a usable numeric area
    skipped_empty: pd.DataFrame      # rows with blank area (optional components left empty)
    notes: list[str] = field(default_factory=list)


def normalize_peak_table(df: pd.DataFrame) -> PeakTable:
    missing = [c for c in REQUIRED_PEAK_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Peak table is missing required column(s): {missing}")
    t = df.copy()
    t.columns = [str(c).strip() for c in t.columns]
    for col in ("sample", "core_level", "component", "role"):
        t[col] = t[col].astype("string").str.strip()
    if "element" not in t.columns:
        t["element"] = pd.NA
    t["element"] = t["element"].astype("string").str.strip()
    blank_el = t["element"].isna() | (t["element"] == "")
    t.loc[blank_el, "element"] = t.loc[blank_el, "core_level"].map(lambda c: infer_element(c) if isinstance(c, str) else pd.NA)
    for col in ("BE", "FWHM"):
        if col in t.columns:
            t[col] = pd.to_numeric(t[col], errors="coerce")

    notes: list[str] = []
    raw_area = t["area"]
    area = pd.to_numeric(raw_area, errors="coerce")
    blank = raw_area.isna() | (raw_area.astype("string").str.strip() == "")
    bad = area.isna() & ~blank
    if bad.any():
        rows = t.loc[bad, ["sample", "component"]].astype(str).agg("/".join, axis=1).tolist()
        raise ValueError(f"Non-numeric area value(s) in rows: {rows}. Fix or blank them explicitly.")
    t["area"] = area
    if (t.loc[~blank, "area"] < 0).any():
        raise ValueError("Negative component areas are not allowed.")

    t["role"] = t["role"].fillna("other")
    unknown = sorted(set(t["role"]) - set(ROLE_LABELS))
    if unknown:
        raise ValueError(f"Unknown role(s): {unknown}. Allowed roles: {list(ROLE_LABELS)}")
    if t["sample"].isna().any():
        raise ValueError("Every row needs a sample name.")

    skipped = t[blank].copy()
    if len(skipped):
        notes.append(f"{len(skipped)} component row(s) with blank area were treated as not fitted: "
                     + ", ".join(skipped["component"].astype(str)))
    return PeakTable(df=t[~blank].reset_index(drop=True), skipped_empty=skipped.reset_index(drop=True), notes=notes)


def read_csv_bytes(data: bytes) -> pd.DataFrame:
    return pd.read_csv(_io.BytesIO(data), encoding="utf-8-sig")
