"""Raw-spectrum and peak-table import (handoff §3, §17.1, §29.7, §33).

Importers are kept in a registry so provider-specific export formats can be
added later without touching the fitting code.
"""
from __future__ import annotations

import io as _io
import re
from dataclasses import dataclass, field
from typing import Callable, Optional

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
    n_segments: int = 1     # monotonic runs separated by a BE discontinuity (in file order)
    max_jump_eV: float = 0.0
    region: Optional[str] = None       # e.g. "C1s" for a region taken from a multi-region file
    core_level: Optional[str] = None   # e.g. "C 1s"


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
    n_seg, max_jump = _segments(x)
    phi = split_phi_blocks(text)
    if n_seg > 1 and len(phi) > 1:
        names = ", ".join(b[0] or "?" for b in phi)
        raise ValueError(
            f"This file contains {len(phi)} separately labelled spectral blocks ({names}) with BE discontinuities "
            f"(largest jump {max_jump:.1f} eV). It looks like a multi-region export (e.g. PHI *_multi.csv); reading it "
            "as one spectrum would concatenate different core levels. Use the 'PHI multi-region' importer and pick "
            "one region.")
    descending = bool(x[0] > x[-1])
    order = np.argsort(x, kind="stable")
    sp = Spectrum(be=x[order], intensity=y[order], n_dropped=dropped, was_descending=descending, had_header=header,
                  n_segments=n_seg, max_jump_eV=max_jump)
    if dropped:
        sp.notes.append(f"{dropped} non-numeric/NaN row(s) were ignored.")
    if len(np.unique(sp.be)) != len(sp.be):
        sp.notes.append("Duplicate BE values present; check the export (e.g. concatenated sweeps).")
    if n_seg > 1:
        sp.notes.append(f"BE sequence has {n_seg} segments separated by discontinuities (largest jump "
                        f"{max_jump:.2f} eV): concatenated sweeps or regions? A single-region fit will be refused "
                        "if the BE span is implausibly wide.")
    return sp


def _segments(x_file_order: np.ndarray) -> tuple[int, float]:
    """Count monotonic runs separated by a BE discontinuity (jump > max(2 eV, 20 x median step))."""
    if len(x_file_order) < 3:
        return 1, 0.0
    d = np.diff(x_file_order)
    step = float(np.median(np.abs(d))) or 1e-9
    thr = max(JUMP_MIN_EV, JUMP_STEP_FACTOR * step)
    sign = np.sign(np.median(d)) or 1.0
    jumps = (np.abs(d) > thr) | (np.sign(d) == -sign) & (np.abs(d) > 2 * step)
    return int(jumps.sum()) + 1, float(np.max(np.abs(d)[jumps])) if jumps.any() else 0.0


JUMP_MIN_EV = 2.0
JUMP_STEP_FACTOR = 20.0

# --- Single-region plausibility (I/O guard, not a scientific parameter) -------------------------------------
MAX_SINGLE_REGION_SPAN_EV = 60.0    # hard stop: wider than any high-resolution core-level window
WARN_SINGLE_REGION_SPAN_EV = 40.0   # strong warning


def check_single_region(be: np.ndarray, core_level: Optional[str] = None,
                        max_span_eV: Optional[float] = MAX_SINGLE_REGION_SPAN_EV,
                        warn_span_eV: float = WARN_SINGLE_REGION_SPAN_EV) -> list[str]:
    """Validate that a spectrum plausibly covers ONE high-resolution region before a single-region fit.

    Raises ValueError (hard stop) when the BE span exceeds ``max_span_eV``; returns warning strings otherwise.
    ``max_span_eV=None`` disables the hard stop (only for deliberate, documented use).
    """
    be = np.asarray(be, float)
    be = be[np.isfinite(be)]
    if len(be) == 0:
        raise ValueError("No finite BE values.")
    span = float(be.max() - be.min())
    tag = f" for {core_level}" if core_level else ""
    if max_span_eV is not None and span > max_span_eV:
        raise ValueError(
            f"BE range {be.min():.2f}–{be.max():.2f} eV (span {span:.1f} eV) is implausibly wide for a single-region "
            f"fit{tag} (limit {max_span_eV:g} eV). The file is probably a survey or several concatenated regions "
            "(e.g. PHI *_multi.csv). Import it with the multi-region importer and select one region.")
    warns = []
    if span > warn_span_eV:
        warns.append(f"BE span {span:.1f} eV is unusually wide for a single high-resolution region{tag}; confirm "
                     "the file holds only one core level.")
    if len(np.unique(be)) != len(be):
        warns.append("Duplicate BE values: several sweeps or regions may be concatenated; the fit treats them as one "
                     "spectrum.")
    return warns


def crop(sp: Spectrum, lo: float, hi: float) -> Spectrum:
    """Return a copy restricted to lo <= BE <= hi (user-chosen fit window; recorded by the caller)."""
    if not lo < hi:
        raise ValueError("Fit window must satisfy low < high.")
    m = (sp.be >= lo) & (sp.be <= hi)
    if m.sum() < MIN_POINTS:
        raise ValueError(f"Only {int(m.sum())} points inside {lo:g}–{hi:g} eV; at least {MIN_POINTS} are required.")
    return Spectrum(be=sp.be[m], intensity=sp.intensity[m], n_dropped=sp.n_dropped, was_descending=sp.was_descending,
                    had_header=sp.had_header, notes=list(sp.notes) + [f"Cropped to {lo:g}–{hi:g} eV for fitting."],
                    n_segments=sp.n_segments, max_jump_eV=sp.max_jump_eV, region=sp.region, core_level=sp.core_level)


# --- PHI MultiPak multi-region ASCII export --------------------------------------------------------------
# Layout (one block per region):  "Full" / blank / "<region>" / "<n>" / "BE,intensity" rows ...
_CORE_RE = re.compile(r"^([A-Z][a-z]?)(\d[spdf])(\d/\d|\d)?$")
SURVEY_NAMES = ("Sur", "Su", "Survey")


def phi_region_to_core(name: str) -> Optional[str]:
    """'C1s' -> 'C 1s', 'Ti2p' -> 'Ti 2p', 'Ti2p3' -> 'Ti 2p'; survey/unknown -> None."""
    m = _CORE_RE.match(str(name).strip())
    return f"{m.group(1)} {m.group(2)}" if m else None


def split_phi_blocks(text: str) -> list[tuple[Optional[str], list[str], list[str]]]:
    """Split text into (region_name, header_lines, data_lines) blocks at header->data transitions."""
    blocks, header, data = [], [], []
    for line in text.splitlines():
        toks = _split(line)
        numeric = len(toks) >= 2 and _is_num(toks[0]) and _is_num(toks[1])
        if numeric:
            data.append(line)
            continue
        if data:                      # text after data: close block
            blocks.append((header, data))
            header, data = [], []
        if line.strip():
            header.append(line.strip())
    if data:
        blocks.append((header, data))
    out = []
    for h, d in blocks:
        name = next((t for t in h if t not in ("Full",) and not re.fullmatch(r"\d+", t)), None)
        out.append((name, h, d))
    return out


def parse_phi_multiregion(data: bytes | str) -> list[Spectrum]:
    """Split a PHI MultiPak multi-region (or survey) ASCII export into one Spectrum per region."""
    text = data.decode("utf-8-sig", errors="replace") if isinstance(data, bytes) else data
    blocks = split_phi_blocks(text)
    if not blocks:
        raise ValueError("No numeric BE/intensity blocks found.")
    out = []
    for name, header, lines in blocks:
        if name is None and len(blocks) > 1:
            raise ValueError("A block without a region label was found; cannot split this file safely.")
        sp = parse_two_column("\n".join(lines))
        sp.region = name
        sp.core_level = phi_region_to_core(name) if name else None
        sp.had_header = bool(header)
        if name and sp.core_level is None:
            sp.notes.append(f"Region '{name}' is not a core level (survey or unrecognized label); not for peak fitting.")
        out.append(sp)
    names = [s.region for s in out]
    if len(set(names)) != len(names):
        raise ValueError(f"Duplicate region labels in file: {names}")
    return out


SPECTRUM_IMPORTERS: dict[str, Callable[[bytes | str], Spectrum]] = {
    "generic two-column (BE, intensity)": parse_two_column,
}
# Importers that return several regions (the UI asks the user to pick one).
MULTIREGION_IMPORTERS: dict[str, Callable[[bytes | str], list[Spectrum]]] = {
    "PHI multi-region / survey (*_multi.csv, *_su.csv)": parse_phi_multiregion,
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
