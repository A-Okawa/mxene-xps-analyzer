"""PHI MultiPak provider-fit import (calibration phase 1, I/O only).

Three layers are kept strictly separate:

1. ``bands``    – the provider's raw fitted peaks exactly as reported (one row per MultiPak band, with the
                  provider's own label, e.g. "Nb3d band 3").  Never modified.
2. ``doublets`` – spin-orbit lines paired into one row per chemical component using region, splitting,
                  area ratio and proximity.  The provider labels of both lines are carried along.
3. ``mapping``  – an optional, tentative link from a doublet/singlet to a knowledge-base component and role.
                  It is a separate table (method + note) and never overwrites the provider label.

The pairing tolerances below are importer heuristics for recognizing which two reported bands belong to one
spin-orbit doublet.  They are NOT the knowledge-base doublet parameters and are not used by fitting or
quantification.
"""
from __future__ import annotations

import io as _io
import re
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from .io import phi_region_to_core
from .models import FitComponent

# --- spin-orbit pairing heuristics -------------------------------------------------------------------------
# (nominal splitting eV, splitting tolerance eV). Ratio = (2j+1) degeneracy ratio of the orbital.
SO_SPLIT: dict[str, tuple[float, float]] = {
    "Nb 3d": (2.72, 0.25), "Mo 3d": (3.13, 0.25), "Zr 3d": (2.40, 0.25), "Y 3d": (2.05, 0.25),
    "Ti 2p": (5.90, 0.45), "V 2p": (7.40, 0.50), "Cr 2p": (9.30, 0.60), "Sc 2p": (4.60, 0.40),
    "Cl 2p": (1.60, 0.20), "S 2p": (1.18, 0.20), "P 2p": (0.87, 0.20),
}
# Spin-orbit components that are normally not resolved and reported as one band per chemical state.
UNRESOLVED_SO = {"Al 2p", "Si 2p", "Ga 3d", "Li 1s"}
SO_RATIO = {"p": 2.0, "d": 1.5, "f": 4.0 / 3.0}
RATIO_REL_TOL = 0.25


def so_expectation(core_level: Optional[str]) -> Optional[tuple[float, float, float]]:
    """(split, split_tol, ratio) when the provider is expected to report two lines per state, else None."""
    if not core_level or core_level in UNRESOLVED_SO:
        return None
    orb = core_level.split()[-1][-1]
    if orb not in SO_RATIO or core_level not in SO_SPLIT:
        return None
    s, tol = SO_SPLIT[core_level]
    return s, tol, SO_RATIO[orb]


# --- text extraction ---------------------------------------------------------------------------------------
_A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"


def pptx_slide_texts(data: bytes) -> list[tuple[int, str]]:
    """Return (slide number, text) for every slide; paragraphs joined by newlines."""
    out = []
    with zipfile.ZipFile(_io.BytesIO(data)) as z:
        names = [n for n in z.namelist() if re.fullmatch(r"ppt/slides/slide\d+\.xml", n)]
        for n in sorted(names, key=lambda k: int(re.search(r"(\d+)", k.rsplit("/", 1)[1]).group(1))):
            info = z.getinfo(n)
            if info.file_size > 20_000_000:
                raise ValueError(f"{n} is implausibly large; refusing to parse.")
            root = ET.fromstring(z.read(n))
            paras = ["".join(t.text or "" for t in p.iter(f"{_A}t")) for p in root.iter(f"{_A}p")]
            out.append((int(re.search(r"slide(\d+)\.xml", n).group(1)), "\n".join(paras)))
    return out


# --- Curve Fit Summary parser ------------------------------------------------------------------------------
_HDR_RE = re.compile(r"^\s*Spec\s+Band\s+Pos", re.I)
_COLMAP = {"Spec": "spec", "Band": "band", "Pos": "BE", "PosSep": "pos_sep", "B_FWHM": "B_FWHM", "FWHM": "FWHM",
           "Height": "height", "%Gauss": "gauss_pct", "Area": "area", "%Area": "area_pct", "ChiSquared": "chi2"}


def _num(t: str) -> Optional[float]:
    try:
        return float(t)
    except ValueError:
        return None


def parse_curve_fit_summaries(text: str, origin: str = "") -> pd.DataFrame:
    """Parse one or more MultiPak 'Curve Fit Summary' blocks into a raw band table (one row per band)."""
    text = text.replace("|", "\n")
    rows: list[dict] = []
    meta: dict = {"summary_index": 0}
    header: Optional[list[str]] = None
    sidx = -1
    for line in text.splitlines():
        s = line.strip()
        if "Curve Fit Summary" in s:
            sidx += 1
            meta, header = {"summary_index": sidx}, None
            continue
        m = re.match(r"Summary File:\s*(\S+)\s*(.*)", s)
        if m:
            meta["summary_date"] = m.group(2).strip() or None
            continue
        m = re.match(r"Setup File:\s*(\S+)", s)
        if m:
            meta["setup_file"] = m.group(1)
            continue
        m = re.match(r"Data File:\s*(\S+).*?Region:\s*(\S+)", s)
        if m:
            meta["data_file"], meta["region"] = m.group(1), m.group(2)
            continue
        m = re.match(r"Data File Comment:\s*(.*)", s)
        if m:
            meta["comment"] = m.group(1).strip()
            continue
        if _HDR_RE.match(s):
            header = [_COLMAP.get(t, t) for t in s.split()]
            continue
        if header is None or not s:
            continue
        toks = s.split()
        nums = [_num(t) for t in toks]
        if any(v is None for v in nums):
            header = None          # text after the band table ends it
            continue
        if len(nums) == len(header):
            rec = dict(zip(header, nums))
        elif len(nums) == len(header) - 2 and header[0] == "spec" and header[-1] == "chi2":
            rec = dict(zip(header[1:-1], nums))
        else:
            continue
        region = meta.get("region")
        band = int(rec["band"])
        rows.append({
            "origin": origin, "summary_index": meta.get("summary_index", 0), "data_file": meta.get("data_file"),
            "comment": meta.get("comment"), "setup_file": meta.get("setup_file"), "summary_date": meta.get("summary_date"),
            "region": region, "core_level": phi_region_to_core(region) if region else None,
            "band": band, "provider_label": f"{region} band {band}",
            **{k: rec.get(k) for k in ("BE", "pos_sep", "B_FWHM", "FWHM", "height", "gauss_pct", "area", "area_pct")},
            "chi2": rec.get("chi2"), "provenance": "provider_raw_fit",
        })
    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError("No MultiPak 'Curve Fit Summary' band rows found.")
    # chi2 is printed once per summary (first row); propagate within the summary.
    df["chi2"] = df.groupby("summary_index")["chi2"].transform(lambda c: c.ffill().bfill())
    return df


def read_multipak_provider(data: bytes | str, name: str = "") -> pd.DataFrame:
    """Read MultiPak Curve Fit Summaries from a .pptx (provider report) or plain text/CSV export."""
    if isinstance(data, bytes) and data[:2] == b"PK":
        frames = []
        for slide, txt in pptx_slide_texts(data):
            if "Curve Fit Summary" in txt:
                d = parse_curve_fit_summaries(txt, origin=f"{name} slide {slide}".strip())
                frames.append(d)
        if not frames:
            raise ValueError("No 'Curve Fit Summary' found in the presentation.")
        df = pd.concat(frames, ignore_index=True)
        df["summary_index"] = df.groupby(["origin", "summary_index"], sort=False).ngroup()
        return df
    text = data.decode("utf-8-sig", errors="replace") if isinstance(data, bytes) else data
    return parse_curve_fit_summaries(text, origin=name)


def summary_keys(bands: pd.DataFrame) -> list[tuple]:
    return list(bands.groupby(["summary_index", "data_file", "comment", "region"], sort=False, dropna=False).groups)


# --- spin-orbit pairing ------------------------------------------------------------------------------------
def pair_spin_orbit(bands: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pair reported spin-orbit lines into one row per component.

    Returns (doublets, bands_with_link).  Pairs are chosen greedily by a score combining splitting deviation and
    area-ratio deviation (normalized by their tolerances); the lower-BE line is the primary (j = l + 1/2) line.
    Lines without an acceptable partner are kept as 'unpaired' rows (area_convention 'primary_line_only' if
    they are the only line type reported, otherwise 'unpaired_check').
    """
    b = bands.copy()
    b["doublet_id"] = None
    b["line"] = None
    out = []
    for sidx, g in b.groupby("summary_index", sort=False):
        core = g["core_level"].iloc[0]
        exp = so_expectation(core)
        idx = list(g.sort_values("BE").index)
        if exp is None:
            for i in idx:
                r = b.loc[i]
                did = f"S{sidx}-{r['provider_label']}"
                b.loc[i, ["doublet_id", "line"]] = [did, "single"]
                out.append(_row(r, None, did, "singlet", "singlet", None, core))
            continue
        split, stol, ratio = exp
        cands = []
        for a in idx:
            for c in idx:
                if a == c or b.at[c, "BE"] <= b.at[a, "BE"]:
                    continue
                sep = b.at[c, "BE"] - b.at[a, "BE"]
                if abs(sep - split) > stol or not (b.at[c, "area"] or 0) > 0:
                    continue
                rr = b.at[a, "area"] / b.at[c, "area"]
                if abs(rr / ratio - 1) > RATIO_REL_TOL:
                    continue
                cands.append((abs(sep - split) / stol + abs(rr / ratio - 1) / RATIO_REL_TOL, a, c))
        used: set = set()
        pairs = []
        for score, a, c in sorted(cands):
            if a in used or c in used:
                continue
            used |= {a, c}
            pairs.append((score, a, c))
        unpaired = [i for i in idx if i not in used]
        conv_unpaired = "primary_line_only" if not pairs else "unpaired_check"
        for score, a, c in pairs:
            ra, rc = b.loc[a], b.loc[c]
            did = f"S{sidx}-{ra['region']} b{int(ra['band'])}+b{int(rc['band'])}"
            b.loc[a, ["doublet_id", "line"]] = [did, "primary"]
            b.loc[c, ["doublet_id", "line"]] = [did, "secondary"]
            out.append(_row(ra, rc, did, "paired", "total_both_lines", round(score, 3), core, ratio))
        for i in unpaired:
            r = b.loc[i]
            did = f"S{sidx}-{r['provider_label']}"
            b.loc[i, ["doublet_id", "line"]] = [did, "unpaired"]
            out.append(_row(r, None, did, "unpaired", conv_unpaired, None, core, ratio))
    d = pd.DataFrame(out).sort_values(["summary_index", "BE"]).reset_index(drop=True)
    d["provenance"] = np.where(d["pairing"] == "paired", "provider_paired_doublet", "provider_raw_fit")
    return d, b


def _row(p: pd.Series, s: Optional[pd.Series], did: str, status: str, conv: str, score, core, ratio=None) -> dict:
    keep = ("origin", "summary_index", "data_file", "comment", "region", "core_level")
    r = {k: p[k] for k in keep}
    r.update({
        "doublet_id": did,
        "provider_label": p["provider_label"] if s is None else f"{p['provider_label']} + {s['provider_label']}",
        "provider_bands": str(int(p["band"])) if s is None else f"{int(p['band'])},{int(s['band'])}",
        "pairing": status, "pairing_score": score, "area_convention": conv,
        "BE": p["BE"], "BE_secondary": None if s is None else s["BE"],
        "split_obs": None if s is None else round(s["BE"] - p["BE"], 3),
        "ratio_obs": None if s is None or not s["area"] else round(p["area"] / s["area"], 3),
        "FWHM": p["FWHM"], "FWHM_secondary": None if s is None else s["FWHM"],
        "gauss_pct": p["gauss_pct"],
        "area_primary": p["area"], "area_secondary": None if s is None else s["area"],
        "area_reported": p["area"] + (0.0 if s is None else s["area"]),
    })
    if conv == "primary_line_only" and ratio:
        r["area_total_estimated"] = p["area"] * (1 + 1 / ratio)
        r["area_note"] = (f"Only one spin-orbit line reported; total estimated with the degeneracy ratio {ratio:g} "
                          "(estimate, flag in reports).")
    else:
        r["area_total_estimated"] = r["area_reported"]
        r["area_note"] = "" if conv != "unpaired_check" else (
            "No partner line found within the pairing tolerances although other lines were paired; check whether "
            "this band is a lone 3/2-type line, a satellite or a separate species.")
    return r


# --- tentative mapping to knowledge-base components (separate layer) ---------------------------------------
@dataclass
class Mapping:
    table: pd.DataFrame
    notes: list[str] = field(default_factory=list)


def map_to_kb(doublets: pd.DataFrame, components: list[FitComponent], method: str = "auto_window_tentative") -> Mapping:
    """Tentatively link provider rows of ONE summary to KB components by fit window and proximity.

    Each KB component is used at most once (nearest first).  Unmatched provider rows keep mapped_component =
    None.  This is a convenience for QC comparison, not a chemical assignment: provider labels are untouched.
    """
    d = doublets.copy()
    d["mapped_component"] = None
    d["mapped_role"] = None
    d["mapping_method"] = None
    d["mapping_note"] = ""
    pairs = []
    for i, r in d.iterrows():
        for c in components:
            if c.bounds[0] <= r["BE"] <= c.bounds[1]:
                pairs.append((abs(r["BE"] - c.center), i, c))
    taken_r, taken_c = set(), set()
    for dist, i, c in sorted(pairs, key=lambda t: t[0]):
        if i in taken_r or c.name in taken_c:
            continue
        taken_r.add(i)
        taken_c.add(c.name)
        d.loc[i, ["mapped_component", "mapped_role", "mapping_method"]] = [c.name, c.role, method]
        if c.kind == "doublet" and d.at[i, "pairing"] == "singlet":
            d.loc[i, "mapping_note"] = "KB component is a doublet but the provider row is a single band."
    notes = []
    um = d[d["mapped_component"].isna()]
    if len(um):
        notes.append("Provider component(s) outside every KB fit window (kept under the provider label, role "
                     "'other' in QC): " + ", ".join(f"{r.provider_label} @ {r.BE:.2f} eV" for r in um.itertuples()))
    return Mapping(table=d, notes=notes)


def apply_user_mapping(mapping: pd.DataFrame, edits: dict[str, Optional[str]], components: list[FitComponent]) -> tuple[pd.DataFrame, list[dict]]:
    """Apply user overrides {doublet_id: KB component name or None}; returns (table, change log)."""
    by = {c.name: c for c in components}
    d = mapping.copy()
    log = []
    for did, comp in edits.items():
        i = d.index[d["doublet_id"] == did]
        if len(i) != 1:
            raise ValueError(f"Unknown provider row {did!r}")
        i = i[0]
        comp = comp or None
        if comp is not None and comp not in by:
            raise ValueError(f"{comp!r} is not a component of the active model")
        old = d.at[i, "mapped_component"]
        if (old or None) == comp:
            continue
        d.loc[i, ["mapped_component", "mapped_role", "mapping_method"]] = [
            comp, by[comp].role if comp else None, "user"]
        log.append({"doublet_id": did, "provider_label": d.at[i, "provider_label"], "auto": old, "user": comp})
    dup = d["mapped_component"].dropna()
    if dup.duplicated().any():
        raise ValueError(f"KB component(s) mapped twice: {sorted(set(dup[dup.duplicated()]))}")
    return d, log


def provider_component_table(mapping: pd.DataFrame) -> pd.DataFrame:
    """Vendor table for validation: component = mapped KB name, else the provider label (role 'other').

    Area is the total over both spin-orbit lines (estimated when only one line was reported; flagged).
    """
    t = pd.DataFrame({
        "component": mapping["mapped_component"].where(mapping["mapped_component"].notna(), mapping["provider_label"]),
        "provider_label": mapping["provider_label"],
        "core_level": mapping["core_level"],
        "BE": mapping["BE"], "FWHM": mapping["FWHM"],
        "area": mapping["area_total_estimated"],
        "area_estimated": mapping["area_convention"].eq("primary_line_only"),
        "role": mapping["mapped_role"].fillna("other"),
        "mapping_method": mapping["mapping_method"].fillna("unmapped"),
    })
    return t.reset_index(drop=True)
