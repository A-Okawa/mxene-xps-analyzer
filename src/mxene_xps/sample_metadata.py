"""Optional per-sample metadata schema (dataset level; never part of the scientific KB)."""
from __future__ import annotations

import io as _io
import re

import pandas as pd

FIELDS = ("sample_id", "provider_file", "nominal_composition", "phase", "synthesis_route", "etchant", "treatment",
          "treatment_temperature_C", "measurement_date", "be_reference_method", "instrument_operator_notes",
          "sputtered", "comments")
PHASES = ("MAX", "MXene")
FIELD_HELP = {
    "sample_id": "free label, e.g. '4_Nb2CTx'; matches the provider 'Data File Comment' when possible",
    "provider_file": "provider data-file stem, e.g. 'sample108' (links provider fits and raw files)",
    "nominal_composition": "e.g. Nb2AlC, (Nb0.9Ti0.1)2AlC, Nb2CTx, (Nb0.6Ti0.4)2CTx",
    "phase": "MAX or MXene",
    "treatment_temperature_C": "number or blank",
    "sputtered": "yes / no / unknown",
}


def template(n_rows: int = 5) -> pd.DataFrame:
    return pd.DataFrame([{f: "" for f in FIELDS} for _ in range(n_rows)])


def read(data: bytes | str) -> pd.DataFrame:
    text = data.decode("utf-8-sig") if isinstance(data, bytes) else data
    df = pd.read_csv(_io.StringIO(text), dtype=str).fillna("")
    return normalize(df)


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d.columns = [str(c).strip() for c in d.columns]
    missing = [f for f in ("sample_id", "phase") if f not in d.columns]
    if missing:
        raise ValueError(f"Sample metadata needs column(s) {missing}; template columns: {list(FIELDS)}")
    for f in FIELDS:
        if f not in d.columns:
            d[f] = ""
    d = d[list(FIELDS) + [c for c in d.columns if c not in FIELDS]].astype(str).apply(lambda s: s.str.strip())
    d = d[d["sample_id"] != ""].reset_index(drop=True)
    bad = sorted(set(d["phase"]) - set(PHASES))
    if bad:
        raise ValueError(f"phase must be one of {PHASES}; got {bad}")
    t = pd.to_numeric(d["treatment_temperature_C"].replace("", None), errors="coerce")
    if (d["treatment_temperature_C"] != "").any() and t[d["treatment_temperature_C"] != ""].isna().any():
        raise ValueError("treatment_temperature_C must be numeric or blank")
    if d["sample_id"].duplicated().any():
        raise ValueError("sample_id values must be unique")
    return d


def metal_fractions(composition: str) -> dict[str, float]:
    """Nominal M-site fractions from strings like '(Nb0.6Ti0.4)2CTx', 'Nb2AlC'. Unknown -> {}."""
    s = composition.replace(" ", "")
    m = re.match(r"^\(([A-Za-z0-9.]+)\)\d", s)
    if m:
        parts = re.findall(r"([A-Z][a-z]?)([0-9.]+)", m.group(1))
        tot = sum(float(x) for _, x in parts)
        return {el: float(x) / tot for el, x in parts} if tot > 0 else {}
    m = re.match(r"^([A-Z][a-z]?)\d", s)
    return {m.group(1): 1.0} if m else {}


def infer_structure(row: dict) -> str | None:
    """Map phase + nominal M-site elements to a KB structure key (user can override). Returns None if unknown."""
    els = set(metal_fractions(str(row.get("nominal_composition", ""))))
    phase = row.get("phase")
    if els == {"Nb"}:
        return "Nb2AlC" if phase == "MAX" else "Nb2CTx"
    if els == {"Nb", "Ti"}:
        return "(Nb,Ti)2AlC" if phase == "MAX" else "(Nb,Ti)2CTx"
    if els == {"Ti"} and phase == "MXene":
        return "Ti3C2Tx"
    return None
