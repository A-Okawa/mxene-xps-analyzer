"""Provider vs Python validation at component and chemical-role level (handoff §18, §19, §20)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .models import FitComponent, worst_status

THRESHOLD_KEYS = ("delta_be_pass_eV", "delta_be_review_eV", "fwhm_rel_pass", "fwhm_rel_review",
                  "area_fraction_abs_pass", "area_fraction_abs_review")

IDENTIFIABILITY_WARNING = ("A statistically excellent fit does not guarantee a unique decomposition among strongly "
                           "overlapping components. Prefer role-grouped conclusions when individual sub-component "
                           "areas are unstable.")


def _grade(value: float, pass_t: float, review_t: float) -> str:
    if value is None or not np.isfinite(value):
        return "CHECK"
    if value <= pass_t:
        return "PASS"
    if value <= review_t:
        return "REVIEW"
    return "FAIL"


def prepare_vendor(vendor: pd.DataFrame, components: list[FitComponent], doublet_area: str = "total") -> pd.DataFrame:
    """Normalize a provider component table for one core level.

    doublet_area: 'total' if the provider area includes both spin-orbit lines,
    'primary' if it is the primary (e.g. 2p3/2, 3d5/2) line only – then it is converted to total
    using the knowledge-base ratio for that component.
    """
    if "component" not in vendor.columns:
        raise ValueError("Provider table needs a 'component' column.")
    v = vendor.copy()
    v["component"] = v["component"].astype(str).str.strip()
    for col in ("BE", "FWHM", "area"):
        v[col] = pd.to_numeric(v[col], errors="coerce") if col in v.columns else np.nan
    kbmap = {c.name: c for c in components}
    if doublet_area == "primary":
        f = v["component"].map(lambda n: (1 + 1 / kbmap[n].ratio_primary_secondary) if n in kbmap and kbmap[n].kind == "doublet" else 1.0)
        v["area"] = v["area"] * f
    elif doublet_area != "total":
        raise ValueError("doublet_area must be 'total' or 'primary'")
    kb_role = v["component"].map(lambda n: kbmap[n].role if n in kbmap else np.nan)
    if "role" in v.columns:
        v["role"] = v["role"].replace("", np.nan).fillna(kb_role)
    else:
        v["role"] = kb_role
    v["role_kb_default"] = kb_role
    v["role"] = v["role"].fillna("other")
    denom = float(v["area"].sum(skipna=True))
    v["area_fraction"] = v["area"] / denom if denom > 0 else np.nan
    return v


def compare_components(vendor: pd.DataFrame, python_fit: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
    v = vendor[["component", "BE", "FWHM", "area", "area_fraction", "role"]].copy()
    p = python_fit[["component", "BE", "FWHM", "area", "area_fraction", "role"]].copy()
    m = pd.merge(v, p, on="component", how="outer", suffixes=("_vendor", "_python"), indicator=True)
    m["presence"] = m["_merge"].map({"both": "both", "left_only": "vendor only", "right_only": "python only"}).astype(str)
    m = m.drop(columns="_merge")
    m["delta_BE_eV"] = m["BE_python"] - m["BE_vendor"]
    m["delta_FWHM_rel"] = (m["FWHM_python"] - m["FWHM_vendor"]).abs() / m["FWHM_vendor"].abs().replace(0, np.nan)
    m["delta_area_fraction_abs"] = (m["area_fraction_python"] - m["area_fraction_vendor"]).abs()
    t = thresholds

    def row_status(r):
        if r["presence"] != "both" or pd.isna(r["BE_vendor"]) or pd.isna(r["BE_python"]):
            return pd.Series({"QC_BE": "CHECK-MISSING", "QC_FWHM": "CHECK-MISSING", "QC_area": "CHECK-MISSING", "QC": "CHECK-MISSING"})
        s_be = _grade(abs(r["delta_BE_eV"]), t["delta_be_pass_eV"], t["delta_be_review_eV"])
        s_fw = _grade(r["delta_FWHM_rel"], t["fwhm_rel_pass"], t["fwhm_rel_review"])
        s_ar = _grade(r["delta_area_fraction_abs"], t["area_fraction_abs_pass"], t["area_fraction_abs_review"])
        return pd.Series({"QC_BE": s_be, "QC_FWHM": s_fw, "QC_area": s_ar, "QC": worst_status([s_be, s_fw, s_ar])})

    m = pd.concat([m, m.apply(row_status, axis=1)], axis=1)
    return m


def compare_grouped(vendor: pd.DataFrame, python_fit: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
    vg = vendor.groupby("role")["area_fraction"].sum(min_count=1).rename("area_fraction_vendor")
    pg = python_fit.groupby("role")["area_fraction"].sum(min_count=1).rename("area_fraction_python")
    out = pd.concat([vg, pg], axis=1).reset_index().rename(columns={"index": "role"})
    out["delta_area_fraction_abs"] = (out["area_fraction_python"] - out["area_fraction_vendor"]).abs()
    out["QC"] = [
        "CHECK-MISSING" if pd.isna(a) or pd.isna(b) else _grade(d, thresholds["area_fraction_abs_pass"], thresholds["area_fraction_abs_review"])
        for a, b, d in zip(out["area_fraction_vendor"], out["area_fraction_python"], out["delta_area_fraction_abs"])
    ]
    return out


@dataclass
class ValidationSummary:
    level: str      # 'error' | 'warning' | 'success'
    status: str     # overall status word
    message: str


def summarize(component_qc: pd.DataFrame, grouped_qc: pd.DataFrame) -> ValidationSummary:
    """UI summary logic of §18.2."""
    g = set(grouped_qc["QC"])
    c = set(component_qc["QC"])
    if "FAIL" in g:
        return ValidationSummary("error", "FAIL", "Role-grouped QC has FAIL: do not interpret chemistry without reviewing the model.")
    if g & {"REVIEW", "CHECK", "CHECK-MISSING"}:
        return ValidationSummary("warning", worst_status(g), "Role-grouped QC requires review (REVIEW/CHECK/CHECK-MISSING present).")
    if c - {"PASS"}:
        return ValidationSummary("warning", "REVIEW",
                                 "Role-grouped partition is consistent, but individual components differ. Internal "
                                 "sub-component fractions are not identifiable/stable; report role-grouped quantities. "
                                 + IDENTIFIABILITY_WARNING)
    return ValidationSummary("success", "PASS", "Component-level and role-grouped comparisons pass the current QC thresholds.")


def role_disagreements(vendor_prepared: pd.DataFrame) -> list[str]:
    d = vendor_prepared[vendor_prepared["role_kb_default"].notna() & (vendor_prepared["role"] != vendor_prepared["role_kb_default"])]
    return [f"{r.component}: provider/user role '{r.role}' differs from knowledge-base role '{r.role_kb_default}'. "
            "Grouped QC uses the provider/user role; confirm this reclassification is justified." for r in d.itertuples()]
