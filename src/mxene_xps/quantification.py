"""Quantification engine (handoff §9, §10, §30).

Two models:
* ``quantify_sample`` – role-based, carbide-C normalized; never forces the metal
  sum to the ideal M count before computing the apparent deficiency.
* ``legacy_quantify`` – exact reproduction of the laboratory Ti3C2Tx workbook.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from .io import normalize_peak_table
from .models import DOPANT_ROLES

INPUT_MODES = ("raw_area", "corrected_area", "atomic_percent")

W_DEFICIENCY = ("Apparent metal deficiency is present after carbide-C normalization. This is a surface-sensitive XPS "
                "index and should not be reported as an absolute vacancy concentration without orthogonal evidence.")
W_EXCESS = ("Framework metal exceeds the ideal M count after carbide-C normalization. Possible causes include "
            "under-counted framework C, RSF/transmission errors, or metal components misclassified as framework; "
            "this is not interpreted structurally.")
W_UNRESOLVED_O = ("O 1s contains unresolved mixed-metal contributions. Low-/mid-BE O components may contain both oxide "
                  "and termination-related species and are not automatically counted as Tx.")
W_DOPANT = ("Dopant amount is reported by chemical role and is not automatically counted as framework M-site occupancy. "
            "Adsorption/intercalation, vacancy trapping, substitution, and secondary phases must be distinguished using "
            "additional evidence where possible.")
W_FCORR = ("Oxide-derived F correction is model-dependent. Always retain/report raw F and justify the correction model "
           "for the actual oxide/oxyfluoride chemistry.")
W_PROVISIONAL = ("apparent_M_deficiency_delta is PROVISIONAL: the Nb 3d ~205.6-205.9 eV component (Nb_suboxide_like, "
                 "calibration candidate) is not validated. It is excluded from framework Nb by default; the alternative "
                 "(termination-modified framework Nb) is reported as *_alt. Validate against O 1s and orthogonal evidence.")
W_CANDIDATE = ("Calibration-candidate role(s) present ({roles}); these components are kept out of framework/Tx sums by "
               "default and are reported separately. Do not use the affected outputs as final values until validated.")
W_MAX = ("MAX-phase structure: termination (z_*) outputs are not meaningful; O/F/Cl are surface oxide/contamination. "
         "M_norm is normalized to MAX carbide C.")
W_FRESIDUE = ("High-BE F residue is not counted as termination by default. It was included here by user choice; "
              "the assignment must be justified.")


@dataclass
class QuantOptions:
    m_target: float
    x_target: float
    input_mode: str = "raw_area"
    include_o_disordered: bool = True
    include_oh_like: bool = True
    include_cl: bool = True
    include_f_residue: bool = False
    oxide_f_correction_factor: float = 0.0
    metal_qc_tolerance: float = 0.35   # software default, user-editable; not a physical constant
    nominal_metal_stoich: dict[str, float] = field(default_factory=dict)
    # Calibration phase 2: False until the Nb 3d ~205.6-205.9 eV candidate is validated (KB calibration_status).
    nb_suboxide_validated: bool = False
    provisional_reasons: tuple = ()  # extra reasons from KB calibration_status (e.g. D4 calibration-sensitive)
    phase: str = "MXene"            # "MXene" or "MAX"
    a_target: Optional[float] = None


@dataclass
class QuantResult:
    sample: str
    values: dict
    table: pd.DataFrame
    warnings: list[str]


def corrected_areas(t: pd.DataFrame, rsf: dict[str, float], input_mode: str) -> pd.Series:
    if input_mode not in INPUT_MODES:
        raise ValueError(f"input_mode must be one of {INPUT_MODES}")
    if input_mode != "raw_area":
        return t["area"].astype(float)
    missing = []
    for core in sorted(set(t["core_level"])):
        v = rsf.get(core)
        if v is None or not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0):
            missing.append(core)
    if missing:
        raise ValueError(f"RSF is required for {missing} in raw_area mode (no default of 1.0 is assumed). "
                         "Enter the provider/instrument RSF or use corrected_area / atomic_percent.")
    return t["area"].astype(float) / t["core_level"].map(lambda c: float(rsf[c]))


def _sum(t: pd.DataFrame, *roles: str) -> float:
    return float(t.loc[t["role"].isin(roles), "corrected_area"].sum())


def quantify_sample(df: pd.DataFrame, rsf: dict[str, float], opt: QuantOptions) -> QuantResult:
    pt = normalize_peak_table(df)
    t = pt.df
    samples = sorted(set(t["sample"]) | set(pt.skipped_empty["sample"]))
    if len(samples) != 1:
        raise ValueError(f"quantify_sample expects exactly one sample, got {samples}")
    sample = samples[0]
    if t.empty:
        raise ValueError(f"{sample}: no component areas entered; fill in the provider areas first.")
    warnings: list[str] = list(pt.notes)
    t = t.copy()
    t["corrected_area"] = corrected_areas(t, rsf, opt.input_mode)

    c_fw = _sum(t, "framework_C")
    if not c_fw > 0:
        raise ValueError(f"{sample}: no positive framework_C component; carbide-C normalization is impossible.")
    scale = float(opt.x_target) / c_fw

    metal = t[t["role"] == "framework_metal"]
    m_fw = float(metal["corrected_area"].sum())
    m_norm = m_fw * scale
    delta = float(opt.m_target) - m_norm
    subox = _sum(t, "metal_suboxide_candidate")       # never folded into framework or oxide
    by_el = {str(k): float(v) * scale for k, v in metal.groupby("element")["corrected_area"].sum().items()}

    o_low = _sum(t, "termination_O")
    o_dis = _sum(t, "termination_O_disordered") if opt.include_o_disordered else 0.0
    oh = _sum(t, "termination_OH_like") if opt.include_oh_like else 0.0
    o_term = o_low + o_dis + oh
    oxide_o = _sum(t, "oxide_O")
    unres_low, unres_mid = _sum(t, "O_unresolved_low"), _sum(t, "O_unresolved_mid")

    f_low = _sum(t, "unresolved_low_BE_F")             # no chemistry assigned; never Tx
    f_term = _sum(t, "termination_F")
    f_res = _sum(t, "residue_F")
    if opt.include_f_residue:
        f_term += f_res
    f_corr = max(0.0, f_term - oxide_o * max(0.0, float(opt.oxide_f_correction_factor)))
    cl = _sum(t, "termination_Cl") if opt.include_cl else 0.0

    v: dict = {
        "C_framework_corrected": c_fw,
        "scale_X_target_per_C": scale,
        "M_target": float(opt.m_target),
        "X_target": float(opt.x_target),
        "M_framework_corrected": m_fw,
        "M_norm": m_norm,
        "apparent_M_deficiency_delta": delta,
        "metal_QC": "PASS" if abs(delta) <= opt.metal_qc_tolerance else "CHECK",
    }
    has_nb = bool((t["element"] == "Nb").any())
    provisional = (has_nb and not opt.nb_suboxide_validated) or subox > 0
    v["apparent_M_deficiency_status"] = "provisional" if provisional else "standard"
    v["M_suboxide_candidate_norm"] = subox * scale
    # Alternative interpretation (termination-modified framework Nb): reported side by side, never substituted.
    v["M_norm_alt_candidate_as_framework"] = (m_fw + subox) * scale
    v["apparent_M_deficiency_delta_alt"] = float(opt.m_target) - (m_fw + subox) * scale
    for el, val in sorted(by_el.items()):
        v[f"{el}_norm"] = val
    z_o = o_term * scale
    z_f_raw, z_f_corr, z_cl = f_term * scale, f_corr * scale, cl * scale
    v.update({
        "z_O_low": o_low * scale, "z_O_disordered": o_dis * scale, "z_OH_like": oh * scale, "z_O_total": z_o,
        "oxide_O_norm": oxide_o * scale,
        "O_unresolved_low_norm": unres_low * scale, "O_unresolved_mid_norm": unres_mid * scale,
        "z_F_raw": z_f_raw, "z_F_corrected": z_f_corr, "F_residue_norm": f_res * scale,
        "F_unresolved_lowBE_norm": f_low * scale,
        "oxide_f_correction_factor": float(opt.oxide_f_correction_factor),
        "z_Cl": z_cl,
        "z_total_raw": z_o + z_f_raw + z_cl,
        "z_total_corrected": z_o + z_f_corr + z_cl,
    })
    zt = v["z_total_corrected"]
    v["O_fraction_corr"] = z_o / zt if zt > 0 else np.nan
    v["F_fraction_corr"] = z_f_corr / zt if zt > 0 else np.nan
    v["Cl_fraction_corr"] = z_cl / zt if zt > 0 else np.nan

    # Nominal vs observed metal fractions (§9.4): selective depletion, separate from total deficiency.
    nom = {k: float(x) for k, x in (opt.nominal_metal_stoich or {}).items() if x is not None and float(x) >= 0}
    nom_total, obs_total = sum(nom.values()), sum(by_el.values())
    if len(nom) > 1 and nom_total > 0 and obs_total > 0:
        for el, n in nom.items():
            nf, of = n / nom_total, by_el.get(el, 0.0) / obs_total
            v[f"{el}_nominal_fraction"] = nf
            v[f"{el}_observed_fraction"] = of
            v[f"{el}_fraction_shift"] = of - nf
            v[f"{el}_retention_ratio_fractional"] = of / nf if nf > 0 else np.nan
            if el not in by_el:
                warnings.append(f"{el} is in the nominal M-site composition but no {el} framework_metal component was supplied.")
        extra = sorted(set(by_el) - set(nom))
        if extra:
            warnings.append(f"framework_metal component(s) for {extra} are not in the nominal M-site composition.")

    a_site = _sum(t, "A_site_metal")
    if a_site > 0 or opt.phase == "MAX":
        v["A_site_metal_norm"] = a_site * scale
        v["A_site_oxide_norm"] = _sum(t, "A_site_oxide") * scale
        if opt.a_target is not None:
            v["A_target"] = float(opt.a_target)

    dop = t[t["role"].isin(DOPANT_ROLES)]
    if len(dop):
        for (el, role), g in dop.groupby(["element", "role"]):
            v[f"{el}_{role}_norm"] = float(g["corrected_area"].sum()) * scale
        for el, g in dop.groupby("element"):
            v[f"{el}_dopant_total_norm"] = float(g["corrected_area"].sum()) * scale
        warnings.append(W_DOPANT)

    if delta > opt.metal_qc_tolerance:
        warnings.append(W_DEFICIENCY)
    elif delta < -opt.metal_qc_tolerance:
        warnings.append(W_EXCESS)
    if unres_low > 0 or unres_mid > 0:
        warnings.append(W_UNRESOLVED_O)
    if opt.oxide_f_correction_factor > 0:
        warnings.append(W_FCORR)
    if opt.include_f_residue and f_res > 0:
        warnings.append(W_FRESIDUE)
    if v["apparent_M_deficiency_status"] == "provisional":
        warnings.append(W_PROVISIONAL)
        for r in opt.provisional_reasons:
            warnings.append(f"delta not finalized: {r}")
    cand_roles = sorted({r for r in ("metal_suboxide_candidate", "unresolved_low_BE_F") if _sum(t, r) > 0})
    if cand_roles:
        warnings.append(W_CANDIDATE.format(roles=", ".join(cand_roles)))
    if opt.phase == "MAX":
        warnings.append(W_MAX)
    if opt.input_mode == "atomic_percent":
        warnings.append("atomic_percent input is treated as already RSF-corrected; confirm the provider's RSF set.")
    return QuantResult(sample=str(sample), values=v, table=t, warnings=warnings)


def quantify_table(df: pd.DataFrame, rsf: dict[str, float], opt: QuantOptions) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    pt = normalize_peak_table(df)
    rows, warns = [], {}
    full = pd.concat([pt.df, pt.skipped_empty], ignore_index=True)
    for sample, sub in full.groupby("sample", sort=False):
        r = quantify_sample(sub, rsf, opt)
        rows.append({"sample": sample, **r.values})
        warns[str(sample)] = r.warnings
    return pd.DataFrame(rows), warns


def role_reclassifications(df: pd.DataFrame, kb_roles: dict[tuple[str, str], str]) -> pd.DataFrame:
    """Audit trail: rows whose role differs from the knowledge-base default.

    ``kb_roles`` maps (core_level, component) -> default role.
    """
    out = []
    for _, r in df.iterrows():
        key = (str(r["core_level"]), str(r["component"]))
        if key in kb_roles and str(r["role"]) != kb_roles[key]:
            out.append({"sample": r["sample"], "core_level": key[0], "component": key[1],
                        "kb_default_role": kb_roles[key], "user_role": r["role"]})
    return pd.DataFrame(out, columns=["sample", "core_level", "component", "kb_default_role", "user_role"])


# --- Legacy Ti3C2Tx workbook model (§10) ------------------------------------------
LEGACY_COMPONENTS = {
    # key: (component name in workbook, core level)
    "C_TiTx": ("C-TiTx", "C 1s"),
    "Ti1": ("C-Ti+1-Tx", "Ti 2p"),
    "Ti2": ("C-Ti+2-Tx", "Ti 2p"),
    "Ti3": ("C-Ti+3-Tx", "Ti 2p"),
    "Ti_oxide": ("TiO2-xF2x (Ti)", "Ti 2p"),
    "O_low": ("C-Ti-O (low-BE O term)", "O 1s"),
    "O_oxide": ("TiO2-xF2x (O)", "O 1s"),
    "O_high": ("high-BE O termination", "O 1s"),
    "F0": ("F0 (C-Ti-F)", "F 1s"),
    "F1": ("F1", "F 1s"),
    "F2": ("F2", "F 1s"),
    "Cl": ("C-Ti-Cl", "Cl 2p"),
}
LEGACY_CONTROLS_DEFAULT = {"F_oxide_per_Ooxide": 0.6667, "Include_O_highBE": 1, "Include_Cl": 1,
                           "Include_F1": 0, "Include_F2": 0}
LEGACY_OUTPUT_COLUMNS = ["C_MX", "Ti_MX", "O_low", "O_high", "O_term_used", "O_oxide", "F0", "F1", "F2",
                         "F_term_raw", "F_oxide_est", "F_term_corr", "Cl_term", "Scale (2/C_MX)", "Ti_norm", "z_O",
                         "z_F_raw", "z_F_corr", "z_Cl", "z_total_raw", "z_total_corr", "O%_corr", "F%_corr",
                         "Cl%_corr", "QC_Ti", "QC_z", "z_O_low", "z_OH_like", "z_O_total", "O_low%_corr",
                         "OH_like%_corr"]
LEGACY_LABEL = ("Legacy/validated laboratory Ti3C2Tx model (historical workbook). Component names such as "
                "C-Ti+1/+2/+3-Tx are historical labels, not oxidation-state assignments; the primary interpretation "
                "model is Brette2024. The 0.6667 F-oxide factor is a legacy, model-dependent correction.")


def legacy_quantify(areas: dict[str, Optional[float]], rsf: dict[str, float],
                    controls: Optional[dict] = None) -> dict:
    """Exact port of the workbook 'Results' sheet formulas.

    Blank inputs behave as in Excel (blank cell / RSF = 0). Operation order follows the workbook.
    """
    c = {**LEGACY_CONTROLS_DEFAULT, **(controls or {})}
    a = {k: (0.0 if areas.get(k) is None or (isinstance(areas.get(k), float) and math.isnan(areas[k])) else float(areas[k]))
         for k in LEGACY_COMPONENTS}
    need = ["C 1s", "Ti 2p", "O 1s", "F 1s", "Cl 2p"]
    missing = [k for k in need if not (rsf.get(k) and rsf[k] > 0)]
    if missing:
        raise ValueError(f"Legacy model requires RSFs for {missing}.")
    rC, rTi, rO, rF, rCl = (float(rsf[k]) for k in need)

    C_MX = a["C_TiTx"] / rC
    Ti_MX = (a["Ti1"] + a["Ti2"] + a["Ti3"]) / rTi
    O_low = a["O_low"] / rO
    O_high = a["O_high"] / rO
    O_term_used = O_low + c["Include_O_highBE"] * O_high
    O_oxide = a["O_oxide"] / rO
    F0, F1, F2 = a["F0"] / rF, a["F1"] / rF, a["F2"] / rF
    F_term_raw = F0 + c["Include_F1"] * F1 + c["Include_F2"] * F2
    F_oxide_est = O_oxide * c["F_oxide_per_Ooxide"]
    F_term_corr = max(0, F_term_raw - F_oxide_est)
    Cl_term = c["Include_Cl"] * (a["Cl"] / rCl)
    if not C_MX > 0:
        raise ValueError("Legacy model: C-TiTx (framework C) must be positive.")
    scale = 2 / C_MX
    Ti_norm = Ti_MX * scale
    z_O = O_term_used * scale
    z_F_raw = F_term_raw * scale
    z_F_corr = F_term_corr * scale
    z_Cl = Cl_term * scale
    z_total_raw = z_O + z_F_raw + z_Cl
    z_total_corr = z_O + z_F_corr + z_Cl
    pos = z_total_corr > 0
    z_O_low = O_low * scale
    z_OH_like = O_high * c["Include_O_highBE"] * scale
    return {
        "C_MX": C_MX, "Ti_MX": Ti_MX, "O_low": O_low, "O_high": O_high, "O_term_used": O_term_used,
        "O_oxide": O_oxide, "F0": F0, "F1": F1, "F2": F2, "F_term_raw": F_term_raw, "F_oxide_est": F_oxide_est,
        "F_term_corr": F_term_corr, "Cl_term": Cl_term, "Scale (2/C_MX)": scale, "Ti_norm": Ti_norm, "z_O": z_O,
        "z_F_raw": z_F_raw, "z_F_corr": z_F_corr, "z_Cl": z_Cl, "z_total_raw": z_total_raw,
        "z_total_corr": z_total_corr,
        "O%_corr": z_O / z_total_corr if pos else np.nan,
        "F%_corr": z_F_corr / z_total_corr if pos else np.nan,
        "Cl%_corr": z_Cl / z_total_corr if pos else np.nan,
        "QC_Ti": "OK" if abs(Ti_norm - 3) <= 0.5 else "CHECK",
        "QC_z": "OK" if 1 <= z_total_corr <= 3 else "CHECK",
        "z_O_low": z_O_low, "z_OH_like": z_OH_like, "z_O_total": z_O_low + z_OH_like,
        "O_low%_corr": z_O_low / z_total_corr if pos else np.nan,
        "OH_like%_corr": z_OH_like / z_total_corr if pos else np.nan,
    }


def legacy_template(sample: str) -> pd.DataFrame:
    return pd.DataFrame([{"sample": sample, "core_level": core, "component": name, "area": np.nan}
                         for name, core in LEGACY_COMPONENTS.values()])


def legacy_quantify_table(df: pd.DataFrame, rsf: dict[str, float], controls: Optional[dict] = None
                          ) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    for col in ("sample", "component", "area"):
        if col not in df.columns:
            raise ValueError(f"Legacy table needs column {col!r}")
    name_to_key = {name: k for k, (name, _) in LEGACY_COMPONENTS.items()}
    rows, warns = [], {}
    for sample, sub in df.groupby("sample", sort=False):
        w: list[str] = [LEGACY_LABEL]
        unknown = sorted(set(sub["component"].astype(str)) - set(name_to_key))
        if unknown:
            w.append(f"Components not used by the legacy model (ignored): {unknown}")
        areas: dict[str, Optional[float]] = {}
        for _, r in sub.iterrows():
            k = name_to_key.get(str(r["component"]))
            if k is None:
                continue
            val = pd.to_numeric(pd.Series([r["area"]]), errors="coerce").iloc[0]
            if k in areas and areas[k] is not None and not pd.isna(val):
                raise ValueError(f"{sample}: duplicate legacy component {r['component']!r}")
            areas[k] = None if pd.isna(val) else float(val)
        blank = [LEGACY_COMPONENTS[k][0] for k in LEGACY_COMPONENTS if areas.get(k) is None]
        if blank:
            w.append(f"Blank legacy inputs treated as 0 (workbook behaviour): {blank}")
        rows.append({"sample": sample, **legacy_quantify(areas, rsf, controls)})
        warns[str(sample)] = w
    return pd.DataFrame(rows), warns
