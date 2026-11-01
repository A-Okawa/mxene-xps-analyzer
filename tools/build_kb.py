"""Build data/knowledge_base.json from the handoff §27 seed JSON.

Numerical values (center, bounds, fwhm, fwhm_bounds, split, ratio, role, source)
are copied verbatim from the seed.  This script only ADDS metadata required by
the handoff (§11 interpretation caveat for every component, §25 source keys,
§32.2 versioning).  Added interpretation text is restricted to wording found in
the handoff (§12-§16); it does not create new chemical assignments.

Usage:  python -I tools/build_kb.py
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "tools" / "kb_seed_from_handoff_s27.json"
OUT = ROOT / "data" / "knowledge_base.json"

TRANSFER_TI = (" Transferred from the Ti3C2Tz reference (Brette2024) as a composite starting model;"
               " not independently validated for (Nb,Ti)2CTx.")
TRANSFER_NB = (" Transferred from the Nb2CTx/Nb1.33CTx reference (Halim2018Vacancy) as a composite"
               " starting model; not independently validated for (Nb,Ti)2CTx.")

ADDED_INTERPRETATION = {
    ("ti3c2_brette", "C 1s", "C_contam_1"): "Contamination carbon (~284.x eV in Brette SI, HF10%). Excluded from framework C.",
    ("ti3c2_brette", "C 1s", "C_contam_2"): "Contamination / oxygenated carbon (~286.x eV in Brette SI). Excluded from framework C.",
    ("ti3c2_brette", "C 1s", "C_contam_3"): "Higher oxygenated / carboxyl contamination (~288.x eV); include only as required. Excluded from framework C.",
    ("nb2c_halim", "C 1s", "C_MXene"): "Literature label C-Nb-Tx (282.5 eV in Halim2018Vacancy). Framework carbide C.",
    ("nb2c_halim", "C 1s", "C_C"): "Literature label C-C. Contamination; excluded from framework C.",
    ("nb2c_halim", "C 1s", "CHx"): "Literature label CHx. Contamination; excluded from framework C.",
    ("nb2c_halim", "C 1s", "C_O"): "Literature label C-O. Contamination; excluded from framework C.",
    ("nb2c_halim", "C 1s", "COO"): "Literature label COO. Contamination; excluded from framework C.",
    ("nb2c_halim", "O 1s", "Nb2O5_O"): "Literature label Nb2O5 (oxide O). Overlaps the C-Nb-Ox window; cross-check against the Nb 3d oxide component.",
    ("nb2c_halim", "O 1s", "Nb_O_term"): "Literature label C-Nb-Ox. O-termination candidate.",
    ("nb2c_halim", "O 1s", "Nb_OH_term"): "Literature label C-Nb-(OH)x. OH-like termination candidate.",
    ("nb2c_halim", "O 1s", "H2O_ads"): "Literature label H2Oads. Adsorbate; excluded from Tx by default.",
    ("nb2c_halim", "F 1s", "Nb_F_term"): "Literature label C-Nb-Fx. F-termination candidate.",
    ("nb2c_halim", "F 1s", "NbF_residue"): "Literature label NbF5 (Nb2CTx) or NbF3 and/or ScF3 (Nb1.33CTx). Fluoride residue candidate; not counted as Tx by default.",
    ("nbti_m2c_defect_aware", "C 1s", "C_C"): "Contamination carbon. Excluded from framework C.",
    ("nbti_m2c_defect_aware", "C 1s", "C_O"): "Contamination / oxygenated carbon. Excluded from framework C.",
    ("nbti_m2c_defect_aware", "C 1s", "COO"): "Higher oxygenated / carboxyl contamination. Excluded from framework C.",
    ("nbti_m2c_defect_aware", "O 1s", "OH_like"): "OH-like termination candidate (C-Nb-(OH)x in Halim2018Vacancy; OH-like in Brette2024). Not attributed to Ti vs Nb.",
    ("nbti_m2c_defect_aware", "O 1s", "O_adsorbate"): "Water/organic surface species candidate. Excluded from Tx by default.",
    ("nbti_m2c_defect_aware", "F 1s", "F_MXene"): "MXene-F termination candidate (C-Nb-Fx in Halim2018Vacancy). Not attributed to Ti- vs Nb-bound F; one F 1s peak can represent several local F environments (Brette2024).",
}

# Components whose role is deliberately left unresolved (handoff §2.6, §16.4).
UNRESOLVED = {("nbti_m2c_defect_aware", "O 1s", "O_low"), ("nbti_m2c_defect_aware", "O 1s", "O_mid")}

EXTRA_REFERENCES = {
    "Chen2024": {
        "title": "Unraveling the Role of Metal Vacancy Sites and Doped Nitrogen in Enhancing Pseudocapacitance Performance of Defective MXene",
        "doi": "10.1002/smll.202307408",
        "journal": "Small",
        "year": 2024,
        "copyright_or_online_year": 2023,
        "use": "Conceptual separation of defect and dopant variables (pristine / defect-only / dopant-only / defect+dopant).",
    },
    "LegacyExcel": {
        "title": "Laboratory workbook mxene_surface_functional_group_quant_template_v2 (private, not distributed)",
        "doi": None,
        "use": "Historical Ti3C2Tx operational quantification logic; reproduced exactly by the legacy model.",
    },
}

# Inconsistencies found inside the handoff document itself.  They are recorded,
# not resolved: §27 JSON is used as the machine-readable source; the user decides.
DECIDED = "2026-10-16"
SPEC_DISCREPANCIES = [
    {
        "id": "D1",
        "topic": "Ti3C2Tx (ti3c2_brette) Ti 2p starting model",
        "detail": "Handoff §12 table (bounds 454.2-455.5, FWHM 1.1 (0.6-2.0), Ti_midBE 455.7, oxide 459.2/1.9; identical to the (Nb,Ti)2CTx profile) vs §27 JSON ti3c2_brette (454.3-455.4, FWHM 1.0 (0.6-1.8), Ti_midBE 455.6, Ti_oxide_oxyfluoride 459.3/1.8).",
        "status": "resolved",
        "decision": "Adopt §27 JSON for Ti3C2Tx. Brette SI experimental centers/FWHMs are kept as literature_reference values; fitting windows (bounds, fwhm_bounds) are editable implementation bounds.",
        "decided_by": "user",
        "decided_at": DECIDED,
    },
    {
        "id": "D2",
        "topic": "Nb2CTx (nb2c_halim) Nb 3d starting model",
        "detail": "Handoff §13 table vs §27 JSON nb2c_halim differed.",
        "status": "resolved",
        "decision": "Reference centers/FWHMs set to the Halim SI Nb2CTx values: Nb 3d5/2 = 204.8/207.3/209.2 eV, FWHM = 1.1/1.7/1.0 eV (C-Nb-Tx/Nb2O5/NbF5). Initial/reference values, not hard-fixed positions; §27 fit windows unchanged.",
        "decided_by": "user",
        "decided_at": DECIDED,
    },
    {
        "id": "D3",
        "topic": "Reference source keys",
        "detail": "§25 recommends Halim2016XPS/Chen2023; §27 JSON uses Halim2016; Chen et al. cited as Small 2024 (§15.4) but filename says 2023 (§26).",
        "status": "resolved",
        "decision": "Keep reference ID Halim2016. Standardize Chen to Chen2024 (year 2024); 2023 retained as copyright_or_online_year metadata.",
        "decided_by": "user",
        "decided_at": DECIDED,
    },
]

# Brette SI experimental values, stored ONLY as literature_reference entries (display/reference metadata).
# They are never used as starting values, constraints or fitting windows.  Values transcribed from
# downloadSupplement.pdf (Brette2024 SI) and checked against the rendered tables.
# D1 (2026-10-16): Ti 2p (Table S3-A) and C 1s (Table S4-A), HF10 % before sputtering.
# 2026-10-16 correction: O 1s / F 1s for HF10 % (Tables S5-B, S8, S4-B) and HF48 % (Table S5-A).
# Mapping to KB components follows the paper's own P1-P4 labelling of O 1s (main text): P1 ~530 eV O in
# ordered Ti3C2T2 -> O_ideal; P2 ~531 eV O on disordered surfaces (mix and/or vacancies on T-groups) ->
# O_disordered; P3 ~532 eV Ti3C2(OH)2-like -> OH_like; P4 adsorbed H2O and/or organic (OR) -> O_adsorbate.
REF_NOTE = "Literature reference value only; not a starting value, constraint or fixed position."
S3A = {"table": "S3-A", "dataset": "HF10", "sample_state": "before Ar sputtering"}
S4A = {"table": "S4-A", "dataset": "HF10", "sample_state": "before Ar sputtering"}
S4B = {"table": "S4-B", "dataset": "HF10", "sample_state": "before Ar sputtering"}
S8 = {"table": "S8", "dataset": "HF10", "sample_state": "before Ar sputtering"}
S5A = {"table": "S5-A", "dataset": "HF48", "sample_state": "no sputtering stated (Part S5)"}
S5B = {"table": "S5-B", "dataset": "HF10", "sample_state": "no sputtering stated (Part S5)"}


def _r(meta, line, center, fwhm, si_label, secondary=None):
    e = {"source": "Brette2024", **meta, "line": line, "center": center, "fwhm": fwhm, "si_label": si_label}
    if secondary:
        e["secondary"] = secondary
    e["note"] = REF_NOTE
    return e


BRETTE_SI = {
    ("Ti 2p", "Ti_lowBE"): [_r(S3A, "Ti 2p3/2", 454.9, 1.0, "Ti1", {"line": "Ti 2p1/2", "center": 460.8, "fwhm": 1.1})],
    ("Ti 2p", "Ti_midBE"): [_r(S3A, "Ti 2p3/2", 455.6, 1.6, "Ti2", {"line": "Ti 2p1/2", "center": 461.5, "fwhm": 1.6})],
    ("Ti 2p", "Ti_highBE"): [_r(S3A, "Ti 2p3/2", 456.8, 1.7, "Ti3", {"line": "Ti 2p1/2", "center": 462.4, "fwhm": 1.6})],
    ("Ti 2p", "Ti_oxide_oxyfluoride"): [_r(S3A, "Ti 2p3/2", 459.3, 1.8, "TiO2-xF2x", {"line": "Ti 2p1/2", "center": 464.8, "fwhm": 3.2})],
    ("C 1s", "C_MXene"): [_r(S4A, "C 1s", 281.8, 0.6, "MXene")],
    ("C 1s", "C_contam_1"): [_r(S4A, "C 1s", 284.1, 1.9, "contamination")],
    ("C 1s", "C_contam_2"): [_r(S4A, "C 1s", 286.4, 1.8, "contamination")],
    ("C 1s", "C_contam_3"): [_r(S4A, "C 1s", 288.9, 2.0, "contamination")],
    ("O 1s", "O_ideal"): [_r(S5B, "O 1s", 529.7, 1.0, "MXene (P1)"), _r(S8, "O 1s", 529.6, 0.9, "MXene (P1)"),
                          _r(S5A, "O 1s", 529.8, 1.0, "MXene (P1)")],
    ("O 1s", "O_disordered"): [_r(S5B, "O 1s", 531.0, 2.0, "MXene (P2)"), _r(S8, "O 1s", 531.1, 2.0, "MXene (P2)"),
                               _r(S5A, "O 1s", 531.0, 2.0, "MXene (P2)")],
    ("O 1s", "OH_like"): [_r(S5B, "O 1s", 532.5, 1.6, "MXene (P3)"), _r(S8, "O 1s", 532.5, 1.5, "MXene (P3)"),
                          _r(S5A, "O 1s", 532.0, 1.7, "MXene (P3)")],
    ("O 1s", "O_adsorbate"): [_r(S5B, "O 1s", 534.0, 1.6, "H2O + OR (P4)"), _r(S8, "O 1s", 534.0, 2.1, "H2O + OR (P4)"),
                              _r(S5A, "O 1s", 533.8, 1.7, "OR (P4)"), _r(S5A, "O 1s", 535.8, 2.0, "H2O (P4)")],
    ("F 1s", "F_MXene"): [_r(S5B, "F 1s", 684.9, 1.3, "MXene"), _r(S4B, "F 1s", 684.8, 1.3, "MXene"),
                          _r(S5A, "F 1s", 685.1, 1.6, "MXene")],
}

# Brette's interpretation of the ~531 eV O 1s contribution (main text and summary table), kept verbatim in spirit.
BRETTE_O_DISORDERED = ("Brette2024: the ~531 eV contribution (P2) is O in Ti3C2Tz with disordered surfaces (mix and/or "
                       "vacancies on T-groups); its magnitude probes surface disorder, especially vacancy-induced "
                       "(it grows at the expense of P1 when F desorbs on Ar sputtering). It is NOT a uniquely assigned "
                       "OH state; the OH-like contribution is the separate ~532 eV P3 (Ti3C2(OH)2 model).")
BRETTE_SI_F_NOTE = "Brette2024 identified a single F 1s contribution in HF10 and HF48 spectra; no F residue value is reported there."

# D2: Halim SI Nb2CTx values (user-supplied) -> starting center/FWHM and literature reference.
HALIM_NB2C_NOTE = "Halim2018Vacancy SI, Nb2CTx (C-Nb-Tx / Nb2O5 / NbF5); initial/reference value, not a hard-fixed position."
HALIM_NB2C = {
    ("Nb 3d", "Nb_MXene"): (204.8, 1.1),
    ("Nb 3d", "Nb2O5"): (207.3, 1.7),
    ("Nb 3d", "NbF5_like"): (209.2, 1.0),
}



# ======================================================================================================================
# Calibration phase 2 (KB 1.3.0, 2026-11-02, user instructions on branch calibration/real-xps).
# Data-derived CANDIDATES from the first real provider fits (PHI MultiPak, Sasaki 2026-10 set; private).  They are
# fitted but not validated: quantification keeps them out of the framework by default and flags dependent outputs
# provisional.  Using them to "validate" the same provider fits is circular.
# ======================================================================================================================
CAL_DATE = "2026-11-02"
CAL_SOURCE = "Sasaki2026MultiPak"
CAL_REFERENCE = {
    CAL_SOURCE: {
        "title": "Provider PHI MultiPak fits of the Sasaki 2026-10 XPS set (Nb2AlC, (Nb,Ti)2AlC, Nb2CTx, (Nb0.6Ti0.4)2CTx); private, not distributed",
        "doi": None,
        "use": "Data-derived calibration candidates and MAX-phase starting values. Not literature; circular if used to validate the same provider fits.",
    }
}
CAL_PROV = (f"Calibration candidate (KB 1.3.0, {CAL_DATE}): starting value/window derived from provider fits ({CAL_SOURCE}); "
            "not literature, not validated.")
NB_SPLIT, NB_RATIO = 2.8, 1.5          # same as the existing Nb 3d components (provider reports 2.72; open item D4)
TI_SPLIT, TI_SPLIT_OX, TI_RATIO = 5.8, 5.7, 2.0

NB_SUBOXIDE = {
    "name": "Nb_suboxide_like", "kind": "doublet", "center": 205.75, "bounds": [205.4, 206.2], "fwhm": 1.5,
    "fwhm_bounds": [0.8, 2.3], "split": NB_SPLIT, "ratio_primary_secondary": NB_RATIO,
    "role": "metal_suboxide_candidate", "source": CAL_SOURCE, "assignment_status": "calibration_candidate",
    "calibration_candidate": True, "validated": False, "confidence": "medium",
    "interpretation": ("Tentative Nb_suboxide_like / Nb(III)-O-like doublet, Nb 3d5/2 ~205.6-205.9 eV (3d3/2 ~208.3-208.6 eV), "
                       "16-25 % of Nb 3d in the provider fits of Nb2CTx and (Nb0.6Ti0.4)2CTx. Non-framework by default. "
                       "Must not be folded into Nb2O5 or framework Nb. Requires cross-check against O 1s."),
    "alternative_interpretation": "Termination-modified framework Nb (framework Nb whose 3d BE is shifted by O-rich/disordered terminations).",
    "alternative_role": "framework_metal", "cross_check": ["O 1s"],
}
F_LOW = {
    "name": "F_lowBE_unresolved", "kind": "single", "center": 682.7, "bounds": [682.0, 683.6], "fwhm": 1.5,
    "fwhm_bounds": [0.7, 2.4], "role": "unresolved_low_BE_F", "source": CAL_SOURCE, "assignment_status": "unresolved",
    "calibration_candidate": True, "validated": False, "confidence": "low",
    "interpretation": ("Unresolved low-BE F 1s component (~682.7 eV, ~10 % of F 1s in the provider fit of Nb2CTx). No chemistry "
                       "assigned; not counted as Tx."),
}
TI_OXIDE_LOW_BOUND = 458.0   # was 458.2; provider TiO2-like 458.17 eV must not be rejected for 0.03 eV
NB3S_WARNING = ("Nb 3s (roughly 466-469 eV) overlaps the Ti 2p1/2 side and the high-BE background of Ti 2p in Nb-containing "
                "samples (also noted by the provider). Do not over-trust Ti 2p1/2 positions/areas; the fixed 2:1 ratio "
                "carries Ti 2p1/2. Compare on Ti 2p3/2 and consider restricting the fit window.")


def _c(name, kind, center, bounds, fwhm, fwhm_bounds, role, interp, status="calibration_candidate", conf="medium",
       split=None, ratio=None, **extra):
    d = {"name": name, "kind": kind, "center": center, "bounds": list(bounds), "fwhm": fwhm, "fwhm_bounds": list(fwhm_bounds),
         "role": role, "source": CAL_SOURCE, "interpretation": interp, "assignment_status": status,
         "calibration_candidate": True, "validated": False, "confidence": conf, "provenance": CAL_PROV}
    if kind == "doublet":
        d["split"], d["ratio_primary_secondary"] = split, ratio
    d.update(extra)
    return d


def _max_profile(material: str, with_ti: bool) -> dict:
    nb = [
        _c("Nb_MAX", "doublet", 202.85, (202.3, 203.4), 1.2, (0.7, 2.0), "framework_metal",
           "Main Nb 3d5/2 line of the MAX phase (provider 202.73-202.95 eV). Metallic MAX lines are asymmetric; a "
           "symmetric pseudo-Voigt needs Nb_MAX_2nd to describe the tail.", split=NB_SPLIT, ratio=NB_RATIO),
        _c("Nb_MAX_2nd", "doublet", 203.9, (203.4, 204.6), 1.8, (0.8, 2.6), "framework_metal",
           "Unresolved: asymmetric tail of the metallic MAX line (symmetric line shape) or a second Nb environment "
           "(provider 203.83-204.42 eV, FWHM 1.7-2.0).", status="unresolved", conf="low", split=NB_SPLIT, ratio=NB_RATIO,
           alternative_interpretation="Near-surface Nb-C-O / sub-oxide environment (non-framework).",
           alternative_role="metal_suboxide_candidate"),
        {**NB_SUBOXIDE, "bounds": [205.4, 206.4], "center": 205.9, "provenance": CAL_PROV,
         "interpretation": NB_SUBOXIDE["interpretation"].replace(
             "of Nb2CTx and (Nb0.6Ti0.4)2CTx", "of Nb2CTx and (Nb0.6Ti0.4)2CTx; MAX provider fits 205.83-206.21 eV")},
        _c("Nb2O5", "doublet", 207.2, (206.7, 207.8), 1.45, (0.8, 2.3), "metal_oxide",
           "Nb2O5-like surface oxide (provider 207.07-207.35 eV).", split=NB_SPLIT, ratio=NB_RATIO),
    ]
    models = {}
    if with_ti:
        models["Ti 2p"] = [
            _c("Ti_MAX", "doublet", 454.5, (453.9, 455.2), 2.0, (0.8, 2.8), "framework_metal",
               "Ti in the MAX lattice (provider 454.44-454.51 eV).", split=TI_SPLIT, ratio=TI_RATIO),
            _c("Ti_oxide", "doublet", 458.4, (TI_OXIDE_LOW_BOUND, 459.5), 1.9, (1.0, 3.0), "metal_oxide",
               "TiO2-like surface oxide (provider 458.36-458.46 eV).", split=TI_SPLIT_OX, ratio=TI_RATIO),
        ]
    models["Nb 3d"] = nb
    models["C 1s"] = [
        _c("C_MAX", "single", 282.4, (281.5, 283.2), 1.3, (0.6, 2.3), "framework_C",
           "Carbide C of the MAX phase (provider 281.90-282.86 eV; wide spread between samples)."),
        _c("C_C", "single", 284.8, (284.2, 285.4), 1.9, (0.8, 2.8), "contaminant", "Adventitious C (charge reference 284.8 eV). Excluded from framework C.", conf="high"),
        _c("C_O", "single", 286.6, (285.8, 287.4), 1.8, (0.8, 2.8), "contaminant", "Oxygenated contamination C. Excluded from framework C.", conf="high"),
        _c("COO", "single", 288.7, (287.8, 289.7), 1.7, (0.8, 2.8), "contaminant", "Carboxyl-type contamination C. Excluded from framework C.", conf="high"),
    ]
    models["O 1s"] = [
        _c("O_lowBE", "single", 529.4, (528.8, 529.9), 1.9, (0.8, 2.5), "O_unresolved_low",
           "Unresolved low-BE O (metal oxide lattice O candidate; provider 529.30-529.39 eV).", status="unresolved", conf="low"),
        _c("O_midBE", "single", 530.3, (529.9, 530.9), 1.6, (0.8, 2.5), "O_unresolved_mid",
           "Unresolved mid-BE O (Nb/Ti oxide O candidate; provider 530.05-530.43 eV).", status="unresolved", conf="low"),
        _c("O_highBE", "single", 531.8, (531.2, 532.6), 2.4, (1.0, 3.0), "other",
           "Unresolved high-BE O (Al2O3/hydroxide/adsorbate candidates; provider 531.66-531.94 eV).", status="unresolved", conf="low"),
    ]
    models["Al 2p"] = [
        _c("Al_MAX", "single", 72.2, (71.6, 72.8), 1.4, (0.7, 2.2), "A_site_metal",
           "Al in the MAX lattice (provider 71.94-72.50 eV). Spin-orbit not resolved; single band."),
        _c("Al_oxide", "single", 74.1, (73.5, 74.8), 1.6, (0.8, 2.4), "A_site_oxide",
           "Al2O3-like surface oxide (provider 73.94-74.32 eV)."),
        _c("Al_highBE", "single", 75.6, (74.9, 76.5), 1.5, (0.8, 2.4), "other",
           "Unresolved high-BE Al (provider 75.26-76.07 eV; e.g. Al-F or hydroxide candidates).", status="unresolved", conf="low"),
    ]
    notes = [
        "MAX-phase profile (KB 1.3.0): ALL starting values are data-derived from the provider fits of the same samples "
        f"({CAL_SOURCE}), not literature. Literature values must be added before these become reference models.",
        "Termination (Tx) outputs are not meaningful for a MAX phase; O/F/Cl are surface oxide/contamination.",
        "Sample-1 provider fit also contains a small Nb 3d doublet at 201.80 eV (3.6 %) that is not modelled here.",
    ]
    prof = {"material": material, "phase": "MAX", "notes": notes, "fit_models": models}
    if with_ti:
        prof["core_level_warnings"] = {"Ti 2p": [NB3S_WARNING]}
    return prof


CAL_STRUCTURES = {
    "Nb2AlC": {"m_target": 2.0, "x_target": 1.0, "a_target": 1.0, "default_profile": "nb2alc_max", "phase": "MAX"},
    "(Nb,Ti)2AlC": {"m_target": 2.0, "x_target": 1.0, "a_target": 1.0, "default_profile": "nbti_211_max", "phase": "MAX"},
}
CAL_D4 = {
    "id": "D4", "topic": "Spin-orbit splitting: KB vs provider",
    "detail": "KB doublets use Nb 3d dSO = 2.8 eV and Ti 2p dSO = 5.6-5.8 eV; the provider (PHI MultiPak) fits use 2.72 eV and 6.17 eV.",
    "status": "open", "decision": None, "decided_by": None, "decided_at": None,
}


def apply_calibration_phase2(kb: dict) -> None:
    kb["kb_version"] = "1.3.0"
    kb["changelog"].append({"version": "1.3.0", "date": CAL_DATE, "note": (
        "Calibration phase 2 (user instructions): Nb_suboxide_like doublet candidate (205.6-205.9 eV, role "
        "metal_suboxide_candidate, non-framework by default, alternative = termination-modified framework Nb, confidence "
        "medium, O 1s cross-check) in nb2c_halim and nbti_m2c_defect_aware; nbti Ti_oxide lower bound 458.2 -> 458.0 eV with "
        "Nb 3s overlap warning; F_lowBE_unresolved (~682.7 eV, role unresolved_low_BE_F, not Tx); separate MAX profiles "
        "nb2alc_max and nbti_211_max (data-derived, not literature); apparent M deficiency provisional until the Nb "
        "candidate is validated; D4 (spin-orbit splitting) recorded as open.")})
    kb["references"].update(CAL_REFERENCE)
    kb["spec_discrepancies"].append(CAL_D4)
    kb["structures"].update(CAL_STRUCTURES)
    for k in ("Ti3C2Tx", "Nb2CTx", "(Nb,Ti)2CTx"):
        kb["structures"][k].setdefault("phase", "MXene")
    kb["calibration_status"] = {
        "phase": 2, "date": CAL_DATE,
        "nb3d_suboxide_like_validated": False,
        "rule": ("Calibration candidates are fitted but not validated. Quantification keeps them out of the framework by "
                 "default, reports the alternative interpretation separately, and flags apparent M-site deficiency as "
                 "provisional until the Nb_suboxide_like component is validated (O 1s cross-check + orthogonal evidence)."),
    }
    for pkey in ("nb2c_halim", "nbti_m2c_defect_aware"):
        prof = kb["profiles"][pkey]
        prof["phase"] = "MXene"
        nb = prof["fit_models"]["Nb 3d"]
        i = next(j for j, c in enumerate(nb) if c["name"] == "Nb2O5")
        nb.insert(i, {**NB_SUBOXIDE, "provenance": CAL_PROV})
        prof["fit_models"]["F 1s"].insert(0, {**F_LOW, "provenance": CAL_PROV})
    ti = kb["profiles"]["nbti_m2c_defect_aware"]["fit_models"]["Ti 2p"]
    tox = next(c for c in ti if c["name"] == "Ti_oxide")
    tox["bounds"] = [TI_OXIDE_LOW_BOUND, tox["bounds"][1]]
    tox["provenance"] += f" Lower BE bound widened 458.2 -> {TI_OXIDE_LOW_BOUND} eV in KB 1.3.0 ({CAL_DATE}) so that the provider TiO2-like 458.17 eV is not rejected solely by a 0.03 eV window edge."
    kb["profiles"]["nbti_m2c_defect_aware"]["core_level_warnings"] = {"Ti 2p": [NB3S_WARNING]}
    kb["profiles"]["ti3c2_brette"]["phase"] = "MXene"
    kb["profiles"]["nb2alc_max"] = _max_profile("Nb2AlC", with_ti=False)
    kb["profiles"]["nbti_211_max"] = _max_profile("(Nb,Ti)2AlC", with_ti=True)


# ======================================================================================================================
# KB 1.3.1 (2026-11-02): metadata only. D4 sensitivity (tools/calibration/d4_nb3d_split_sensitivity.py) and the O 1s
# cross-check of the Nb candidate (tools/calibration/o1s_crosscheck_nb_candidate.py). No fit value is changed.
# ======================================================================================================================
D4_RESULT = {
    "status": "calibration_sensitive (Nb 3d); decided (Ti 2p)",
    "decision": ("Nb 3d: KB dSO 2.8 eV kept (not changed). Fixed-grid refit (2.70/2.72/2.75/2.80 eV) of the Nb2CTx and "
                 "(Nb0.6Ti0.4)2CTx spectra: sample (Nb0.6Ti0.4)2CTx stable, Nb2CTx exceeds the declared materiality "
                 "criteria (component fraction span >= 0.05 or delta span >= 0.05) in 3 of 4 window/background settings "
                 "(delta span up to 0.095; at 2.80 eV Nb_suboxide_like sits on its lower window bound 205.40 eV). AIC "
                 "favours 2.75-2.80 eV in the provider window and is nearly flat on the full raw region. D4 is therefore "
                 "calibration-sensitive and delta is not finalized. Ti 2p: the provider's 6.17 eV is NOT adopted (Nb 3s "
                 "overlap compromises the 2p1/2 side); literature-based Ti splitting kept; Ti 2p3/2 prioritized for "
                 "interpretation."),
    "decided_by": "user (Ti 2p decision, Nb grid design and criteria); sensitivity analysis by the app maintainer",
    "decided_at": "2026-11-02",
}
O1S_CROSSCHECK = {
    "date": "2026-11-02", "candidate": "Nb_suboxide_like", "result": "inconclusive",
    "method": ("Stoichiometric O budget in the low-BE O 1s window (<= 530.8 eV): R = O_low - 2.5 Nb(Nb2O5) - 2 Ti(TiO2), "
               "compared with x * Nb(candidate) for x = 1 (NbO), 1.5 (Nb(III)), 2 (NbO2), over provider and Python "
               "decompositions of Nb 3d (incl. the D4 grid) and O 1s (Shirley/linear). Supports if min rho_1 >= 1; "
               "against if max rho_1 < 0.5."),
    "summary": ("Nb2CTx: rho_1 0.18-1.06 (median 0.78; provider decompositions 0.55-0.74, rho_2 0.28-0.37): the oxide-O "
                "budget does not cover a stoichiometric NbO2/Nb(III) oxide but cannot exclude a partial sub-oxide. "
                "(Nb0.6Ti0.4)2CTx: rho_1 -0.82-1.82 because the O 1s decomposition itself is unstable (low-BE O fraction "
                "0.38-0.65 between backgrounds); no discriminating power."),
    "biases": ["provider RSF/transmission set", "overlayer attenuation (O 1s vs Nb 3d kinetic energy)",
               "ordered O terminations near 530 eV inflate R", "possible O deficiency of surface Nb2O5"],
}


def apply_calibration_d4_o1s(kb: dict) -> None:
    kb["kb_version"] = "1.3.1"
    kb["changelog"].append({"version": "1.3.1", "date": "2026-11-02", "note": (
        "Metadata only: D4 marked calibration-sensitive for Nb 3d (KB 2.8 eV kept; fixed-grid sensitivity) and decided "
        "for Ti 2p (literature splitting kept, provider 6.17 eV not adopted, Ti 2p3/2 prioritized); O 1s cross-check of "
        "Nb_suboxide_like recorded as inconclusive; candidate remains unvalidated; delta not finalized. MAX profiles unchanged.")})
    d4 = next(d for d in kb["spec_discrepancies"] if d["id"] == "D4")
    d4.update(D4_RESULT)
    cs = kb["calibration_status"]
    cs["d4_nb3d_dSO_calibration_sensitive"] = True
    cs["o1s_crosscheck_nb_suboxide_like"] = O1S_CROSSCHECK
    cs["delta_finalization_blocked_by"] = [
        "Nb_suboxide_like (Nb 3d ~205.6-205.9 eV) not validated; O 1s cross-check inconclusive",
        "D4: Nb 3d spin-orbit splitting is calibration-sensitive for Nb2CTx",
    ]
    cs["ti2p_interpretation"] = "Prioritize Ti 2p3/2; do not over-trust Ti 2p1/2 (Nb 3s overlap)."


# ======================================================================================================================
# KB 1.4.0 (2026-11-02): MAX-phase literature calibration (user instructions). The provider-derived MAX starting
# values of 1.3.0 are replaced. Three reference layers are kept apart on every MAX component:
#   literature_reference      values reported for the same material (citation, or citation explicitly pending)
#   transferred_reference     values reported for related materials (not independently validated for the target)
#   calibration_observations  provider (PHI MultiPak) fit values of the Sasaki set: observations only
# Nb_suboxide_like is NOT modified (unvalidated, non-framework, delta provisional, D4 calibration-sensitive).
# Spin-orbit splittings are not changed (D4); literature splittings are recorded in the reference entries only.
# ======================================================================================================================
MAX_DATE = "2026-11-02"
MAX_REFERENCES = {
    "Naslund2020": {
        "title": "X-ray Photoelectron Spectroscopy of Ti3AlC2, Ti3C2Tz, and TiC Provides Evidence for the Electrostatic Interaction between Laminated Layers in MAX-Phase Materials",
        "authors": "L.-A. Naslund, P. O. A. Persson, J. Rosen", "journal": "J. Phys. Chem. C", "volume": 124,
        "pages": "27732-27742", "year": 2020, "doi": "10.1021/acs.jpcc.0c07413",
        "be_reference": "Fermi edge of each sample (not adventitious C)",
        "use": "Ti3AlC2/TiC Ti 2p, C 1s, Al 2p, O 1s fit tables; Al 2p core-level shifts of Ti2AlC, V2AlC, Nb2AlC vs Al metal.",
    },
    "Greczynski2018Carbides": {
        "title": "Reference binding energies of transition metal carbides by core-level x-ray photoelectron spectroscopy free from Ar+ etching artefacts",
        "authors": "G. Greczynski, D. Primetzhofer, L. Hultman", "journal": "Appl. Surf. Sci.", "volume": 436,
        "pages": "102-110", "year": 2018, "doi": "10.1016/j.apsusc.2017.11.264",
        "be_reference": "ISO-calibrated spectrometer (Au 4f7/2 83.96, Ag 3d5/2 368.21, Cu 2p3/2 932.62 eV)",
        "use": "Binary carbide reference values (TiC, NbC) and their native-oxide lines; supporting/transferred only.",
    },
    "UserNb2AlC2026": {
        "title": "Nb2AlC XPS reference values supplied by the user on 2026-11-02 (Nb-C 3d5/2 202.9 / 3d3/2 205.5 eV; Nb oxide 206.8 / 209.6 eV; C-Nb 282.5 eV)",
        "doi": None, "citation_status": "pending",
        "use": "Nb2AlC starting values. The original publication has not been identified yet; the citation must be added before these are treated as fully cited literature.",
    },
}
REF_FRAME_NOTE = ("Literature BE scales differ (Fermi edge / ISO calibration) from the provider data (adventitious C 1s = 284.8 eV); "
                  "offsets of a few 0.1 eV are expected. Values are starting points; fit windows are editable implementation bounds.")


def _lit(source, line, center, material, **kw):
    return {"reference_type": "literature_reference", "source": source, "material": material, "line": line,
            "center": center, "note": "Literature reference value only; not a constraint or fixed position.", **kw}


def _tr(source, line, center, material, confidence, **kw):
    return {"reference_type": "transferred_reference", "source": source, "material": material, "line": line,
            "center": center, "confidence": confidence, **kw}


def _obs(samples, be_range, n, **kw):
    return {"reference_type": "provider_calibration_observation", "source": CAL_SOURCE, "samples": samples,
            "BE_range": list(be_range), "n_fits": n,
            "note": "Provider PHI MultiPak fit of the real data; observation only, never a KB reference or starting value.", **kw}


def _m(name, kind, center, bounds, fwhm, fwhm_bounds, role, source, source_type, interp, conf, status="candidate",
       split=None, ratio=None, lit=None, tr=None, obs=None, fwhm_basis=None, **extra):
    d = {"name": name, "kind": kind, "center": center, "bounds": list(bounds), "fwhm": fwhm, "fwhm_bounds": list(fwhm_bounds),
         "role": role, "source": source, "source_type": source_type, "interpretation": interp, "assignment_status": status,
         "confidence": conf, "calibration_candidate": False, "validated": True,
         "provenance": (f"KB 1.4.0 ({MAX_DATE}): starting center from {source_type} reference ({source}); bounds/fwhm_bounds are "
                        "editable implementation windows, not literature constants."),
         "fwhm_basis": fwhm_basis or "implementation default (resolution/pass-energy dependent), not a literature constant"}
    if kind == "doublet":
        d["split"], d["ratio_primary_secondary"] = split, ratio
    if lit:
        d["literature_reference"] = lit
    if tr:
        d["transferred_reference"] = tr
    if obs:
        d["calibration_observations"] = obs
    d.update(extra)
    return d


TR_NOTE_NBTI = "not independently validated for (Nb,Ti)2AlC"
S1, S23 = ["1_Nb2AlC"], ["2_(Nb0.9Ti0.1)2AlC", "3_(Nb0.6Ti0.4)2AlC"]


def _nb_block(target: str, kb: dict, old: dict) -> list:
    """Nb 3d for nb2alc_max (literature, citation pending) or nbti_211_max (transferred from Nb2AlC)."""
    nbti = target == "nbti"
    lit_c = [_lit("UserNb2AlC2026", "Nb 3d5/2", 202.9, "Nb2AlC", secondary={"line": "Nb 3d3/2", "center": 205.5},
                  implied_split=2.6, citation_status="pending")]
    lit_o = [_lit("UserNb2AlC2026", "Nb 3d5/2", 206.8, "Nb2AlC (surface Nb oxide)", secondary={"line": "Nb 3d3/2", "center": 209.6},
                  implied_split=2.8, citation_status="pending")]
    obs_c = [_obs(S1 if not nbti else S23, (202.95, 202.95) if not nbti else (202.73, 202.81), 1 if not nbti else 2)]
    obs_o = [_obs(S1 if not nbti else S23, (207.35, 207.35) if not nbti else (207.07, 207.13), 1 if not nbti else 2)]
    obs_t = [_obs(S1 if not nbti else S23, (204.42, 204.42) if not nbti else (203.83, 203.86), 1 if not nbti else 2)]
    tr_note = {"note": TR_NOTE_NBTI} if nbti else {}
    nb_max = _m("Nb_MAX", "doublet", 202.9, (202.3, 203.5), 1.1, (0.6, 2.0), "framework_metal",
                "UserNb2AlC2026", "transferred" if nbti else "literature",
                ("Nb-C(-Al) 3d5/2 of the MAX phase (Nb2AlC reference 202.9 eV; 3d3/2 205.5 eV). Metallic MAX lines are "
                 "asymmetric; with symmetric pseudo-Voigt profiles the tail is carried by Nb_MAX_2nd. NOTE: Nb_MAX 3d3/2 "
                 "(~205.5-205.7 eV) overlaps the Nb_suboxide_like 3d5/2 window.")
                + (f" Transferred from Nb2AlC; {TR_NOTE_NBTI}." if nbti else ""),
                "low" if nbti else "medium", split=NB_SPLIT, ratio=NB_RATIO,
                lit=None if nbti else lit_c,
                tr=[{**{k: v for k, v in lit_c[0].items() if k not in ("reference_type", "note")}, "reference_type": "transferred_reference",
                     "confidence": "low", **tr_note}] if nbti else None,
                obs=obs_c)
    tail = _m("Nb_MAX_2nd", "doublet", 203.9, (203.4, 204.6), 1.8, (0.8, 2.6), "framework_metal",
              "Naslund2020", "implementation_surrogate",
              ("Line-shape surrogate: high-BE tail of the asymmetric metallic MAX Nb line (Naslund2020 fit MAX metal-C lines "
               "with asymmetric profiles toward higher BE; this engine uses symmetric pseudo-Voigts). Start = Nb_MAX + 1.0 eV "
               "(implementation rule). Alternative: a second Nb environment (non-framework)."),
              "low", status="unresolved", split=NB_SPLIT, ratio=NB_RATIO, obs=obs_t,
              alternative_interpretation="Near-surface Nb-C-O / sub-oxide environment (non-framework).",
              alternative_role="metal_suboxide_candidate")
    subox = next(c for c in old["fit_models"]["Nb 3d"] if c["name"] == "Nb_suboxide_like")   # unchanged
    nb_ox = _m("Nb2O5", "doublet", 206.8, (206.3, 207.8), 1.5, (0.8, 2.3), "metal_oxide",
               "UserNb2AlC2026", "transferred" if nbti else "literature",
               "Surface Nb oxide (Nb2O5-like) on the MAX phase (Nb2AlC reference 206.8 / 209.6 eV)."
               + (f" Transferred from Nb2AlC; {TR_NOTE_NBTI}." if nbti else ""),
               "low" if nbti else "medium", split=NB_SPLIT, ratio=NB_RATIO,
               lit=None if nbti else lit_o,
               tr=([{**{k: v for k, v in lit_o[0].items() if k not in ("reference_type", "note")}, "reference_type": "transferred_reference",
                     "confidence": "low", **tr_note}] if nbti else [])
               + [_tr("Greczynski2018Carbides", "Nb 3d3/2", 209.8, "native oxide on NbC film", "medium",
                      note="supporting value from a related material (binary NbC), not used as the start value")],
               obs=obs_o)
    return [nb_max, tail, subox, nb_ox]


def _c_block(target: str, kb: dict) -> list:
    nbti = target == "nbti"
    halim = {c["name"]: c for c in kb["profiles"]["nb2c_halim"]["fit_models"]["C 1s"]}
    if nbti:
        c_max = _m("C_MAX", "single", 282.2, (281.5, 283.0), 0.9, (0.5, 2.3), "framework_C", "Naslund2020", "transferred",
                   ("Carbide C of (Nb,Ti)2AlC. Start 282.2 eV = implementation choice between the transferred references "
                    f"Nb2AlC 282.5 eV and Ti3AlC2 281.88 eV; {TR_NOTE_NBTI}."), "low",
                   tr=[_tr("UserNb2AlC2026", "C 1s", 282.5, "Nb2AlC", "low", citation_status="pending", note=TR_NOTE_NBTI),
                       _tr("Naslund2020", "C 1s", 281.88, "Ti3AlC2 (312 MAX; Ti2AlC value not reported)", "medium", fwhm=0.57,
                           note=TR_NOTE_NBTI)],
                   obs=[_obs(S23, (281.90, 282.34), 2)])
    else:
        c_max = _m("C_MAX", "single", 282.5, (281.6, 283.2), 0.9, (0.5, 2.3), "framework_C", "UserNb2AlC2026", "literature",
                   "Carbide C-Nb of Nb2AlC (reference 282.5 eV).", "medium",
                   lit=[_lit("UserNb2AlC2026", "C 1s", 282.5, "Nb2AlC", citation_status="pending")],
                   obs=[_obs(S1, (282.86, 282.86), 1)])
    out = [c_max]
    for name in ("C_C", "C_O", "COO"):
        h = halim[name]
        out.append(_m(name, "single", h["center"], h["bounds"], h["fwhm"], h["fwhm_bounds"], "contaminant", h["source"],
                      "transferred", f"{h['interpretation']} Transferred from the Nb2CTx literature model (nb2c_halim).", "medium",
                      tr=[_tr(h["source"], "C 1s", h["center"], "Nb2CTx (contamination C)", "medium")]))
    return out


def _al_block(target: str) -> list:
    nbti = target == "nbti"
    nb_al = _lit("Naslund2020", "Al 2p3/2", 72.40, "Nb2AlC", derived="Al metal 72.76 eV + reported Nb2AlC shift -0.36 eV")
    ti_al = _tr("Naslund2020", "Al 2p3/2", 71.49, "Ti2AlC", "medium", derived="Al metal 72.76 eV + reported Ti2AlC shift -1.27 eV",
                note=TR_NOTE_NBTI)
    al2o3 = _tr("Naslund2020", "Al 2p", 74.5, "Al2O3 impurity on Ti3AlC2 (figure caption value, ~74.5 eV)", "medium")
    if nbti:
        al = _m("Al_MAX", "single", 72.0, (71.2, 72.8), 1.2, (0.6, 2.2), "A_site_metal", "Naslund2020", "transferred",
                ("Al in the (Nb,Ti)2AlC lattice. Start 72.0 eV = implementation choice between transferred Nb2AlC (72.40) and "
                 f"Ti2AlC (71.49) values; {TR_NOTE_NBTI}. Spin-orbit (0.40 eV) unresolved: single band."), "low",
                tr=[{**nb_al, "reference_type": "transferred_reference", "confidence": "medium", "note": TR_NOTE_NBTI}, ti_al],
                obs=[_obs(S23, (71.94, 72.04), 2)])
    else:
        al = _m("Al_MAX", "single", 72.4, (71.7, 72.9), 1.2, (0.6, 2.2), "A_site_metal", "Naslund2020", "literature",
                "Al in the Nb2AlC lattice (Al 2p3/2 72.40 eV derived from the reported -0.36 eV shift vs Al metal). Single band.",
                "medium", lit=[nb_al], obs=[_obs(S1, (72.50, 72.50), 1)])
    ox = _m("Al_oxide", "single", 74.5, (73.6, 75.2), 1.6, (0.8, 2.4), "A_site_oxide", "Naslund2020", "transferred",
            "Al2O3-like surface oxide (transferred from the Al2O3 impurity value reported on Ti3AlC2)."
            + (f" {TR_NOTE_NBTI}." if nbti else ""), "medium", tr=[al2o3],
            obs=[_obs(S1 if not nbti else S23, (74.32, 74.32) if not nbti else (73.94, 73.97), 1 if not nbti else 2)])
    return [al, ox]


def _o_block(target: str, kb: dict) -> list:
    nbti = target == "nbti"
    halim_o = next(c for c in kb["profiles"]["nb2c_halim"]["fit_models"]["O 1s"] if c["name"] == "Nb2O5_O")
    trn = {"note": TR_NOTE_NBTI} if nbti else {}
    return [
        _m("O_oxide", "single", 530.5, (529.6, 531.0), 1.4, (0.7, 2.4), "oxide_O", "Naslund2020", "transferred",
           "Metal-oxide lattice O (Nb2O5/TiO2-like) on the MAX surface.", "medium",
           tr=[_tr(halim_o["source"], "O 1s", halim_o["center"], "Nb2CTx (Nb2O5 O)", "medium", **trn),
               _tr("Naslund2020", "O 1s", 530.55, "TiO2 on Ti3AlC2", "medium", fwhm=1.4, **trn)],
           obs=[_obs(S1 + S23, (529.30, 530.43), 5, note_extra="provider low-BE O at 529.30-529.39 eV (samples 1-2) has no literature counterpart here")]),
        _m("O_OH_ads", "single", 531.9, (531.2, 532.4), 1.6, (0.8, 2.6), "adsorbate", "Naslund2020", "transferred",
           "OH adsorbed on surface oxide (OH-TiO2 value on Ti3AlC2).", "medium",
           tr=[_tr("Naslund2020", "O 1s", 531.9, "OH on TiO2 (Ti3AlC2)", "medium", fwhm=1.6, **trn)],
           obs=[_obs(S1 + S23, (531.66, 531.94), 3)]),
        _m("O_Al2O3", "single", 532.8, (532.3, 533.6), 1.6, (0.8, 2.6), "oxide_O", "Naslund2020", "transferred",
           "Al2O3 O (value on Ti3AlC2).", "medium",
           tr=[_tr("Naslund2020", "O 1s", 532.8, "Al2O3 on Ti3AlC2", "medium", fwhm=1.6, **trn)]),
    ]


def _ti_block() -> list:
    return [
        _m("Ti_MAX", "doublet", 454.75, (454.0, 455.5), 1.0, (0.6, 2.6), "framework_metal", "Naslund2020", "transferred",
           (f"Ti-C of the MAX phase, transferred from Ti3AlC2 (no direct Ti2AlC / (Nb,Ti)2AlC peak table found); {TR_NOTE_NBTI}. "
            "Prioritize Ti 2p3/2 (Nb 3s overlaps Ti 2p1/2)."), "medium", split=TI_SPLIT, ratio=TI_RATIO,
           tr=[_tr("Naslund2020", "Ti 2p3/2", 454.75, "Ti3AlC2", "medium", fwhm=0.88, secondary={"line": "Ti 2p1/2", "center": 460.75},
                   literature_split=6.00, note=TR_NOTE_NBTI),
               _tr("Greczynski2018Carbides", "Ti 2p3/2", 454.7, "TiC (binary carbide)", "low", note=TR_NOTE_NBTI)],
           obs=[_obs(S23, (454.44, 454.51), 2)]),
        _m("Ti_oxide", "doublet", 458.8, (458.0, 459.6), 1.6, (1.0, 3.0), "metal_oxide", "Naslund2020", "transferred",
           f"TiO2-like surface oxide (TiO2 values on TiC / Ti3AlC2); {TR_NOTE_NBTI}.", "medium", split=TI_SPLIT_OX, ratio=TI_RATIO,
           tr=[_tr("Naslund2020", "Ti 2p3/2", 458.8, "TiO2 on TiC", "medium", fwhm=1.2, literature_split=5.8),
               _tr("Naslund2020", "Ti 2p3/2", 459.0, "TiO2 on Ti3AlC2", "medium", fwhm=1.2, literature_split=5.8),
               _tr("Greczynski2018Carbides", "Ti 2p3/2", 458.30, "native oxide on TiC film", "low")],
           obs=[_obs(S23, (458.36, 458.46), 2)]),
    ]


def apply_max_literature(kb: dict) -> None:
    kb["kb_version"] = "1.4.0"
    kb["references"].update(MAX_REFERENCES)
    kb["changelog"].append({"version": "1.4.0", "date": MAX_DATE, "note": (
        "MAX-phase literature calibration: nb2alc_max starting values from Nb2AlC references (user-supplied values, citation "
        "pending; Al 2p from Naslund2020); nbti_211_max = transferred model (Nb2AlC + Ti3AlC2/TiC/Ti2AlC references; no direct "
        "(Nb,Ti)2AlC peak table found), labelled source_type=transferred. Provider fit values kept only as "
        "calibration_observations. Provider-only components dropped (Al ~75.3-76.1 eV; O ~529.3 eV; documented). "
        "Nb_suboxide_like, D4 status and spin-orbit splittings unchanged.")})
    for key, target, material in (("nb2alc_max", "nb2alc", "Nb2AlC"), ("nbti_211_max", "nbti", "(Nb,Ti)2AlC")):
        old = kb["profiles"][key]
        models = {}
        if target == "nbti":
            models["Ti 2p"] = _ti_block()
        models["Nb 3d"] = _nb_block(target, kb, old)
        models["C 1s"] = _c_block(target, kb)
        models["O 1s"] = _o_block(target, kb)
        models["Al 2p"] = _al_block(target)
        notes = [
            ("Literature-calibrated MAX profile (KB 1.4.0). Nb 3d / C 1s Nb2AlC values are user-supplied with the citation "
             "PENDING; Al 2p from Naslund2020." if target == "nb2alc" else
             "TRANSFERRED MAX model (KB 1.4.0): no direct XPS peak table for pristine (Nb,Ti)2AlC was found. Starting values are "
             "transferred from Nb2AlC (citation pending), Ti3AlC2/TiC (Naslund2020, Greczynski2018Carbides) and Ti2AlC Al 2p "
             f"(Naslund2020); {TR_NOTE_NBTI}."),
            REF_FRAME_NOTE,
            "Provider fit values are stored only as calibration_observations (never reference or starting values).",
            "Not modelled (provider-only, no literature counterpart found): Al 2p ~75.3-76.1 eV, O 1s ~529.3-529.4 eV"
            + (", Nb 3d 201.80 eV (sample 1)." if target == "nb2alc" else "."),
            "Nb_MAX 3d3/2 (~205.5-205.7 eV) overlaps the Nb_suboxide_like 3d5/2 window: expect strong area correlation.",
            "Spin-orbit splittings follow the KB convention (D4 calibration-sensitive); literature splittings are recorded "
            "in the reference entries (e.g. Nb2AlC pair implies 2.6 eV; Ti3AlC2 Ti-C 6.00 eV).",
            "Termination (Tx) outputs are not meaningful for a MAX phase; O/F/Cl are surface oxide/contamination.",
        ]
        prof = {"material": material, "phase": "MAX", "notes": notes, "fit_models": models}
        if "core_level_warnings" in old:
            prof["core_level_warnings"] = old["core_level_warnings"]
        kb["profiles"][key] = prof
    kb["calibration_status"]["max_profiles"] = {
        "nb2alc_max": "literature (Nb 3d / C 1s citation pending; Al 2p Naslund2020); O 1s, Al oxide, contamination C transferred",
        "nbti_211_max": "transferred (no direct (Nb,Ti)2AlC peak table found); not independently validated for (Nb,Ti)2AlC",
    }


# ======================================================================================================================
# KB 1.5.0 (2026-11-02): citation for the Nb2AlC values + splitting decisions (user instructions).
#  - UserNb2AlC2026 (citation pending) -> Su2024JMCC (verified; SI Fig. S7). confidence medium: the SI is not primarily an
#    Nb2AlC XPS study and does not clearly document the BE reference.
#  - Nb 3d fit constraint stays 2.8 eV; the implied 2.6 eV is stored only as a literature-reported splitting.
#  - Provider 2.72 eV (Nb 3d) and 6.17 eV (Ti 2p) are stored only as calibration observations.
#  - nbti_211_max ONLY: transferred Ti-C dSO 5.8 -> 6.00 eV (Naslund2020 Ti3AlC2), not independently validated.
#  - No MXene profile is touched.
# ======================================================================================================================
SU2024 = {
    "Su2024JMCC": {
        "title": "MXene quantum dot functionalized bimetallic MOFs as a label-free aptasensor for highly sensitive detection of exosomes",
        "authors": "X. Su, Q. You, P. Wang, L. Li, M. Ge, L. Yang, W.-F. Dong, Z. Chang",
        "journal": "J. Mater. Chem. C", "volume": 12, "pages": "8561-8568", "year": 2024, "doi": "10.1039/D4TC01594K",
        "location": "Supporting Information, Fig. S7 (Nb2AlC XPS)",
        "citation_status": "verified",
        "verification": ("Citation and SI Fig. S7 values supplied by the user (2026-11-02); bibliographic record cross-checked "
                         "against the RSC article page. The SI PDF itself was not machine-read in this session."),
        "be_reference": "not clearly documented in the SI",
        "confidence": "medium",
        "use": ("Nb2AlC: C-Nb C 1s 282.5 eV; Nb-C Nb 3d5/2 202.9 / 3d3/2 205.5 eV; Nb(V) 206.8 / 209.6 eV. Confidence medium: "
                "the SI is not primarily an Nb2AlC XPS study and the BE-reference method is not clearly documented."),
    }
}
SPLIT_IMPLIED = "implied by the reported 3d5/2 and 3d3/2 positions (literature-reported, NOT used as a fit constraint)"


def _walk_refs(prof):
    for comps in prof["fit_models"].values():
        for c in comps:
            for key in ("literature_reference", "transferred_reference"):
                for r in c.get(key) or []:
                    yield c, r


def apply_su2024_and_splits(kb: dict) -> None:
    kb["kb_version"] = "1.5.0"
    kb["references"].pop("UserNb2AlC2026", None)
    kb["references"].update(SU2024)
    for key in ("nb2alc_max", "nbti_211_max"):
        prof = kb["profiles"][key]
        for c, r in _walk_refs(prof):
            if r.get("source") == "UserNb2AlC2026":
                r["source"] = "Su2024JMCC"
                r["citation_status"] = "verified"
                r["location"] = "SI Fig. S7"
                if "implied_split" in r:
                    r["literature_split"] = r.pop("implied_split")
                    r["literature_split_basis"] = SPLIT_IMPLIED
        for comps in prof["fit_models"].values():
            for c in comps:
                if c.get("source") == "UserNb2AlC2026":
                    c["source"] = "Su2024JMCC"
                    c["provenance"] = c["provenance"].replace("UserNb2AlC2026", "Su2024JMCC")
                if c["name"] in ("Nb_MAX", "Nb2O5", "Nb_MAX_2nd"):
                    c["split_basis"] = ("KB Nb 3d convention 2.8 eV (D4 calibration-sensitive); literature-reported / provider "
                                        "values are stored separately and not used as the constraint")
                    for o in c.get("calibration_observations") or []:
                        o["split_obs"] = 2.72
                c["interpretation"] = c["interpretation"].replace("(Nb2AlC reference ", "(Nb2AlC reference, Su2024JMCC SI Fig. S7: ")
                c["interpretation"] = c["interpretation"].replace("Carbide C-Nb of Nb2AlC (reference 282.5 eV).",
                                                                  "Carbide C-Nb of Nb2AlC (Su2024JMCC SI Fig. S7: 282.5 eV).")
        prof["notes"] = [n.replace("Nb 3d / C 1s Nb2AlC values are user-supplied with the citation PENDING",
                                   "Nb 3d / C 1s Nb2AlC values from Su2024JMCC SI Fig. S7 (verified citation; confidence medium, BE reference not documented)")
                          .replace("Nb2AlC (citation pending)", "Nb2AlC (Su2024JMCC)") for n in prof["notes"]]
    ti = next(c for c in kb["profiles"]["nbti_211_max"]["fit_models"]["Ti 2p"] if c["name"] == "Ti_MAX")
    ti["split"] = 6.00
    ti["split_basis"] = ("transferred: Naslund2020 Ti3AlC2 Ti-C dSO 6.00 eV (direct MAX literature); "
                         "not independently validated for (Nb,Ti)2AlC")
    ti["provenance"] += " dSO 5.8 -> 6.00 eV in KB 1.5.0 (transferred from Naslund2020 Ti3AlC2)."
    for o in ti.get("calibration_observations") or []:
        o["split_obs"] = 6.17
        o["split_obs_note"] = "provider dSO from the sample whose Ti 2p1/2 was fitted (Nb 3s overlap); observation only"
    tox = next(c for c in kb["profiles"]["nbti_211_max"]["fit_models"]["Ti 2p"] if c["name"] == "Ti_oxide")
    tox["split_basis"] = "KB Ti oxide convention 5.7 eV (unchanged); Naslund2020 TiO2 5.8 eV recorded in transferred_reference"
    for c in kb["profiles"]["nbti_211_max"]["fit_models"]["Nb 3d"]:
        pass
    kb["profiles"]["nbti_211_max"]["notes"] = [
        n.replace("Ti3AlC2 Ti-C 6.00 eV)", "Ti3AlC2 Ti-C 6.00 eV, which IS used as the transferred Ti_MAX constraint in nbti_211_max)")
        for n in kb["profiles"]["nbti_211_max"]["notes"]]
    d4 = next(d for d in kb["spec_discrepancies"] if d["id"] == "D4")
    d4["decision"] += (" 2026-11-02 (KB 1.5.0): Nb 3d constraint stays 2.8 eV in all profiles; the Nb2AlC 2.6 eV implied by "
                       "Su2024JMCC is a literature-reported splitting only; provider 2.72 eV and 6.17 eV are calibration "
                       "observations; nbti_211_max Ti_MAX uses the transferred Naslund2020 Ti3AlC2 value 6.00 eV; MXene Ti 2p "
                       "splittings unchanged.")
    kb["calibration_status"]["max_profiles"]["nb2alc_max"] = ("literature (Nb 3d / C 1s Su2024JMCC SI Fig. S7, verified, confidence "
                                                             "medium; Al 2p Naslund2020); O 1s, Al oxide, contamination C transferred")
    kb["changelog"].append({"version": "1.5.0", "date": "2026-11-02", "note": (
        "Nb2AlC citation verified: UserNb2AlC2026 -> Su2024JMCC (J. Mater. Chem. C 12 (2024) 8561, SI Fig. S7; confidence medium). "
        "Splitting layers: literature-reported (Nb2AlC 2.6 eV, stored only), fit constraint (Nb 3d 2.8 eV unchanged; "
        "nbti_211_max Ti_MAX 5.8 -> 6.00 eV transferred from Naslund2020), provider observations (2.72, 6.17 eV). "
        "MXene profiles unchanged.")})


def main() -> None:
    seed = json.loads(SEED.read_text(encoding="utf-8"))
    kb = copy.deepcopy(seed)

    kb["kb_version"] = "1.2.0"
    kb["seed_version"] = seed["version"]
    kb["created_at"] = "2026-10-15"
    kb["reviewed_by"] = None  # not yet reviewed by a domain expert
    kb["changelog"] = [
        {"version": "0.1.0", "note": "Prototype seed (handoff §27)."},
        {"version": "1.0.0", "note": "Numerical values unchanged from 0.1.0. Added per-component interpretation caveats (wording from handoff §12-§16), assignment_status, provenance, Chen2023/LegacyExcel references, spec_discrepancies."},
        {"version": "1.1.0", "date": DECIDED, "note": "User decisions D1-D3: D1 §27 kept for Ti3C2Tx, Brette SI values added as literature_reference; D2 nb2c_halim Nb 3d center/FWHM -> Halim SI Nb2CTx 204.8/207.3/209.2 eV, 1.1/1.7/1.0 eV (windows unchanged); D3 Chen2023 -> Chen2024, Halim2016 kept."},
        {"version": "1.2.0", "date": DECIDED, "note": "literature_reference is now a list of SI-table-labelled entries. Added Brette SI O 1s/F 1s values for HF10 (Tables S5-B, S8, S4-B; before sputtering / no sputtering stated) and HF48 (Table S5-A) as reference metadata only; §27 starting values and fitting windows unchanged. Added Brette interpretation of ~531 eV O 1s (disordered/mixed terminations and T-group vacancies, not unique OH) to O_disordered."},
    ]
    kb["references"].update(EXTRA_REFERENCES)
    kb["spec_discrepancies"] = SPEC_DISCREPANCIES
    kb["legacy_excel"]["source"] = "LegacyExcel"
    kb["profiles"]["ti3c2_brette"].setdefault("notes", []).append(BRETTE_SI_F_NOTE)

    bri = {c["name"]: c.get("interpretation") for c in seed["profiles"]["ti3c2_brette"]["fit_models"]["Ti 2p"]}
    bri["Ti_oxide"] = bri["Ti_oxide_oxyfluoride"]
    nbi = {c["name"]: c.get("interpretation") for c in seed["profiles"]["nb2c_halim"]["fit_models"]["Nb 3d"]}

    for pkey, prof in kb["profiles"].items():
        for core, comps in prof["fit_models"].items():
            for c in comps:
                key = (pkey, core, c["name"])
                if not c.get("interpretation"):
                    if key in ADDED_INTERPRETATION:
                        c["interpretation"] = ADDED_INTERPRETATION[key]
                    elif pkey == "nbti_m2c_defect_aware" and core == "Ti 2p":
                        c["interpretation"] = bri[c["name"]] + TRANSFER_TI
                    elif pkey == "nbti_m2c_defect_aware" and core == "Nb 3d":
                        c["interpretation"] = nbi[c["name"]] + TRANSFER_NB
                    else:
                        raise SystemExit(f"No interpretation caveat for {key}")
                c["assignment_status"] = "unresolved" if key in UNRESOLVED else "candidate"
                c["provenance"] = (
                    f"Starting center informed by {c['source']}; bounds and fwhm_bounds are software fit windows, "
                    "not literature constants."
                )
                if pkey == "ti3c2_brette" and (core, c["name"]) in BRETTE_SI:
                    c["literature_reference"] = BRETTE_SI[(core, c["name"])]
                if pkey == "ti3c2_brette" and (core, c["name"]) == ("O 1s", "O_disordered"):
                    c["literature_interpretation"] = BRETTE_O_DISORDERED
                if pkey == "nb2c_halim" and (core, c["name"]) in HALIM_NB2C:
                    center, fwhm = HALIM_NB2C[(core, c["name"])]
                    lo, hi = c["bounds"]
                    flo, fhi = c["fwhm_bounds"]
                    if not (lo <= center <= hi and flo <= fwhm <= fhi):
                        raise SystemExit(f"D2 value outside unchanged window for {c['name']}")
                    c["center"], c["fwhm"] = center, fwhm
                    c["literature_reference"] = [{"source": "Halim2018Vacancy", "dataset": "Nb2CTx", "line": "Nb 3d5/2",
                                                  "center": center, "fwhm": fwhm, "note": HALIM_NB2C_NOTE}]
                    c["provenance"] = ("Starting center/FWHM = Halim SI Nb2CTx value (user decision D2); bounds and "
                                       "fwhm_bounds are editable software fit windows (§27), not literature constants.")
    apply_calibration_phase2(kb)
    apply_calibration_d4_o1s(kb)
    apply_max_literature(kb)
    apply_su2024_and_splits(kb)
    OUT.write_text(json.dumps(kb, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
