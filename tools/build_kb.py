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
    OUT.write_text(json.dumps(kb, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
