"""Measurement-region selection and provider fitting request (handoff §7, §8, §28)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .knowledge_base import KnowledgeBase
from .models import A_SITE_TO_CORE, ELEMENT_TO_CORE, FitComponent, SampleSetup

FLUORIDE_KEYS = ("hf", "lif", "naf", "kf", "nh4f", "fluor")
CHLORIDE_KEYS = ("hcl", "chlor")
N_KEYS = ("nh4", "tbaoh", "amine", "ionic liquid")


def metal_core(element: str) -> tuple[str, bool]:
    """Return (core level, known).  Unknown elements must be chosen manually."""
    if element in ELEMENT_TO_CORE:
        return ELEMENT_TO_CORE[element], True
    return f"{element} core level (select manually)", False


def required_regions(setup: SampleSetup) -> list[dict]:
    et = " ".join(setup.etchants).lower()
    regions: list[dict] = [
        {"region": "Survey", "reason": "Global composition / unexpected residues", "trigger": "always"},
    ]
    for m in setup.m_elements:
        core, known = metal_core(m)
        regions.append({"region": core, "reason": f"MXene framework {m}, oxidation, depletion",
                        "trigger": "M-site element" + ("" if known else " (no seed core-level mapping)")})
    regions += [
        {"region": "C 1s", "reason": "Framework carbide C (normalization route) and contamination", "trigger": "always"},
        {"region": "O 1s", "reason": "Termination / oxide / adsorbate discrimination", "trigger": "always"},
    ]
    if any(k in et for k in FLUORIDE_KEYS):
        regions.append({"region": "F 1s", "reason": "F termination vs fluoride/fluoroaluminate/salt residue", "trigger": "F-containing etchant"})
    if any(k in et for k in CHLORIDE_KEYS) or setup.cl_interest:
        regions.append({"region": "Cl 2p", "reason": "Cl termination / residual chloride", "trigger": "HCl/chloride chemistry or Cl interest"})
    if setup.a_element:
        a_core = A_SITE_TO_CORE.get(setup.a_element, f"{setup.a_element} core level (select manually)")
        regions.append({"region": a_core, "reason": f"Residual A-site element ({setup.a_element}) / reaction products", "trigger": "A-site element"})
    if "lif" in et:
        regions.append({"region": "Li 1s", "reason": "Li-containing residue / intercalant", "trigger": "LiF"})
    if "naf" in et:
        regions.append({"region": "Na 1s", "reason": "Na-containing residue", "trigger": "NaF"})
    if "kf" in et:
        regions.append({"region": "K 2p", "reason": "K-containing residue", "trigger": "KF"})
    if any(k in et for k in N_KEYS):
        regions.append({"region": "N 1s", "reason": "NH4+/TBA+/N-containing reagent residue", "trigger": "N-containing reagent"})
    for d in setup.dopants:
        core, _ = metal_core(d)
        regions.append({"region": core, "reason": f"Dopant {d}: lattice/vacancy incorporation vs adsorbed/secondary phase",
                        "trigger": "post-synthesis dopant"})
    seen, out = set(), []
    for r in regions:
        if r["region"] not in seen:
            seen.add(r["region"]); out.append(r)
    return out


def _component_line(c: FitComponent) -> str:
    s = (f"  - {c.name}: {c.kind}; start BE {c.center:.2f} eV; allowed BE {c.bounds[0]:.2f}–{c.bounds[1]:.2f} eV; "
         f"start FWHM {c.fwhm:.2f} eV (window {c.fwhm_bounds[0]:.2f}–{c.fwhm_bounds[1]:.2f})")
    if c.kind == "doublet":
        s += f"; ΔSO {c.split:.2f} eV; primary:secondary area ratio {c.ratio_primary_secondary:g}:1"
        if c.split_basis:
            s += f" (ΔSO basis: {c.split_basis})"
    s += f"; role = {c.role}"
    if c.assignment_status == "unresolved":
        s += " (UNRESOLVED – please do not assign by BE alone)"
    if c.calibration_candidate and not c.validated:
        s += f" (CALIBRATION CANDIDATE, confidence {c.confidence or 'n/a'} – fit it as a separate component; do not fold it into another"
        s += " component)" if not c.cross_check else f" component; cross-check: {', '.join(c.cross_check)})"
        if c.alternative_interpretation:
            s += f"\n      alternative interpretation: {c.alternative_interpretation}"
    s += f"; source = {c.source}"
    if c.source_type:
        s += f"; start-value basis = {c.source_type}" + (f" (confidence {c.confidence})" if c.confidence else "")
    for r in c.literature_reference or []:
        tag = ", ".join(x for x in (r["source"], r.get("material"), r.get("dataset"), f"Table {r['table']}" if r.get("table") else None,
                                    r.get("sample_state"), r["line"]) if x)
        s += f"\n      literature reference ({tag}): {r['center']:.2f} eV"
        if r.get("fwhm") is not None:
            s += f", FWHM {r['fwhm']:.2f} eV"
        if r.get("si_label"):
            s += f" [{r['si_label']}]"
        if r.get("secondary"):
            sec = r["secondary"]
            s += f"; {sec['line']}: {sec['center']:.2f} eV" + (f", FWHM {sec['fwhm']:.2f} eV" if sec.get("fwhm") is not None else "")
        if r.get("citation_status") == "pending":
            s += " [citation pending]"
        s += " (reference value, not a fixed position)"
    for r in c.transferred_reference or []:
        s += (f"\n      transferred reference ({r['source']}, {r.get('material')}, {r['line']}): {r['center']:.2f} eV "
              f"(confidence {r.get('confidence')}; {r.get('note') or 'related material'})")
    if c.literature_interpretation:
        s += f"\n      note: {c.literature_interpretation}"
    return s


def generate_request(setup: SampleSetup, kb: KnowledgeBase) -> tuple[str, pd.DataFrame]:
    regions = required_regions(setup)
    requested = {r["region"] for r in regions}
    L: list[str] = []
    L += ["XPS measurement / fitting request", "",
          f"Sample: {setup.sample_name}",
          f"Nominal structure: {setup.structure}",
          f"M-site elements: {', '.join(setup.m_elements)}",
          f"Etching / processing chemistry: {', '.join(setup.etchants) if setup.etchants else 'not specified'}",
          f"A-site element to check: {setup.a_element or 'not specified'}"]
    if setup.post_treatment:
        L.append(f"Post-synthesis treatment: {setup.post_treatment}")
    if setup.dopants:
        L.append(f"Post-synthesis dopant(s): {', '.join(setup.dopants)}")
    L += ["", "Purpose:",
          "Evaluate MXene framework chemistry, surface terminations, oxidation/oxyfluoride products,",
          "residual salts, selective metal depletion, and—where relevant—defect/dopant-related changes.",
          "", "Requested spectra:"]
    L += [f"- {r['region']}: {r['reason']}" for r in regions]

    if setup.expected_defects:
        L += ["", "Defect-aware request:",
              "- M-site deficiency/vacancies are expected or possible. Please do NOT force a vacancy-specific",
              "  component solely from a binding-energy shift.",
              "- Please retain raw spectra and provide a residual so vacancy/depletion hypotheses can be",
              "  cross-checked against elemental ratios and other characterization."]
    if setup.dopants:
        L += ["", "Dopant-aware request:"]
        L += [f"- {d}: where data permit, distinguish lattice/vacancy incorporation from adsorbed/intercalated species,"
              f" oxide/fluoride, or other secondary phases. Please do not assign substitution from BE alone." for d in setup.dopants]
    if len(setup.m_elements) > 1:
        L += ["", "Mixed-metal MXene note:",
              "- Fit each metal core level independently.",
              "- Do not impose identical chemical-state fractions on different metals.",
              "- Report integrated MXene-framework contribution for each metal separately so selective",
              "  depletion can be assessed without renormalizing total M-site occupancy to the ideal value."]

    L += ["", "General fitting requirements:",
          "- Provide raw spectrum, background, individual fitted components, total envelope, and residual.",
          "- For every component provide BE, integrated area, area %, FWHM, line shape, and constraints.",
          "- For doublets, state whether the reported area is the total of both spin-orbit components or the primary line only.",
          "- State binding-energy referencing method explicitly.",
          "- State X-ray source, pass energy, step size, spot size, charge-neutralization settings, and sputtering.",
          "- Do not add components solely to improve visual agreement; optional components need statistical/",
          "  chemical justification.",
          "", "Metal-core guidance:",
          "- For spin-orbit doublets constrain splitting, area ratio, and paired FWHM physically."]
    if "Ti" in setup.m_elements or "Ti" in setup.dopants:
        L.append("- Ti 2p: do not force one-to-one Ti2+/Ti3+/Ti4+ interpretation; termination/local structure/electronic"
                 " effects matter strongly.")
    if "Nb" in setup.m_elements or "Nb" in setup.dopants:
        L.append("- Nb 3d: begin with framework Nb + oxide + fluoride-like contribution when required; do not force"
                 " a unique vacancy-only Nb peak.")
        if any(c.name == "Nb_suboxide_like" for c in kb.fit_models(setup.profile_key).get("Nb 3d", [])):
            L.append("- Nb 3d: report any ~205.6–205.9 eV (3d5/2) doublet as a SEPARATE component; do not merge it into"
                     " Nb2O5 or framework Nb (calibration candidate, under evaluation).")
    L += ["", "C/O/F guidance:",
          "- C 1s: separate framework carbide C from adventitious/oxygenated carbon.",
          "- O 1s: keep oxide O, termination-related O, OH-like O, adsorbed H2O/organic O conceptually separate."]
    if len(setup.m_elements) > 1:
        L.append("  In mixed-metal samples, low/mid-BE O components are not assigned by BE alone.")
    L += ["- F 1s: separate main MXene-F from higher-BE fluoride/fluoroaluminate/salt residue if needed.",
          "  Do not count all F automatically as Tx."]

    models = kb.fit_models(setup.profile_key)
    L += ["", f"Recommended starting fit model (knowledge base v{kb.version}, profile '{setup.profile_key}'):",
          "- Positions are literature-informed starting values and windows, not fixed assignments.",
          "- If a candidate component is not statistically/chemically required, please report that rather than forcing it."]
    for core, comps in models.items():
        if core not in requested:
            continue
        n_d = sum(c.kind == "doublet" for c in comps)
        label = f"{len(comps)} candidate doublet(s)" if n_d == len(comps) else f"{len(comps)} candidate component(s)"
        L.append(f"\n{core}: {label}")
        L += [_component_line(c) for c in comps]
    no_model = [r["region"] for r in regions if r["region"] not in models and r["region"] != "Survey"]
    if no_model:
        L += ["", "Regions without an automated literature model in this app (fit manually; report components and constraints):",
              "  " + ", ".join(no_model)]

    L += ["", "Requested deliverables:",
          "1. Native fitting project file if possible",
          "2. Raw BE-intensity data for every high-resolution region",
          "3. Component table: name, BE, FWHM, area, area %, line shape, constraints",
          "4. Survey-derived concentrations and RSFs used",
          "5. Exported fitted/component traces if available", ""]
    return "\n".join(L), pd.DataFrame(regions)


def make_peak_template(setup: SampleSetup, kb: KnowledgeBase) -> pd.DataFrame:
    """Provider component template (§28), restricted to requested regions with a KB model.

    As in the §28 example, 'contaminant' components (adventitious/oxygenated C) are omitted because
    they never enter quantification; they remain in the fit model and the request text.
    """
    requested = {r["region"] for r in required_regions(setup)}
    rows = []
    for core, comps in kb.fit_models(setup.profile_key).items():
        if core not in requested:
            continue
        for c in comps:
            if c.role == "contaminant":
                continue
            rows.append({"sample": setup.sample_name, "core_level": core, "element": core.split()[0],
                         "component": c.name, "BE": c.center, "FWHM": c.fwhm, "area": np.nan,
                         "role": c.role, "notes": ""})
    for d in setup.dopants:
        core, _ = metal_core(d)
        for suffix, role in (("lattice_or_vacancy", "dopant_lattice_candidate"),
                             ("adsorbed", "dopant_adsorbed"), ("oxide", "dopant_oxide")):
            rows.append({"sample": setup.sample_name, "core_level": core, "element": d,
                         "component": f"{d}_{suffix}", "BE": np.nan, "FWHM": np.nan, "area": np.nan,
                         "role": role, "notes": ""})
    return pd.DataFrame(rows, columns=["sample", "core_level", "element", "component", "BE", "FWHM", "area", "role", "notes"])
