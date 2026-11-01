"""Typed schemas and the role ontology (handoff §6, §24)."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Literal, Optional

# --- Role ontology (§6) -------------------------------------------------------
ROLE_LABELS: dict[str, str] = {
    "framework_C": "MXene framework C",
    "framework_metal": "MXene framework metal",
    "metal_oxide": "metal oxide / oxyfluoride",
    "metal_fluoride": "metal fluoride residue",
    "metal_suboxide_candidate": "sub-oxide / Nb(III)-O-like candidate (non-framework by default; calibration candidate)",
    "termination_O": "O termination",
    "termination_O_disordered": "disordered O termination candidate",
    "termination_OH_like": "OH-like termination candidate",
    "oxide_O": "oxide O",
    "O_unresolved_low": "unresolved low-BE O (not counted as Tx)",
    "O_unresolved_mid": "unresolved mid-BE O (not counted as Tx)",
    "termination_F": "F termination",
    "residue_F": "F-containing residue (not Tx by default)",
    "unresolved_low_BE_F": "unresolved low-BE F (no chemistry assigned; not counted as Tx)",
    "termination_Cl": "Cl termination",
    "A_site_metal": "A-site metal in the MAX lattice (MAX profiles only)",
    "A_site_oxide": "A-site oxide / reaction product (e.g. Al2O3-like)",
    "dopant_lattice_candidate": "dopant lattice/vacancy-incorporation candidate",
    "dopant_adsorbed": "adsorbed/intercalated dopant",
    "dopant_oxide": "dopant oxide / secondary phase",
    "contaminant": "contamination",
    "adsorbate": "adsorbate",
    "other": "other / unclassified",
}
ROLES = tuple(ROLE_LABELS)
SOURCE_TYPES = ("literature", "transferred", "implementation_surrogate", "provider_derived")
DOPANT_ROLES = ("dopant_lattice_candidate", "dopant_adsorbed", "dopant_oxide")

# --- Status vocabulary (§20) -------------------------------------------------
STATUSES = ("PASS", "CHECK", "REVIEW", "CHECK-MISSING", "FAIL")
_SEVERITY = {s: i for i, s in enumerate(STATUSES)}


def worst_status(statuses) -> str:
    """Return the most severe status (FAIL > CHECK-MISSING > REVIEW > CHECK > PASS)."""
    items = [s for s in statuses if s in _SEVERITY]
    return max(items, key=_SEVERITY.__getitem__) if items else "CHECK"


# --- Element / core-level maps (§7) -------------------------------------------
ELEMENT_TO_CORE: dict[str, str] = {
    "Ti": "Ti 2p", "Nb": "Nb 3d", "V": "V 2p", "Mo": "Mo 3d", "Zr": "Zr 3d",
    "Cr": "Cr 2p", "Mn": "Mn 2p", "La": "La 3d", "Ru": "Ru 3p", "Fe": "Fe 2p",
    "Co": "Co 2p", "Ni": "Ni 2p",
}
A_SITE_TO_CORE: dict[str, str] = {"Al": "Al 2p", "Si": "Si 2p", "Ga": "Ga 3d"}


def infer_element(core_level: str) -> str:
    """'Ti 2p' -> 'Ti'.  Core-level strings are '<Element> <orbital>'."""
    return str(core_level).strip().split()[0] if str(core_level).strip() else ""


# --- Schemas (§24) ------------------------------------------------------------
@dataclass(frozen=True)
class FitComponent:
    name: str
    core_level: str
    kind: Literal["single", "doublet"]
    center: float
    bounds: tuple[float, float]
    fwhm: float
    fwhm_bounds: tuple[float, float]
    role: str
    source: str
    interpretation: Optional[str] = None
    split: Optional[float] = None
    ratio_primary_secondary: Optional[float] = None
    assignment_status: str = "candidate"
    provenance: Optional[str] = None
    # Literature-reported values (e.g. Brette SI, Halim SI), one entry per dataset/SI table. Display and
    # reference metadata only: kept separate from the starting value and the editable fitting windows and
    # never used as a constraint by fitting, validation or quantification.
    literature_reference: Optional[list] = None
    literature_interpretation: Optional[str] = None
    # Calibration-phase metadata (KB >= 1.3.0). A calibration candidate is fitted, but its chemical role is not
    # validated: quantification keeps it out of the framework by default and flags dependent outputs provisional.
    calibration_candidate: bool = False
    validated: bool = True
    confidence: Optional[str] = None                  # e.g. "medium"
    alternative_interpretation: Optional[str] = None
    alternative_role: Optional[str] = None            # role under the alternative interpretation (reported, not used)
    cross_check: Optional[list] = None                # core levels that must be checked, e.g. ["O 1s"]
    # KB >= 1.4.0: provenance class of the STARTING value, and three separate reference layers (display/metadata
    # only; never constraints):
    #   literature_reference          – values reported for this material (citation required, or citation pending)
    #   transferred_reference         – values reported for related materials, transferred (not independently validated)
    #   calibration_observations      – provider fit values of the real data (observations, never reference values)
    source_type: Optional[str] = None                 # literature | transferred | implementation_surrogate | provider_derived
    transferred_reference: Optional[list] = None
    calibration_observations: Optional[list] = None
    fwhm_basis: Optional[str] = None
    split_basis: Optional[str] = None                 # provenance of the doublet splitting used as a FIT CONSTRAINT

    def __post_init__(self) -> None:
        if self.kind not in ("single", "doublet"):
            raise ValueError(f"{self.name}: kind must be single/doublet")
        if self.role not in ROLE_LABELS:
            raise ValueError(f"{self.name}: unknown role {self.role!r}")
        lo, hi = self.bounds
        if not lo <= self.center <= hi:
            raise ValueError(f"{self.name}: center {self.center} outside bounds {self.bounds}")
        flo, fhi = self.fwhm_bounds
        if not 0 < flo <= self.fwhm <= fhi:
            raise ValueError(f"{self.name}: fwhm {self.fwhm} outside fwhm_bounds {self.fwhm_bounds}")
        if self.kind == "doublet":
            if self.split is None or self.ratio_primary_secondary is None or self.ratio_primary_secondary <= 0:
                raise ValueError(f"{self.name}: doublet requires split and a positive ratio_primary_secondary")
        if not self.source:
            raise ValueError(f"{self.name}: every component needs a source key")
        if self.source_type is not None and self.source_type not in SOURCE_TYPES:
            raise ValueError(f"{self.name}: unknown source_type {self.source_type!r}")
        for lst, kind in ((self.literature_reference, "literature_reference"),
                          (self.transferred_reference, "transferred_reference"),
                          (self.calibration_observations, "provider_calibration_observation")):
            for e in lst or []:
                if e.get("reference_type", kind) != kind:
                    raise ValueError(f"{self.name}: entry of type {e.get('reference_type')!r} stored under {kind}")
        if not self.interpretation:
            raise ValueError(f"{self.name}: every component needs an interpretation caveat (§11)")
        if self.alternative_role is not None and self.alternative_role not in ROLE_LABELS:
            raise ValueError(f"{self.name}: unknown alternative_role {self.alternative_role!r}")
        if self.calibration_candidate and self.validated and self.assignment_status == "calibration_candidate":
            raise ValueError(f"{self.name}: a calibration candidate marked validated must change assignment_status")

    @classmethod
    def from_dict(cls, d: dict, core_level: str) -> "FitComponent":
        return cls(
            name=d["name"], core_level=core_level, kind=d.get("kind", "single"),
            center=float(d["center"]), bounds=tuple(float(v) for v in d["bounds"]),
            fwhm=float(d["fwhm"]), fwhm_bounds=tuple(float(v) for v in d["fwhm_bounds"]),
            role=d["role"], source=d["source"], interpretation=d.get("interpretation"),
            split=None if d.get("split") is None else float(d["split"]),
            ratio_primary_secondary=None if d.get("ratio_primary_secondary") is None else float(d["ratio_primary_secondary"]),
            assignment_status=d.get("assignment_status", "candidate"), provenance=d.get("provenance"),
            literature_reference=d.get("literature_reference"),
            literature_interpretation=d.get("literature_interpretation"),
            calibration_candidate=bool(d.get("calibration_candidate", False)),
            validated=bool(d.get("validated", not d.get("calibration_candidate", False))),
            confidence=d.get("confidence"), alternative_interpretation=d.get("alternative_interpretation"),
            alternative_role=d.get("alternative_role"), cross_check=d.get("cross_check"),
            source_type=d.get("source_type"), transferred_reference=d.get("transferred_reference"),
            calibration_observations=d.get("calibration_observations"), fwhm_basis=d.get("fwhm_basis"),
            split_basis=d.get("split_basis"),
        )

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def _ref_text(r: dict) -> str:
        tag = "/".join(str(x) for x in (r.get("material"), r.get("dataset"),
                                         f"Table {r['table']}" if r.get("table") else None) if x)
        t = f"{r['center']} eV" + (f" (FWHM {r['fwhm']})" if r.get("fwhm") is not None else "")
        if r.get("si_label"):
            t += f" {r['si_label']}"
        if r.get("secondary"):
            sec = r["secondary"]
            t += f"; {sec['line']} {sec['center']} eV" + (f" (FWHM {sec['fwhm']})" if sec.get("fwhm") is not None else "")
        if r.get("citation_status") == "pending":
            t += " [citation pending]"
        return f"{r['source']} {tag}: {t}" if tag else f"{r['source']}: {t}"

    def literature_summary(self) -> str:
        """One-line, human-readable list of literature reference values (display only)."""
        return " | ".join(self._ref_text(r) for r in self.literature_reference or [])

    def transferred_summary(self) -> str:
        return " | ".join(self._ref_text(r) + f" (transferred, {r.get('confidence')})" for r in self.transferred_reference or [])

    def observation_summary(self) -> str:
        return " | ".join(f"provider {o['BE_range'][0]:.2f}-{o['BE_range'][1]:.2f} eV ({', '.join(o['samples'])})"
                          for o in self.calibration_observations or [])


@dataclass
class PeakRow:
    sample: str
    core_level: str
    element: str
    component: str
    area: float
    role: str
    BE: Optional[float] = None
    FWHM: Optional[float] = None
    notes: Optional[str] = None


BE_REFERENCE_METHODS = (
    "unknown/other", "Fermi edge", "metallic reference (specify)", "adventitious C 1s",
    "framework carbide C internal reference", "provider correction",
)
SURFACE_STATES = ("as-received/as-prepared", "after Ar+ sputtering", "depth-profiled")


@dataclass
class MeasurementMetadata:
    xray_source: Optional[str] = None
    pass_energy_eV: Optional[float] = None
    step_eV: Optional[float] = None
    spot_size: Optional[str] = None
    charge_neutralization: Optional[str] = None
    be_reference_method: Optional[str] = None
    be_reference_detail: Optional[str] = None
    surface_state: str = "as-received/as-prepared"
    sputter_species: Optional[str] = None
    sputter_energy_eV: Optional[float] = None
    sputter_time_s: Optional[float] = None
    provider: Optional[str] = None
    instrument: Optional[str] = None

    @property
    def sputtered(self) -> bool:
        return self.surface_state != "as-received/as-prepared"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["sputtered"] = self.sputtered
        return d


ORTHOGONAL_METHODS = ("ICP", "bulk/particle EDS", "HAADF-STEM", "EELS", "ToF-SIMS", "EPR", "XAS/XANES/EXAFS", "XRD/Rietveld")


@dataclass
class OrthogonalEvidence:
    method: str
    available: bool = False
    result_summary: Optional[str] = None
    supports: list[str] = field(default_factory=list)


@dataclass
class SampleSetup:
    sample_name: str
    structure: str
    m_target: float
    x_target: float
    profile_key: str
    m_elements: list[str]
    etchants: list[str] = field(default_factory=list)
    a_element: str = "Al"
    expected_defects: bool = False
    dopants: list[str] = field(default_factory=list)
    nominal_metal_stoich: dict[str, float] = field(default_factory=dict)
    parent_composition: Optional[str] = None
    post_treatment: Optional[str] = None
    termination_family: Optional[str] = None
    cl_interest: bool = False
    notes: Optional[str] = None
