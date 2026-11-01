"""Shared calibration-phase constants (single source for the app and tools/calibration scripts).

These are QC / sensitivity-analysis settings, not scientific constants, and they do not change the KB.
"""
from __future__ import annotations

# D4 fixed spin-orbit grid for Nb 3d (tools/calibration/d4_nb3d_split_sensitivity.py; CALIBRATION_LOG "D4").
D4_NB3D_SPLIT_GRID = (2.70, 2.72, 2.75, 2.80)
# Materiality criteria declared for D4 and reused for "model-sensitive" flags.
MATERIAL_AREA_FRACTION_SPAN = 0.05
MATERIAL_DELTA_SPAN = 0.05

# Identifiability thresholds (software QC defaults, editable in the UI; not physical constants).
IDENT_THRESHOLDS = {
    "rel_stderr_weak": 0.30,      # area relative standard error above which a component is weakly identifiable
    "rel_stderr_non": 1.00,       # ... non-identifiable
    "corr_weak": 0.80,            # |r| between areas (matches the existing correlation diagnostic)
    "corr_non": 0.95,
    "alt_span_weak": MATERIAL_AREA_FRACTION_SPAN,   # area-fraction span over fitting alternatives
    "alt_span_non": 0.10,
    "min_area_fraction": 0.01,    # matches the existing "< 1 %" diagnostic
    "near_bound_frac": 0.05,      # within 5 % of the window span counts as "near a bound"
}

# Core levels that are recognised/imported but have no validated KB fit model yet.
UNSUPPORTED_REGIONS = {
    "Cl 2p": "Spin-orbit doublet (2p3/2, 2p1/2); provider pairs are imported and kept as doublets.",
    "S 2p": "Spin-orbit doublet (2p3/2, 2p1/2); provider pairs are imported and kept as doublets.",
    "Li 1s": ("May be weak or absent; the provider annotated a possible overlap with an Nb line ('Li1s?'). "
              "Do not interpret without a validated model."),
}
NO_MODEL_TEXT = "no validated KB model"

# Calibration observations already studied (CALIBRATION_LOG). Used only to label provider components that have no
# KB counterpart; never a chemical assignment and never a fit constraint.
KNOWN_OBSERVATIONS = (
    {"label": "unresolved_low_BE_Nb", "core_level": "Nb 3d", "be_range": (201.3, 202.3), "phase": "MAX",
     "identifiability": "non-identifiable",
     "note": ("Sample-1 low-BE Nb 3d intensity: present in every fitting setting but position, width, area and 3d3/2 partner are "
              "model-dependent (tools/calibration/s1_lowbe_nb_sensitivity.py). Not assigned to any species."),
     "study": "CALIBRATION_LOG: Sample-1 low-BE Nb 3d5/2 feature"},
)
