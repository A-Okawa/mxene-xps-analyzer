"""Release / reproducibility metadata embedded in every export and in release/release_metadata_*.json."""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone
from importlib import metadata as _md
from pathlib import Path

from . import __version__

REPO = Path(__file__).resolve().parents[2]
PACKAGES = ("numpy", "scipy", "pandas", "streamlit", "matplotlib")

# Known scientific open issues at this release (see docs/CALIBRATION_LOG.md); unchanged by software releases.
OPEN_ISSUES = (
    "Nb 3d ~205.7-205.9 eV component (Nb_suboxide_like) is not chemically validated; O 1s cross-check inconclusive.",
    "Apparent M-site deficiency delta is model-dependent and provisional (candidate excluded vs counted as framework).",
    "D4: Nb 3d spin-orbit splitting is calibration-sensitive (KB constraint 2.8 eV kept).",
    "Ti 2p in mixed Nb/Ti samples is affected by Nb 3s overlap; interpret Ti 2p3/2 first.",
    "O1: MAX Nb_MAX and Nb_MAX_2nd are strongly correlated; interpret their sum.",
    "O2: provider-fraction vs Python-fraction decompositions give different delta.",
    "O3: sample 2 provider Ti quantification is biased low (only Ti 2p3/2 fitted).",
    "O4: sample 1 MAX M_norm = 1.00 is a numerical coincidence; not interpreted.",
    "Sample-1 low-BE Nb 3d feature (unresolved_low_BE_Nb) is real extra intensity but non-identifiable; unassigned.",
    "Cl 2p, S 2p and Li 1s have no validated KB fit model.",
    "Transferred models (e.g. nbti_211_max) are not equivalent to direct literature validation.",
)


def package_versions() -> dict:
    out = {}
    for p in PACKAGES:
        try:
            out[p] = _md.version(p)
        except _md.PackageNotFoundError:
            out[p] = None
    return out


def release_info(kb=None, commit: str | None = "auto") -> dict:
    if commit == "auto":
        from .reporting import git_commit
        commit = git_commit(REPO)
    return {
        "app": "mxene-xps-analyzer", "app_version": __version__,
        "kb_version": None if kb is None else kb.version,
        "git_commit": commit, "python": sys.version.split()[0], "platform": platform.platform(),
        "packages": package_versions(),
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "scientific_open_issues": list(OPEN_ISSUES),
        "interpretation_levels": {
            "fit_succeeded": "optimizer converged and residuals are reported; says nothing about uniqueness or chemistry",
            "statistically_identifiable": "component parameters are determined by the data within this model (see identifiability status)",
            "chemically_validated": "assignment confirmed by orthogonal evidence; never claimed by this app",
        },
    }
