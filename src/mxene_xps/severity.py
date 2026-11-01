"""Consistent message severity for user-facing warnings.

  INFO     neutral information (import notes, inferred metadata, how a value was computed)
  CAUTION  the fit/result is usable but its chemical meaning is uncertain or provisional (calibration candidates,
           provisional delta, unresolved O, transferred models, referencing caveats, estimated areas, user edits)
  WARNING  a statistical/model problem that limits what can be read from this fit (strong (anti)correlation,
           boundary-limited or non-identifiable components, model-sensitive quantities, QC FAIL, BE discontinuities)
  INVALID  the operation is not valid and was refused or failed (wrong core-level model, implausibly wide BE span,
           failed import/fit, invalid input)

Messages keep their original text; the severity is attached, never used to change a result.
"""
from __future__ import annotations

import re

LEVELS = ("INFO", "CAUTION", "WARNING", "INVALID")

_RULES = (
    ("INVALID", r"implausibly wide|refus|unrelated model|fit failed|import failed|comparison failed|assembly failed|"
                r"invalid model edit|non-numeric area|missing required column|rsf is required|rsf required|"
                r"no positive framework_c|no component areas entered|at least \d+ (finite )?(data )?points"),
    ("WARNING", r"non-identifiable|boundary-limited|anticorrelated|strongly correlated|edge of its allowed window|"
                r"outside the measured be range|area < 1 %|model-sensitive|\bfail\b|role-grouped qc has fail|"
                r"unpaired band|no provider summary matches|segments separated by discontinuities|duplicate be|"
                r"covariance could not be estimated"),
    ("CAUTION", r"provisional|calibration.candidate|not validated|apparent metal deficiency|exceeds the ideal|unresolved|"
                r"max.phase|atomic_percent input|nb 3s|transferred|reclassif|roles changed|edit\(s\) relative|"
                r"mapping change|only one spin-orbit line|estimated|model-dependent|weakly identifiable|\breview\b|"
                r"dopant|unusually wide|f residue|referencing|reference is|sputter|not yet recorded|"
                r"no validated kb model|delta not finalized|not in the nominal|no .* framework_metal component|"
                r"role-grouped qc requires review|individual components differ|outside every kb fit window"),
)


def classify(text: str, default: str = "INFO") -> str:
    t = str(text).lower()
    for level, pat in _RULES:
        if re.search(pat, t):
            return level
    return default


def tag(text: str, level: str | None = None) -> str:
    lv = level or classify(text)
    return f"[{lv}] {text}"


def records(messages, default: str = "INFO") -> list[dict]:
    return [{"severity": classify(m, default), "text": str(m)} for m in messages or []]


QC_STATUS_SEVERITY = {"PASS": "INFO", "CHECK": "CAUTION", "REVIEW": "CAUTION", "CHECK-MISSING": "CAUTION", "FAIL": "WARNING"}
