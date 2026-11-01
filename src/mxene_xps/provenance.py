"""Provenance vocabulary for every reported value, and guards that keep provider observations out of references."""
from __future__ import annotations

PROVENANCE = (
    "raw_measurement",          # imported BE/intensity spectrum
    "provider_raw_fit",         # provider band exactly as reported (incl. provider atomic-concentration table)
    "provider_paired_doublet",  # provider bands paired into one spin-orbit component by the importer
    "python_fit",               # independent Python refit
    "literature_reference",     # value reported for the same material (KB literature_reference)
    "transferred_reference",    # value reported for a related material (KB transferred_reference)
    "calibration_observation",  # provider value of the real data stored as an observation (never a reference)
    "user_override",            # any user edit (model window/start, mapping, role reclassification)
    "derived_quantity",         # computed from other values (quantification, ratios, spans)
)

SOURCE_TYPE_TO_PROVENANCE = {
    "literature": "literature_reference",
    "transferred": "transferred_reference",
    "implementation_surrogate": "derived_quantity",
    "provider_derived": "calibration_observation",
    None: "literature_reference",   # pre-1.4 KB components: §27 / literature-seeded starting values
}


def start_value_provenance(component) -> str:
    """Provenance class of a KB component's starting value."""
    if component.calibration_candidate and not component.validated:
        return "calibration_observation"
    return SOURCE_TYPE_TO_PROVENANCE.get(component.source_type, "literature_reference")


def check_kb_provenance(kb) -> list[str]:
    """Return violations: provider/data-derived sources used as literature or as constraints of non-candidates."""
    bad = []
    for pkey in kb.profile_keys:
        for core, comps in kb.fit_models(pkey).items():
            for c in comps:
                ref = kb.references.get(c.source, {})
                cited = bool(ref.get("doi")) or ref.get("citation_status") == "verified"
                if not cited and not c.calibration_candidate:
                    bad.append(f"{pkey}/{core}/{c.name}: starting value from uncited source {c.source!r} without being a calibration candidate")
                if c.source_type == "literature" and not cited:
                    bad.append(f"{pkey}/{core}/{c.name}: source_type literature but source {c.source!r} is not cited")
                for r in (c.literature_reference or []) + (c.transferred_reference or []):
                    rr = kb.references.get(r.get("source"), {})
                    if not (rr.get("doi") or rr.get("citation_status") == "verified"):
                        bad.append(f"{pkey}/{core}/{c.name}: reference entry from uncited source {r.get('source')!r}")
    return bad
