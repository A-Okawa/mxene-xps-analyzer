"""Load and validate the versioned knowledge base (handoff §11, §25, §27)."""
from __future__ import annotations

import dataclasses
import json
import math
from functools import lru_cache
from pathlib import Path

from .models import FitComponent

DEFAULT_KB_PATH = Path(__file__).resolve().parents[2] / "data" / "knowledge_base.json"


class KnowledgeBase:
    def __init__(self, raw: dict):
        self.raw = raw
        self.references: dict = raw["references"]
        self.structures: dict = raw["structures"]
        self.validation_defaults: dict = raw["validation_defaults"]
        self.legacy: dict = raw["legacy_excel"]
        self._models: dict[str, dict[str, list[FitComponent]]] = {}
        for pkey, prof in raw["profiles"].items():
            self._models[pkey] = {
                core: [FitComponent.from_dict(c, core) for c in comps]
                for core, comps in prof.get("fit_models", {}).items()
            }
        self._check_sources()

    def _check_sources(self) -> None:
        for pkey, models in self._models.items():
            for core, comps in models.items():
                names = [c.name for c in comps]
                if len(names) != len(set(names)):
                    raise ValueError(f"{pkey}/{core}: duplicate component names")
                for c in comps:
                    if c.source not in self.references:
                        raise ValueError(f"{pkey}/{core}/{c.name}: source {c.source!r} not in references")
        from .provenance import check_kb_provenance   # provider observations must never act as references
        bad = check_kb_provenance(self)
        if bad:
            raise ValueError("KB provenance violation(s): " + "; ".join(bad))

    @property
    def version(self) -> str:
        return str(self.raw.get("kb_version", self.raw.get("version")))

    @property
    def profile_keys(self) -> list[str]:
        return list(self._models)

    def profile(self, key: str) -> dict:
        return self.raw["profiles"][key]

    def fit_models(self, key: str) -> dict[str, list[FitComponent]]:
        return self._models[key]

    def components(self, key: str, core: str) -> list[FitComponent]:
        return self._models[key][core]

    def reference_text(self, source_key: str) -> str:
        ref = self.references.get(source_key, {})
        doi = ref.get("doi")
        return f"{ref.get('title', source_key)}" + (f" — DOI {doi}" if doi else "")

    def role_of(self, profile_key: str, component: str, core: str | None = None) -> str | None:
        for c_core, comps in self._models[profile_key].items():
            if core and c_core != core:
                continue
            for c in comps:
                if c.name == component:
                    return c.role
        return None

    def to_json(self) -> str:
        return json.dumps(self.raw, indent=2, ensure_ascii=False)


EDITABLE_FIELDS = ("center", "bounds", "fwhm", "fwhm_bounds")


def with_overrides(components: list[FitComponent], overrides: dict[str, dict]) -> tuple[list[FitComponent], list[dict]]:
    """Apply user edits of starting values / fitting windows (implementation bounds) to a fit model.

    ``overrides`` maps component name -> {field: value} for fields in EDITABLE_FIELDS.  Roles, kind,
    spin-orbit split/ratio, sources and literature_reference cannot be overridden.  The knowledge base
    itself is not modified.  Returns (new components, change log rows).  Invalid edits (e.g. a starting
    value outside its window) raise ValueError via FitComponent validation.
    """
    by_name = {c.name: c for c in components}
    unknown = sorted(set(overrides) - set(by_name))
    if unknown:
        raise ValueError(f"Overrides for unknown component(s): {unknown}")
    out, log = [], []
    for c in components:
        edits = overrides.get(c.name, {})
        bad = sorted(set(edits) - set(EDITABLE_FIELDS))
        if bad:
            raise ValueError(f"{c.name}: field(s) {bad} are not editable (editable: {list(EDITABLE_FIELDS)})")
        changes = {}
        for f, v in edits.items():
            new = tuple(float(x) for x in v) if f in ("bounds", "fwhm_bounds") else float(v)
            if isinstance(new, tuple):
                if len(new) != 2 or not new[0] < new[1]:
                    raise ValueError(f"{c.name}: {f} must be (low, high) with low < high")
            if any(not math.isfinite(x) for x in (new if isinstance(new, tuple) else (new,))):
                raise ValueError(f"{c.name}: {f} must be finite")
            old = getattr(c, f)
            if new != old:
                changes[f] = new
                log.append({"component": c.name, "field": f, "kb_value": old, "user_value": new})
        out.append(dataclasses.replace(c, **changes) if changes else c)
    return out, log


def load_kb(path: str | Path = DEFAULT_KB_PATH) -> KnowledgeBase:
    with open(path, "r", encoding="utf-8") as f:
        return KnowledgeBase(json.load(f))


@lru_cache(maxsize=1)
def default_kb() -> KnowledgeBase:
    return load_kb(DEFAULT_KB_PATH)
