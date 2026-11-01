"""Dataset-level assembly for cross-sample comparison (provider report + optional raw spectra + sample metadata).

Nothing here assigns chemistry: provider components are mapped tentatively to the KB model of each sample's
structure (same logic as the validation tab), unmapped components keep their provider label plus a neutral
'unresolved_*' tag, and every value carries a provenance class.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from .calibration import KNOWN_OBSERVATIONS, UNSUPPORTED_REGIONS
from .fitting import fit_core_level
from .identifiability import assess
from .io import phi_region_to_core
from .knowledge_base import KnowledgeBase
from .multipak import map_to_kb, pair_spin_orbit, pptx_slide_texts, read_multipak_provider
from .quantification import QuantOptions, provisional_report, quantify_sample
from .sample_metadata import infer_structure, metal_fractions


# --- provider atomic concentration table (survey quantification) ---------------------------------------------
def parse_atomic_concentration(text: str, origin: str = "") -> pd.DataFrame:
    """Parse MultiPak 'Atomic Concentration Table' blocks -> long table (data_file, comment, core_level, RSF, at%)."""
    text = text.replace("|", "\n")
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    out, meta = [], {}
    for i, l in enumerate(lines):
        m = re.match(r"File Name:\s*(\S+)", l)
        if m:
            meta = {"data_file": m.group(1)}
        m = re.match(r"Comment:\s*(.*)", l)
        if m and "Data File" not in l:
            meta["comment"] = m.group(1).strip()
        if l.startswith("Abscissa"):
            regions = l.split()[1:]
            rsf, conc = None, None
            for nxt in lines[i + 1:i + 6]:
                toks = nxt.split()
                if toks and toks[-1] == "RSF" and rsf is None:
                    rsf = [float(x) for x in toks[:-1]]
                elif toks and toks[-1] == "CorrectedRSF":
                    continue
                elif rsf is not None and len(toks) == len(regions) + 1:
                    try:
                        conc = [float(x) for x in toks[1:]]
                    except ValueError:
                        conc = None
                    break
            if rsf and conc and len(rsf) == len(regions):
                for reg, r, c in zip(regions, rsf, conc):
                    out.append({"origin": origin, "data_file": meta.get("data_file"), "comment": meta.get("comment"),
                                "region": reg, "core_level": phi_region_to_core(reg), "RSF": r, "at_percent": c,
                                "provenance": "provider_raw_fit", "source_detail": "provider atomic concentration table"})
    return pd.DataFrame(out)


def read_provider_report(data: bytes, name: str = "") -> tuple[pd.DataFrame, pd.DataFrame]:
    bands = read_multipak_provider(data, name)
    acts = []
    if isinstance(data, bytes) and data[:2] == b"PK":
        for slide, txt in pptx_slide_texts(data):
            if "Atomic Concentration Table" in txt:
                acts.append(parse_atomic_concentration(txt, origin=f"{name} slide {slide}"))
    else:
        text = data.decode("utf-8-sig", errors="replace") if isinstance(data, bytes) else data
        for block in re.split(r"(?=File Information)", text):
            if "Atomic Concentration Table" in block:
                acts.append(parse_atomic_concentration(block, origin=name))
    act = pd.concat(acts, ignore_index=True) if acts else pd.DataFrame(
        columns=["origin", "data_file", "comment", "region", "core_level", "RSF", "at_percent", "provenance", "source_detail"])
    return bands, act


# --- per-sample assembly -------------------------------------------------------------------------------------
def composition_from_comment(comment: str) -> str:
    return comment.split("_", 1)[1] if "_" in comment else comment


def phase_from_composition(comp: str) -> Optional[str]:
    c = comp.replace(" ", "")
    if re.search(r"AlC", c):
        return "MAX"
    if re.search(r"CT[xz]", c):
        return "MXene"
    return None


def unresolved_tag(core: str, be: float, comps) -> str:
    el = core.split()[0]
    if comps:
        lo = min(c.bounds[0] for c in comps)
        hi = max(c.bounds[1] for c in comps)
        if be < lo:
            return f"unresolved_low_BE_{el}"
        if be > hi:
            return f"unresolved_high_BE_{el}"
    return f"unresolved_{el}_{be:.2f}"


def known_observation(core: str, be: float, phase: Optional[str]) -> Optional[dict]:
    for o in KNOWN_OBSERVATIONS:
        if o["core_level"] == core and o["be_range"][0] <= be <= o["be_range"][1] and (o.get("phase") in (None, phase)):
            return o
    return None


@dataclass
class SampleResult:
    sample: str
    data_file: str
    meta: dict
    structure: Optional[str]
    profile: Optional[str]
    mappings: dict = field(default_factory=dict)        # core -> mapping table (provider rows)
    act: pd.DataFrame = field(default_factory=pd.DataFrame)
    quant_values: dict = field(default_factory=dict)
    quant_warnings: list = field(default_factory=list)
    provisional: pd.DataFrame = field(default_factory=pd.DataFrame)
    python: dict = field(default_factory=dict)          # core -> {"fit": FitResult, "ident": DataFrame}
    notes: list = field(default_factory=list)


def build_sample(data_file: str, bands: pd.DataFrame, act: pd.DataFrame, kb: KnowledgeBase, meta: Optional[dict] = None,
                 raw_regions: Optional[list] = None) -> SampleResult:
    b = bands[bands["data_file"] == data_file]
    comment = str(b["comment"].iloc[0]) if len(b) else data_file
    meta = dict(meta or {})
    notes = []
    if not meta.get("nominal_composition"):
        meta["nominal_composition"] = composition_from_comment(comment)
        notes.append("composition inferred from the provider comment")
    if not meta.get("phase"):
        meta["phase"] = phase_from_composition(meta["nominal_composition"]) or ""
        notes.append("phase inferred from the composition string")
    meta.setdefault("sample_id", comment)
    structure = meta.get("structure") or infer_structure(meta)
    profile = kb.structures[structure]["default_profile"] if structure in kb.structures else None
    res = SampleResult(sample=str(meta["sample_id"]), data_file=data_file, meta=meta, structure=structure, profile=profile,
                       notes=notes)
    models = kb.fit_models(profile) if profile else {}
    for sidx, one in b.groupby("summary_index", sort=False):
        core = one["core_level"].iloc[0]
        doublets, _ = pair_spin_orbit(one)
        comps = models.get(core, [])
        m = map_to_kb(doublets, comps).table if comps else doublets.assign(mapped_component=None, mapped_role=None,
                                                                           mapping_method=None, mapping_note="")
        m["provenance"] = np.where(m["pairing"] == "paired", "provider_paired_doublet", "provider_raw_fit")
        tags, status, notes_o = [], [], []
        for r in m.itertuples():
            if isinstance(r.mapped_component, str):
                tags.append(r.mapped_component); status.append(""); notes_o.append("")
                continue
            if core in UNSUPPORTED_REGIONS or not comps:
                tags.append(f"{core} (no validated KB model)"); status.append("not modelled"); notes_o.append(UNSUPPORTED_REGIONS.get(core, ""))
                continue
            o = known_observation(core, r.BE, meta.get("phase"))
            tags.append(o["label"] if o else unresolved_tag(core, r.BE, comps))
            status.append(o["identifiability"] if o else "not modelled")
            notes_o.append(o["note"] if o else "provider component outside every KB window")
        m["comparison_label"] = tags
        m["identifiability"] = status
        m["observation_note"] = notes_o
        res.mappings[core] = m
    a = act[act["data_file"] == data_file] if len(act) else act
    res.act = a
    if len(a) and profile:
        at = dict(zip(a["core_level"], a["at_percent"]))
        rows = []
        for core, m in res.mappings.items():
            if core not in at or at[core] <= 0:
                continue
            tot = m["area_total_estimated"].sum()
            for r in m.itertuples():
                role = r.mapped_role if isinstance(r.mapped_role, str) else ("termination_Cl" if core == "Cl 2p" else "other")
                rows.append({"sample": res.sample, "core_level": core, "component": r.comparison_label,
                             "area": at[core] * r.area_total_estimated / tot, "role": role})
        if rows and any(r["role"] == "framework_C" for r in rows):
            nom = metal_fractions(meta["nominal_composition"])
            extra = []
            for core, m in res.mappings.items():
                if (m["area_convention"] == "primary_line_only").any():
                    extra.append(f"{core}: provider reported only one spin-orbit line; component totals are estimated and the "
                                 "provider survey at.% may be underestimated")
            res.notes += extra
            spec = kb.structures.get(structure, {})
            opt = QuantOptions(m_target=float(spec.get("m_target", 2.0)), x_target=float(spec.get("x_target", 1.0)),
                               input_mode="atomic_percent", nominal_metal_stoich=nom if len(nom) > 1 else {},
                               nb_suboxide_validated=bool(kb.raw.get("calibration_status", {}).get("nb3d_suboxide_like_validated", False)),
                               phase=spec.get("phase", meta.get("phase") or "MXene"), a_target=spec.get("a_target"),
                               provisional_reasons=tuple(kb.raw.get("calibration_status", {}).get("delta_finalization_blocked_by", [])) + tuple(extra))
            q = quantify_sample(pd.DataFrame(rows), {}, opt)
            res.quant_values, res.quant_warnings = q.values, q.warnings
            res.provisional = provisional_report(q.values, opt)
            res.notes.append("quantification: provider at.% x provider (paired) area fractions, atomic_percent mode")
    for sp in raw_regions or []:
        core = sp.core_level
        if core in models and core not in UNSUPPORTED_REGIONS:
            fr = fit_core_level(sp.be, sp.intensity, models[core], spectrum_core_level=core)
            res.python[core] = {"fit": fr, "ident": assess(fr)}
    return res


# --- comparison tables ----------------------------------------------------------------------------------------
def long_table(samples: list[SampleResult]) -> pd.DataFrame:
    rows = []
    for s in samples:
        g = {"sample": s.sample, "phase": s.meta.get("phase"), "composition": s.meta.get("nominal_composition"),
             "treatment": " ".join(x for x in (str(s.meta.get("treatment") or ""), f"{s.meta.get('treatment_temperature_C')} C"
                                              if s.meta.get("treatment_temperature_C") else "") if x) or "none/unknown",
             "structure": s.structure, "profile": s.profile}
        for core, m in s.mappings.items():
            tot = m["area_total_estimated"].sum()
            for r in m.itertuples():
                base = {**g, "core_level": core, "component": r.comparison_label, "provider_label": r.provider_label,
                        "identifiability": r.identifiability}
                rows.append({**base, "quantity": "BE (eV)", "value": r.BE, "provenance": r.provenance})
                rows.append({**base, "quantity": "area fraction", "value": r.area_total_estimated / tot if tot else np.nan,
                             "provenance": "derived_quantity"})
        for core, py in s.python.items():
            for r in py["ident"].itertuples():
                base = {**g, "core_level": core, "component": r.component, "provider_label": None, "identifiability": r.status}
                rows.append({**base, "quantity": "BE (eV) [python]", "value": r.BE, "provenance": "python_fit"})
                rows.append({**base, "quantity": "area fraction [python]", "value": r.area_fraction, "provenance": "python_fit"})
        for r in s.act.itertuples() if len(s.act) else []:
            rows.append({**g, "core_level": r.core_level, "component": "(survey)", "provider_label": None, "identifiability": "",
                         "quantity": "at.%", "value": r.at_percent, "provenance": "provider_raw_fit"})
    return pd.DataFrame(rows)


def quant_summary(samples: list[SampleResult]) -> pd.DataFrame:
    rows = []
    for s in samples:
        base = {"sample": s.sample, "phase": s.meta.get("phase"), "composition": s.meta.get("nominal_composition")}
        if s.provisional is None or not len(s.provisional):
            rows.append({**base, "quantity": "(no quantification)", "nominal (default)": None,
                         "alternative (candidate as framework)": None, "provisional": None, "why provisional": "",
                         "provenance": ""})
            continue
        for r in s.provisional.to_dict("records"):
            rows.append({**base, **r})
    return pd.DataFrame(rows)


def wide(long: pd.DataFrame, quantity: str, core: Optional[str] = None) -> pd.DataFrame:
    d = long[long["quantity"] == quantity]
    if core:
        d = d[d["core_level"] == core]
    if d.empty:
        return pd.DataFrame()
    return d.pivot_table(index=["core_level", "component"], columns="sample", values="value", aggfunc="first")
