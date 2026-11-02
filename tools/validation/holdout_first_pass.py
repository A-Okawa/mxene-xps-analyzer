"""Frozen first-pass validation run (calibration or hold-out set). Does NOT change the KB or any model.

Usage:
  python -I tools/validation/holdout_first_pass.py --data DIR --label holdout --out OUT_DIR [--metadata META.csv]

DIR must contain one provider report (MultiPak .pptx or Curve-Fit-Summary text) and the raw exports
(*_multi.csv; optional *curve.csv). OUT_DIR must be outside the repository (private data).

Frozen-state guard: refuses to run unless app == 1.6.0, KB == 1.5.0 and data/knowledge_base.json is byte-identical to
the pre-validation baseline commit f466b94 (and src/ has no uncommitted changes).

Per sample and core level it records: import status + messages, fit success, RMSE, provider-vs-Python QC, fitted
BE/FWHM/area fraction, identifiability, bound hits, strongest correlations, sensitivity spans (model-sensitive),
provisional/alternative quantities, unresolved provider components, unsupported regions, and every message with its
severity. Each issue gets an AUTO-SUGGESTED category (A software/import, B numerical, C expected limitation of the frozen
model, D potentially new chemistry / unsupported component) that must be reviewed by a human before acting on it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from mxene_xps import __version__  # noqa: E402
from mxene_xps import dataset as ds  # noqa: E402
from mxene_xps import sample_metadata as sm  # noqa: E402
from mxene_xps import sensitivity as sens  # noqa: E402
from mxene_xps.calibration import KNOWN_OBSERVATIONS, UNSUPPORTED_REGIONS  # noqa: E402
from mxene_xps.fitting import fit_core_level  # noqa: E402
from mxene_xps.identifiability import assess, warnings_from  # noqa: E402
from mxene_xps.io import parse_phi_multiregion, parse_two_column  # noqa: E402
from mxene_xps.knowledge_base import default_kb  # noqa: E402
from mxene_xps.multipak import provider_component_table  # noqa: E402
from mxene_xps.quantification import QuantOptions, quantify_sample  # noqa: E402
from mxene_xps.severity import QC_STATUS_SEVERITY, classify  # noqa: E402
from mxene_xps.validation import compare_components, compare_grouped, prepare_vendor, summarize  # noqa: E402

FROZEN = {"app": "1.6.0", "kb": "1.5.0", "commit": "f466b94"}
# Known, already documented patterns (CALIBRATION_LOG / README limitations) -> category C when reproduced.
KNOWN_PAIRS = ({"Nb_MXene", "Nb_suboxide_like"}, {"Nb_MAX", "Nb_MAX_2nd"}, {"Nb_suboxide_like", "Nb2O5"})


def _git(*args):
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True).stdout.strip()


def frozen_guard(kb) -> dict:
    head = _git("rev-parse", "--short", "HEAD")
    kb_diff = _git("diff", FROZEN["commit"], "--", "data/knowledge_base.json")
    src_dirty = _git("status", "--porcelain", "--", "src", "data")
    ok = __version__ == FROZEN["app"] and kb.version == FROZEN["kb"] and not kb_diff and not src_dirty
    if not ok:
        raise SystemExit(f"Frozen-state guard failed: app {__version__}, KB {kb.version}, KB diff vs {FROZEN['commit']}: "
                         f"{bool(kb_diff)}, uncommitted src/data changes: {bool(src_dirty)}")
    return {"head": head, "baseline": FROZEN["commit"], "app_version": __version__, "kb_version": kb.version,
            "kb_sha256": hashlib.sha256((ROOT / "data" / "knowledge_base.json").read_bytes()).hexdigest()}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def msg(text, source, level=None):
    return {"severity": level or classify(text), "source": source, "text": str(text)}


def issue(sample, core, category, text, evidence=""):
    return {"sample": sample, "core_level": core, "category_auto": category, "text": text, "evidence": evidence,
            "review_status": "auto-suggested; requires human review"}


def run(data: Path, label: str, meta_path: Path | None) -> dict:
    kb = default_kb()
    guard = frozen_guard(kb)
    reports = sorted(list(data.glob("*.pptx")) + [p for p in data.glob("*.txt") if "Curve Fit Summary" in p.read_text(errors="ignore")])
    out = {"label": label, "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "frozen": guard,
           "inputs": {p.name: sha(p) for p in sorted(data.iterdir()) if p.is_file()}, "samples": [], "issues": [],
           "global_messages": []}
    if len(reports) != 1:
        out["issues"].append(issue("*", "*", "A", f"expected exactly one provider report, found {len(reports)}"))
        return out
    try:
        bands, act = ds.read_provider_report(reports[0].read_bytes(), reports[0].name)
    except Exception as e:
        out["issues"].append(issue("*", "*", "A", f"provider report import failed: {e}"))
        return out
    meta = sm.read(meta_path.read_bytes()) if meta_path else None
    for f in dict.fromkeys(bands["data_file"]):
        stem = str(f).rsplit(".", 1)[0]
        comment = str(bands.loc[bands["data_file"] == f, "comment"].iloc[0])
        rec = {"data_file": f, "comment": comment, "regions": [], "messages": []}
        m_row = None
        if meta is not None:
            hit = meta[(meta["provider_file"] == stem) | (meta["sample_id"] == comment)]
            m_row = hit.iloc[0].to_dict() if len(hit) else None
        # --- raw imports (existing importers only)
        raw = {}
        for mf in sorted(data.glob(f"{stem}_*multi.csv")):
            try:
                for r in parse_phi_multiregion(mf.read_bytes()):
                    raw[r.core_level] = (r, mf.name)
                    for n in r.notes:
                        rec["messages"].append(msg(n, f"import {mf.name}:{r.region}"))
            except Exception as e:
                rec["messages"].append(msg(f"import failed {mf.name}: {e}", "import", "INVALID"))
                out["issues"].append(issue(comment, "*", "A", f"multi-region import failed: {e}", mf.name))
        curves = {}
        for cf in sorted(data.glob(f"{stem}_*curve.csv")):
            try:
                sp = parse_two_column(cf.read_bytes())
                reg = cf.name[len(stem) + 1:].rsplit("curve", 1)[0].split("_")[-1]
                from mxene_xps.io import phi_region_to_core
                core = phi_region_to_core(reg)
                if core:
                    curves[core] = (sp, cf.name)
            except Exception as e:
                rec["messages"].append(msg(f"import failed {cf.name}: {e}", "import", "INVALID"))
                out["issues"].append(issue(comment, "*", "A", f"curve CSV import failed: {e}", cf.name))
        try:
            s = ds.build_sample(f, bands, act, kb, meta=m_row)
        except Exception as e:
            out["issues"].append(issue(comment, "*", "A", f"sample assembly failed: {e}", traceback.format_exc(limit=2)))
            out["samples"].append(rec)
            continue
        rec.update({"sample": s.sample, "structure": s.structure, "profile": s.profile, "meta": s.meta, "notes": s.notes})
        rec["messages"] += [msg(n, "dataset") for n in s.notes]
        rec["messages"] += [msg(w, "quantification") for w in s.quant_warnings]
        rec["quant_values"] = {k: v for k, v in s.quant_values.items() if isinstance(v, (int, float, str))}
        rec["provisional"] = json.loads(s.provisional.to_json(orient="records")) if len(s.provisional) else []
        rec["survey_at_percent"] = json.loads(s.act[["core_level", "at_percent"]].to_json(orient="records")) if len(s.act) else []
        if s.structure is None:
            out["issues"].append(issue(s.sample, "*", "A", "structure could not be inferred (metadata needed)"))
        models = kb.fit_models(s.profile) if s.profile else {}
        # provider-based quant table for derived-quantity sensitivity
        qrows = []
        at = dict(zip(s.act["core_level"], s.act["at_percent"])) if len(s.act) else {}
        for core, mp in s.mappings.items():
            if at.get(core, 0) > 0:
                tot = mp["area_total_estimated"].sum()
                for r in mp.itertuples():
                    role = r.mapped_role if isinstance(r.mapped_role, str) else ("termination_Cl" if core == "Cl 2p" else "other")
                    qrows.append({"sample": s.sample, "core_level": core, "component": r.comparison_label,
                                  "area": at[core] * r.area_total_estimated / tot, "role": role})
        qtab = pd.DataFrame(qrows)
        spec = kb.structures.get(s.structure, {})
        qopt = QuantOptions(m_target=float(spec.get("m_target", 2)), x_target=float(spec.get("x_target", 1)), input_mode="atomic_percent",
                            phase=spec.get("phase", "MXene"), nominal_metal_stoich={})
        for core in sorted(set(s.mappings) | set(raw) | set(curves)):
            reg = {"core_level": core, "supported": core in models, "fits": []}
            mp = s.mappings.get(core)
            if mp is not None:
                reg["provider"] = json.loads(mp[["provider_label", "pairing", "area_convention", "BE", "FWHM", "split_obs", "ratio_obs",
                                                 "area_total_estimated", "mapped_component", "mapped_role", "comparison_label",
                                                 "identifiability", "provenance"]].to_json(orient="records"))
                for r in mp.itertuples():
                    if isinstance(r.mapped_component, str):
                        continue
                    known = any(o["label"] == r.comparison_label for o in KNOWN_OBSERVATIONS)
                    if core in UNSUPPORTED_REGIONS or core not in models:
                        continue
                    out["issues"].append(issue(s.sample, core, "C" if known else "D",
                                               f"provider component outside every KB window: {r.comparison_label} @ {r.BE:.2f} eV "
                                               f"({r.provider_label})", "known calibration observation" if known else "not seen in KB windows"))
                if (mp["area_convention"] == "primary_line_only").any():
                    out["issues"].append(issue(s.sample, core, "C", "provider reported only one spin-orbit line (totals estimated)"))
            if core in UNSUPPORTED_REGIONS or core not in models:
                reg["status"] = "unsupported (no validated KB model)"
                rec["messages"].append(msg(f"{core}: no validated KB model", "region", "CAUTION"))
                out["issues"].append(issue(s.sample, core, "C", f"{core}: no validated KB model (unsupported region)"))
                rec["regions"].append(reg)
                continue
            comps = models[core]
            for window, src in (("raw_full_region", raw.get(core)), ("provider_curve_window", curves.get(core))):
                if src is None:
                    continue
                sp, fname = src
                fitrec = {"window": window, "file": fname, "BE_range": [float(sp.be.min()), float(sp.be.max())]}
                try:
                    fr = fit_core_level(sp.be, sp.intensity, comps, spectrum_core_level=core)
                except Exception as e:
                    fitrec.update({"fit_success": False, "error": str(e)})
                    rec["messages"].append(msg(f"Fit failed ({core}, {window}): {e}", "fit", "INVALID"))
                    out["issues"].append(issue(s.sample, core, "A" if "refus" in str(e) or "wide" in str(e) else "B",
                                               f"fit failed ({window}): {e}"))
                    reg["fits"].append(fitrec)
                    continue
                try:
                    derived_fn = sens.derived_from_fractions(qtab, s.sample, core, lambda t: quantify_sample(t, {}, qopt).values) \
                        if len(qtab) and (qtab["core_level"] == core).any() else None
                    sres = sens.run(sp.be, sp.intensity, comps, sens.standard_alternatives(core, comps, sp.be),
                                    derived_fn=derived_fn, spectrum_core_level=core)
                    spans = sres.spans
                except Exception as e:
                    sres, spans = None, pd.DataFrame()
                    out["issues"].append(issue(s.sample, core, "B", f"sensitivity run failed ({window}): {e}"))
                ident = assess(fr, alt_spans=sres.component_spans() if sres is not None else None)
                fitrec.update({
                    "fit_success": bool(fr.metrics["success"]), "RMSE": fr.metrics["RMSE"], "R2": fr.metrics["R2"],
                    "AIC": fr.metrics["AIC"], "BIC": fr.metrics["BIC"],
                    "components": json.loads(fr.components[["component", "BE", "FWHM", "area_fraction", "area_rel_stderr",
                                                            "at_bound", "role"]].to_json(orient="records")),
                    "identifiability": json.loads(ident[["component", "status", "status_reasons", "near_bound", "at_bound",
                                                         "strongest_area_corr_with", "strongest_area_r", "strongest_param_corr",
                                                         "strongest_param_r", "strongly_anticorrelated_with"]].to_json(orient="records")),
                    "model_sensitive": json.loads(spans[spans["model_sensitive"]].to_json(orient="records")) if len(spans) else [],
                    "sensitivity_spans": json.loads(spans.to_json(orient="records")) if len(spans) else [],
                })
                for d in fr.diagnostics:
                    rec["messages"].append(msg(d, f"fit {core} {window}"))
                for w in warnings_from(ident):
                    rec["messages"].append(msg(w, f"identifiability {core} {window}"))
                if mp is not None:
                    vp = prepare_vendor(provider_component_table(mp), comps, "total")
                    cqc = compare_components(vp, fr.components, kb.validation_defaults)
                    gqc = compare_grouped(vp, fr.components, kb.validation_defaults)
                    summ = summarize(cqc, gqc)
                    fitrec["qc"] = {"status": summ.status, "message": summ.message,
                                    "grouped": json.loads(gqc.to_json(orient="records")),
                                    "components": json.loads(cqc[["component", "presence", "delta_BE_eV", "delta_FWHM_rel",
                                                                  "delta_area_fraction_abs", "QC"]].to_json(orient="records"))}
                    rec["messages"].append(msg(f"provider-vs-Python QC {summ.status} ({core}, {window}): {summ.message}", "qc",
                                               QC_STATUS_SEVERITY.get(summ.status, "CAUTION")))
                # auto-classification of fit-level issues
                if not fr.metrics["success"]:
                    out["issues"].append(issue(s.sample, core, "B", f"optimizer did not report success ({window})"))
                # Statistical non-identifiability / boundary hits / model sensitivity are properties of the frozen model on
                # these data -> C (expected limitation).  B is reserved for numerical failures (non-convergence,
                # non-finite results, covariance or sensitivity failures).  'documented' marks patterns already recorded.
                for r in ident.itertuples():
                    if r.status == "identifiable":
                        continue
                    pair = {r.component, r.strongest_area_corr_with} if r.strongest_area_corr_with else set()
                    documented = pair in KNOWN_PAIRS or (core == "Ti 2p" and "Nb 3d" in models)
                    out["issues"].append(issue(s.sample, core, "C", f"{r.component}: {r.status} ({window})",
                                               r.status_reasons + ("; documented pattern" if documented else "")))
                    if not np.isfinite(r.BE) or not np.isfinite(r.area_fraction):
                        out["issues"].append(issue(s.sample, core, "B", f"{r.component}: non-finite fit result ({window})"))
                if any("covariance could not be estimated" in d.lower() for d in fr.diagnostics):
                    out["issues"].append(issue(s.sample, core, "B", f"parameter covariance not estimable ({window})"))
                for r in (spans[spans["model_sensitive"]].itertuples() if len(spans) else []):
                    out["issues"].append(issue(s.sample, core, "C", f"model-sensitive {r.quantity} (span {r.span:.3f}, {window})",
                                               "D4 documented" if core == "Nb 3d" else "background/window sensitivity"))
                reg["fits"].append(fitrec)
            rec["regions"].append(reg)
        out["samples"].append(rec)
    return out


def summarize_run(r: dict) -> dict:
    fits = [f for s in r["samples"] for g in s.get("regions", []) for f in g.get("fits", [])]
    ok = [f for f in fits if f.get("fit_success")]
    idents = [i["status"] for f in ok for i in f["identifiability"]]
    msgs = [m["severity"] for s in r["samples"] for m in s.get("messages", [])]
    qc = [f["qc"]["status"] for f in ok if f.get("qc")]
    return {
        "label": r["label"], "n_samples": len(r["samples"]),
        "n_regions": sum(len(s.get("regions", [])) for s in r["samples"]),
        "n_unsupported_regions": sum(1 for s in r["samples"] for g in s.get("regions", []) if not g["supported"]),
        "import_failures": sum(1 for i in r["issues"] if i["category_auto"] == "A"),
        "n_fits": len(fits), "fit_success_rate": len(ok) / len(fits) if fits else None,
        "RMSE": {k: float(v) for k, v in pd.Series([f["RMSE"] for f in ok]).describe().items()} if ok else {},
        "RMSE_by_window": {w: float(np.median([f["RMSE"] for f in ok if f["window"] == w])) for w in {f["window"] for f in ok}},
        "qc_status_counts": pd.Series(qc).value_counts().to_dict(),
        "identifiability_counts": pd.Series(idents).value_counts().to_dict(),
        "boundary_hits": sum(1 for f in ok for c in f["components"] if c["at_bound"]),
        "unresolved_provider_components": sum(1 for s in r["samples"] for g in s.get("regions", []) for p in g.get("provider", [])
                                              if not p["mapped_component"] and g["supported"]),
        "provisional_quantities": sum(1 for s in r["samples"] for p in s.get("provisional", []) if p.get("provisional")),
        "model_sensitive_quantities": sum(len(f["model_sensitive"]) for f in ok),
        "message_severity_counts": pd.Series(msgs).value_counts().to_dict(),
        "issue_category_counts": pd.Series([i["category_auto"] for i in r["issues"]]).value_counts().to_dict(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--metadata")
    a = ap.parse_args()
    out = Path(a.out).expanduser().resolve()
    private = (ROOT / "private").resolve()
    if (ROOT in out.parents or out == ROOT) and not (out == private or private in out.parents):
        raise SystemExit("Write outputs outside the repository or under the git-ignored private/ folder (private data).")
    r = run(Path(a.data).expanduser(), a.label, Path(a.metadata) if a.metadata else None)
    r["summary"] = summarize_run(r)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"first_pass_{a.label}.json").write_text(json.dumps(r, indent=1, default=float))
    pd.DataFrame(r["issues"]).to_csv(out / f"issues_{a.label}.csv", index=False)
    print(json.dumps(r["summary"], indent=1, default=float))


if __name__ == "__main__":
    main()
