"""Compare a frozen first-pass hold-out run against the calibration-set run (no model change).

Usage:
  python -I tools/validation/compare_sets.py CALIBRATION_JSON HOLDOUT_JSON OUT_DIR

Writes comparison.json and comparison.md (OUT_DIR outside the repo or under private/).
Patterns are compared by issue *signature* (core level + category + normalized description, sample-independent):
  reproduced  = signature present in both sets
  new         = present only in the hold-out set
  not seen    = present only in the calibration set
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def signature(i: dict) -> str:
    t = re.sub(r"\(.*?\)", "", i["text"])                  # drop window / numbers in brackets
    t = re.sub(r"@ [0-9.]+ eV", "@ BE", t)
    t = re.sub(r"span [0-9.]+", "span", t)
    t = re.sub(r"\s+", " ", t).strip()
    return f"{i['core_level']} | {i['category_auto']} | {t}"


def compare(cal: dict, hold: dict) -> dict:
    sc, sh = cal["summary"], hold["summary"]
    sig_c = {signature(i) for i in cal["issues"]}
    sig_h = {signature(i) for i in hold["issues"]}
    new = sorted(sig_h - sig_c)
    new_d = [s for s in new if " | D | " in s]
    new_ab = [s for s in new if " | A | " in s or " | B | " in s]
    metrics = ["n_samples", "n_regions", "n_unsupported_regions", "import_failures", "n_fits", "fit_success_rate",
               "boundary_hits", "unresolved_provider_components", "provisional_quantities", "model_sensitive_quantities"]
    table = [{"metric": m, "calibration": sc.get(m), "holdout": sh.get(m)} for m in metrics]
    for key in ("RMSE", "qc_status_counts", "identifiability_counts", "message_severity_counts", "issue_category_counts",
                "RMSE_by_window"):
        for k in sorted(set(sc.get(key, {})) | set(sh.get(key, {}))):
            table.append({"metric": f"{key}.{k}", "calibration": sc.get(key, {}).get(k), "holdout": sh.get(key, {}).get(k)})
    fail_criteria = {
        "import_failures_in_holdout": sh.get("import_failures", 0) > 0,
        "fit_success_rate_below_1": (sh.get("fit_success_rate") or 0) < 1.0,
        "new_A_or_B_patterns": bool(new_ab),
    }
    return {
        "calibration_label": cal["label"], "holdout_label": hold["label"],
        "frozen": {"calibration": cal["frozen"], "holdout": hold["frozen"]},
        "same_frozen_state": cal["frozen"]["kb_sha256"] == hold["frozen"]["kb_sha256"]
        and cal["frozen"]["app_version"] == hold["frozen"]["app_version"],
        "metrics": table,
        "reproduced_patterns": sorted(sig_c & sig_h), "new_patterns": new, "not_seen_in_holdout": sorted(sig_c - sig_h),
        "new_D_patterns": new_d, "new_A_or_B_patterns": new_ab,
        "software_pass": not any(fail_criteria.values()), "fail_criteria": fail_criteria,
        "kb_revision_possibly_needed": bool(new_d),
        "note": ("software_pass concerns import/fit robustness only. Scientific conclusions are unchanged by construction "
                 "(frozen KB); new D patterns are candidates for review, not KB changes."),
    }


def to_md(c: dict) -> str:
    L = [f"# Hold-out comparison: {c['calibration_label']} vs {c['holdout_label']}", "",
         f"Same frozen state (app/KB sha): **{c['same_frozen_state']}**  ", f"Software pass: **{c['software_pass']}** {c['fail_criteria']}  ",
         f"KB revision possibly needed (new D patterns): **{c['kb_revision_possibly_needed']}**", "", c["note"], "",
         "## Metrics", "", "| metric | calibration | holdout |", "|---|---|---|"]
    L += [f"| {r['metric']} | {r['calibration']} | {r['holdout']} |" for r in c["metrics"]] + [""]
    for key, title in (("reproduced_patterns", "Reproduced patterns"), ("new_patterns", "New patterns (hold-out only)"),
                       ("not_seen_in_holdout", "Calibration patterns not seen in the hold-out set")):
        L += [f"## {title} ({len(c[key])})", ""] + [f"- {s}" for s in c[key]] + [""]
    return "\n".join(L)


def main():
    cal = json.loads(Path(sys.argv[1]).read_text())
    hold = json.loads(Path(sys.argv[2]).read_text())
    out = Path(sys.argv[3]).expanduser().resolve()
    private = (ROOT / "private").resolve()
    if (ROOT in out.parents or out == ROOT) and not (out == private or private in out.parents):
        raise SystemExit("Write outputs outside the repository or under private/.")
    c = compare(cal, hold)
    out.mkdir(parents=True, exist_ok=True)
    (out / "comparison.json").write_text(json.dumps(c, indent=1, default=float))
    (out / "comparison.md").write_text(to_md(c))
    print(json.dumps({k: c[k] for k in ("same_frozen_state", "software_pass", "fail_criteria", "kb_revision_possibly_needed")},
                     indent=1), "\nreproduced", len(c["reproduced_patterns"]), "new", len(c["new_patterns"]))


if __name__ == "__main__":
    main()
