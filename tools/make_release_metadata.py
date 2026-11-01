"""Write release/release_metadata_v<version>.json (machine-readable reproducibility record).

Usage:
  python -I tools/make_release_metadata.py [--e2e KEY=RESULT ...] [--numerical-outputs-unchanged yes|no|not-checked]

Runs the test suite without private data (environment variables removed) and, if XPS_SASAKI_DIR is set, again with
private data. Records app/KB version, git commit (the commit the file was generated from; the release commit that
adds this file is its child), Python and package versions, test counts, date and the known scientific open issues.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mxene_xps.knowledge_base import default_kb  # noqa: E402
from mxene_xps.release import release_info  # noqa: E402
from mxene_xps.reporting import git_commit  # noqa: E402


def run_pytest(env: dict) -> dict:
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"], cwd=ROOT, env=env,
                       capture_output=True, text=True)
    last = [l for l in p.stdout.splitlines() if re.search(r"\d+ (passed|failed)", l)]
    line = last[-1] if last else p.stdout[-300:]
    counts = {k: int(v) for v, k in re.findall(r"(\d+) (passed|failed|skipped|errors?)", line)}
    return {"summary": line.strip("= ").strip(), "returncode": p.returncode, **counts}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--e2e", nargs="*", default=[])
    ap.add_argument("--numerical-outputs-unchanged", default="not-checked")
    a = ap.parse_args()
    kb = default_kb()
    info = release_info(kb, commit=git_commit(ROOT))
    clean = {k: v for k, v in os.environ.items() if k not in ("XPS_SASAKI_DIR", "XPS_LEGACY_WORKBOOK")}
    tests = {"without_private_data": run_pytest(clean)}
    private = bool(os.environ.get("XPS_SASAKI_DIR"))
    if private:
        tests["with_private_data"] = run_pytest(dict(os.environ))
    rec = {
        "release": f"v{info['app_version']}-rc", "app_version": info["app_version"], "kb_version": info["kb_version"],
        "generated_from_commit": info["git_commit"],
        "commit_note": "The release commit that adds this file is the child of generated_from_commit.",
        "python": info["python"], "platform": info["platform"], "packages": info["packages"],
        "date": info["generated_at"], "tests": tests,
        "private_data_regression_tests_run": private,
        "legacy_workbook_test_run": bool(os.environ.get("XPS_LEGACY_WORKBOOK")),
        "e2e": dict(x.split("=", 1) for x in a.e2e),
        "numerical_outputs_unchanged_vs_kb150_baseline": a.numerical_outputs_unchanged,
        "frozen_scientific_state": kb.raw.get("calibration_status", {}).get("delta_finalization_blocked_by", []),
        "scientific_open_issues": info["scientific_open_issues"],
        "interpretation_levels": info["interpretation_levels"],
    }
    out = ROOT / "release" / f"release_metadata_v{info['app_version']}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(rec, indent=2) + "\n")
    print(out, json.dumps(tests))


if __name__ == "__main__":
    main()
