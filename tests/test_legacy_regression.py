"""Legacy Excel regression (§10.4, §29.8).

The historical workbook is private and NOT stored in the repository. Point the
test at a local copy with:   export XPS_LEGACY_WORKBOOK=/path/to/mxene_surface_functional_group_quant_template_v2.xlsx
(or place it at private/legacy_workbook.xlsx, which is git-ignored). Without it, the test is skipped.
"""
import math
import os
from pathlib import Path

import pytest

from mxene_xps.quantification import LEGACY_OUTPUT_COLUMNS, legacy_quantify

ROOT = Path(__file__).resolve().parents[1]
WB = Path(os.environ.get("XPS_LEGACY_WORKBOOK", ROOT / "private" / "legacy_workbook.xlsx"))
TOL = 1e-10

# Peak_Areas rows of the workbook -> legacy keys
ROWS = {"C_TiTx": 5, "Ti1": 11, "Ti2": 12, "Ti3": 13, "Ti_oxide": 14, "O_low": 17, "O_oxide": 18, "O_high": 19,
        "F0": 23, "F1": 24, "F2": 25, "Cl": 28}


def load_workbook_cases():
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.load_workbook(WB, data_only=True)  # cached values computed by Excel
    ctl = wb["Controls_RSF"]
    controls = {ctl[f"A{r}"].value: ctl[f"B{r}"].value for r in range(4, 9)}
    rsf = {ctl[f"A{r}"].value: ctl[f"B{r}"].value for r in range(13, 18)}
    pa, res = wb["Peak_Areas"], wb["Results"]
    header = [c.value for c in res[3]]
    cases = []
    for j, col in enumerate("DEFGHI"):
        areas = {k: pa[f"{col}{r}"].value for k, r in ROWS.items()}
        expected = dict(zip(header, [c.value for c in res[4 + j]]))
        cases.append((expected["Sample"], areas, expected))
    return controls, rsf, cases


@pytest.mark.skipif(not WB.exists(), reason="private legacy workbook not available (set XPS_LEGACY_WORKBOOK)")
def test_legacy_workbook_regression():
    controls, rsf, cases = load_workbook_cases()
    assert len(cases) == 6
    n, max_err = 0, 0.0
    for sample, areas, expected in cases:
        got = legacy_quantify(areas, rsf, controls)
        for k in LEGACY_OUTPUT_COLUMNS:
            e, g = expected[k], got[k]
            if isinstance(e, str):
                assert g == e, (sample, k)
            else:
                err = abs(g - e)
                max_err = max(max_err, err)
                assert err <= TOL, (sample, k, g, e)
            n += 1
    print(f"legacy regression: {n} values compared over {len(cases)} samples, max |error| = {max_err:.3e}")
    assert not math.isnan(max_err)
