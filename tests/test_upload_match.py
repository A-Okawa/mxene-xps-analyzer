"""Matching uploaded raw exports to provider-report samples (Import & compare tab). Synthetic files only."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from upload_match import file_kind, match_stem, match_uploads  # noqa: E402

DEMO = ROOT / "examples" / "synthetic_demo"


class F:
    def __init__(self, name, data):
        self.name, self._d = name, data if isinstance(data, bytes) else data.encode()

    def getvalue(self):
        return self._d


def _curve(lo, hi):
    x = np.linspace(hi, lo, 120)
    return "\n".join(f"{a:.3f},{1000 + 50 * np.exp(-((a - (lo + hi) / 2) ** 2)):.2f}" for a in x)


def test_file_kinds():
    assert file_kind("okawa121_Nb2CTx_multi.csv") == "multi"
    assert file_kind("okawa121_Nb2CTx_Ti2pcurve.csv") == "curve"
    assert file_kind("okawa121_Nb2CTx_C1scurve.csv") == "curve"
    assert file_kind("okawa121_Nb2CTx_su.csv") == "survey"
    assert file_kind("okawa121_notes.csv") == "other"


def test_stem_is_the_provider_data_file_not_the_first_underscore():
    stems = ["okawa121", "okawa1211", "my_sample_3"]
    assert match_stem("okawa1211_x_multi.csv", stems) == "okawa1211"
    assert match_stem("okawa121_x_multi.csv", stems) == "okawa121"
    assert match_stem("my_sample_3_Ti2pcurve.csv", stems) == "my_sample_3"     # stem with underscores
    assert match_stem("other_x.csv", stems) is None


def test_multi_curve_and_survey_files_of_one_sample_do_not_overwrite_each_other():
    multi = (DEMO / "SYNdemo203_3_SYNTHETIC_multi.csv").read_bytes()
    files = [F("SYNdemo203_3_SYNTHETIC_su.csv", multi),                    # survey: never fitted
             F("SYNdemo203_3_SYNTHETIC_multi.csv", multi),
             F("SYNdemo203_3_SYNTHETIC_Ti2pcurve.csv", _curve(452, 468)),  # superseded by the multi Ti 2p
             F("SYNdemo203_3_SYNTHETIC_Cl2pcurve.csv", _curve(195, 205)),  # core level not in the multi file -> used
             F("unrelated_file.csv", "1,2\n3,4\n")]
    m = match_uploads(files, ["SYNdemo203.spe", "SYNdemo201.spe"])
    st = {r["file"]: r["status"] for r in m.rows}
    assert st == {"SYNdemo203_3_SYNTHETIC_su.csv": "survey", "SYNdemo203_3_SYNTHETIC_multi.csv": "used",
                  "SYNdemo203_3_SYNTHETIC_Ti2pcurve.csv": "superseded", "SYNdemo203_3_SYNTHETIC_Cl2pcurve.csv": "used",
                  "unrelated_file.csv": "no_sample"}
    cores = {sp.core_level: sp for sp in m.raw_regions["SYNdemo203"]}
    assert {"C 1s", "O 1s", "Nb 3d", "Ti 2p", "F 1s", "Cl 2p"} <= set(cores)
    assert len(cores["Ti 2p"].be) == 161 and cores["Cl 2p"].region == "Cl2p"   # Ti 2p from the multi file (161 pts), not the curve (120)
    assert m.raw_regions["SYNdemo201"] == []


def test_curve_only_upload_is_used_for_refit():
    m = match_uploads([F("SYNdemo201_1_C1scurve.csv", _curve(280, 292))], ["SYNdemo201.spe"])
    assert [r["status"] for r in m.rows] == ["used"]
    assert [sp.core_level for sp in m.raw_regions["SYNdemo201"]] == ["C 1s"]


def test_upload_order_does_not_matter():
    multi = (DEMO / "SYNdemo201_1_SYNTHETIC_multi.csv").read_bytes()
    a = [F("SYNdemo201_1_multi.csv", multi), F("SYNdemo201_1_su.csv", multi), F("SYNdemo201_1_C1scurve.csv", _curve(280, 292))]
    r1 = {sp.core_level: len(sp.be) for sp in match_uploads(a, ["SYNdemo201.spe"]).raw_regions["SYNdemo201"]}
    r2 = {sp.core_level: len(sp.be) for sp in match_uploads(a[::-1], ["SYNdemo201.spe"]).raw_regions["SYNdemo201"]}
    assert r1 == r2 and r1["C 1s"] == 161
