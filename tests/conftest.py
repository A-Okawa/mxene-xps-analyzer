import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mxene_xps.knowledge_base import default_kb  # noqa: E402

FIX = ROOT / "tests" / "fixtures"


@pytest.fixture(scope="session")
def kb():
    return default_kb()


@pytest.fixture(scope="session")
def truth():
    return json.loads((FIX / "synthetic_truth.json").read_text())
