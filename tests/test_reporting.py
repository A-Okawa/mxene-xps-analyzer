import json

import pandas as pd

from mxene_xps.models import MeasurementMetadata, OrthogonalEvidence, SampleSetup
from mxene_xps.reporting import evidence_summary, html_report, metadata_warnings, project_dict


def test_metadata_warnings():
    w = metadata_warnings(MeasurementMetadata())
    assert any("unknown" in x for x in w)
    w = metadata_warnings(MeasurementMetadata(be_reference_method="adventitious C 1s", surface_state="after Ar+ sputtering"))
    assert any("Fermi-edge" in x for x in w) and any("sputtering" in x.lower() for x in w) and any("284.8" in x for x in w)
    w = metadata_warnings(MeasurementMetadata(be_reference_method="Fermi edge", xray_source="Al Ka", pass_energy_eV=20,
                                              step_eV=0.05, charge_neutralization="on"))
    assert w == []


def test_evidence_summary():
    _, s = evidence_summary([OrthogonalEvidence("ICP")])
    assert "XPS-only" in s
    lines, s = evidence_summary([OrthogonalEvidence("ICP", True, "Ti/Nb = 0.3", ["Ti depletion"])])
    assert "1 orthogonal" in s and "Ti depletion" in lines[0]


def test_project_and_html_export_serializable():
    setup = SampleSetup("s", "(Nb,Ti)2CTx", 2, 1, "nbti_m2c_defect_aware", ["Nb", "Ti"])
    res = pd.DataFrame([{"sample": "s", "M_norm": 1.6}])
    d = project_dict(setup, MeasurementMetadata(), [], "1.2.0", "nbti_m2c_defect_aware", {"C 1s": 0.3}, {"m_target": 2},
                     res, res, {"s": ["w"]}, None, None, ["Brette2024"])
    json.dumps(d)
    assert d["kb_version"] == "1.2.0" and d["quantification_results"][0]["M_norm"] == 1.6
    h = html_report("t <x>", [("a", "text"), ("b", ["w1"]), ("c", res)])
    assert "<table" in h and "t &lt;x&gt;" in h
