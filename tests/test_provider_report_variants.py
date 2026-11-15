"""MultiPak export variants seen in real provider reports (synthetic numbers only)."""
from mxene_xps import dataset as ds
from mxene_xps.sample_metadata import infer_structure

ACT_BARE_HEADER = """
--------------------------
     File Information
--------------------------
File Name:  sample77.spe
File Type:  XPS SPECTRUM
Comment:  Ti3C2Tx

--------------------------
Atomic Concentration Table
--------------------------
     C1s        O1s        F1s        Ti2p
    0.314      0.733      1.000      2.077    RSF
   10.000     20.000     30.000     40.000    CorrectedRSF
    30.00      20.00      25.00      25.00
"""

ACT_ABSCISSA = """File Name:  sample78.spe
Comment:  2_Nb2CTx
Abscissa  C1s  O1s  Nb3d
  0.314  0.733  2.921  RSF
  12.5  30.0  57.5  12.5
"""


def test_act_without_abscissa_label():
    a = ds.parse_atomic_concentration(ACT_BARE_HEADER, "s")
    assert list(a["core_level"]) == ["C 1s", "O 1s", "F 1s", "Ti 2p"]
    assert list(a["at_percent"]) == [30.0, 20.0, 25.0, 25.0]
    assert list(a["RSF"]) == [0.314, 0.733, 1.0, 2.077]
    assert set(a["data_file"]) == {"sample77.spe"} and set(a["comment"]) == {"Ti3C2Tx"}


def test_act_with_abscissa_label_still_parsed():
    a = ds.parse_atomic_concentration(ACT_ABSCISSA, "s")
    assert list(a["core_level"]) == ["C 1s", "O 1s", "Nb 3d"] and list(a["at_percent"]) == [30.0, 57.5, 12.5]


def test_ti3c2tx_comment_is_an_mxene():
    assert ds.phase_from_composition("Ti3C2Tx") == "MXene"
    assert ds.phase_from_composition("Nb4C3Tz") == "MXene"
    assert ds.phase_from_composition("Nb2CTx") == "MXene"
    assert ds.phase_from_composition("Ti3AlC2") == "MAX"
    assert infer_structure({"nominal_composition": "Ti3C2Tx", "phase": "MXene"}) == "Ti3C2Tx"
