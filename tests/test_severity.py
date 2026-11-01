from mxene_xps.severity import LEVELS, QC_STATUS_SEVERITY, classify, records


def test_levels_and_examples():
    assert LEVELS == ("INFO", "CAUTION", "WARNING", "INVALID")
    cases = {
        "BE range -9.60–1200.00 eV (span 1209.6 eV) is implausibly wide for a single-region fit": "INVALID",
        "Spectrum is Cl 2p but the fit model is for ['Nb 3d']; refusing to fit": "INVALID",
        "Areas of Nb_MXene and Nb_suboxide_like are strongly correlated (r = -0.95)": "WARNING",
        "Nb_suboxide_like: BE at the edge of its allowed window": "WARNING",
        "Nb_MXene: non-identifiable (|r| 0.98 with Nb_suboxide_like)": "WARNING",
        "model-sensitive: apparent_M_deficiency_delta": "WARNING",
        "apparent_M_deficiency_delta is PROVISIONAL: the Nb 3d candidate is not validated": "CAUTION",
        "Calibration-candidate role(s) present (metal_suboxide_candidate)": "CAUTION",
        "Nb 3s (roughly 466-469 eV) overlaps the Ti 2p1/2 side": "CAUTION",
        "Cl 2p: no validated KB model": "CAUTION",
        "Cropped to 202–213 eV for fitting.": "INFO",
        "composition inferred from the provider comment": "INFO",
    }
    for text, lv in cases.items():
        assert classify(text) == lv, text
    assert QC_STATUS_SEVERITY["FAIL"] == "WARNING" and QC_STATUS_SEVERITY["PASS"] == "INFO"
    assert records(["x"]) == [{"severity": "INFO", "text": "x"}]


def test_successful_fit_with_uncertain_chemistry_is_not_an_error():
    assert classify("Calibration candidate(s) in this model (fitted, not validated)") == "CAUTION"
