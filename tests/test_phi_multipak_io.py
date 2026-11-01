"""Calibration phase 1 (I/O only): PHI multi-region split, single-region range guard, MultiPak provider import.

Synthetic text only; the real provider data are private.  Set XPS_SASAKI_DIR to the private provider-data folder to
run the optional real-data checks.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mxene_xps.fitting import fit_core_level
from mxene_xps.io import (check_single_region, crop, parse_phi_multiregion, parse_two_column,
                          phi_region_to_core)
from mxene_xps.multipak import (apply_user_mapping, map_to_kb, pair_spin_orbit, parse_curve_fit_summaries,
                                provider_component_table, read_multipak_provider)


def _phi_block(region, lo, hi, n=81, peak=None):
    x = np.linspace(hi, lo, n)                       # PHI exports descending BE
    c = peak if peak is not None else (lo + hi) / 2
    y = 1000 + 5000 * np.exp(-0.5 * ((x - c) / 0.6) ** 2)
    return "Full\r\n \r\n" + f"{region}\r\n1\r\n" + "".join(f"{a:.4f},{b:.4f}\r\n" for a, b in zip(x, y))


MULTI = _phi_block("C1s", 280, 300, peak=284.8) + _phi_block("O1s", 525, 545, peak=531) + \
    _phi_block("Nb3d", 199, 219, peak=204.5)


def test_region_names():
    assert phi_region_to_core("C1s") == "C 1s"
    assert phi_region_to_core("Ti2p") == "Ti 2p"
    assert phi_region_to_core("Ti2p3") == "Ti 2p"
    assert phi_region_to_core("Nb3d") == "Nb 3d"
    assert phi_region_to_core("Sur1") is None


def test_generic_importer_refuses_multi_region_file():
    with pytest.raises(ValueError, match="multi-region export"):
        parse_two_column(MULTI)


def test_phi_multiregion_split():
    regs = parse_phi_multiregion(MULTI.encode())
    assert [r.region for r in regs] == ["C1s", "O1s", "Nb3d"]
    assert [r.core_level for r in regs] == ["C 1s", "O 1s", "Nb 3d"]
    for r in regs:
        assert len(r.be) == 81 and r.be.max() - r.be.min() == pytest.approx(20.0)
        assert (np.diff(r.be) > 0).all() and r.was_descending and r.n_segments == 1


def test_phi_survey_is_single_non_core_region():
    regs = parse_phi_multiregion(_phi_block("Sur1", -10, 1200, n=1513).encode())
    assert len(regs) == 1 and regs[0].core_level is None
    assert any("not a core level" in n for n in regs[0].notes)


def test_junk_row_inside_single_region_is_not_a_split():
    lines = [f"{480 - i * 0.1:.2f},{100 + i}" for i in range(60)]
    lines.insert(30, "not a number row")
    sp = parse_two_column("\n".join(lines))
    assert sp.n_segments == 1 and sp.n_dropped == 1


def test_span_guard_hard_stop_and_warning():
    with pytest.raises(ValueError, match="implausibly wide"):
        check_single_region(np.linspace(-10, 1200, 500), "Nb 3d")
    assert check_single_region(np.linspace(450, 470, 200)) == []
    w = check_single_region(np.linspace(440, 485, 200))
    assert w and "unusually wide" in w[0]


def test_fit_refuses_concatenated_spectrum(kb):
    x = np.linspace(-5, 700, 1449)
    y = 1000 + np.exp(-0.5 * ((x - 204.5) / 0.6) ** 2)
    with pytest.raises(ValueError, match="implausibly wide"):
        fit_core_level(x, y, kb.components("nb2c_halim", "Nb 3d"))


def test_crop():
    sp = parse_phi_multiregion(MULTI)[2]
    c = crop(sp, 202.0, 213.0)
    assert c.be.min() >= 202.0 and c.be.max() <= 213.0 and "Cropped" in c.notes[-1]
    with pytest.raises(ValueError):
        crop(sp, 210, 205)


SUMMARY = """ Curve Fit Summary
==============================================================================
Summary File: CfSpecSummary.csv 30-Oct-2026
Setup File: Nb3d.fit
==============================================================================
Data File: demo110.spe (XPS SPECTRUM)   Region: Nb3d
Data File Comment: 5_demo
==============================================================================
Spec Band    Pos PosSep B_FWHM  FWHM Height %Gauss  Area %Area ChiSquared
   1    1 204.49   0.00   1.31  1.31   7213     97 10181 30.57       9.21
        2 205.65   1.16   1.84  1.84   2521     98  4981 14.96
        3 207.21   2.72   1.38  1.38   4630    100  6821 20.48
        4 207.49   3.00   1.69  1.69   2439     81  4778 14.35
        5 208.37   3.88   1.55  1.55   2016    100  3337 10.02
        6 210.21   5.72   1.82  1.82   1656    100  3201  9.61
 Curve Fit Summary
Data File: demo104.spe (XPS SPECTRUM)   Region: Ti2p
Data File Comment: 2_demo
Spec Band    Pos PosSep B_FWHM  FWHM Height %Gauss  Area %Area ChiSquared
   1    1 454.44   0.00   1.89  1.89    170     80   386 38.70       1.10
        2 458.46   4.02   1.67  1.67    300     80   611 61.30
 Curve Fit Summary
Data File: demo104.spe (XPS SPECTRUM)   Region: O1s
Data File Comment: 2_demo
Spec Band    Pos PosSep B_FWHM  FWHM Height %Gauss  Area %Area ChiSquared
   1    1 529.39   0.00   1.91  1.91   1500     73  3690 14.00       2.00
        2 530.23   0.84   1.47  1.47   7000     83 12404 47.05
"""


@pytest.fixture()
def bands():
    return parse_curve_fit_summaries(SUMMARY, origin="demo")


def test_summary_parser_keeps_raw_bands(bands):
    assert len(bands) == 10 and bands["summary_index"].nunique() == 3
    nb = bands[bands.region == "Nb3d"]
    assert list(nb["provider_label"]) == [f"Nb3d band {i}" for i in range(1, 7)]
    assert nb["area"].sum() == pytest.approx(33299)
    assert nb["chi2"].eq(9.21).all() and nb["core_level"].eq("Nb 3d").all()
    assert nb.iloc[0]["setup_file"] == "Nb3d.fit" and nb.iloc[0]["comment"] == "5_demo"


def test_pairing_by_split_ratio_proximity(bands):
    d, b = pair_spin_orbit(bands)
    nb = d[d.region == "Nb3d"].reset_index(drop=True)
    assert list(nb["provider_bands"]) == ["1,3", "2,5", "4,6"]
    assert (nb["pairing"] == "paired").all() and nb["split_obs"].tolist() == pytest.approx([2.72, 2.72, 2.72])
    assert nb["ratio_obs"].tolist() == pytest.approx([1.493, 1.493, 1.493], abs=2e-3)
    assert nb["area_reported"].sum() == pytest.approx(33299)       # nothing lost or double counted
    # raw bands are untouched and linked to their doublet
    assert b["provider_label"].equals(bands["provider_label"]) and b["BE"].equals(bands["BE"])
    assert set(b.loc[b.region == "Nb3d", "line"]) == {"primary", "secondary"}


def test_lone_primary_lines_flagged_and_estimated(bands):
    d, _ = pair_spin_orbit(bands)
    ti = d[d.region == "Ti2p"]
    assert (ti["pairing"] == "unpaired").all() and (ti["area_convention"] == "primary_line_only").all()
    assert ti["area_total_estimated"].tolist() == pytest.approx([386 * 1.5, 611 * 1.5])
    o = d[d.region == "O1s"]
    assert (o["pairing"] == "singlet").all()


def test_mapping_is_separate_from_provider_labels(bands, kb):
    d, _ = pair_spin_orbit(bands)
    nb = d[d.region == "Nb3d"].reset_index(drop=True)
    comps = kb.components("nb2c_halim", "Nb 3d")
    m = map_to_kb(nb, comps)
    t = m.table
    assert t["provider_label"].equals(nb["provider_label"])            # never overwritten
    assert t.set_index("provider_bands").loc["1,3", "mapped_component"] == "Nb_MXene"
    assert (t["mapping_method"].dropna() == "auto_window_tentative").all()
    vt = provider_component_table(t)
    assert set(vt.columns) >= {"component", "provider_label", "role", "mapping_method", "area"}
    for r in vt.itertuples():
        if r.mapping_method == "unmapped":
            assert r.component == r.provider_label and r.role == "other"


def test_user_mapping_logged_and_unique(bands, kb):
    d, _ = pair_spin_orbit(bands)
    nb = d[d.region == "Nb3d"].reset_index(drop=True)
    comps = kb.components("nb2c_halim", "Nb 3d")
    t = map_to_kb(nb, comps).table
    did = t.loc[t["provider_bands"] == "2,5", "doublet_id"].item()
    t2, log = apply_user_mapping(t, {did: None}, comps)
    assert t2.loc[t2.doublet_id == did, "mapped_component"].isna().all()
    with pytest.raises(ValueError, match="mapped twice"):
        apply_user_mapping(t, {did: "Nb_MXene"}, comps)


SASAKI = os.environ.get("XPS_SASAKI_DIR")


@pytest.mark.skipif(not SASAKI, reason="private provider data (set XPS_SASAKI_DIR)")
def test_real_provider_report_and_multi_files():
    root = Path(SASAKI)
    pptx = next(root.glob("*.pptx"))
    b = read_multipak_provider(pptx.read_bytes(), pptx.name)
    assert b["summary_index"].nunique() == 27 and len(b) == 105
    d, _ = pair_spin_orbit(b)
    nb = d[d.core_level == "Nb 3d"]
    assert (nb["pairing"] == "paired").all()
    assert nb["split_obs"].between(2.70, 2.74).all()
    for f in root.glob("*_multi.csv"):
        with pytest.raises(ValueError, match="multi-region"):
            parse_two_column(f.read_bytes())
        regs = parse_phi_multiregion(f.read_bytes())
        assert len(regs) >= 5 and all(r.be.max() - r.be.min() < 30 for r in regs)
