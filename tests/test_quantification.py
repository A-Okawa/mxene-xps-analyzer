import numpy as np
import pandas as pd
import pytest

from mxene_xps.quantification import (W_DEFICIENCY, W_DOPANT, W_UNRESOLVED_O, QuantOptions, legacy_quantify,
                                      quantify_sample, quantify_table, role_reclassifications)


def row(comp, core, area, role, el=None, sample="s1"):
    return {"sample": sample, "core_level": core, "component": comp, "area": area, "role": role,
            "element": el or core.split()[0]}


BASE = [
    row("C_MXene", "C 1s", 1.0, "framework_C"),
    row("Ti_lowBE", "Ti 2p", 0.2, "framework_metal"),
    row("Ti_midBE", "Ti 2p", 0.2, "framework_metal"),
    row("Ti_oxide", "Ti 2p", 0.3, "metal_oxide"),
    row("Nb_MXene", "Nb 3d", 1.2, "framework_metal"),
    row("Nb2O5", "Nb 3d", 0.5, "metal_oxide"),
    row("O_low", "O 1s", 0.3, "O_unresolved_low"),
    row("O_mid", "O 1s", 0.25, "O_unresolved_mid"),
    row("OH_like", "O 1s", 0.4, "termination_OH_like"),
    row("O_adsorbate", "O 1s", 0.2, "adsorbate"),
    row("F_MXene", "F 1s", 0.6, "termination_F"),
    row("F_residue", "F 1s", 0.3, "residue_F"),
]
OPT = QuantOptions(m_target=2, x_target=1, input_mode="corrected_area", nominal_metal_stoich={"Ti": 0.8, "Nb": 1.2})


def test_mixed_metal_deficiency_29_2():
    r = quantify_sample(pd.DataFrame(BASE), {}, OPT)
    v = r.values
    assert v["M_norm"] == pytest.approx(1.6, abs=1e-12)
    assert v["apparent_M_deficiency_delta"] == pytest.approx(0.4, abs=1e-12)
    assert v["Ti_norm"] == pytest.approx(0.4, abs=1e-12)
    assert v["Nb_norm"] == pytest.approx(1.2, abs=1e-12)
    assert W_DEFICIENCY in r.warnings
    # selective depletion is reported separately from total deficiency (§9.4)
    assert v["Ti_observed_fraction"] == pytest.approx(0.25)
    assert v["Ti_fraction_shift"] == pytest.approx(-0.15)
    assert v["Ti_retention_ratio_fractional"] == pytest.approx(0.625)
    assert v["Nb_fraction_shift"] == pytest.approx(0.15)


def test_metal_sum_not_forced_to_ideal():
    v = quantify_sample(pd.DataFrame(BASE), {}, OPT).values
    assert v["Ti_norm"] + v["Nb_norm"] != pytest.approx(2.0)


def test_unresolved_O_and_F_residue_excluded_by_default():
    r = quantify_sample(pd.DataFrame(BASE), {}, OPT)
    v = r.values
    assert v["z_O_total"] == pytest.approx(0.4)          # OH_like only
    assert v["O_unresolved_low_norm"] == pytest.approx(0.3)
    assert v["O_unresolved_mid_norm"] == pytest.approx(0.25)
    assert v["z_F_raw"] == pytest.approx(0.6)            # residue excluded
    assert v["F_residue_norm"] == pytest.approx(0.3)
    assert W_UNRESOLVED_O in r.warnings
    v2 = quantify_sample(pd.DataFrame(BASE), {}, QuantOptions(**{**OPT.__dict__, "include_f_residue": True})).values
    assert v2["z_F_raw"] == pytest.approx(0.9)


def test_metal_oxide_and_fluoride_not_framework():
    v = quantify_sample(pd.DataFrame(BASE), {}, OPT).values
    assert v["M_framework_corrected"] == pytest.approx(1.6)


def test_dopant_separation_29_3():
    doped = pd.DataFrame(BASE + [row("V_lattice_or_vacancy", "V 2p", 0.10, "dopant_lattice_candidate"),
                                 row("V_oxide", "V 2p", 0.05, "dopant_oxide")])
    r = quantify_sample(doped, {}, OPT)
    assert r.values["M_norm"] == pytest.approx(1.6, abs=1e-12)
    assert r.values["V_dopant_total_norm"] == pytest.approx(0.15, abs=1e-12)
    assert r.values["V_dopant_lattice_candidate_norm"] == pytest.approx(0.10, abs=1e-12)
    assert r.values["V_dopant_oxide_norm"] == pytest.approx(0.05, abs=1e-12)
    assert W_DOPANT in r.warnings
    # explicit user reclassification is the only way into framework M
    re = doped.copy()
    re.loc[re["component"] == "V_lattice_or_vacancy", "role"] = "framework_metal"
    v = quantify_sample(re, {}, OPT).values
    assert v["M_norm"] == pytest.approx(1.7) and v["V_norm"] == pytest.approx(0.1)
    audit = role_reclassifications(re, {("V 2p", "V_lattice_or_vacancy"): "dopant_lattice_candidate"})
    assert len(audit) == 1 and audit.iloc[0]["user_role"] == "framework_metal"


def test_scale_uses_x_target():
    rows = [row("C_MXene", "C 1s", 0.5, "framework_C"), row("Ti", "Ti 2p", 0.9, "framework_metal")]
    v = quantify_sample(pd.DataFrame(rows), {}, QuantOptions(m_target=3, x_target=2, input_mode="corrected_area")).values
    assert v["scale_X_target_per_C"] == pytest.approx(4.0)
    assert v["M_norm"] == pytest.approx(3.6)


def test_raw_area_rsf_division_and_missing_rsf_blocks():
    rows = [row("C", "C 1s", 0.314, "framework_C"), row("Ti", "Ti 2p", 2.077 * 3, "framework_metal")]
    v = quantify_sample(pd.DataFrame(rows), {"C 1s": 0.314, "Ti 2p": 2.077},
                        QuantOptions(m_target=3, x_target=2)).values
    assert v["M_norm"] == pytest.approx(6.0)
    with pytest.raises(ValueError, match="RSF is required"):
        quantify_sample(pd.DataFrame(rows), {"C 1s": 0.314}, QuantOptions(m_target=3, x_target=2))
    with pytest.raises(ValueError, match="RSF is required"):
        quantify_sample(pd.DataFrame(rows), {"C 1s": 0.314, "Ti 2p": float("nan")}, QuantOptions(m_target=3, x_target=2))


def test_missing_framework_C():
    with pytest.raises(ValueError, match="framework_C"):
        quantify_sample(pd.DataFrame([row("Ti", "Ti 2p", 1, "framework_metal")]), {}, OPT)


def test_empty_optional_components_and_bad_input():
    rows = BASE + [row("Ti_highBE", "Ti 2p", np.nan, "framework_metal"), row("X", "Ti 2p", "", "framework_metal")]
    r = quantify_sample(pd.DataFrame(rows), {}, OPT)
    assert r.values["M_norm"] == pytest.approx(1.6)
    assert any("blank area" in w for w in r.warnings)
    with pytest.raises(ValueError, match="Non-numeric"):
        quantify_sample(pd.DataFrame(BASE + [row("Y", "Ti 2p", "abc", "framework_metal")]), {}, OPT)
    with pytest.raises(ValueError, match="Unknown role"):
        quantify_sample(pd.DataFrame(BASE + [row("Y", "Ti 2p", 1, "Ti3+")]), {}, OPT)
    with pytest.raises(ValueError, match="Negative"):
        quantify_sample(pd.DataFrame(BASE + [row("Y", "Ti 2p", -1, "framework_metal")]), {}, OPT)


def test_quantify_table_multiple_samples():
    df = pd.DataFrame(BASE + [dict(r, sample="s2") for r in BASE])
    res, warns = quantify_table(df, {}, OPT)
    assert list(res["sample"]) == ["s1", "s2"] and set(warns) == {"s1", "s2"}


def test_role_model_matches_legacy_formulas():
    """Role-based engine with the legacy role mapping == legacy model (F1/F2 off)."""
    a = {"C_TiTx": 2891, "Ti1": 9497, "Ti2": 14535, "Ti3": 9558, "Ti_oxide": 8599, "O_low": 4734, "O_oxide": 2465,
         "O_high": 1724, "F0": 11079, "F1": 3385, "F2": 2608, "Cl": 753}
    rsf = {"C 1s": 0.314, "Ti 2p": 2.077, "O 1s": 0.733, "F 1s": 1.0, "Cl 2p": 0.954}
    leg = legacy_quantify(a, rsf)
    rows = [row("C", "C 1s", a["C_TiTx"], "framework_C"),
            row("T1", "Ti 2p", a["Ti1"], "framework_metal"), row("T2", "Ti 2p", a["Ti2"], "framework_metal"),
            row("T3", "Ti 2p", a["Ti3"], "framework_metal"), row("TO", "Ti 2p", a["Ti_oxide"], "metal_oxide"),
            row("OL", "O 1s", a["O_low"], "termination_O"), row("OX", "O 1s", a["O_oxide"], "oxide_O"),
            row("OH", "O 1s", a["O_high"], "termination_OH_like"), row("F0", "F 1s", a["F0"], "termination_F"),
            row("F1", "F 1s", a["F1"], "residue_F"), row("F2", "F 1s", a["F2"], "residue_F"),
            row("Cl", "Cl 2p", a["Cl"], "termination_Cl")]
    v = quantify_sample(pd.DataFrame(rows), rsf, QuantOptions(m_target=3, x_target=2, oxide_f_correction_factor=0.6667)).values
    pairs = {"Ti_norm": "M_norm", "z_O": "z_O_total", "z_F_raw": "z_F_raw", "z_F_corr": "z_F_corrected", "z_Cl": "z_Cl",
             "z_total_corr": "z_total_corrected", "O%_corr": "O_fraction_corr", "F%_corr": "F_fraction_corr"}
    for lk, rk in pairs.items():
        assert v[rk] == pytest.approx(leg[lk], rel=1e-12, abs=1e-12), lk


def test_all_blank_table_gives_clear_error():
    rows = [dict(r, area=np.nan) for r in BASE]
    with pytest.raises(ValueError, match="no component areas entered"):
        quantify_table(pd.DataFrame(rows), {}, OPT)
