from mxene_xps.models import SampleSetup
from mxene_xps.request_generator import generate_request, make_peak_template, required_regions


def nbti_setup(**kw):
    base = dict(sample_name="NbTi_MXene_01", structure="(Nb,Ti)2CTx", m_target=2, x_target=1,
                profile_key="nbti_m2c_defect_aware", m_elements=["Nb", "Ti"], etchants=["LiF/HCl"],
                a_element="Al", expected_defects=True)
    base.update(kw)
    return SampleSetup(**base)


def test_regions_29_1():
    regions = {r["region"] for r in required_regions(nbti_setup())}
    assert regions == {"Survey", "Ti 2p", "Nb 3d", "C 1s", "O 1s", "F 1s", "Cl 2p", "Al 2p", "Li 1s"}


def test_fit_candidate_counts_and_warnings_29_1(kb):
    txt, _ = generate_request(nbti_setup(), kb)
    assert "Ti 2p: 4 candidate doublet(s)" in txt
    assert "Nb 3d: 4 candidate doublet(s)" in txt     # KB 1.3.0: + Nb_suboxide_like calibration candidate
    assert "C 1s: 4 candidate component(s)" in txt
    assert "O 1s: 4 candidate component(s)" in txt
    assert "F 1s: 3 candidate component(s)" in txt     # KB 1.3.0: + F_lowBE_unresolved
    assert "CALIBRATION CANDIDATE" in txt and "cross-check: O 1s" in txt
    assert "do not merge it into Nb2O5 or framework Nb" in txt
    assert "do NOT force a vacancy-specific" in txt
    assert "Fit each metal core level independently" in txt
    assert "UNRESOLVED" in txt  # O_low / O_mid flagged
    assert "ΔSO 5.80 eV" in txt and "ΔSO 2.80 eV" in txt


def test_cl_not_requested_without_chloride():
    regions = {r["region"] for r in required_regions(nbti_setup(etchants=["HF"]))}
    assert "Cl 2p" not in regions and "F 1s" in regions and "Li 1s" not in regions
    regions = {r["region"] for r in required_regions(nbti_setup(etchants=["HF"], cl_interest=True))}
    assert "Cl 2p" in regions


def test_conditional_regions():
    s = nbti_setup(etchants=["NaF/HCl", "KF/HCl", "TBAOH"], dopants=["V", "Ru"], a_element="Si")
    regions = [r["region"] for r in required_regions(s)]
    for r in ("Na 1s", "K 2p", "N 1s", "V 2p", "Ru 3p", "Si 2p"):
        assert r in regions
    assert len(regions) == len(set(regions))


def test_unknown_element_not_invented():
    regions = [r["region"] for r in required_regions(nbti_setup(m_elements=["Hf"]))]
    assert "Hf core level (select manually)" in regions


def test_peak_template_28_with_dopants(kb):
    t = make_peak_template(nbti_setup(dopants=["V"]), kb)
    assert len(t) == 16 + 3
    spec28 = ["C_MXene", "Ti_lowBE", "Ti_midBE", "Ti_highBE", "Ti_oxide", "Nb_MXene", "Nb2O5", "NbF5_like",
              "O_low", "O_mid", "OH_like", "O_adsorbate", "F_MXene", "F_residue"]
    cal = ["Nb_suboxide_like", "F_lowBE_unresolved"]                 # KB 1.3.0 calibration candidates
    assert sorted(t["component"][:16]) == sorted(spec28 + cal)
    assert list(t.columns) == ["sample", "core_level", "element", "component", "BE", "FWHM", "area", "role", "notes"]
    v = t[t["element"] == "V"]
    assert set(v["role"]) == {"dopant_lattice_candidate", "dopant_adsorbed", "dopant_oxide"}
    assert set(v["component"]) == {"V_lattice_or_vacancy", "V_adsorbed", "V_oxide"}
    assert t["area"].isna().all()
