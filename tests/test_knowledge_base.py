import json

from conftest import ROOT


# The only numerical deviations from the §27 seed: user decision D2 (2026-10-16).
D2_EXCEPTIONS = {("nb2c_halim", "Nb 3d", n, f) for n in ("Nb_MXene", "Nb2O5", "NbF5_like") for f in ("center", "fwhm")}


def test_numerical_values_identical_to_handoff_seed(kb):
    """§27: values must be preserved exactly (except the explicit D2 decision); only metadata may be added."""
    seed = json.loads((ROOT / "tools" / "kb_seed_from_handoff_s27.json").read_text())
    keys = ("name", "kind", "center", "bounds", "fwhm", "fwhm_bounds", "split", "ratio_primary_secondary", "role", "source")
    for pkey, prof in seed["profiles"].items():
        built = kb.raw["profiles"][pkey]["fit_models"]
        assert list(built) == list(prof["fit_models"])
        for core, comps in prof["fit_models"].items():
            assert len(built[core]) == len(comps)
            for a, b in zip(comps, built[core]):
                for k in keys:
                    if (pkey, core, a["name"], k) in D2_EXCEPTIONS:
                        continue
                    assert a.get(k) == b.get(k), (pkey, core, a["name"], k)
                if "interpretation" in a:
                    assert a["interpretation"] == b["interpretation"]
    for k in ("structures", "legacy_excel", "solid_solution_reference", "principles"):
        for kk, vv in seed[k].items() if isinstance(seed[k], dict) else enumerate(seed[k]):
            assert kb.raw[k][kk] == vv


def test_every_component_has_source_doi_and_caveat(kb):
    for pkey in kb.profile_keys:
        for core, comps in kb.fit_models(pkey).items():
            for c in comps:
                assert c.interpretation and c.source in kb.references
                assert kb.references[c.source].get("doi"), c.source


def test_mixed_metal_O_low_mid_remain_unresolved(kb):
    o = {c.name: c for c in kb.components("nbti_m2c_defect_aware", "O 1s")}
    assert o["O_low"].role == "O_unresolved_low" and o["O_low"].assignment_status == "unresolved"
    assert o["O_mid"].role == "O_unresolved_mid" and o["O_mid"].assignment_status == "unresolved"


def test_no_vacancy_specific_nb_peak(kb):
    for pkey in kb.profile_keys:
        for c in kb.fit_models(pkey).get("Nb 3d", []):
            assert "vacan" not in c.name.lower()
            assert c.role in ("framework_metal", "metal_oxide", "metal_fluoride")


def test_ti_components_not_named_by_oxidation_state(kb):
    for pkey in kb.profile_keys:
        for c in kb.fit_models(pkey).get("Ti 2p", []):
            assert "+" not in c.name and "Ti2" not in c.name and "Ti3" not in c.name and "Ti4" not in c.name
