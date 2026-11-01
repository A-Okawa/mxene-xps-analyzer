"""KB 1.4.0 MAX-phase literature calibration: literature reference vs transferred reference vs provider observation."""
import pytest

from mxene_xps.models import FitComponent

MAX = ("nb2alc_max", "nbti_211_max")
PROVIDER = "Sasaki2026MultiPak"


def _c(kb, p, core, name):
    return next(c for c in kb.components(p, core) if c.name == name)


def _all(kb, p):
    return [c for cs in kb.fit_models(p).values() for c in cs]


# ------------------------------------------------------------------ 1. literature_reference
def test_nb2alc_literature_values_as_specified(kb):
    nb = _c(kb, "nb2alc_max", "Nb 3d", "Nb_MAX")
    ox = _c(kb, "nb2alc_max", "Nb 3d", "Nb2O5")
    cm = _c(kb, "nb2alc_max", "C 1s", "C_MAX")
    for comp, center, sec in ((nb, 202.9, 205.5), (ox, 206.8, 209.6), (cm, 282.5, None)):
        assert comp.source_type == "literature" and comp.center == center
        (ref,) = comp.literature_reference
        assert ref["reference_type"] == "literature_reference" and ref["material"].startswith("Nb2AlC")
        assert ref["center"] == center and ref["source"] == comp.source
        if sec is not None:
            assert ref["secondary"]["center"] == sec


def test_literature_positions_not_hard_fixed(kb):
    """Reference values live apart from the editable fit windows; positions are free within the window."""
    for p in MAX:
        for c in _all(kb, p):
            assert c.bounds[0] < c.center < c.bounds[1], (p, c.name)
            assert c.bounds[1] - c.bounds[0] >= 0.8, (p, c.name)
            for r in c.literature_reference or []:
                assert "constraint" in r["note"] or "not" in r["note"]


def test_literature_citation_metadata_preserved(kb):
    for key in ("Naslund2020", "Greczynski2018Carbides"):
        ref = kb.references[key]
        assert ref["doi"] and ref["year"] and ref["journal"] and ref["be_reference"]
    assert kb.references["Naslund2020"]["doi"] == "10.1021/acs.jpcc.0c07413"
    assert kb.references["Greczynski2018Carbides"]["doi"] == "10.1016/j.apsusc.2017.11.264"
    user = kb.references["UserNb2AlC2026"]
    assert user["doi"] is None and user["citation_status"] == "pending"     # never a fabricated citation
    for p in MAX:
        for c in _all(kb, p):
            for r in (c.literature_reference or []) + (c.transferred_reference or []):
                assert r["source"] in kb.references
                if r["source"] == "UserNb2AlC2026":
                    assert r.get("citation_status") == "pending"


def test_al2p_nb2alc_from_naslund_shift(kb):
    al = _c(kb, "nb2alc_max", "Al 2p", "Al_MAX")
    (ref,) = al.literature_reference
    assert ref["source"] == "Naslund2020" and ref["center"] == pytest.approx(72.76 - 0.36)
    assert "derived" in ref


# ------------------------------------------------------------------ 2. transferred reference
def test_nbti_211_is_transferred_not_literature(kb):
    """No direct (Nb,Ti)2AlC peak table: nothing in nbti_211_max may claim to be literature for that material."""
    for c in _all(kb, "nbti_211_max"):
        assert not c.literature_reference, c.name
        if c.name == "Nb_suboxide_like":
            continue
        assert c.source_type in ("transferred", "implementation_surrogate"), c.name
        if c.source_type == "transferred":
            assert c.confidence in ("low", "medium"), c.name
            assert c.transferred_reference, c.name
            for r in c.transferred_reference:
                assert r["reference_type"] == "transferred_reference"
                assert "(Nb,Ti)2AlC" not in r["material"]
                assert r["confidence"] in ("low", "medium")
    for name in ("Ti_MAX", "Nb_MAX", "C_MAX", "Al_MAX"):
        c = next(x for x in _all(kb, "nbti_211_max") if x.name == name)
        assert "not independently validated for (Nb,Ti)2AlC" in c.interpretation
    assert "TRANSFERRED" in kb.profile("nbti_211_max")["notes"][0]


def test_transferred_sources(kb):
    ti = _c(kb, "nbti_211_max", "Ti 2p", "Ti_MAX")
    mats = {r["material"] for r in ti.transferred_reference}
    assert any("Ti3AlC2" in m for m in mats) and ti.center == 454.75
    al = _c(kb, "nbti_211_max", "Al 2p", "Al_MAX")
    assert {r["material"] for r in al.transferred_reference} == {"Nb2AlC", "Ti2AlC"}
    nb = _c(kb, "nbti_211_max", "Nb 3d", "Nb_MAX")
    assert nb.transferred_reference[0]["material"] == "Nb2AlC" and nb.confidence == "low"


# ------------------------------------------------------------------ 3. provider calibration observation
def test_provider_values_only_as_observations(kb):
    for p in MAX:
        for c in _all(kb, p):
            if c.name == "Nb_suboxide_like":            # unvalidated candidate, deliberately unchanged
                continue
            assert c.source != PROVIDER, (p, c.name)
            assert c.source_type != "provider_derived", (p, c.name)
            for r in (c.literature_reference or []) + (c.transferred_reference or []):
                assert r["source"] != PROVIDER
            for o in c.calibration_observations or []:
                assert o["reference_type"] == "provider_calibration_observation" and o["source"] == PROVIDER
                assert "never a KB reference" in o["note"]


def test_reference_layers_cannot_be_mixed():
    base = dict(name="x", core_level="C 1s", kind="single", center=282.5, bounds=(282.0, 283.0), fwhm=1.0,
                fwhm_bounds=(0.5, 2.0), role="framework_C", source="s", interpretation="i")
    obs = {"reference_type": "provider_calibration_observation", "center": 282.5}
    with pytest.raises(ValueError, match="stored under literature_reference"):
        FitComponent(**base, literature_reference=[obs])
    with pytest.raises(ValueError, match="source_type"):
        FitComponent(**base, source_type="provider_fit")


# ------------------------------------------------------------------ unchanged items
def test_candidate_d4_and_splittings_unchanged(kb):
    for p in MAX + ("nb2c_halim", "nbti_m2c_defect_aware"):
        c = _c(kb, p, "Nb 3d", "Nb_suboxide_like")
        assert c.calibration_candidate and not c.validated and c.role == "metal_suboxide_candidate"
    cs = kb.raw["calibration_status"]
    assert cs["nb3d_suboxide_like_validated"] is False and cs["d4_nb3d_dSO_calibration_sensitive"] is True
    for p in MAX:
        assert {c.split for c in kb.components(p, "Nb 3d")} == {2.8}
    assert {c.split for c in kb.components("nbti_211_max", "Ti 2p")} == {5.8, 5.7}


def test_mxene_profiles_untouched_by_max_step(kb):
    nb = _c(kb, "nb2c_halim", "Nb 3d", "Nb_MXene")
    assert (nb.center, nb.bounds, nb.source_type) == (204.8, (204.1, 205.6), None)
