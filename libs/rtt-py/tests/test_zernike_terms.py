"""Zernike sag terms in Python (#316, ADR 0051, PR 2): the model classes, the layout field terms,
the result format 0.1.23, the paraxial rule for a rotationally symmetric and an astigmatic
term, and the sag of Fringe against Noll."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import layout, model, results

MIRROR = REFERENCE_DIR / "m8" / "zernike_mirror.rtt.json"


def test_model_classes_of_the_reference_file() -> None:
    s = rt.load(MIRROR)
    mirror = s.root.children[1]
    assert isinstance(mirror, model.Element)
    noll, fringe = mirror.surfaces[0].shape.terms
    assert isinstance(noll, model.ZernikeSag) and isinstance(fringe, model.ZernikeFringeSag)
    assert noll.normalization_radius.value == 22.0
    assert noll.coefficients[10].value == 2e-05 and noll.coefficients[10].variable
    assert [c.value for c in fringe.coefficients] == [0.0, 0.0, 0.0, -1e-05]


def test_layout_lists_the_terms_and_reports_the_base() -> None:
    stop, mirror, detector = layout.surfaces(rt.compile(rt.load(MIRROR)))
    assert mirror.shape == "conic"  # the base of the stack
    assert mirror.terms == [("zernike_sag", {"normalization_radius": 22.0, "terms": 11.0}),
                            ("zernike_fringe_sag", {"normalization_radius": 22.0, "terms": 4.0})]
    assert stop.terms == [] and detector.terms == []


def test_result_format_carries_the_terms() -> None:
    surfaces = layout.surfaces(rt.compile(rt.load(MIRROR)))
    text = surfaces[1].to_json()
    doc = json.loads(text)
    assert doc["schema_version"] == results.SCHEMA_VERSION == "0.1.25"
    assert doc["data"]["terms"][1] == ["zernike_fringe_sag",
                                      {"normalization_radius": 22.0, "terms": 4.0}]
    old = dict(doc)
    old["schema_version"] = "0.1.22"
    old["data"] = {k: v for k, v in doc["data"].items() if k != "terms"}
    assert results.load_json(json.dumps(old)).data["terms"] == []


def test_a_symmetric_term_changes_the_paraxial_curvature() -> None:
    """Noll Z11 (spherical) and Fringe term 4 (defocus 2 rho^2 - 1) add to the vertex curvature
    (ADR 0041, point 2): c = 1/R + 2 sqrt(5) a11 (-6) / rho^2 + 2 a4 (2) / rho^2. The mirror
    focal length is -1 / (2 c) for the light along +z (f' = R/2 of a mirror with R < 0)."""
    cs = rt.compile(rt.load(MIRROR))
    efl = rt.paraxial.first_order(cs, path=0).efl
    rho = 22.0
    c = 1.0 / -200.0 + 2.0 * math.sqrt(5.0) * 2e-05 * -6.0 / rho**2 + 2.0 * -1e-05 * 2.0 / rho**2
    plain = json.loads(MIRROR.read_text(encoding="utf-8"))
    surface = plain["root"]["children"][1]["surfaces"][0]
    surface["shape"] = {"base": {"type": "conic", "radius": 1.0 / c, "conic": 0.0}}
    efl_plain = rt.paraxial.first_order(rt.compile(rt.System.from_json(json.dumps(plain))),
                                        path=0).efl
    assert efl is not None and efl_plain is not None
    assert math.isclose(efl, efl_plain, rel_tol=1e-13)
    # content guard: the terms move the focus away from that of the base sphere alone
    surface["shape"] = {"base": {"type": "conic", "radius": -200.0, "conic": 0.0}}
    efl_base = rt.paraxial.first_order(rt.compile(rt.System.from_json(json.dumps(plain))),
                                       path=0).efl
    assert efl_base is not None and abs(efl - efl_base) > 0.01  # 0.024 mm, against 1e-11


def test_an_astigmatic_term_has_no_paraxial_model() -> None:
    data = json.loads(MIRROR.read_text(encoding="utf-8"))
    terms = data["root"]["children"][1]["surfaces"][0]["shape"]["terms"]
    terms[1]["coefficients"] = [0.0, 0.0, 0.0, 0.0, 1e-05]  # Fringe 5: rho^2 cos 2 theta
    cs = rt.compile(rt.System.from_json(json.dumps(data)))
    with pytest.raises(rt.errors.ParaxialError):
        rt.paraxial.first_order(cs, path=0)


def test_fringe_equals_noll_through_the_layout_sag() -> None:
    """Fringe term 9 (6 rho^4 - 6 rho^2 + 1) equals Noll Z11 / sqrt(5) (ADR 0051, point 3)."""
    def surface_sag(term: dict[str, object], x: float, y: float) -> float:
        data = json.loads(MIRROR.read_text(encoding="utf-8"))
        data["root"]["children"][1]["surfaces"][0]["shape"] = {
            "base": {"type": "plane"}, "terms": [term]}
        cs = rt.compile(rt.System.from_json(json.dumps(data)))
        return float(rt.layout.sag(cs, 1, np.array([x]), np.array([y]))[0])

    a = 3e-5
    fringe = {"type": "zernike_fringe_sag", "normalization_radius": 22.0,
              "coefficients": [0.0] * 8 + [a]}
    noll = {"type": "zernike_sag", "normalization_radius": 22.0,
            "coefficients": [0.0] * 10 + [a / math.sqrt(5.0)]}
    for x, y in ((0.0, 0.0), (5.0, -3.0), (15.0, 10.0)):
        assert math.isclose(surface_sag(fringe, x, y), surface_sag(noll, x, y), rel_tol=1e-14,
                            abs_tol=1e-20)
    assert surface_sag(fringe, 0.0, 0.0) == a  # content guard: 1 at the vertex
