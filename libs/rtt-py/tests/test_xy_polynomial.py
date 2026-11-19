"""The XY polynomial term in Python (#316, ADR 0051, point 4, PR 3): the model class, the layout
field terms, the result format, the paraxial rule and the sag through the layout."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import layout, model, results

MIRROR = REFERENCE_DIR / "m8" / "xy_mirror.rtt.json"


def _with_term(term: dict[str, object]) -> rt.CompiledSystem:
    data = json.loads(MIRROR.read_text(encoding="utf-8"))
    data["root"]["children"][1]["surfaces"][0]["shape"] = {"base": {"type": "plane"},
                                                           "terms": [term]}
    return rt.compile(rt.System.from_json(json.dumps(data)))


def test_model_class_of_the_reference_file() -> None:
    s = rt.load(MIRROR)
    mirror = s.root.children[1]
    assert isinstance(mirror, model.Element)
    (term,) = mirror.surfaces[0].shape.terms
    assert isinstance(term, model.XyPolynomial)
    assert term.normalization_radius.value == 22.0
    assert len(term.coefficients) == 21
    assert term.coefficients[20].value == 5e-06 and term.coefficients[20].variable


def test_layout_lists_the_term_and_the_result_format_carries_it() -> None:
    stop, mirror, detector = layout.surfaces(rt.compile(rt.load(MIRROR)))
    assert mirror.shape == "conic"  # the base of the stack
    assert mirror.terms == [("xy_polynomial", {"normalization_radius": 22.0, "terms": 21.0})]
    assert stop.terms == [] and detector.terms == []
    doc = json.loads(mirror.to_json())
    assert doc["schema_version"] == results.SCHEMA_VERSION
    assert doc["data"]["terms"] == [["xy_polynomial",
                                     {"normalization_radius": 22.0, "terms": 21.0}]]


def test_sag_is_the_sum_of_the_monomials() -> None:
    """Terms 1 (x), 4 (xy), 12 (x^2 y^2) and 13 (x y^3) after [UG] "Extended Polynomial"."""
    r = 22.0
    c = [0.0] * 13
    c[0], c[3], c[11], c[12] = 1e-3, -2e-3, 5e-4, 3e-4
    cs = _with_term({"type": "xy_polynomial", "normalization_radius": r, "coefficients": c})
    for x, y in ((0.0, 0.0), (5.0, -3.0), (-15.0, 10.0)):
        u, v = x / r, y / r
        expected = 1e-3 * u - 2e-3 * u * v + 5e-4 * u * u * v * v + 3e-4 * u * v**3
        got = float(rt.layout.sag(cs, 1, np.array([x]), np.array([y]))[0])
        assert math.isclose(got, expected, rel_tol=1e-14, abs_tol=1e-20)
    assert float(rt.layout.sag(cs, 1, np.array([5.0]), np.array([-3.0]))[0]) != 0.0


def test_a_rotational_term_changes_the_paraxial_curvature() -> None:
    """c_x2 = c_y2 = 1e-3 adds 2e-3 / 22^2 to the vertex curvature; the fourth-order part and the
    x^6 term do not enter first order (ADR 0041, point 2)."""
    cs = rt.compile(rt.load(MIRROR))
    efl = rt.paraxial.first_order(cs, path=0).efl
    c = 1.0 / -200.0 + 2.0 * 1e-3 / (22.0 * 22.0)
    plain = json.loads(MIRROR.read_text(encoding="utf-8"))
    plain["root"]["children"][1]["surfaces"][0]["shape"] = {
        "base": {"type": "conic", "radius": 1.0 / c}}
    efl_plain = rt.paraxial.first_order(rt.compile(rt.System.from_json(json.dumps(plain))),
                                        path=0).efl
    assert efl is not None and efl_plain is not None
    assert math.isclose(efl, efl_plain, rel_tol=1e-13)
    # content guard: the term moves the focus away from that of the base sphere (0.08 mm)
    plain["root"]["children"][1]["surfaces"][0]["shape"] = {
        "base": {"type": "conic", "radius": -200.0}}
    efl_base = rt.paraxial.first_order(rt.compile(rt.System.from_json(json.dumps(plain))),
                                       path=0).efl
    assert efl_base is not None and abs(efl - efl_base) > 0.01


def test_an_astigmatic_term_has_no_paraxial_model() -> None:
    cs = _with_term({"type": "xy_polynomial", "normalization_radius": 22.0,
                     "coefficients": [0.0, 0.0, 1e-3, 0.0, -1e-3]})
    with pytest.raises(rt.errors.ParaxialError):
        rt.paraxial.first_order(cs, path=0)
