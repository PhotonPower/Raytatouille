"""Aperture dimensions as Params in Python (#320, ADR 0050, point 1)."""

from __future__ import annotations

import json

from conftest import REFERENCE_DIR

import raytatouille as rt

STOP = "/root/children/0/surfaces/0/aperture"


def test_aperture_dimensions_are_params() -> None:
    s = rt.load(REFERENCE_DIR / "m0/singlet.rtt.json")
    element = s.root.children[0]
    assert isinstance(element, rt.model.Element)
    aperture = element.surfaces[0].aperture
    assert isinstance(aperture, rt.model.CircularAperture)
    assert isinstance(aperture.radius, rt.model.Param)
    radius = aperture.radius.value
    assert radius is not None and radius > 0.0 and not aperture.radius.variable
    assert aperture.inner_radius.value == 0.0
    # The edit form writes the Param as an object, the canonical form as a number.
    assert s.json_at(STOP + "/radius") == {"value": radius, "variable": False}
    canonical = json.loads(s.to_json())
    assert canonical["root"]["children"][0]["surfaces"][0]["aperture"]["radius"] == radius


def test_a_variable_stop_radius_round_trips() -> None:
    ed = rt.Editor(rt.load(REFERENCE_DIR / "m0/singlet.rtt.json"))
    ed.set(STOP + "/radius", {"value": 8.0, "variable": True, "min": 2.0})
    element = ed.system.root.children[0]
    assert isinstance(element, rt.model.Element)
    aperture = element.surfaces[0].aperture
    assert isinstance(aperture, rt.model.CircularAperture)
    assert (aperture.radius.value, aperture.radius.variable, aperture.radius.min) == (8.0, True, 2.0)
    back = rt.System.from_json(ed.system.to_json())
    assert back.to_json() == ed.system.to_json()
    # Undo restores the plain radius.
    ed.undo()
    assert ed.system.json_at(STOP + "/radius")["variable"] is False
