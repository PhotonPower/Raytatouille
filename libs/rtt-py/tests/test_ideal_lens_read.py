"""The ideal lens and the ideal cylinder lens from Python (#178, ADR 0031): read types in
raytatouille.model, the edit form, and the interim errors until the tracer and rtt-paraxial
know them."""

from __future__ import annotations

import json

import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import model

FILE = REFERENCE_DIR / "r2" / "ideal_lens.rtt.json"


def lens_interactions(system: rt.System) -> tuple[model.IdealLens, model.IdealCylinderLens]:
    il = system.root.children[1]
    cl = system.root.children[2]
    assert isinstance(il, model.Element) and isinstance(cl, model.Element)
    assert il.kind == model.ElementKind.THIN_ELEMENT
    a, b = il.surfaces[0].interaction, cl.surfaces[0].interaction
    assert isinstance(a, model.IdealLens) and isinstance(b, model.IdealCylinderLens)
    return a, b


def test_read_types() -> None:
    lens, cylinder = lens_interactions(rt.load(FILE))
    assert lens.focal_length.value == 50.0 and not lens.focal_length.variable
    assert lens.object_distance is not None and lens.object_distance.value == 75.0
    assert cylinder.focal_length.value == 50.0
    assert cylinder.axis_deg == 0.0
    assert cylinder.object_distance is not None and cylinder.object_distance.value == 75.0


def test_object_at_infinity_has_no_object_distance() -> None:
    data = json.loads(FILE.read_text(encoding="utf-8"))
    data["root"]["children"][1]["surfaces"][0]["interaction"].pop("object_distance")
    lens, _ = lens_interactions(rt.System.from_json(json.dumps(data)))
    assert lens.object_distance is None


def test_edit_form_and_patch() -> None:
    system = rt.load(FILE)
    pointer = "/root/children/2/surfaces/0/interaction"
    assert system.json_at(pointer + "/axis_deg") == 0.0  # always in the edit form
    changed = rt.apply_patch(system, [
        {"op": "replace", "path": pointer + "/axis_deg", "value": 90.0},
        {"op": "replace", "path": "/root/children/1/surfaces/0/interaction/focal_length/value",
         "value": -40.0},
    ])
    lens, cylinder = lens_interactions(changed)
    assert cylinder.axis_deg == 90.0 and lens.focal_length.value == -40.0
    with pytest.raises(rt.EditError):  # f = 0 is interaction.focal_length_invalid
        rt.apply_patch(system, [{"op": "replace",
                                 "path": "/root/children/1/surfaces/0/interaction/focal_length/value",
                                 "value": 0.0}])


def test_interim_errors_until_the_tracer_and_paraxial_know_the_lens() -> None:
    # ADR 0031, point 8: never silently wrong in between.
    cs = rt.compile(rt.load(FILE))
    rays = rt.trace.RayBatch(1)
    rays.pos_z[0] = -75.0  # the object point, in front of the stop plane z = -1
    rt.trace.trace(cs, rays, path="ideal lens")
    assert rays.status[0] == int(rt.trace.RayStatus.EVENT_IMPOSSIBLE)
    with pytest.raises(rt.ParaxialError):
        rt.paraxial.first_order(cs, path="ideal lens")
