"""The ideal lens and the ideal cylinder lens from Python (#178, ADR 0031): read types in
raytatouille.model, the edit form, the tracer and the paraxial data."""

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


def test_the_tracer_images_the_object_point_and_paraxial_knows_the_lens() -> None:
    # The axial ray of the object point reaches its image at z = 150 on the axis (PR 2); the
    # paraxial data have phi = |n| / f (PR 3): rear focal length f, image at 150, m = -2. The
    # cylinder lens has no paraxial model (ADR 0031, point 6).
    cs = rt.compile(rt.load(FILE))
    rays = rt.trace.RayBatch(1)
    rays.pos_z[0] = -75.0  # the object point, in front of the stop plane z = -1
    rt.trace.trace(cs, rays, path="ideal lens")
    assert rays.status[0] == int(rt.trace.RayStatus.ALIVE)
    assert (rays.pos_x[0], rays.pos_y[0], rays.pos_z[0]) == (0.0, 0.0, 150.0)
    fo = rt.paraxial.first_order(cs, path="ideal lens")
    assert fo.rear_focal_length is not None and fo.image_z is not None
    assert fo.lateral_magnification is not None
    assert abs(fo.rear_focal_length - 50.0) <= 1e-12
    assert abs(fo.image_z - 150.0) <= 1e-12 and abs(fo.lateral_magnification + 2.0) <= 1e-14
    with pytest.raises(rt.ParaxialError):
        rt.paraxial.first_order(cs, path="cylinder lens")
