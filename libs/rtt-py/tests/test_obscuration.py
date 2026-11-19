"""Obscurations in Python (#320, ADR 0050, point 3): the model classes, the edit form, the
layout field obscurations, the result format 0.1.22 and an end-to-end trace through
m8/spider_telescope."""

from __future__ import annotations

import json

from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import layout, model, results

TELESCOPE = REFERENCE_DIR / "m8" / "spider_telescope.rtt.json"
STOP = "/root/children/0/surfaces/0"


def test_model_classes_of_the_reference_file() -> None:
    s = rt.load(TELESCOPE)
    stop = s.root.children[0]
    mirror = s.root.children[1]
    assert isinstance(stop, model.Element) and isinstance(mirror, model.Element)
    disc, bar_x, bar_y = stop.surfaces[0].obscurations
    assert isinstance(disc, model.CircularObscuration)
    assert (disc.radius.value, disc.decenter_x.value, disc.decenter_y.value) == (10.0, 0.0, 0.0)
    assert isinstance(bar_x, model.RectangularObscuration)
    assert (bar_x.half_width_x.value, bar_x.half_width_y.value, bar_x.rotation_deg) == (
        31.0, 0.5, 0.0)
    assert isinstance(bar_y, model.RectangularObscuration) and bar_y.rotation_deg == 90.0
    (clip,) = mirror.surfaces[0].obscurations
    assert isinstance(clip, model.EllipticalObscuration)
    assert (clip.semi_axis_x.value, clip.semi_axis_y.value, clip.rotation_deg) == (3.0, 1.5, 30.0)
    assert (clip.decenter_x.value, clip.decenter_y.value) == (20.0, 5.0)
    assert clip.decenter_y.variable
    # A surface without obscurations has an empty list.
    detector = s.root.children[2]
    assert isinstance(detector, model.Element) and detector.surfaces[0].obscurations == []


def test_edit_form_adds_an_obscuration() -> None:
    ed = rt.Editor(rt.load(REFERENCE_DIR / "m0" / "singlet.rtt.json"))
    assert ed.system.json_at(STOP + "/obscurations") == []
    ed.insert(STOP + "/obscurations/0", {"type": "circular", "radius": 1.5})
    element = ed.system.root.children[0]
    assert isinstance(element, model.Element)
    (o,) = element.surfaces[0].obscurations
    assert isinstance(o, model.CircularObscuration) and o.radius.value == 1.5
    # The edit form writes every Param as an object; the canonical form only what differs.
    assert ed.system.json_at(STOP + "/obscurations/0/decenter_x") == {"value": 0.0,
                                                                      "variable": False}
    canonical = json.loads(ed.system.to_json())
    assert canonical["root"]["children"][0]["surfaces"][0]["obscurations"] == [
        {"type": "circular", "radius": 1.5}]
    ed.undo()
    assert ed.system.json_at(STOP + "/obscurations") == []


def test_layout_lists_the_obscurations() -> None:
    stop, mirror, detector = layout.surfaces(rt.compile(rt.load(TELESCOPE)))
    assert stop.obscurations == [
        ("circular", {"radius": 10.0, "decenter_x": 0.0, "decenter_y": 0.0}),
        ("rectangular", {"half_width_x": 31.0, "half_width_y": 0.5, "decenter_x": 0.0,
                         "decenter_y": 0.0, "rotation_deg": 0.0}),
        ("rectangular", {"half_width_x": 31.0, "half_width_y": 0.5, "decenter_x": 0.0,
                         "decenter_y": 0.0, "rotation_deg": 90.0})]
    assert mirror.obscurations == [
        ("elliptical", {"semi_axis_x": 3.0, "semi_axis_y": 1.5, "decenter_x": 20.0,
                        "decenter_y": 5.0, "rotation_deg": 30.0})]
    assert detector.obscurations == []
    # Obscurations do not cut the profile: the aperture stays the circle.
    assert mirror.aperture == ("circular", {"radius": 34.0, "inner_radius": 0.0})


def test_result_format_carries_the_obscurations() -> None:
    surfaces = layout.surfaces(rt.compile(rt.load(TELESCOPE)))
    text = surfaces[1].to_json()
    doc = json.loads(text)
    assert doc["schema_version"] == results.SCHEMA_VERSION == "0.1.22"
    assert doc["data"]["obscurations"] == [
        ["elliptical", {"semi_axis_x": 3.0, "semi_axis_y": 1.5, "decenter_x": 20.0,
                        "decenter_y": 5.0, "rotation_deg": 30.0}]]
    loaded = results.load_json(text).data
    assert loaded["obscurations"] == doc["data"]["obscurations"]
    # An older file has no key: no obscurations.
    old = dict(doc)
    old["schema_version"] = "0.1.21"
    old["data"] = {k: v for k, v in doc["data"].items() if k != "obscurations"}
    assert results.load_json(json.dumps(old)).data["obscurations"] == []


def test_the_spider_blocks_its_arms_and_the_centre_only() -> None:
    """End to end through m8/spider_telescope (ADR 0050, point 3).

    EPD 60 mm, object at infinity: pupil coordinate p maps to the height 30 p mm at the stop.
    - (0.5, 0): (15, 0) mm lies on the arm along x (|y| < 0.5): VIGNETTED at the stop.
    - (0, -0.6): (0, -18) mm lies on the arm along y (rotated by 90 deg): VIGNETTED at the stop.
    - (0.2, 0.2): (6, 6) mm lies in the centred circle (radius 10, r = 8.49) and on no arm
      (|x|, |y| = 6 > 0.5): VIGNETTED at the stop.
    - (0.5, 0.5): (15, 15) mm is free; it reaches the focus inside the +-2 mm square: ALIVE."""
    cs = rt.compile(rt.load(TELESCOPE))
    rays = rt.trace.make_rays(
        cs, rt.trace.PointListPupil([0.5, 0.0, 0.2, 0.5], [0.0, -0.6, 0.2, 0.5]))
    rt.trace.trace(cs, rays)
    v = rt.trace.RayStatus.VIGNETTED
    assert list(rays.status) == [v, v, v, rt.trace.RayStatus.ALIVE]
    ids = cs.surface_ids
    assert [ids[int(s)] for s in rays.last_surface] == ["STO", "STO", "STO", "F"]
    # The spot uses the virtual chief ray through the centre of the spider; the lost rays are
    # obscured, not clipped by the stop.
    spot = rt.analysis.spot(cs, field=0)
    assert not any(w.code == "stop.clips_beam" for w in spot.warnings)
    assert spot.rays_arrived > 0
