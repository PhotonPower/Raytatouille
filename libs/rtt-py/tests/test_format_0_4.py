"""Schema 0.4 from Python (#162): relative poses (ADR 0028), the parameter table and
configurations, bound Params with bounds (ADR 0029), and the warning io.pickup_dropped when a
0.3 file with pickups is loaded."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import model

ZOOM = """{
  "schema_version": "0.4.0",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [{"um": 0.5876, "reference": true}],
  "aperture": {"type": "epd", "value": 10.0},
  "fields": {"points": [{}]},
  "configurations": [{"name": "wide"}, {"name": "tele"}],
  "parameters": [
    {"name": "TOTAL", "value": 40.0},
    {"name": "G", "values": [20.0, 5.0], "variable": true, "min": 2.0, "max": 30.0},
    {"name": "B", "expression": "TOTAL - G"}
  ],
  "root": {"type": "assembly", "name": "zoom", "children": [
    {"type": "stop", "name": "stop",
     "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]},
    {"type": "detector", "name": "image",
     "pose": {"reference": "relative_to_sibling", "order": "rotate_first",
              "position": [0.0, 0.0, {"param": "B"}],
              "rotation_deg": [{"value": 1.0, "variable": true, "min": -2.0, "max": 2.0}, 0.0, 0.0]},
     "surfaces": [{"id": "IMG"}]}
  ]},
  "paths": [{"name": "main", "events": "auto"}]
}"""


def zoom() -> rt.System:
    return rt.System.from_json(ZOOM)


def test_configurations_and_rows() -> None:
    s = zoom()
    assert [c.name for c in s.configurations] == ["wide", "tele"]
    total, g, b = s.parameters
    assert (total.name, total.value, total.values, total.expression) == ("TOTAL", 40.0, None, None)
    assert (total.variable, total.min, total.max) == (False, None, None)
    assert (g.name, g.value, g.values, g.expression) == ("G", None, [20.0, 5.0], None)
    assert (g.variable, g.min, g.max) == (True, 2.0, 30.0)
    assert (b.name, b.value, b.values, b.expression) == ("B", None, None, "TOTAL - G")
    assert not b.variable


def test_pose_reference_and_order() -> None:
    pose = zoom().root.children[1].pose
    assert pose.reference == model.PoseReference.RELATIVE_TO_SIBLING
    assert pose.order == model.PoseOrder.ROTATE_FIRST
    stop = zoom().root.children[0].pose
    assert stop.reference == model.PoseReference.ABSOLUTE
    assert stop.order == model.PoseOrder.TRANSLATE_FIRST


def test_bound_param_has_no_value() -> None:
    # ADR 0029, point 3: the value of a bound Param is None in Python.
    z = zoom().root.children[1].pose.position[2]
    assert (z.param, z.value, z.variable, z.min, z.max) == ("B", None, False, None, None)
    assert "param='B'" in repr(z)
    tilt = zoom().root.children[1].pose.rotation_deg[0]
    assert (tilt.param, tilt.value, tilt.variable, tilt.min, tilt.max) == (None, 1.0, True,
                                                                            -2.0, 2.0)
    assert repr(tilt) == "Param(value=1.0, variable=True, min=-2.0, max=2.0)"


def test_pickup_is_gone_from_param() -> None:
    assert not hasattr(zoom().root.children[1].pose.position[0], "pickup")


def test_round_trip_and_edit_form() -> None:
    s = zoom()
    assert rt.System.from_json(s.to_json()) == s
    data = s.to_dict()
    assert data["root"]["children"][1]["pose"]["position"][2] == {"param": "B"}
    assert data["parameters"][2] == {"name": "B", "expression": "TOTAL - G", "variable": False}
    assert data["root"]["children"][0]["pose"]["reference"] == "absolute"


def test_loading_a_0_3_file_warns_about_dropped_pickups(tmp_path: Path) -> None:
    # Migration 0.3 -> 0.4 (ADR 0029, point 6): the pickup goes, the value stays, and each
    # dropped pickup is a RaytatouilleWarning io.pickup_dropped at its pointer.
    data = json.loads(ZOOM)
    data["schema_version"] = "0.3.0"
    for key in ("configurations", "parameters"):
        del data[key]
    image = data["root"]["children"][1]
    image["pose"] = {"position": [0.0, 0.0, {"value": 6.0, "pickup": "2 * 3"}]}
    text = json.dumps(data)
    for load in (lambda: rt.System.from_json(text), lambda: rt.load(write(tmp_path, text))):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            s = load()
        found = [w.message for w in caught if isinstance(w.message, rt.RaytatouilleWarning)]
        assert len(found) == 1
        assert found[0].code == "io.pickup_dropped"
        assert found[0].location == "/root/children/1/pose/position/2/pickup"
        assert "2 * 3" in found[0].message
        z = s.root.children[1].pose.position[2]
        assert (z.value, z.param) == (6.0, None)
        assert s.schema_version == "0.4.0"


def test_a_0_4_file_loads_without_warnings() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")
        zoom()


def write(directory: Path, text: str) -> Path:
    path = directory / "pickup_0_3.rtt.json"
    path.write_text(text, encoding="utf-8")
    return path


def test_pickup_in_a_0_4_file_is_an_error() -> None:
    data = json.loads(ZOOM)
    data["root"]["children"][1]["pose"]["position"][0] = {"value": 1.0, "pickup": "x"}
    with pytest.raises(rt.ParseError) as info:
        rt.System.from_json(json.dumps(data))
    assert info.value.pointer == "/root/children/1/pose/position/0/pickup"
