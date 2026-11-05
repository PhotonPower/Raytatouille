"""Anamorphic base shapes in Python (ADR 0042, points 1 and 5; #315): the read-only classes
Biconic, Toroid and Cylinder, their edit form, the layout fields curvature_x and
conic_constant_x with the result format 0.1.17, and the examples under examples/anamorphic.

Oracle: the parameters written into the test systems. The layout gives the curvature c = 1/R
of the y section (biconic) or of the y-z profile (toroid, cylinder) and, in curvature_x, the
curvature of the x section through the vertex: c_x = 1/R_x (biconic), c_R = 1/R of the
rotation (toroid) and 0 (cylinder); for rotationally symmetric shapes it equals curvature.
conic_constant_x is k_x of the biconic, 0 for the toroid (its x section through the vertex is
the circle of rotation), 0 for the cylinder (a straight line) and conic_constant otherwise.

No test traces these systems: aiming needs the paraxial data per section, which follow with
the second part of #315 (ADR 0042, points 2 and 3).
"""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]  # no py.typed; test tool (ADR 0018)
import numpy as np
import pytest
from conftest import REPO_ROOT

import raytatouille as rt
from raytatouille import model

SYSTEM_SCHEMA = json.loads(
    (REPO_ROOT / "schema" / "raytatouille.schema.json").read_text(encoding="utf-8"))
RESULT_SCHEMA = json.loads(
    (REPO_ROOT / "schema" / "raytatouille-result.schema.json").read_text(encoding="utf-8"))
EXAMPLE_DIR = REPO_ROOT / "examples" / "anamorphic"

#: The first surface of lens Li has BASES[i]; the second surface is plane.
BASES: list[dict[str, Any]] = [
    {"type": "plane"},
    {"type": "conic", "radius": 40.0, "conic": -0.5},
    {"type": "even_asphere", "radius": -30.0, "conic": 0.25, "coefficients": [1e-5, -2e-8]},
    {"type": "biconic", "radius_x": -60.0, "radius_y": 40.0, "conic_x": -0.5, "conic_y": 0.5},
    {"type": "toroid", "radius": 50.0, "conic": -0.3, "coefficients": [1e-6],
     "radius_of_rotation": -80.0},
    {"type": "cylinder", "radius": 25.0, "conic": -2.0, "coefficients": [3e-6, 1e-9]},
]


def system_text(bases: list[dict[str, Any]]) -> str:
    """VACUUM, a stop at z = 0 and one CONST:1.5 lens per base at z = 10 (i + 1), 2 mm thick,
    the image at z = 10 (len + 1) + 20."""
    lenses = [{
        "type": "lens", "name": f"L{i}", "pose": {"position": [0.0, 0.0, 10.0 * (i + 1)]},
        "material": "CONST:1.5",
        "surfaces": [
            {"id": f"L{i}.S1", "shape": {"base": base},
             "aperture": {"type": "circular", "radius": 5.0}},
            {"id": f"L{i}.S2", "pose": {"position": [0.0, 0.0, 2.0]},
             "aperture": {"type": "circular", "radius": 5.0}},
        ]} for i, base in enumerate(bases)]
    data = {
        "schema_version": "0.9.0",
        "units": {"length": "mm", "wavelength": "um"},
        "environment": {"medium": "VACUUM"},
        "wavelengths": [{"um": 0.5876, "reference": True}],
        "aperture": {"type": "epd", "value": 4.0},
        "fields": {"points": [{}]},
        "root": {"type": "assembly", "name": "root", "children": [
            {"type": "stop", "name": "stop",
             "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 5.0}}]},
            *lenses,
            {"type": "detector", "name": "image",
             "pose": {"position": [0.0, 0.0, 10.0 * (len(bases) + 1) + 20.0]},
             "surfaces": [{"id": "IMG"}]}]},
        "paths": [{"name": "main", "events": "auto"}],
    }
    return json.dumps(data)


def base_of(system: rt.System, i: int) -> Any:
    lens = system.root.children[i + 1]
    assert isinstance(lens, model.Element) and lens.name == f"L{i}"
    return lens.surfaces[0].shape.base


def attributes(obj: Any) -> set[str]:
    return {name for name in dir(type(obj))
            if not name.startswith("_") and isinstance(getattr(type(obj), name), property)}


def test_the_new_classes_are_read() -> None:
    system = rt.System.from_json(system_text(BASES))
    biconic, toroid, cylinder = (base_of(system, i) for i in (3, 4, 5))
    assert isinstance(biconic, model.Biconic)
    assert (biconic.radius_x.value, biconic.radius_y.value, biconic.conic_x.value,
            biconic.conic_y.value) == (-60.0, 40.0, -0.5, 0.5)
    assert isinstance(toroid, model.Toroid)
    assert (toroid.radius.value, toroid.conic.value, [a.value for a in toroid.coefficients],
            toroid.radius_of_rotation.value) == (50.0, -0.3, [1e-6], -80.0)
    assert isinstance(cylinder, model.Cylinder)
    assert (cylinder.radius.value, cylinder.conic.value,
            [a.value for a in cylinder.coefficients]) == (25.0, -2.0, [3e-6, 1e-9])
    # Exactly the parameters of the file format (ADR 0042, point 1), as read-only copies.
    assert attributes(biconic) == {"radius_x", "radius_y", "conic_x", "conic_y"}
    assert attributes(toroid) == {"radius", "conic", "coefficients", "radius_of_rotation"}
    assert attributes(cylinder) == {"radius", "conic", "coefficients"}
    for obj in (biconic, toroid, cylinder):
        for name in attributes(obj):
            with pytest.raises(AttributeError):
                setattr(obj, name, getattr(obj, name))
        with pytest.raises(TypeError):
            type(obj)()
    assert biconic == base_of(rt.System.from_json(system_text(BASES)), 3)
    assert biconic != base_of(system, 1)
    # Exported by raytatouille.model.
    assert {"Biconic", "Toroid", "Cylinder"} <= set(model.__all__)


def test_missing_optional_keys_have_the_parser_defaults() -> None:
    system = rt.System.from_json(system_text([
        {"type": "biconic", "radius_x": 10.0, "radius_y": 20.0},
        {"type": "toroid", "radius": 10.0, "radius_of_rotation": 30.0},
        {"type": "cylinder", "radius": 10.0}]))
    biconic, toroid, cylinder = (base_of(system, i) for i in range(3))
    assert (biconic.conic_x.value, biconic.conic_y.value) == (0.0, 0.0)
    assert (toroid.conic.value, toroid.coefficients) == (0.0, [])
    assert (cylinder.conic.value, cylinder.coefficients) == (0.0, [])


def test_the_edit_form_has_every_key_and_round_trips() -> None:
    system = rt.System.from_json(system_text(BASES))
    data = system.to_dict()
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(data)
    assert rt.System.from_json(json.dumps(data)) == system
    lenses = data["root"]["children"][1:-1]
    keys = [set(lens["surfaces"][0]["shape"]["base"]) for lens in lenses]
    assert keys[3] == {"type", "radius_x", "radius_y", "conic_x", "conic_y"}
    assert keys[4] == {"type", "radius", "conic", "coefficients", "radius_of_rotation"}
    assert keys[5] == {"type", "radius", "conic", "coefficients"}
    toroid = lenses[4]["surfaces"][0]["shape"]["base"]
    assert toroid["radius_of_rotation"] == {"value": -80.0, "variable": False}
    assert toroid["coefficients"] == [{"value": 1e-6, "variable": False}]


def test_a_patch_changes_the_new_parameters() -> None:
    system = rt.System.from_json(system_text(BASES))
    toroid = system.locate_surface("L4.S1")
    biconic = system.locate_surface("L3.S1")
    assert toroid is not None and biconic is not None
    changed = rt.apply_patch(system, [
        {"op": "replace", "path": toroid + "/shape/base/radius_of_rotation/value",
         "value": 100.0},
        {"op": "replace", "path": biconic + "/shape/base/conic_x/variable", "value": True}])
    assert base_of(changed, 4).radius_of_rotation.value == 100.0
    assert base_of(changed, 3).conic_x.variable
    assert base_of(system, 4).radius_of_rotation.value == -80.0  # the original is unchanged
    layout = {s.id: s for s in rt.layout.surfaces(rt.compile(changed))}
    assert layout["L4.S1"].curvature_x == pytest.approx(0.01, rel=1e-15)


def layouts() -> dict[str, rt.layout.SurfaceLayout]:
    cs = rt.compile(rt.System.from_json(system_text(BASES)))
    return {s.id: s for s in rt.layout.surfaces(cs)}


def test_layout_of_the_rotationally_symmetric_shapes_has_x_equal_y() -> None:
    by_id = layouts()
    for sid, shape in (("STO", "plane"), ("L0.S1", "plane"), ("L1.S1", "conic"),
                       ("L2.S1", "even_asphere"), ("IMG", "plane")):
        s = by_id[sid]
        assert s.shape == shape, sid
        assert s.curvature_x == s.curvature and s.conic_constant_x == s.conic_constant, sid
    conic, asphere = by_id["L1.S1"], by_id["L2.S1"]
    assert conic.curvature == pytest.approx(1.0 / 40.0, rel=1e-15)
    assert conic.conic_constant == -0.5 and conic.coefficients == []
    assert asphere.curvature == pytest.approx(-1.0 / 30.0, rel=1e-15)
    assert (asphere.conic_constant, asphere.coefficients) == (0.25, [1e-5, -2e-8])
    assert by_id["L0.S1"].curvature == 0.0 and by_id["L0.S1"].conic_constant == 0.0


def test_layout_of_the_biconic() -> None:
    s = layouts()["L3.S1"]
    assert s.shape == "biconic"
    # The y section in curvature and conic_constant, the x section in the new fields.
    assert s.curvature == pytest.approx(1.0 / 40.0, rel=1e-15)
    assert s.conic_constant == 0.5
    assert s.curvature_x == pytest.approx(-1.0 / 60.0, rel=1e-15)
    assert s.conic_constant_x == -0.5
    assert s.coefficients == []
    # Domain a x^2 + b y^2 <= 1 with a = (1 + k_x) c_x^2, b = (1 + k_y) c_y^2; max_radius is
    # the outer extent 1/sqrt(min(a, b)) (ADR 0042, point 1).
    a = 0.5 / 60.0**2
    b = 1.5 / 40.0**2
    assert s.max_radius is not None
    assert s.max_radius == pytest.approx(1.0 / math.sqrt(min(a, b)), rel=1e-12)


def test_layout_of_the_toroid() -> None:
    s = layouts()["L4.S1"]
    assert s.shape == "toroid"
    assert s.curvature == pytest.approx(1.0 / 50.0, rel=1e-15)
    assert (s.conic_constant, s.coefficients) == (-0.3, [1e-6])
    assert s.curvature_x == pytest.approx(-1.0 / 80.0, rel=1e-15)  # c_R = 1/radius_of_rotation
    assert s.conic_constant_x == 0.0  # the circle of rotation
    # Every point of the domain has |x| <= |d| with d = R - z_y(y), and d = R at the vertex:
    # the outer extent is at least |R| = 80 mm.
    assert s.max_radius is not None and s.max_radius >= 80.0


def test_layout_of_the_cylinder() -> None:
    s = layouts()["L5.S1"]
    assert s.shape == "cylinder"
    assert s.curvature == pytest.approx(1.0 / 25.0, rel=1e-15)
    assert (s.conic_constant, s.coefficients) == (-2.0, [3e-6, 1e-9])
    assert s.curvature_x == 0.0
    assert s.conic_constant_x == 0.0  # a straight line
    assert s.max_radius is None  # a strip, unbounded in x (ADR 0042, point 1)


def test_the_profiles_follow_the_sections() -> None:
    # Content guard of the layout: the x-z and y-z profiles of the cylinder differ (flat in x),
    # those of the conic agree.
    cs = rt.compile(rt.System.from_json(system_text(BASES)))
    for sid, flat_in_x in (("L5.S1", True), ("L1.S1", False)):
        (yz,) = rt.layout.profile(cs, sid, "yz", samples=11)
        (xz,) = rt.layout.profile(cs, sid, "xz", samples=11)
        z0 = float(yz[5, 2])
        assert float(np.max(yz[:, 2])) - z0 > 0.1, sid  # curved in y
        rise_x = float(np.max(np.abs(xz[:, 2] - z0)))
        assert (rise_x == 0.0) if flat_in_x else (rise_x > 0.1), sid


def test_the_surface_layout_result_has_the_new_keys_at_the_end() -> None:
    validator = jsonschema.Draft202012Validator(RESULT_SCHEMA)
    for s in layouts().values():
        text = s.to_json()
        envelope = json.loads(text)
        assert envelope["schema_version"] == rt.results.SCHEMA_VERSION == "0.1.17"
        assert list(envelope["data"])[-2:] == ["curvature_x", "conic_constant_x"]
        validator.validate(envelope)
        loaded = rt.results.load_json(text)
        assert loaded.data["curvature_x"] == s.curvature_x
        assert loaded.data["conic_constant_x"] == s.conic_constant_x
        assert loaded.data["shape"] == s.shape


def test_an_older_surface_layout_reads_with_x_equal_y() -> None:
    # Files before 0.1.17 only have plane, conic and even_asphere: their x section is the y
    # section (ADR 0023, addendum #315).
    validator = jsonschema.Draft202012Validator(RESULT_SCHEMA)
    envelope = json.loads(layouts()["L2.S1"].to_json())
    for key in ("curvature_x", "conic_constant_x"):
        del envelope["data"][key]
    for version in ("0.1.0", "0.1.16"):
        old = copy.deepcopy(envelope)
        old["schema_version"] = version
        validator.validate(old)
        loaded = rt.results.load_json(json.dumps(old))
        assert loaded.data["curvature_x"] == loaded.data["curvature"]
        assert loaded.data["conic_constant_x"] == loaded.data["conic_constant"] == 0.25
        assert list(loaded.data)[-2:] == ["curvature_x", "conic_constant_x"]
    for version in ("0.1.17", "0.1.20"):
        new = copy.deepcopy(envelope)
        new["schema_version"] = version
        with pytest.raises(jsonschema.ValidationError):
            validator.validate(new)
        with pytest.raises(ValueError, match="misses curvature_x, conic_constant_x"):
            rt.results.load_json(json.dumps(new))


EXAMPLES = sorted(EXAMPLE_DIR.glob("*.rtt.json"))


def test_the_examples_are_there() -> None:
    assert [f.name for f in EXAMPLES] == ["biconic_singlet.rtt.json",
                                          "crossed_cylinders.rtt.json",
                                          "cylinder_planoconvex.rtt.json"]


@pytest.mark.parametrize("file", EXAMPLES, ids=lambda f: f.name)
def test_an_example_is_canonical_and_valid(file: Path) -> None:
    system = rt.load(file)
    assert system.to_json().encode("utf-8") == file.read_bytes()  # as `rtt format` writes it
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(
        json.loads(file.read_text(encoding="utf-8")))
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(system.to_dict())
    assert [d for d in rt.validate(system) if d.severity == rt.Severity.ERROR] == []
    assert rt.System.from_json(json.dumps(system.to_dict())) == system


def test_the_example_shapes() -> None:
    def layout(name: str) -> dict[str, rt.layout.SurfaceLayout]:
        return {s.id: s for s in rt.layout.surfaces(rt.compile(rt.load(EXAMPLE_DIR / name)))}

    n = 1.5168
    cylinder = layout("cylinder_planoconvex.rtt.json")["L1.S1"]
    assert (cylinder.shape, cylinder.curvature_x) == ("cylinder", 0.0)
    assert 1.0 / cylinder.curvature / (n - 1.0) == pytest.approx(100.0, rel=1e-12)  # f_y

    crossed = layout("crossed_cylinders.rtt.json")
    first, second = crossed["L1.S1"], crossed["L2.S1"]
    assert first.shape == second.shape == "cylinder"
    assert 1.0 / second.curvature / (n - 1.0) == pytest.approx(75.0, rel=1e-12)
    # L2 is turned by 90 degree about z: its local y axis (the profile) is the global -x axis.
    assert np.max(np.abs(first.rotation - np.eye(3))) == 0.0
    assert np.max(np.abs(second.rotation[:, 1] - np.array([-1.0, 0.0, 0.0]))) <= 1e-15
    # Thin-lens check of the design: both sections focus at the image (z = 106.363 mm).
    image = crossed["IMG"].translation[2]
    assert 9.0 + 100.0 - 4.0 / n == pytest.approx(image, abs=1e-12)
    assert 34.0 + 75.0 - 4.0 / n == pytest.approx(image, abs=1e-12)

    biconic = layout("biconic_singlet.rtt.json")["L1.S1"]
    assert biconic.shape == "biconic"
    assert (biconic.conic_constant, biconic.conic_constant_x) == (0.0, -0.5)
    assert 1.0 / biconic.curvature / (n - 1.0) == pytest.approx(100.0, rel=1e-12)
    assert 1.0 / biconic.curvature_x / (n - 1.0) == pytest.approx(150.0, rel=1e-12)
