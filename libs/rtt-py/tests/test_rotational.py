"""Odd asphere and axicon in Python (ADR 0043, points 1 and 4; #314): the read-only classes
OddAsphere and Axicon, their edit form, the layout fields normalization_radius and
axicon_angle_deg with the result format 0.1.18, and the examples under examples/rotational.

Oracle: the parameters written into the test systems. The layout gives the curvature c and the
conic constant k of the base (for the odd asphere the Param ``curvature`` itself, not 1/R; for
the axicon 0 and 0), the coefficients A1, A2, ... of the odd asphere, and the base angle of the
axicon in degrees; normalization_radius is None for both (it belongs to Q-con and Q-bfs).
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
EXAMPLE_DIR = REPO_ROOT / "examples" / "rotational"

#: The first surface of lens Li has BASES[i]; the second surface is plane.
BASES: list[dict[str, Any]] = [
    {"type": "conic", "radius": 40.0, "conic": -0.5},
    {"type": "odd_asphere", "curvature": 0.02, "conic": -0.5,
     "coefficients": [1e-3, 1e-4, -2e-6, 1e-7]},
    {"type": "axicon", "angle_deg": -5.0},
    {"type": "odd_asphere", "curvature": 0.0, "coefficients": [0.05]},
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
    odd, axicon = base_of(system, 1), base_of(system, 2)
    assert isinstance(odd, model.OddAsphere)
    assert (odd.curvature.value, odd.conic.value, [a.value for a in odd.coefficients]) == (
        0.02, -0.5, [1e-3, 1e-4, -2e-6, 1e-7])
    assert isinstance(axicon, model.Axicon)
    assert axicon.angle_deg.value == -5.0
    # Exactly the parameters of the file format (ADR 0043, point 4), as read-only copies.
    assert attributes(odd) == {"curvature", "conic", "coefficients"}
    assert attributes(axicon) == {"angle_deg"}
    for obj in (odd, axicon):
        for name in attributes(obj):
            with pytest.raises(AttributeError):
                setattr(obj, name, getattr(obj, name))
        with pytest.raises(TypeError):
            type(obj)()
    assert odd == base_of(rt.System.from_json(system_text(BASES)), 1)
    assert odd != base_of(system, 3)
    assert {"OddAsphere", "Axicon"} <= set(model.__all__)


def test_missing_optional_keys_have_the_parser_defaults() -> None:
    system = rt.System.from_json(system_text([{"type": "odd_asphere", "curvature": 0.01}]))
    odd = base_of(system, 0)
    assert (odd.curvature.value, odd.conic.value, odd.coefficients) == (0.01, 0.0, [])


def test_the_edit_form_has_every_key_and_round_trips() -> None:
    system = rt.System.from_json(system_text(BASES))
    data = system.to_dict()
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(data)
    assert rt.System.from_json(json.dumps(data)) == system
    lenses = data["root"]["children"][1:-1]
    keys = [set(lens["surfaces"][0]["shape"]["base"]) for lens in lenses]
    assert keys[1] == {"type", "curvature", "conic", "coefficients"}
    assert keys[2] == {"type", "angle_deg"}
    odd = lenses[3]["surfaces"][0]["shape"]["base"]
    assert odd["curvature"] == {"value": 0.0, "variable": False}
    assert odd["conic"] == {"value": 0.0, "variable": False}


def test_a_patch_changes_the_new_parameters() -> None:
    system = rt.System.from_json(system_text(BASES))
    odd = system.locate_surface("L1.S1")
    axicon = system.locate_surface("L2.S1")
    assert odd is not None and axicon is not None
    changed = rt.apply_patch(system, [
        {"op": "replace", "path": odd + "/shape/base/curvature/value", "value": -0.01},
        {"op": "replace", "path": axicon + "/shape/base/angle_deg/variable", "value": True}])
    assert base_of(changed, 1).curvature.value == -0.01
    assert base_of(changed, 2).angle_deg.variable
    assert base_of(system, 1).curvature.value == 0.02  # the original is unchanged
    layout = {s.id: s for s in rt.layout.surfaces(rt.compile(changed))}
    assert layout["L1.S1"].curvature == -0.01


def layouts() -> dict[str, rt.layout.SurfaceLayout]:
    cs = rt.compile(rt.System.from_json(system_text(BASES)))
    return {s.id: s for s in rt.layout.surfaces(cs)}


def test_layout_of_the_odd_asphere() -> None:
    s = layouts()["L1.S1"]
    assert s.shape == "odd_asphere"
    # The curvature of the base, not the vertex curvature c + 2 A2 (ADR 0043, point 4).
    assert (s.curvature, s.conic_constant) == (0.02, -0.5)
    assert s.coefficients == [1e-3, 1e-4, -2e-6, 1e-7]
    assert (s.curvature_x, s.conic_constant_x) == (s.curvature, s.conic_constant)
    assert s.normalization_radius is None and s.axicon_angle_deg is None
    # The domain of the base conic: 1 / (|c| sqrt(1 + k)).
    assert s.max_radius == pytest.approx(1.0 / (0.02 * math.sqrt(0.5)), rel=1e-15)
    plane_base = layouts()["L3.S1"]
    assert (plane_base.shape, plane_base.curvature, plane_base.max_radius) == (
        "odd_asphere", 0.0, None)


def test_layout_of_the_axicon() -> None:
    s = layouts()["L2.S1"]
    assert s.shape == "axicon"
    assert (s.curvature, s.conic_constant, s.coefficients) == (0.0, 0.0, [])
    assert (s.curvature_x, s.conic_constant_x) == (0.0, 0.0)
    assert s.axicon_angle_deg == -5.0
    assert s.normalization_radius is None
    assert s.max_radius is None
    # The sag is r tan(alpha): the profile is a cone (ADR 0043, point 1).
    cs = rt.compile(rt.System.from_json(system_text(BASES)))
    (yz,) = rt.layout.profile(cs, "L2.S1", "yz", samples=11)
    local_z = yz[:, 2] - s.translation[2]
    assert np.max(np.abs(local_z - math.tan(math.radians(-5.0)) * np.abs(yz[:, 1]))) <= 1e-14


def test_the_other_shapes_have_no_new_values() -> None:
    by_id = layouts()
    for sid in ("STO", "L0.S1", "L0.S2", "IMG"):
        assert by_id[sid].normalization_radius is None, sid
        assert by_id[sid].axicon_angle_deg is None, sid


def test_the_surface_layout_result_has_the_new_keys() -> None:
    validator = jsonschema.Draft202012Validator(RESULT_SCHEMA)
    for s in layouts().values():
        text = s.to_json()
        envelope = json.loads(text)
        assert envelope["schema_version"] == rt.results.SCHEMA_VERSION == "0.1.18"
        assert {"normalization_radius", "axicon_angle_deg"} <= set(envelope["data"])
        validator.validate(envelope)
        loaded = rt.results.load_json(text)
        assert loaded.data["shape"] == s.shape
        assert loaded.data["axicon_angle_deg"] == s.axicon_angle_deg
        assert loaded.data["normalization_radius"] is None


def test_an_older_surface_layout_reads_without_the_new_values() -> None:
    # Files before 0.1.18 have no odd asphere, Q type or axicon: both keys read as None (ADR
    # 0023, addendum #314); from 0.1.18 on they are required.
    validator = jsonschema.Draft202012Validator(RESULT_SCHEMA)
    envelope = json.loads(layouts()["L0.S1"].to_json())
    for key in ("normalization_radius", "axicon_angle_deg"):
        del envelope["data"][key]
    for version in ("0.1.0", "0.1.16", "0.1.17"):
        old = copy.deepcopy(envelope)
        old["schema_version"] = version
        if version != "0.1.17":
            for key in ("curvature_x", "conic_constant_x"):
                del old["data"][key]
        validator.validate(old)
        loaded = rt.results.load_json(json.dumps(old))
        assert loaded.data["normalization_radius"] is None
        assert loaded.data["axicon_angle_deg"] is None
        assert loaded.data["curvature_x"] == loaded.data["curvature"]
    for version in ("0.1.18", "0.1.20"):
        new = copy.deepcopy(envelope)
        new["schema_version"] = version
        with pytest.raises(jsonschema.ValidationError):
            validator.validate(new)
        with pytest.raises(ValueError, match="misses normalization_radius, axicon_angle_deg"):
            rt.results.load_json(json.dumps(new))


def test_the_axicon_and_the_cone_term_have_no_paraxial_model() -> None:
    # ADR 0043, point 3 (F6): first_order raises at the surface.
    for bases in ([BASES[2]], [{"type": "odd_asphere", "curvature": 0.02,
                                "coefficients": [1e-3]}]):
        cs = rt.compile(rt.System.from_json(system_text(bases)))
        with pytest.raises(rt.ParaxialError, match="no paraxial model at the vertex"):
            rt.paraxial.first_order(cs)


EXAMPLES = sorted(EXAMPLE_DIR.glob("*.rtt.json"))


def test_the_examples_are_there() -> None:
    assert [f.name for f in EXAMPLES] == ["axicon_planoconvex.rtt.json",
                                          "odd_asphere_singlet.rtt.json"]


@pytest.mark.parametrize("file", EXAMPLES, ids=lambda f: f.name)
def test_an_example_is_canonical_and_valid(file: Path) -> None:
    system = rt.load(file)
    assert system.to_json().encode("utf-8") == file.read_bytes()  # as `rtt format` writes it
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(
        json.loads(file.read_text(encoding="utf-8")))
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(system.to_dict())
    assert [d for d in rt.validate(system) if d.severity == rt.Severity.ERROR] == []
    assert rt.System.from_json(json.dumps(system.to_dict())) == system


def test_the_axicon_example_bends_by_the_cone_angle() -> None:
    # ADR 0043, acceptance 8, from Python: an axis-parallel ray at h = 2 mm leaves the cone at
    # beta = asin(1.5 sin 5 deg) - 5 deg towards the axis.
    cs = rt.compile(rt.load(EXAMPLE_DIR / "axicon_planoconvex.rtt.json"))
    rays = rt.trace.RayBatch(1)
    rays.pos_y[0] = 2.0
    rays.pos_z[0] = -1.0
    rays.dir_z[0] = 1.0
    _, paths = rt.trace.trace(cs, rays, record_path=True)
    assert int(rays.status[0]) == int(rt.trace.RayStatus.ALIVE)
    d = paths.direction[0, 3]  # after L1.S2
    beta = math.asin(1.5 * math.sin(math.radians(5.0))) - math.radians(5.0)
    assert math.atan2(-d[1], d[2]) == pytest.approx(beta, abs=1e-13)


def test_the_odd_asphere_example_has_the_vertex_curvature_c_plus_2_a2() -> None:
    # ADR 0043, point 3: c_v = c + 2 A2 = 0.0202 / mm; the thick-lens power of a plano-convex
    # lens with the curved side first is (n - 1) c_v.
    cs = rt.compile(rt.load(EXAMPLE_DIR / "odd_asphere_singlet.rtt.json"))
    fo = rt.paraxial.first_order(cs)
    assert fo.efl is not None
    assert 1.0 / fo.efl == pytest.approx(0.5168 * (0.02 + 2e-4), rel=1e-12)
    with pytest.raises(rt.ParaxialError, match="no expansion to fourth order"):
        rt.paraxial.seidel(cs)
