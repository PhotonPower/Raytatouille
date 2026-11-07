"""Forbes Q-con and Q-bfs in Python (ADR 0043, points 1 and 4; #314): the read-only classes QCon
and QBfs, their edit form, the layout with normalization_radius in the result format 0.1.18, and
the examples qcon_singlet and qbfs_singlet under examples/forbes.

Oracle: the parameters written into the test systems. The layout gives the curvature c and the
conic constant k of the base (not the vertex curvature of Q-bfs), the coefficients a_0, a_1, ...
in mm and the normalization radius rho; axicon_angle_deg is None.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]  # no py.typed; test tool (ADR 0018)
import pytest
from conftest import REPO_ROOT

import raytatouille as rt
from raytatouille import model

SYSTEM_SCHEMA = json.loads(
    (REPO_ROOT / "schema" / "raytatouille.schema.json").read_text(encoding="utf-8"))
RESULT_SCHEMA = json.loads(
    (REPO_ROOT / "schema" / "raytatouille-result.schema.json").read_text(encoding="utf-8"))
EXAMPLE_DIR = REPO_ROOT / "examples" / "forbes"

#: The first surface of lens Li has BASES[i]; the second surface is plane.
BASES: list[dict[str, Any]] = [
    {"type": "qcon", "curvature": 0.02, "conic": -0.5, "normalization_radius": 5.0,
     "coefficients": [1e-3, -2e-4, 5e-5]},
    {"type": "qbfs", "curvature": -0.01, "conic": 0.25, "normalization_radius": 4.0,
     "coefficients": [3e-3, 1e-3]},
    {"type": "qcon", "curvature": 0.0, "normalization_radius": 5.0, "coefficients": [0.0]},
]


def system_text(bases: list[dict[str, Any]]) -> str:
    """VACUUM, a stop at z = 0 and one CONST:1.5 lens per base at z = 10 (i + 1), 2 mm thick,
    aperture radius 4 mm, the image at z = 10 (len + 1) + 20."""
    lenses = [{
        "type": "lens", "name": f"L{i}", "pose": {"position": [0.0, 0.0, 10.0 * (i + 1)]},
        "material": "CONST:1.5",
        "surfaces": [
            {"id": f"L{i}.S1", "shape": {"base": base},
             "aperture": {"type": "circular", "radius": 4.0}},
            {"id": f"L{i}.S2", "pose": {"position": [0.0, 0.0, 2.0]},
             "aperture": {"type": "circular", "radius": 4.0}},
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
    qcon, qbfs = base_of(system, 0), base_of(system, 1)
    assert isinstance(qcon, model.QCon)
    assert (qcon.curvature.value, qcon.conic.value, qcon.normalization_radius.value,
            [a.value for a in qcon.coefficients]) == (0.02, -0.5, 5.0, [1e-3, -2e-4, 5e-5])
    assert isinstance(qbfs, model.QBfs)
    assert (qbfs.curvature.value, qbfs.conic.value, qbfs.normalization_radius.value,
            [a.value for a in qbfs.coefficients]) == (-0.01, 0.25, 4.0, [3e-3, 1e-3])
    for obj in (qcon, qbfs):
        assert attributes(obj) == {"curvature", "conic", "normalization_radius", "coefficients"}
        for name in attributes(obj):
            with pytest.raises(AttributeError):
                setattr(obj, name, getattr(obj, name))
        with pytest.raises(TypeError):
            type(obj)()
    assert {"QCon", "QBfs"} <= set(model.__all__)
    assert base_of(system, 2).conic.value == 0.0  # the parser default


def test_the_edit_form_has_every_key_and_round_trips() -> None:
    system = rt.System.from_json(system_text(BASES))
    data = system.to_dict()
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(data)
    assert rt.System.from_json(json.dumps(data)) == system
    for lens in data["root"]["children"][1:-1]:
        base = lens["surfaces"][0]["shape"]["base"]
        assert list(base) == ["type", "curvature", "conic", "normalization_radius", "coefficients"]


def test_a_patch_changes_the_normalization_radius() -> None:
    system = rt.System.from_json(system_text(BASES))
    where = system.locate_surface("L1.S1")
    assert where is not None
    changed = rt.apply_patch(system, [
        {"op": "replace", "path": where + "/shape/base/normalization_radius/value", "value": 6.0}])
    assert base_of(changed, 1).normalization_radius.value == 6.0
    layout = {s.id: s for s in rt.layout.surfaces(rt.compile(changed))}
    assert layout["L1.S1"].normalization_radius == 6.0


def layouts() -> dict[str, rt.layout.SurfaceLayout]:
    cs = rt.compile(rt.System.from_json(system_text(BASES)))
    return {s.id: s for s in rt.layout.surfaces(cs)}


def test_layout_of_qcon_and_qbfs() -> None:
    by_id = layouts()
    qcon, qbfs = by_id["L0.S1"], by_id["L1.S1"]
    assert (qcon.shape, qcon.curvature, qcon.conic_constant) == ("qcon", 0.02, -0.5)
    assert qcon.coefficients == [1e-3, -2e-4, 5e-5]
    assert qcon.normalization_radius == 5.0
    assert (qbfs.shape, qbfs.curvature, qbfs.conic_constant) == ("qbfs", -0.01, 0.25)
    assert qbfs.coefficients == [3e-3, 1e-3]  # a_m, not the auxiliary b_m
    assert qbfs.normalization_radius == 4.0
    for s in (qcon, qbfs):
        assert (s.curvature_x, s.conic_constant_x) == (s.curvature, s.conic_constant)
        assert s.axicon_angle_deg is None
    # The domain of the base conic: 1 / (|c| sqrt(1 + k)).
    assert qbfs.max_radius == pytest.approx(1.0 / (0.01 * math.sqrt(1.25)), rel=1e-15)
    assert by_id["L2.S1"].max_radius is None


def test_the_surface_layout_result_carries_the_normalization_radius() -> None:
    validator = jsonschema.Draft202012Validator(RESULT_SCHEMA)
    s = layouts()["L1.S1"]
    envelope = json.loads(s.to_json())
    assert envelope["schema_version"] == rt.results.SCHEMA_VERSION == "0.1.18"
    assert envelope["data"]["normalization_radius"] == 4.0
    validator.validate(envelope)
    assert rt.results.load_json(s.to_json()).data["normalization_radius"] == 4.0


def test_an_aperture_beyond_the_normalization_radius_warns() -> None:
    # ADR 0043, point 2 (F3): aperture radius 4 mm, rho 5 and 4 mm: no warning; rho 3 mm: one.
    assert not [d for d in rt.validate(rt.System.from_json(system_text(BASES)))
                if d.code == "shape.normalization_radius_exceeded"]
    bases = [{**BASES[0], "normalization_radius": 3.0}]
    found = [d for d in rt.validate(rt.System.from_json(system_text(bases)))
             if d.code == "shape.normalization_radius_exceeded"]
    assert len(found) == 1
    assert found[0].location.endswith("/shape/base/normalization_radius")


EXAMPLES = sorted(EXAMPLE_DIR.glob("*.rtt.json"))


def test_the_forbes_examples_are_there() -> None:
    assert [f.name for f in EXAMPLES] == ["qbfs_singlet.rtt.json", "qcon_singlet.rtt.json"]


@pytest.mark.parametrize("file", EXAMPLES, ids=lambda f: f.name)
def test_a_forbes_example_is_canonical_and_valid(file: Path) -> None:
    system = rt.load(file)
    assert system.to_json().encode("utf-8") == file.read_bytes()  # as `rtt format` writes it
    jsonschema.Draft202012Validator(SYSTEM_SCHEMA).validate(
        json.loads(file.read_text(encoding="utf-8")))
    assert [d for d in rt.validate(system) if d.severity == rt.Severity.ERROR] == []


def test_the_qbfs_example_has_the_vertex_curvature_of_eq_3_16() -> None:
    # [F10b], Eq. (3.16): c_v = c + 2 S_0 / rho^2 with S_0 = a_0 + 13 a_1 / sqrt(19) ([F07],
    # Eq. (14)); the power of the plano-convex singlet is (n - 1) c_v.
    cs = rt.compile(rt.load(EXAMPLE_DIR / "qbfs_singlet.rtt.json"))
    fo = rt.paraxial.first_order(cs)
    assert fo.efl is not None
    cv = 1.0 / 60.0 + 2.0 * (0.1 - 13.0 * 0.02 / math.sqrt(19.0)) / 12.7**2
    assert 1.0 / fo.efl == pytest.approx(0.5168 * cv, rel=1e-12)
    rt.paraxial.seidel(cs)  # order 4: the Seidel sums exist
