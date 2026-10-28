"""The medium of the object space, object.medium (#282, ADR 0037; schema 0.7)."""

from __future__ import annotations

import json
import math
from typing import Any

import numpy as np
import pytest

import raytatouille as rt


def document(medium: str | None, *, at_infinity: bool = False) -> dict[str, Any]:
    """A plate of VACUUM behind the object space, a detector behind it; VACUUM around."""
    obj: dict[str, Any] = {"at_infinity": at_infinity}
    if not at_infinity:
        obj["distance"] = 2.0
    if medium is not None:
        obj["medium"] = medium
    return {
        "schema_version": "0.7.0",
        "units": {"length": "mm", "wavelength": "um"},
        "environment": {"medium": "VACUUM"},
        "object": obj,
        "wavelengths": [{"um": 0.5876, "reference": True}],
        "aperture": {"type": "epd", "value": 10.0},
        "fields": {"type": "object_height", "points": [{}]},
        "root": {"type": "assembly", "name": "system", "children": [
            {"type": "plate", "name": "P", "pose": {"position": [0, 0, 0]}, "material": "VACUUM",
             "surfaces": [{"id": "P.S1"}, {"id": "P.S2", "pose": {"position": [0, 0, 5.0]}}]},
            {"type": "detector", "name": "D", "pose": {"position": [0, 0, 40.0]},
             "surfaces": [{"id": "D"}]},
        ]},
        "paths": [{"name": "main", "events": "auto"}],
    }


def load(doc: dict[str, Any]) -> rt.System:
    return rt.System.from_json(json.dumps(doc))


def test_the_object_medium_is_read_and_none_without_one() -> None:
    assert load(document("CONST:1.53")).object_space.medium == "CONST:1.53"
    assert load(document(None)).object_space.medium is None


def test_the_object_medium_survives_a_json_round_trip() -> None:
    system = load(document("CONST:1.53"))
    assert rt.System.from_json(system.to_json()) == system
    assert json.loads(system.to_json())["object"]["medium"] == "CONST:1.53"
    assert "medium" not in json.loads(load(document(None)).to_json())["object"]


def test_the_object_medium_needs_a_finite_object() -> None:
    errors = [d for d in rt.validate(load(document("CONST:1.5", at_infinity=True)))
              if d.severity == rt.Severity.ERROR]
    assert [(e.code, e.location) for e in errors] == [("object.medium_infinite", "/object/medium")]


def test_an_unknown_object_medium_is_material_unknown() -> None:
    with pytest.raises(rt.CompileError) as error:
        rt.compile(load(document("NOPE:X")))
    assert [(d.code, d.location) for d in error.value.diagnostics] == [
        ("material.unknown", "/object/medium")]


def test_rays_from_the_object_medium_refract_by_snell_at_the_first_surface() -> None:
    # sin(theta') = 1.5 sin(theta) for glass of index 1.5 into the VACUUM plate (Snell's law, as
    # in rtt-trace/tests/test_sequential.cpp: de Greve 2006, Eq. (23)); 45 degrees is beyond the
    # critical angle asin(1 / 1.5) = 41.81 degrees.
    comp = rt.compile(load(document("CONST:1.5")))
    degrees = np.array([10.0, 30.0, 45.0])
    theta = np.radians(degrees)
    rays = rt.trace.RayBatch(len(theta))
    rays.pos_x[:] = 0.0
    rays.pos_y[:] = 0.0
    rays.pos_z[:] = -2.0
    rays.dir_x[:] = 0.0
    rays.dir_y[:] = np.sin(theta)
    rays.dir_z[:] = np.cos(theta)
    rays.wl[:] = 0
    rt.trace.trace(comp, rays, path="main")
    status = np.asarray(rays.status)
    assert (status[:2] == int(rt.trace.RayStatus.ALIVE)).all()
    assert status[2] == int(rt.trace.RayStatus.TIR)
    out = np.degrees(np.arcsin(np.asarray(rays.dir_y)[:2]))
    assert out == pytest.approx([15.0980866, 48.5903779], abs=1e-6)
    assert math.isclose(float(np.asarray(rays.dir_y)[0]), 1.5 * math.sin(theta[0]), abs_tol=1e-12)
