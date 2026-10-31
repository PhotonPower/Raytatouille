"""object.reference from Python (#288, ADR 0038): read, compile, edit."""

from __future__ import annotations

import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt

FILE = REFERENCE_DIR / "r3" / "object_reference.rtt.json"


def test_the_reference_is_read() -> None:
    s = rt.load(FILE)
    assert s.object_space.reference == "IL"
    assert s.object_space.distance.value == 75.0


def test_object_z_is_in_front_of_the_vertex() -> None:
    # The ideal lens at z = 20, the object 75 mm in front of it: z = -55; from the origin -75.
    s = rt.load(FILE)
    assert rt.compile(s).object_z == -55.0
    origin = rt.apply_patch(s, [{"op": "remove", "path": "/object/reference"}])
    assert origin.object_space.reference is None
    assert rt.compile(origin).object_z == -75.0


def test_the_reference_is_added_by_a_patch() -> None:
    s = rt.load(FILE)
    origin = rt.apply_patch(s, [{"op": "remove", "path": "/object/reference"}])
    back = rt.apply_patch(origin, [
        {"op": "add", "path": "/object/reference", "value": {"surface": "IL"}}])
    assert back.object_space.reference == "IL"
    assert rt.compile(back).object_z == -55.0


def test_an_unknown_surface_is_an_error() -> None:
    s = rt.load(FILE)
    bad = rt.apply_patch(s, [
        {"op": "replace", "path": "/object/reference", "value": {"surface": "X"}}],
        check="structure_only")
    with pytest.raises(rt.CompileError) as error:
        rt.compile(bad)
    assert "object.reference_surface_unknown" in error.value.codes
