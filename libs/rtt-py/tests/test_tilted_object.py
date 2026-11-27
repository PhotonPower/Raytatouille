"""A tilted object plane in Python (#358, ADR 0055): ObjectSpace.slope as two Params, the field
points on the plane z = z_obj + sx x + sy y, and Seidel rejecting the tilt (the numbers are
checked in C++, test_tilted_object.cpp)."""
from __future__ import annotations

import pytest

from conftest import REFERENCE_DIR

import raytatouille as rt

TILTED = REFERENCE_DIR / "m1" / "tilted_object_ideal_lens.rtt.json"


def test_object_space_has_the_slopes() -> None:
    system = rt.load(TILTED)
    slope = system.object_space.slope
    assert len(slope) == 2
    assert (slope[0].value, slope[1].value) == (0.0, -1.0)
    # the edit form writes every Param as an object (ADR 0024)
    assert [p["value"] for p in system.to_dict()["object"]["slope"]] == [0.0, -1.0]


def test_rays_start_on_the_tilted_plane() -> None:
    cs = rt.compile(rt.load(TILTED), rt.MaterialLibrary())
    rays = rt.trace.make_rays(cs, rt.trace.SinglePupilPoint(), fields=[2])
    assert rays.pos_y[0] == 5.0
    assert rays.pos_z[0] == -105.0  # z = -100 - 1 * 5
    spot = rt.analysis.spot(cs, field=2)
    assert spot.rays_arrived == 127
    assert spot.stats.rms_centroid > 0.1  # defocused on the detector normal to the axis


def test_seidel_rejects_the_tilt() -> None:
    cs = rt.compile(rt.load(TILTED), rt.MaterialLibrary())
    with pytest.raises(rt.ParaxialError, match="tilted"):
        rt.paraxial.seidel(cs)
    assert rt.paraxial.first_order(cs).efl == pytest.approx(50.0)
