"""aperture.reading from Python (#296, ADR 0040): read, aim, the ApertureError with .code."""

from __future__ import annotations

import math

import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt

FILE = REFERENCE_DIR / "r3" / "immersion_real_na.rtt.json"


def test_the_reading_is_read() -> None:
    s = rt.load(FILE)
    assert s.aperture.reading == rt.model.ApertureReading.REAL
    paraxial = rt.apply_patch(s, [{"op": "remove", "path": "/aperture/reading"}])
    assert paraxial.aperture.reading == rt.model.ApertureReading.PARAXIAL


def test_the_marginal_ray_has_the_real_na() -> None:
    # n sin U = NA for the ray at p = (0, 1), up to rounding (test_real_na.cpp).
    cs = rt.compile(rt.load(FILE))
    rays = rt.trace.make_rays(cs, rt.trace.SinglePupilPoint(px=0.0, py=1.0), fields=[0])
    x, y, z = float(rays.dir_x[0]), float(rays.dir_y[0]), float(rays.dir_z[0])
    assert 1.52 * math.hypot(x, y) / math.sqrt(x * x + y * y + z * z) == pytest.approx(
        1.45, abs=1e-12)


def test_an_unreachable_na_is_an_aperture_error() -> None:
    s = rt.apply_patch(rt.load(FILE), [{"op": "replace", "path": "/aperture/value", "value": 1.6}])
    cs = rt.compile(s)
    with pytest.raises(rt.ApertureError) as error:
        rt.trace.make_rays(cs, rt.trace.SinglePupilPoint(px=0.0, py=1.0), fields=[0])
    assert error.value.code == "aperture.na_unreachable"
    assert error.value.location == "/aperture/value"
    assert isinstance(error.value, ValueError)
