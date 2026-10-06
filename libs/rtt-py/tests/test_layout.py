"""Geometry export from Python (#81, raytatouille.layout) against analytic values.

The C++ reference tests (libs/rtt-compile/tests/test_layout.cpp) cover the formulas in detail;
here the Python API: broadcasting, shapes, surface ids, plane shorthands and errors.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt


def singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def sphere_sag(radius: float, r: np.ndarray) -> np.ndarray:
    # Conic sag with k = 0 (Forbes 2011, Eq. (2.1)): c r^2 / (1 + sqrt(1 - c^2 r^2)).
    c = 1.0 / radius
    result: np.ndarray = c * r * r / (1.0 + np.sqrt(1.0 - c * c * r * r))
    return result


def test_surfaces_and_elements(reference_dir: Path) -> None:
    # m1/singlet_const: STO at z = 0, L1.S1 (conic R = 51.68) at z = 5, L1.S2 plane at z = 9,
    # IMG at z = 106.363; media front/back of the lens AIR | glass and glass | AIR.
    cs = singlet(reference_dir)
    surfaces = rt.layout.surfaces(cs)
    assert [s.id for s in surfaces] == ["STO", "L1.S1", "L1.S2", "IMG"]
    assert [s.kind for s in surfaces] == ["stop", "lens", "lens", "detector"]
    s1 = surfaces[1]
    assert s1.shape == "conic" and abs(s1.curvature - 1.0 / 51.68) <= 1e-15
    assert s1.conic_constant == 0.0 and s1.coefficients == []
    assert s1.aperture == ("circular", {"radius": 12.7, "inner_radius": 0.0})
    assert np.max(np.abs(s1.translation - np.array([0.0, 0.0, 5.0]))) <= 1e-12
    assert np.max(np.abs(s1.rotation - np.eye(3))) <= 1e-15
    assert s1.to_global.shape == (4, 4) and s1.to_global[2, 3] == s1.translation[2]
    media = [m.reference for m in cs.media]
    assert (media[s1.medium_front], media[s1.medium_back]) == ("AIR", "CONST:1.5168")
    s2 = surfaces[2]
    assert (media[s2.medium_front], media[s2.medium_back]) == ("CONST:1.5168", "AIR")
    assert s2.shape == "plane" and s2.max_radius is None
    elements = rt.layout.elements(cs)
    assert [e.name for e in elements] == ["stop", "L1", "image"]
    assert elements[1].kind == "lens" and elements[1].first_surface == 1
    assert elements[1].surface_count == 2 and elements[1].segmented
    assert [cs.media[i].reference for i in elements[1].media] == ["CONST:1.5168"]


def test_sag_and_normal_broadcast(reference_dir: Path) -> None:
    # Sag of L1.S1 against the closed formula, broadcasting a column of x against a row of y;
    # normal (-dz/dx, -dz/dy, 1)/|...| with dz/dr = c r / sqrt(1 - c^2 r^2). Tolerance 1e-12.
    cs = singlet(reference_dir)
    x = np.linspace(-8.0, 8.0, 9)[:, None]
    y = np.linspace(-6.0, 6.0, 7)[None, :]
    z = rt.layout.sag(cs, "L1.S1", x, y)
    assert z.shape == (9, 7)
    r = np.hypot(x, y)
    assert np.max(np.abs(z - sphere_sag(51.68, r))) <= 1e-12
    n = rt.layout.normal(cs, 1, x, y)
    assert n.shape == (9, 7, 3)
    c = 1.0 / 51.68
    slope = c * r / np.sqrt(1.0 - c * c * r * r)
    with np.errstate(invalid="ignore", divide="ignore"):
        gx = np.where(r > 0, slope * x / r, 0.0)
        gy = np.where(r > 0, slope * y / r, 0.0)
    expected = np.stack(np.broadcast_arrays(-gx, -gy, np.ones_like(r)), axis=-1)
    expected /= np.linalg.norm(expected, axis=-1, keepdims=True)
    assert np.max(np.abs(n - expected)) <= 1e-12
    assert np.max(np.abs(rt.layout.normal(cs, 1, 0.0, 0.0, frame="global") -
                         np.array([0.0, 0.0, 1.0]))) <= 1e-15


def test_profile_and_outline_of_the_singlet(reference_dir: Path) -> None:
    # Profile of L1.S1 in the y-z plane: (0, y, 5 + sag(y)), y from -12.7 to 12.7; the outline
    # of L1 is closed and its area is 8 a - I(R, a) with I the integral of the sphere sag
    # (derivation and tolerance 1e-5 mm^2 as in the C++ test).
    cs = singlet(reference_dir)
    pieces = rt.layout.profile(cs, "L1.S1", "yz", samples=51)
    assert len(pieces) == 1 and pieces[0].shape == (51, 3)
    y = np.linspace(-12.7, 12.7, 51)
    expected = np.stack([np.zeros(51), y, 5.0 + sphere_sag(51.68, np.abs(y))], axis=1)
    assert np.max(np.abs(pieces[0] - expected)) <= 1e-12
    outline = rt.layout.outlines(cs, 1, samples=2001)
    assert len(outline) == 1
    p = outline[0]
    assert np.array_equal(p[0], p[-1])
    area = 0.5 * abs(np.sum(p[:-1, 1] * p[1:, 2] - p[1:, 1] * p[:-1, 2]))
    a, radius = 12.7, 51.68
    integral = 2 * a * radius - a * math.sqrt(radius**2 - a**2) - radius**2 * math.asin(a / radius)
    assert abs(area - (8.0 * a - integral)) <= 1e-5
    assert rt.layout.outlines(cs, 0) == []  # the stop has no body


def test_doublet_outlines_and_plane_tuple(reference_dir: Path, catalog_dir: Path) -> None:
    # m2/achromat: two segments, two closed outlines; a plane given as (point, normal) equal to
    # "yz" gives the same result.
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    cs = rt.compile(rt.load(reference_dir / "m2" / "achromat.rtt.json"), lib)
    outlines = rt.layout.outlines(cs, 1, samples=31)
    assert len(outlines) == 2
    assert all(np.array_equal(o[0], o[-1]) for o in outlines)
    same = rt.layout.outlines(cs, 1, ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0)), samples=31)
    assert all(np.array_equal(a, b) for a, b in zip(outlines, same))


def test_layout_errors(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    with pytest.raises(ValueError, match="M10"):
        rt.layout.profile(cs, "L1.S1", ((0.0, 0.0, 5.0), (0.0, 0.0, 1.0)))
    with pytest.raises(ValueError, match="unknown surface"):
        rt.layout.sag(cs, "nope", 0.0, 0.0)
    with pytest.raises(ValueError, match="section plane"):
        rt.layout.profile(cs, 1, "xy")
    with pytest.raises(ValueError, match="frame"):
        rt.layout.normal(cs, 1, 0.0, 0.0, frame="world")
    with pytest.raises(ValueError, match="samples"):
        rt.layout.profile(cs, 1, samples=1)
    for bad in (99, -1):
        with pytest.raises(IndexError):
            rt.layout.sag(cs, bad, 0.0, 0.0)
        with pytest.raises(IndexError):
            rt.layout.outlines(cs, bad)
    with pytest.raises(ValueError, match="samples"):
        rt.layout.outlines(cs, 1, samples=-3)
    with pytest.raises(TypeError, match="integer"):
        rt.layout.profile(cs, 1, samples=2.5)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="integer"):
        rt.layout.outlines(cs, 1, samples=True)


def test_plane_shorthands_run_along_plus_y_and_plus_x(reference_dir: Path) -> None:
    # "yz" runs along +y, "xz" along +x (normal (0, -1, 0), t = z x n = +x): the profile of the
    # plane L1.S2 (r = 12.7 at z = 9) starts at -12.7 in both.
    cs = singlet(reference_dir)
    yz = rt.layout.profile(cs, "L1.S2", "yz", samples=3)[0]
    xz = rt.layout.profile(cs, "L1.S2", "xz", samples=3)[0]
    assert np.max(np.abs(yz - np.array([[0, -12.7, 9], [0, 0, 9], [0, 12.7, 9]]))) <= 1e-12
    assert np.max(np.abs(xz - np.array([[-12.7, 0, 9], [0, 0, 9], [12.7, 0, 9]]))) <= 1e-12
