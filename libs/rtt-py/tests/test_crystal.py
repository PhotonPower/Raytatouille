"""Uniaxial crystals from Python against closed forms (#134, ADR 0026).

Calcite plate tests/reference/m4/calcite_walkoff.rtt.json: n_O = 1.6584, n_E = 1.4864 (CONST),
optic axis (1, 0, 1) at 45 degree in the x-z plane, thickness t = 2 mm from z = 10 to z = 12,
detector at z = 20, vacuum. At normal incidence the wave normal stays +z for both modes:
- index of the e-mode along k at 45 degree to the axis, Lam, Eq. (2.39):
  1 / n_e^2 = cos^2(45) / n_O^2 + sin^2(45) / n_E^2;
- walk-off of S from the normal of the K-surface (Lam, p. 107; derivation in docs/quellen.md):
  tan(rho) = (n_O^2 - n_E^2) / (n_O^2 + n_E^2), away from the axis for a negative crystal;
- OPL in the crystal n l (k . S) with l = t / cos(rho) and k . S = cos(rho), Lam, Eq. (2.17):
  OPL = 18 + n t at the detector for rays from z = 0;
- each mode carries half of unpolarized light (projection model, ADR 0026, point 5).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille.trace import RayBatch, RayStatus

N_O = 1.6584
N_E = 1.4864
THICKNESS = 2.0
N_E45 = 1.0 / math.sqrt(0.5 / N_O**2 + 0.5 / N_E**2)
TAN_RHO = (N_O**2 - N_E**2) / (N_O**2 + N_E**2)


def plate() -> rt.CompiledSystem:
    """The calcite plate with the paths "o entry" and "e entry", which end at P.S1."""
    s = rt.load(REFERENCE_DIR / "m4" / "calcite_walkoff.rtt.json")
    s = rt.apply_patch(s, [
        {"op": "add", "path": "/paths/-", "value": {"name": name, "events": [
            {"surface": "P.S1", "kind": kind}]}}
        for name, kind in (("o entry", "ordinary"), ("e entry", "extraordinary"))])
    return rt.compile(s, rt.MaterialLibrary())


def start() -> RayBatch:
    """Three rays along +z from z = 0 at x = -0.5, 0, 0.5 mm and y = 0.3 mm; only dir is set
    (reading rule: wave := dir at mode_index 0)."""
    rays = RayBatch(3)
    for k in range(3):
        rays.pos_x[k] = 0.5 * (k - 1)
        rays.pos_y[k] = 0.3
    return rays


def test_closed_forms() -> None:
    # The values of the M4 acceptance (test_crystal_trace.cpp), computed independently.
    assert N_E45 == pytest.approx(1.565356826060665, abs=1e-12)
    assert THICKNESS * TAN_RHO == pytest.approx(0.218121366133243, abs=1e-12)
    assert math.degrees(math.atan(TAN_RHO)) == pytest.approx(6.224117602966, abs=1e-9)


def test_compiled_medium_of_the_crystal() -> None:
    media = plate().media
    assert not media[0].is_crystal and media[0].optic_axis is None
    assert media[0].reference_extraordinary == "" and media[0].index_extraordinary.size == 0
    crystal = [m for m in media if m.is_crystal]
    assert len(crystal) == 1
    c = crystal[0]
    assert c.reference == "CONST:1.6584" and c.reference_extraordinary == "CONST:1.4864"
    assert c.index == [complex(N_O, 0.0)]
    np.testing.assert_array_equal(c.index_extraordinary, [N_E])
    assert c.optic_axis is not None
    np.testing.assert_allclose(c.optic_axis, [math.sqrt(0.5), 0.0, math.sqrt(0.5)], rtol=0,
                               atol=1e-15)
    assert not c.optic_axis.flags.writeable and not c.index_extraordinary.flags.writeable


@pytest.mark.parametrize("mode", ["o", "e"])
def test_inside_the_crystal(mode: str) -> None:
    # After the entry: k = +z, n = n_O or n_e(45); S along k (o) or tilted by rho away from
    # the axis in the x-z plane (e).
    rays = start()
    rt.trace.trace(plate(), rays, path=f"{mode} entry")
    assert (np.asarray(rays.status) == int(RayStatus.ALIVE)).all()
    n = N_O if mode == "o" else N_E45
    np.testing.assert_allclose(rays.mode_index, n, rtol=0, atol=1e-12)
    for column, value in (("wave_x", 0.0), ("wave_y", 0.0), ("wave_z", 1.0)):
        np.testing.assert_allclose(getattr(rays, column), value, rtol=0, atol=1e-12)
    rho = 0.0 if mode == "o" else math.atan(TAN_RHO)
    np.testing.assert_allclose(rays.dir_x, -math.sin(rho), rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.dir_y, 0.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.dir_z, math.cos(rho), rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.weight, 0.5, rtol=0, atol=1e-12)


def test_double_image() -> None:
    # At the detector: both images parallel to the incidence, the e-image shifted by
    # -t tan(rho) in x, OPL 18 + n t; the o- and e-power of unpolarized light sum to 1.
    cs = plate()
    images = {}
    for mode in ("o", "e"):
        rays = start()
        rt.trace.trace(cs, rays, path=mode)
        assert (np.asarray(rays.status) == int(RayStatus.ALIVE)).all()
        assert not np.asarray(rays.mode_index).any()
        np.testing.assert_array_equal(rays.wave_z, rays.dir_z)
        np.testing.assert_allclose(rays.dir_z, 1.0, rtol=0, atol=1e-12)
        images[mode] = rays
    o, e = images["o"], images["e"]
    x0 = 0.5 * (np.arange(3) - 1)
    np.testing.assert_allclose(o.pos_x, x0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(e.pos_x, x0 - THICKNESS * TAN_RHO, rtol=0, atol=1e-10)
    np.testing.assert_allclose(o.opl, 18.0 + N_O * THICKNESS, rtol=0, atol=1e-10)
    np.testing.assert_allclose(e.opl, 18.0 + N_E45 * THICKNESS, rtol=0, atol=1e-10)
    np.testing.assert_allclose(np.asarray(o.weight) + np.asarray(e.weight), 1.0, rtol=0,
                               atol=1e-12)

