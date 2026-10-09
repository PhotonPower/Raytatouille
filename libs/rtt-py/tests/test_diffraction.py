"""Diffraction orders from Python against closed forms (#134, ADR 0025).

Grating bench tests/reference/m4/grating_transmission.rtt.json: stop at z = 0, a linear grating
of G = 300 lines/mm (lines along y) on a thin element at z = 10, screen at z = 60, vacuum,
lambda0 = 0.5876 um. In direction cosines the grating equation (Palmer, Diffraction Grating
Handbook, Eq. (2-1), with the sign mapping of ADR 0025) reads t'_x - t_x = m lambda0 G, for the
transmitted and for the reflected order; the order adds m lambda0 G x to the OPL at the hit
point x (ADR 0025, point 3).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille.trace import RayBatch, RayStatus

LAMBDA_G = 0.5876e-3 * 300.0  # lambda0 G, lambda0 in mm
SCREEN = 60.0
GRATING = 10.0
START = -5.0


def bench() -> rt.System:
    return rt.load(REFERENCE_DIR / "m4" / "grating_transmission.rtt.json")


def grating_pointer(system: rt.System) -> str:
    """JSON pointer of the grating surface G in the edit form."""
    pointer = system.locate_surface("G")
    assert pointer is not None
    return pointer


def start(dx: float = 0.0) -> RayBatch:
    """Five rays at x = -2 ... 2 mm, y = 0.5 mm, z = -5 (in front of the stop), direction (dx, 0, sqrt(1 - dx^2))."""
    rays = RayBatch(5)
    for k in range(5):
        rays.pos_x[k] = k - 2.0
        rays.pos_y[k] = 0.5
        rays.pos_z[k] = -5.0
        rays.dir_x[k] = dx
        rays.dir_z[k] = math.sqrt(1.0 - dx * dx)
    return rays


def trace(system: rt.System, path: str, rays: RayBatch) -> RayBatch:
    cs = rt.compile(system, rt.MaterialLibrary())
    rt.trace.trace(cs, rays, path=path)
    return rays


@pytest.mark.parametrize("dx", [0.0, -0.1], ids=["normal", "oblique"])
@pytest.mark.parametrize("order", [-1, 0, 1])
def test_transmitted_orders_follow_the_grating_equation(order: int, dx: float) -> None:
    name = "order 0" if order == 0 else f"order {order:+d}"
    rays = trace(bench(), name, start(dx))
    assert (np.asarray(rays.status) == int(RayStatus.ALIVE)).all()
    sx = dx + order * LAMBDA_G
    sz = math.sqrt(1.0 - sx * sx)
    np.testing.assert_allclose(rays.dir_x, sx, rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.dir_y, 0.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.dir_z, sz, rtol=0, atol=1e-12)
    # OPL: path to the grating, m lambda0 G x there, path to the screen.
    cz = math.sqrt(1.0 - dx * dx)
    x_grating = np.arange(5) - 2.0 + (GRATING - START) * dx / cz
    expected = (GRATING - START) / cz + order * LAMBDA_G * x_grating + (SCREEN - GRATING) / sz
    np.testing.assert_allclose(rays.opl, expected, rtol=0, atol=1e-10)
    np.testing.assert_allclose(rays.pos_x, x_grating + (SCREEN - GRATING) * sx / sz,
                               rtol=0, atol=1e-10)


def test_order_6_is_evanescent_at_normal_incidence_and_propagates_obliquely() -> None:
    # 6 lambda0 G = 1.0577 > 1: no real direction at normal incidence (status Evanescent at the
    # grating, surface 1). With t_x = -0.2 the order exists, t'_x = 0.8577, and leaves the
    # screen (half width 40 mm).
    normal = trace(bench(), "order +6", start())
    assert (np.asarray(normal.status) == int(RayStatus.EVANESCENT)).all()
    assert (np.asarray(normal.last_surface) == 1).all()
    oblique = trace(bench(), "order +6", start(-0.2))
    assert (np.asarray(oblique.status) == int(RayStatus.VIGNETTED)).all()
    assert (np.asarray(oblique.last_surface) == 2).all()


def test_reflection_grating() -> None:
    # Ideal mirror with the grating: the reflected order +1 runs back with t'_x = lambda0 G.
    s = bench()
    s = rt.apply_patch(s, [
        {"op": "replace", "path": grating_pointer(s) + "/interaction",
         "value": {"type": "ideal_mirror"}},
        {"op": "add", "path": "/paths/-", "value": {"name": "back", "events": [
            {"surface": "STO", "kind": "transmit"},
            {"surface": "G", "kind": "reflect", "order": 1},
            {"surface": "STO", "kind": "transmit"}]}},
    ])
    rays = trace(s, "back", start())
    assert (np.asarray(rays.status) == int(RayStatus.ALIVE)).all()
    sz = math.sqrt(1.0 - LAMBDA_G**2)
    np.testing.assert_allclose(rays.dir_x, LAMBDA_G, rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.dir_z, -sz, rtol=0, atol=1e-12)
    np.testing.assert_allclose(rays.pos_z, 0.0, rtol=0, atol=1e-12)
    x = np.arange(5) - 2.0
    np.testing.assert_allclose(rays.opl, GRATING - START + LAMBDA_G * x + GRATING / sz, rtol=0,
                               atol=1e-10)


@pytest.mark.parametrize(("order", "efficiency"), [(-1, 0.0), (0, 0.3), (1, 0.4)])
def test_diffraction_efficiency_scales_the_weight(order: int, efficiency: float) -> None:
    # Orders not listed have efficiency 0: weight 0, status Alive (ADR 0025, point 5). The
    # interface in vacuum transmits completely, so the transmission is the efficiency.
    s = bench()
    s = rt.apply_patch(s, [{
        "op": "add", "path": grating_pointer(s) + "/diffraction_efficiency",
        "value": [{"order": 0, "efficiency": 0.3}, {"order": 1, "efficiency": 0.4}]}])
    cs = rt.compile(s, rt.MaterialLibrary())
    name = "order 0" if order == 0 else f"order {order:+d}"
    t = an.path_transmission(cs, name, start=start())
    assert t.rays_arrived == 5
    np.testing.assert_allclose(t.rays.weight, efficiency, rtol=0, atol=1e-12)
    assert t.mean == pytest.approx(efficiency, abs=1e-12)
