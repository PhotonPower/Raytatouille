"""Irradiance map on a surface from traced rays (#379, ADR 0048): both forms, the incident
weight, the event choice and the errors, against exact values."""

from __future__ import annotations

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt


def cell_centred_bundle() -> rt.trace.RayBatch:
    """64 collimated rays along +z from z = -1 at the centres of 2 x 2 sub-cells of the grid
    [-1, 1]^2 with 4 x 4 cells (4 rays per cell)."""
    rays = rt.trace.RayBatch(64)
    c = -1.0 + 0.125 + 0.25 * np.arange(8)
    x, y = np.meshgrid(c, c)
    rays.pos_x[:] = x.ravel()
    rays.pos_y[:] = y.ravel()
    rays.pos_z[:] = -1.0
    return rays


def test_both_forms_give_the_same_constant_map() -> None:
    # m4/ghost_plates: plane plates keep x and y at normal incidence; every ray reaches IMG with
    # the same weight w, so each cell of 0.5 x 0.5 mm holds 4 w / 0.25 = 16 w.
    cs = rt.compile(rt.load(REFERENCE_DIR / "m4" / "ghost_plates.rtt.json"))
    rays = cell_centred_bundle()
    _, paths = rt.trace.trace(cs, rays, record_path=True)
    w = float(rays.weight[0])
    assert 0.5 < w < 1.0 and np.all(rays.weight == w)  # content guard: Fresnel losses
    extent = (-1.0, 1.0, -1.0, 1.0)
    by_paths = rt.analysis.irradiance_map(cs, paths, surface="IMG", nx=4, ny=4, extent=extent)
    by_rays = rt.analysis.irradiance_map(cs, rays, nx=4, ny=4, extent=extent)
    for m in (by_paths, by_rays):
        assert isinstance(m, rt.analysis.IrradianceMap)
        assert m.values.shape == (4, 4)
        assert np.all(m.values == 16.0 * w)
        assert list(m.x_edges) == [-1.0, -0.5, 0.0, 0.5, 1.0]
        assert m.surface_id == "IMG" and m.rays_counted == 64
        assert m.to_global.shape == (4, 4) and m.to_global[3, 3] == 1.0
        assert m.not_reached.count == 0 and m.vignetted.count == 0 and m.outside_grid.count == 0
    assert by_paths.values.tobytes() == by_rays.values.tobytes()  # bitwise
    # The result format (ADR 0023) carries the grid and the counts.
    data = by_rays.to_dict()
    assert data["surface_id"] == "IMG" and data["not_reached"] == {"count": 0, "weight": 0.0}


def test_incident_weight_and_event_choice_on_the_michelson() -> None:
    # Path "reference arm": BS (transmit, T = 0.5), M1, BS (reflect), DET. The rays meet BS
    # with weight 1 and then with 0.5; the incident weight counts.
    cs = rt.compile(rt.load(REFERENCE_DIR / "m0" / "michelson.rtt.json"))
    rays = rt.trace.RayBatch(9)
    rays.pos_x[:] = np.repeat([-1.0, 0.0, 1.0], 3)
    rays.pos_y[:] = np.tile([-1.0, 0.0, 1.0], 3)
    _, paths = rt.trace.trace(cs, rays, path="reference arm", record_path=True)
    extent = (-20.0, 20.0, -20.0, 20.0)
    arm = "reference arm"
    first = rt.analysis.irradiance_map(cs, paths, surface="BS", path=arm, nx=8, ny=8,
                                       extent=extent, event=0)
    both = rt.analysis.irradiance_map(cs, paths, surface="BS", path=arm, nx=8, ny=8,
                                      extent=extent)
    assert first.weight_counted == 9.0 and list(first.events) == [0]
    assert list(both.events) == [0, 2] and both.rays_counted == 18
    assert both.weight_counted == pytest.approx(13.5, rel=1e-12)
    with pytest.raises(ValueError):
        rt.analysis.irradiance_map(cs, paths, surface="BS", path=arm, extent=extent,
                                   event=1)  # M1


def test_errors() -> None:
    cs = rt.compile(rt.load(REFERENCE_DIR / "m4" / "ghost_plates.rtt.json"))
    rays = cell_centred_bundle()
    _, paths = rt.trace.trace(cs, rays, record_path=True)
    with pytest.raises(ValueError, match="no aperture"):
        rt.analysis.irradiance_map(cs, rays)  # IMG has no aperture and there is no extent
    with pytest.raises(ValueError):
        rt.analysis.irradiance_map(cs, paths)  # a RayPaths needs the surface
    with pytest.raises(ValueError):
        rt.analysis.irradiance_map(cs, paths, surface="NOPE", extent=(-1.0, 1.0, -1.0, 1.0))
    with pytest.raises(ValueError):
        rt.analysis.irradiance_map(cs, rays, extent=(1.0, 1.0, -1.0, 1.0))
    with pytest.raises(ValueError):
        rt.analysis.irradiance_map(cs, rays, surface="STO", extent=(-1.0, 1.0, -1.0, 1.0))
    with pytest.raises(ValueError):
        rt.analysis.irradiance_map(cs, rays, extent=(-1.0, 1.0, -1.0, 1.0), event=0)
    # The default extent of the stop: its circular aperture.
    m = rt.analysis.irradiance_map(cs, paths, surface="STO", nx=2, ny=2)
    r = m.x_edges[-1]
    assert m.x_edges[0] == -r and m.y_edges[0] == -r and m.y_edges[-1] == r
    assert m.weight_counted == 64.0  # incident on the first surface: weight 1


def test_energy_and_centroid_on_the_singlet_image() -> None:
    # #379, reference cases "energy" and "spot": m1/singlet_const at field 1 (3.5 degree). The
    # integral of the map is the counted weight; binning moves a ray by at most half a cell, so
    # the centroid of the map lies within half a cell of the ray centroid.
    cs = rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json"))
    rays = rt.trace.make_rays(cs, rt.trace.GridPupil(21), fields=[1])
    rt.trace.trace(cs, rays)
    alive = rays.status == rt.trace.RayStatus.ALIVE
    w = rays.weight[alive]
    cx = float(np.average(rays.pos_x[alive], weights=w))
    cy = float(np.average(rays.pos_y[alive], weights=w))
    m = rt.analysis.irradiance_map(cs, rays, nx=40, ny=40,
                                   extent=(cx - 4.0, cx + 4.0, cy - 4.0, cy + 4.0))
    assert m.rays_counted == int(alive.sum()) > 300 and m.outside_grid.count == 0
    dx = m.x_edges[1] - m.x_edges[0]
    dy = m.y_edges[1] - m.y_edges[0]
    assert abs(m.values.sum() * dx * dy - m.weight_counted) <= 1e-12 * m.weight_counted
    assert abs(m.weight_counted - float(w.sum())) <= 1e-12 * float(w.sum())
    xc = 0.5 * (m.x_edges[:-1] + m.x_edges[1:])
    yc = 0.5 * (m.y_edges[:-1] + m.y_edges[1:])
    total = m.values.sum()
    assert abs(float((m.values.sum(axis=0) * xc).sum() / total) - cx) <= 0.5 * dx + 1e-12
    assert abs(float((m.values.sum(axis=1) * yc).sum() / total) - cy) <= 0.5 * dy + 1e-12


def test_grid_size_errors_are_value_errors() -> None:
    cs = rt.compile(rt.load(REFERENCE_DIR / "m4" / "ghost_plates.rtt.json"))
    rays = cell_centred_bundle()
    rt.trace.trace(cs, rays)
    for nx in (0, -1):
        with pytest.raises(ValueError, match="nx and ny"):
            rt.analysis.irradiance_map(cs, rays, nx=nx, extent=(-1.0, 1.0, -1.0, 1.0))

