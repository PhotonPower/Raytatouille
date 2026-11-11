"""Irradiance on a surface from traced rays, as a texture for a GUI (#379, ADR 0048).

Traces a 41 x 41 grid pupil through the singlet at field 3.5 degree, bins the hits on the image
surface IMG into 24 x 24 cells around the spot and prints the map as text; the same for the
first lens surface L1.S1 from recorded paths, where the incident weight counts (before the
Fresnel loss of that surface). ``to_global`` and the cell edges place the map on the geometry.

Run from the repository root after building the package (pip install .):

    python examples/python/irradiance.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
SINGLET = REPO_ROOT / "tests" / "reference" / "m1" / "singlet_const.rtt.json"


def main() -> int:
    cs = rt.compile(rt.load(SINGLET))
    rays = rt.trace.make_rays(cs, rt.trace.GridPupil(41), fields=[1])
    _, paths = rt.trace.trace(cs, rays, record_path=True, max_recorded_rays=len(rays))

    # On the image: from the batch after the trace (any number of rays), around the spot. The
    # extent is in local coordinates of IMG; here they equal the global x and y (IMG is not
    # tilted), so the centroid of the global positions centres the grid.
    alive = rays.status == rt.trace.RayStatus.ALIVE
    cy = float(np.average(rays.pos_y[alive], weights=rays.weight[alive]))
    image = rt.analysis.irradiance_map(cs, rays, nx=24, ny=24,
                                       extent=(-0.25, 0.25, cy - 0.25, cy + 0.25))
    print(f"IMG: {image.rays_counted} rays, weight {image.weight_counted:.4f}, "
          f"peak {image.values.max():.1f} per mm^2")
    shades = " .:-=+*#%@"
    peak = image.values.max()
    for row in image.values[::-1]:  # +y up
        print("".join(shades[int(v / peak * (len(shades) - 1))] for v in row))

    # On the first lens surface: from the recorded paths, default extent = its aperture.
    lens = rt.analysis.irradiance_map(cs, paths, surface="L1.S1", nx=8, ny=8)
    print(f"L1.S1: {lens.rays_counted} rays, incident weight {lens.weight_counted:.1f} "
          f"(1 per ray), grid {lens.x_edges[0]:.1f} .. {lens.x_edges[-1]:.1f} mm")
    dx = image.x_edges[1] - image.x_edges[0]
    dy = image.y_edges[1] - image.y_edges[0]
    total = image.values.sum() * dx * dy
    conserved = abs(total - image.weight_counted) <= 1e-12 * image.weight_counted
    incident = lens.weight_counted == float(lens.rays_counted)
    print(f"energy conserved: {conserved}, incident weight on L1.S1: {incident}")
    return 0 if conserved and incident and image.outside_grid.count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
