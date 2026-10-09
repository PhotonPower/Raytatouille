"""Diffraction orders of a transmission grating (#127, #134, ADR 0025).

The grating bench of tests/reference/m4 has a linear grating of 300 lines/mm on a thin element in
vacuum and one path per order (-1, 0, +1, +6). For each order the example traces a collimated
bundle and compares the direction with the grating equation t'_x = t_x + m lambda0 G (Palmer,
Diffraction Grating Handbook, Eq. (2-1), in direction cosines). Order +6 has no real direction
at normal incidence (6 lambda0 G > 1): its rays end with status EVANESCENT at the grating.

Then the same bench with diffraction efficiencies (order 0: 0.3, order +1: 0.4, all others 0):
the path transmission of each order is its efficiency.

Run from the repository root after building the package (pip install .):

    python examples/python/grating_orders.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH = REPO_ROOT / "tests" / "reference" / "m4" / "grating_transmission.rtt.json"
LINES_PER_MM = 300.0


def collimated(n: int = 5) -> rt.trace.RayBatch:
    """n rays along +z at x = -2 ... 2 mm, at z = -5 mm, in front of the stop."""
    rays = rt.trace.RayBatch(n)
    for k in range(n):
        rays.pos_x[k] = -2.0 + 4.0 * k / (n - 1)
        rays.pos_z[k] = -5.0
    return rays


def orders(cs: rt.CompiledSystem) -> float:
    """Prints the direction of every order; returns the largest deviation from the grating
    equation."""
    lambda_g = cs.wavelengths_um[cs.reference_wavelength] * 1e-3 * LINES_PER_MM
    print(f"grating 300 lines/mm, lambda0 = {lambda_g / LINES_PER_MM * 1e3:.4f} um, "
          f"lambda0 G = {lambda_g:.5f}")
    worst = 0.0
    for name in cs.path_names:
        m = int(name.split()[1])
        rays = collimated()
        rt.trace.trace(cs, rays, path=name)
        status = rt.trace.RayStatus(rays.status[0])
        if status != rt.trace.RayStatus.ALIVE:
            print(f"  {name:9s} {status.name} (m lambda0 G = {m * lambda_g:.4f})")
            continue
        sx = m * lambda_g
        worst = max(worst, abs(rays.dir_x[0] - sx), abs(rays.dir_z[0] - math.sqrt(1 - sx * sx)))
        print(f"  {name:9s} angle {math.degrees(math.asin(rays.dir_x[0])):+8.4f} deg, "
              f"spot at x = {rays.pos_x[2]:+8.4f} mm")
    return worst


def grating_pointer(system: rt.System) -> str:
    """JSON pointer of the grating surface G in the edit form."""
    pointer = system.locate_surface("G")
    if pointer is None:
        raise ValueError("no surface G")
    return pointer


def efficiencies(system: rt.System) -> float:
    """Path transmission per order with efficiencies; returns the largest deviation."""
    table = {0: 0.3, 1: 0.4}
    system = rt.apply_patch(system, [{
        "op": "add", "path": grating_pointer(system) + "/diffraction_efficiency",
        "value": [{"order": m, "efficiency": e} for m, e in table.items()]}])
    cs = rt.compile(system, rt.MaterialLibrary())
    print("with efficiencies (order 0: 0.3, order +1: 0.4, all others 0):")
    worst = 0.0
    for name in ("order -1", "order 0", "order +1"):
        t = rt.analysis.path_transmission(cs, name, start=collimated())
        worst = max(worst, abs(t.mean - table.get(int(name.split()[1]), 0.0)))
        print(f"  {name:9s} transmission {t.mean:.3f}")
    return worst


def main() -> int:
    system = rt.load(BENCH)
    worst = max(orders(rt.compile(system, rt.MaterialLibrary())), efficiencies(system))
    return 0 if worst <= 1e-12 else 1


if __name__ == "__main__":
    sys.exit(main())
