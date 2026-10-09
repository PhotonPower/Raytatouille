"""Double image behind a calcite plate (#132, #134, ADR 0026).

The plate of tests/reference/m4/calcite_walkoff.rtt.json is 2 mm thick, n_O = 1.6584 and
n_E = 1.4864, with the optic axis at 45 degree in the x-z plane. A collimated bundle at normal
incidence splits into the ordinary ray (path "o", straight through) and the extraordinary ray
(path "e"): inside the plate its wave normal stays along +z, but the energy runs at the walk-off
angle rho to it, so the e-image is shifted by t tan(rho). The example traces both paths and
compares the shift, the index of the e-mode and the optical path difference with the closed
forms (Lam, Anisotropic Ray Trace, Eqs. (2.39) and (2.17); walk-off from the normal of the
K-surface, docs/quellen.md). Each mode carries half of unpolarized light.

Run from the repository root after building the package (pip install .):

    python examples/python/calcite_double_image.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
PLATE = REPO_ROOT / "tests" / "reference" / "m4" / "calcite_walkoff.rtt.json"
THICKNESS = 2.0  # mm


def main() -> int:
    system = rt.load(PLATE)
    # Inside the plate the rays of the path "e" are in the e-mode; a path that ends at the
    # entry face shows the state there (wave normal, mode index, energy direction).
    entry = rt.apply_patch(system, [{"op": "add", "path": "/paths/-", "value": {
        "name": "e entry", "events": [{"surface": "P.S1", "kind": "extraordinary"}]}}])
    cs = rt.compile(entry, rt.MaterialLibrary())
    crystal = next(m for m in cs.media if m.is_crystal)
    n_o = crystal.index[cs.reference_wavelength].real
    n_e = float(crystal.index_extraordinary[cs.reference_wavelength])
    axis = crystal.optic_axis
    assert axis is not None  # is_crystal
    print(f"calcite: n_O = {n_o}, n_E = {n_e}, optic axis {axis.round(4)}")

    rays = rt.trace.RayBatch(1)  # along +z from the origin; only dir is set (reading rule)
    rt.trace.trace(cs, rays, path="e entry")
    rho = math.atan2(-rays.dir_x[0], rays.dir_z[0])
    print(f"inside: wave normal ({rays.wave_x[0]:.4f}, {rays.wave_y[0]:.4f}, "
          f"{rays.wave_z[0]:.4f}), mode index {rays.mode_index[0]:.6f}, "
          f"walk-off {math.degrees(rho):.4f} deg")

    images = {}
    for mode in ("o", "e"):
        rays = rt.trace.RayBatch(1)
        rt.trace.trace(cs, rays, path=mode)
        images[mode] = (rays.pos_x[0], rays.opl[0], rays.weight[0])
        print(f"{mode}-image at x = {rays.pos_x[0]:+.6f} mm, OPL {rays.opl[0]:.6f} mm, "
              f"power {rays.weight[0]:.3f}")
    shift = images["o"][0] - images["e"][0]
    delta = images["o"][1] - images["e"][1]
    print(f"image separation {shift:.6f} mm, OPL(o) - OPL(e) = {delta:.6f} mm")

    # Closed forms at 45 degree between k and the axis.
    n_e45 = 1.0 / math.sqrt(0.5 / n_o**2 + 0.5 / n_e**2)
    tan_rho = (n_o**2 - n_e**2) / (n_o**2 + n_e**2)
    worst = max(abs(shift - THICKNESS * tan_rho), abs(delta - (n_o - n_e45) * THICKNESS),
                abs(math.tan(rho) - tan_rho))
    print(f"closed form: separation {THICKNESS * tan_rho:.6f} mm, "
          f"OPL difference {(n_o - n_e45) * THICKNESS:.6f} mm")
    return 0 if worst <= 1e-10 else 1


if __name__ == "__main__":
    sys.exit(main())
