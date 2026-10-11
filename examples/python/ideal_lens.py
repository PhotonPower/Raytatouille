"""Ideal lens and ideal cylinder lens as placeholders in a system sketch (#178, ADR 0031).

tests/reference/r2/ideal_lens.rtt.json has an ideal lens f = 50 mm (path "ideal lens") and an
ideal cylinder lens f = 50 mm (path "cylinder lens"), both at z = 0, an object at 75 mm and the
detector at 150 mm (1/s' = 1/s + 1/f, magnification -2).

- The ideal lens images every ray of an object point into its image point: the spot is a point
  for every field, and its paraxial data are those of a thin lens (EFL = f / n_air in air).
- Its optical path is exact for the design conjugate object_distance = 75 mm: the OPD map is 0.
  Without it (object at infinity as design), the OPD shows a rest of some waves, although the
  rays still meet in one point.
- The cylinder lens focuses only along y (axis_deg = 0): an object point becomes a line. It has
  no paraxial model, so its rays are set by hand here.

main() returns 1 if a spot or an OPD is not as described.

Run from the repository root after building the package (pip install .):

    python examples/python/ideal_lens.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
FILE = REPO_ROOT / "tests" / "reference" / "r2" / "ideal_lens.rtt.json"


def main() -> int:
    system = rt.load(FILE)
    cs = rt.compile(system)
    ok = True

    fo = rt.paraxial.first_order(cs, path="ideal lens")
    print(f"ideal lens: EFL {fo.efl:.6f} mm, rear focal length {fo.rear_focal_length:.6f} mm, "
          f"image at z = {fo.image_z:.6f} mm, magnification {fo.lateral_magnification:.6f}")

    for field in range(3):
        spot = rt.analysis.spot(cs, path="ideal lens", field=field, wavelength=0)
        opd = rt.analysis.opd_map(cs, path="ideal lens", field=field)
        print(f"field {field}: RMS spot {spot.stats.rms_centroid:.2e} mm, "
              f"OPD PV {opd.pv:.2e} waves")
        ok = ok and spot.stats.rms_centroid < 1e-9 and opd.pv < 1e-6

    # The same lens with the design conjugate at infinity: the OPD of the object at 75 mm.
    ed = rt.Editor(system)
    lens = ed.system.locate_surface("IL")
    assert lens is not None
    ed.apply([{"op": "remove", "path": lens + "/interaction/object_distance"}])
    other = rt.compile(ed.system)
    pv = rt.analysis.opd_map(other, path="ideal lens", field=0).pv
    print(f"design conjugate at infinity, object at 75 mm: OPD PV {pv:.2f} waves")
    ok = ok and pv > 1.0

    # Cylinder lens: rays from (3, 2, -75) through the stop end on the line y = -4 at z = 150.
    rays = rt.trace.RayBatch(9)
    for i, (x, y) in enumerate([(a, b) for a in (-6.0, 0.0, 6.0) for b in (-6.0, 0.0, 6.0)]):
        dx, dy, dz = x - 3.0, y - 2.0, 74.0
        n = math.sqrt(dx * dx + dy * dy + dz * dz)
        rays.pos_x[i], rays.pos_y[i], rays.pos_z[i] = 3.0, 2.0, -75.0
        rays.dir_x[i], rays.dir_y[i], rays.dir_z[i] = dx / n, dy / n, dz / n
    rt.trace.trace(cs, rays, path="cylinder lens")
    ys = [float(y) for y in rays.pos_y]
    xs = [float(x) for x in rays.pos_x]
    print(f"cylinder lens: y from {min(ys):.9f} to {max(ys):.9f} mm, "
          f"x from {min(xs):.3f} to {max(xs):.3f} mm (a line)")
    ok = ok and all(abs(y + 4.0) < 1e-9 for y in ys) and max(xs) - min(xs) > 1.0
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
