"""Load a system, print its first-order data and trace a hexapolar bundle.

Run from the repository root after building the package (pip install .):

    python examples/python/singlet.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]


def main(file: Path = REPO_ROOT / "tests" / "reference" / "m1" / "singlet_const.rtt.json") -> int:
    system = rt.load(file)
    problems = [d for d in rt.validate(system) if d.severity == rt.Severity.ERROR]
    if problems:
        for d in problems:
            print(d)
        return 1

    compiled = rt.compile(system)  # VACUUM, AIR and CONST: need no catalogue
    fo = rt.paraxial.first_order(compiled, path="main")
    print(f"{system.name}")
    print(f"  EFL {fo.efl:.6f} mm, BFL {fo.bfl:.6f} mm (absolute indices, AIR = Ciddor air)")

    for field in range(compiled.field_count):
        rays = rt.trace.make_rays(compiled, rt.trace.HexapolarPupil(rings=6), fields=[field])
        rt.trace.trace(compiled, rays)
        alive = rays.status == rt.trace.RayStatus.ALIVE
        x, y = rays.pos_x[alive], rays.pos_y[alive]
        rms = float(np.sqrt(np.mean((x - x.mean()) ** 2 + (y - y.mean()) ** 2)))
        print(f"  field {field}: {alive.sum()} rays, RMS spot radius {1000.0 * rms:.3f} um")
    return 0


if __name__ == "__main__":
    sys.exit(main())
