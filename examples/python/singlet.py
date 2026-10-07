"""Load a system, print its first-order and prescription data and trace a hexapolar bundle.

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

    # Paraxial marginal and chief ray per surface, as in a prescription report.
    p = rt.paraxial.prescription(compiled, path="main")
    q = p.surfaces
    print("  surface        z          y          u        i     y_bar      u_bar    i_bar")
    for k in range(len(q)):
        name = compiled.surface_ids[int(q.surface[k])]
        print(f"  {name:8s} {q.z[k]:8.3f} {q.y[k]:10.5f} {q.u[k]:10.6f} {q.i[k]:8.5f}"
              f" {q.y_bar[k]:9.5f} {q.u_bar[k]:10.6f} {q.i_bar[k]:8.5f}")
    print(f"  total track {p.total_track:.3f} mm, paraxial working F/# "
          f"{p.paraxial_working_f_number:.4f}, paraxial image NA {p.paraxial_image_na:.5f}, "
          f"Lagrange invariant {p.lagrange_invariant:.6f} mm")

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
