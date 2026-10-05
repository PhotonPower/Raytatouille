"""Polarization from Python: Malus's law, a quarter-wave plate and crossed polarizers.

Uses tests/reference/m3/polarizer_qwp.rtt.json (x polarizer, quarter-wave plate at 45 deg,
crossed analyzer with extinction ratio 1e-4). P and weight follow ADR 0021: P is
power-normalised, weight is the power for an unpolarized source.

Run from the repository root after building the package (pip install .):

    python examples/python/polarization.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
FILE = REPO_ROOT / "tests" / "reference" / "m3" / "polarizer_qwp.rtt.json"


def main(file: Path = FILE) -> int:
    compiled = rt.compile(rt.load(file))
    worst = 0.0

    # Malus: linear light at angle theta through the x polarizer (path "circular": polarizer,
    # then the lossless quarter-wave plate) transmits cos^2(theta).
    rays = rt.trace.make_rays(compiled, rt.trace.SinglePupilPoint(0.0, 0.0), path="circular",
                              fields=[0])
    rt.trace.trace(compiled, rays, path="circular")
    print("Malus's law (on axis)")
    for deg in (0, 30, 45, 60, 90):
        theta = math.radians(deg)
        state = np.array([math.cos(theta), math.sin(theta), 0.0], dtype=np.complex128)
        t = float(rt.polar.transmission(rays, state)[0])
        worst = max(worst, abs(t - math.cos(theta) ** 2))
        print(f"  theta {deg:2d} deg: transmission {t:.12f}, cos^2 {math.cos(theta) ** 2:.12f}")

    # After the quarter-wave plate at 45 deg the x light is circular: |S3| = S0.
    s = rt.polar.stokes(rays, np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]))[0]
    print(f"Stokes after the quarter-wave plate: {np.array2string(s, precision=12)}")
    worst = max(worst, abs(abs(s[3]) - s[0]))

    # Unpolarized through polarizer, plate and analyzer: (1 + eps) / 4.
    main_rays = rt.trace.make_rays(compiled, rt.trace.SinglePupilPoint(0.0, 0.0), path="main",
                                   fields=[0])
    rt.trace.trace(compiled, main_rays, path="main")
    weight = float(main_rays.weight[0])
    d = rt.polar.diattenuation(main_rays)
    print(f"Unpolarized weight through all three: {weight:.12f} (expected {1.0001 / 4.0:.12f}), "
          f"diattenuation {d.value[0]:.12f}")
    worst = max(worst, abs(weight - 1.0001 / 4.0))
    print(f"Largest deviation from the analytic values: {worst:.1e}")
    return 0 if worst <= 1e-12 else 1


if __name__ == "__main__":
    sys.exit(main())
