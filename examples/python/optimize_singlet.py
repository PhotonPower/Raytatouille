"""Optimize a singlet for EFL 100 mm and a small spot, then take the result back with undo (#169).

The singlet of tests/reference/m5/singlet_optim.rtt.json has three variables: both radii and
the image distance. Its merit function is the section "optimization" of the file: the operand
EFL = 100 mm (weight 1e4) and the generator rms_spot (Gaussian pupil quadrature, ADR 0030).
rt.optim.optimize returns a new system and an RFC 6902 patch from the input to it; an Editor
applies the patch as one step, so undo restores the input bit for bit.

Run from the repository root after building the package (pip install .):

    python examples/python/optimize_singlet.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
SINGLET = REPO_ROOT / "tests" / "reference" / "m5" / "singlet_optim.rtt.json"


def main() -> int:
    system = rt.load(SINGLET)

    result = rt.optim.optimize(system)
    h = result.history
    print(f"Status {result.status.name} after {result.iterations} solves, "
          f"{result.evaluations} merit evaluations")
    print("solve  phi           mu")
    for k in range(len(h)):
        print(f"{h.k[k]:5}  {h.phi[k]:.6e}  {h.mu[k]:.3e}")
    for v in result.variables:
        print(f"{v.pointer:<45} {v.start:10.4f} -> {v.end:10.4f}")
    [efl] = result.operands
    [spot] = result.generators
    print(f"EFL {efl.value:.6f} mm, RMS spot {spot.rms * 1e3:.3f} um")

    # The result as one undo step of an Editor.
    editor = rt.Editor(system)
    editor.apply(result.patch)
    optimized = editor.system.to_json() == result.system.to_json()
    editor.undo()
    restored = editor.system.to_json() == system.to_json()
    print(f"patch applied: {optimized}, undo restores the input: {restored}")
    return 0 if optimized and restored else 1


if __name__ == "__main__":
    sys.exit(main())
