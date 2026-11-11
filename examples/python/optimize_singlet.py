"""Optimize a singlet for EFL 100 mm and a small spot, then take the result back with undo (#169).

The singlet of tests/reference/m5/singlet_optim.rtt.json has three variables: both radii and
the image distance. Its merit function is the section "optimization" of the file: the operand
EFL = 100 mm (weight 1e4) and the generator rms_spot (Gaussian pupil quadrature, ADR 0030).
rt.optim.optimize returns a new system and an RFC 6902 patch from the input to it; an Editor
applies the patch as one step, so undo restores the input bit for bit. ``on_iteration`` shows
the state after every accepted solve, as a live display in a GUI would (#378): it follows the run
by applying each patch to a copy of the input.

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

    # Live state (#378): one call per accepted solve; the patch leads from the input to that state.
    live: list[str] = []

    def show(state: rt.optim.IterationState) -> None:
        shown = rt.apply_patch(system, state.patch)  # what a display would draw
        print(f"  accepted solve {state.k:3}: phi {state.phi:.6e}, "
              f"{len(shown.to_json())} bytes of system")
        live.append(state.patch)

    result = rt.optim.optimize(system, on_iteration=show)
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
    followed = bool(live) and live[-1] == result.patch
    print(f"the live display ended at the result: {followed}")
    return 0 if optimized and restored and followed else 1


if __name__ == "__main__":
    sys.exit(main())
