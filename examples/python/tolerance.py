"""Tolerance an achromat: the items of its tolerance set, sensitivity and Monte Carlo (#261).

The achromat of tests/reference/m7/achromat_tolerances.rtt.json carries a tolerance set (ADR
0035): the class "normal" (the objective "precise"), three explicit items, the image position as
compensator and three criteria (the RMS spot on axis and in the zone, at most 0.02 and 0.03 mm,
and the EFL within 0.5 mm of nominal). rt.tolerance.expand lists the items the classes generate.
Sensitivity and Monte Carlo then run on the three explicit items alone (class none), which keeps
the example short:

- sensitivity: every item at +T and -T, the RSS per criterion, a limit T' per item from the
  allowance of the criteria, and a confirmation run at +-T' (ADR 0036, point 6);
- Monte Carlo: N trials with the seed of the file, the yield with its Wilson interval and the
  statistics per criterion (ADR 0036, point 7).

Run from the repository root after building the package (pip install .):

    python examples/python/tolerance.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
ACHROMAT = REPO_ROOT / "tests" / "reference" / "m7" / "achromat_tolerances.rtt.json"
SCHOTT = REPO_ROOT / "tests" / "catalogs" / "schott.agf"
TRIALS = 20


def main() -> int:
    materials = rt.MaterialLibrary()
    materials.add_catalog(SCHOTT)
    system = rt.load(ACHROMAT)

    items = rt.tolerance.expand(system, materials=materials)
    print(f"{len(items.items)} items, the first five:")
    for item in items.items[:5]:
        print(f"  {item.id:<38} {item.kind.name:<12} +-{item.limit:g}")

    explicit = rt.apply_patch(system, [
        {"op": "replace", "path": "/tolerances/class", "value": "none"},
        {"op": "remove", "path": "/tolerances/overrides/4"},
        {"op": "remove", "path": "/tolerances/overrides/0"},
    ])

    s = rt.tolerance.sensitivity(explicit, materials=materials)
    print("Sensitivity (unit of the criterion):")
    for c in s.criteria:
        print(f"  {c.id:<10} {c.kind.name:<9} nominal {c.nominal:.6g}  RSS {c.rss:.6g}"
              f"  {c.inverse.note}")
    print("Limits T -> T' (unit of the item):")
    for trials, limit in zip(s.items, s.limits):
        new = "unlimited" if limit is None else f"{limit:.4g}"
        print(f"  {trials.item.id:<10} {trials.item.limit:g} -> {new}")

    mc = rt.tolerance.monte_carlo(explicit, materials=materials, trials=TRIALS)
    y = mc.yield_
    print(f"Monte Carlo, seed {mc.seed}, {y.total} trials: yield {y.fraction:.2f}"
          f" ({y.lower:.2f} ... {y.upper:.2f}, 95 %)")
    for c in mc.criteria:
        median = c.percentiles[0].value
        print(f"  {c.id:<10} mean {c.mean:.6g}  sd {c.standard_deviation:.3g}  median {median:.6g}")

    ok = (s.status == rt.tolerance.SensitivityStatus.COMPLETED
          and mc.status == rt.tolerance.MonteCarloStatus.COMPLETED and y.total == TRIALS)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
