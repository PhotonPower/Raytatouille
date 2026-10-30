"""Writes example files of the result format (ADR 0023) to tests/reference/results.

Run from the repository root with the Python package of a build on the path, e.g.:

    PYTHONPATH=build/py/python python tools/result_examples.py ToleranceItems ToleranceSensitivity

Each name writes <name>.result.json with to_json(indent=1) at the current SCHEMA_VERSION. Only
the types listed here have a generator; the older examples were written by hand in their PRs
and keep their version (ADR 0023: an example is not rewritten when the format version rises).
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "tests" / "reference"
CATALOGS = ROOT / "tests" / "catalogs"


def tolerance_items() -> Any:
    """expand of the M7 reference achromat with the SCHOTT catalogue (#261)."""
    import raytatouille as rt

    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOGS / "schott.agf")
    system = rt.load(REFERENCE / "m7" / "achromat_tolerances.rtt.json")
    return rt.tolerance.expand(system, materials=lib)


def tolerance_sensitivity() -> Any:
    """sensitivity of the three explicit items of the M7 reference achromat (class none), with
    the confirmation, so that the example holds TrialResult lists at both levels (#261)."""
    import raytatouille as rt

    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOGS / "schott.agf")
    system = rt.apply_patch(rt.load(REFERENCE / "m7" / "achromat_tolerances.rtt.json"), [
        {"op": "replace", "path": "/tolerances/class", "value": "none"},
        {"op": "remove", "path": "/tolerances/overrides/4"},
        {"op": "remove", "path": "/tolerances/overrides/0"},
    ])
    return rt.tolerance.sensitivity(system, materials=lib)


GENERATORS: dict[str, Callable[[], Any]] = {
    "ToleranceItems": tolerance_items,
    "ToleranceSensitivity": tolerance_sensitivity,
}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("names", nargs="+", choices=sorted(GENERATORS),
                        help="result types to write")
    args = parser.parse_args(argv)
    for name in args.names:
        out = REFERENCE / "results" / f"{name}.result.json"
        text = GENERATORS[name]().to_json(indent=1)
        out.write_text(text if text.endswith("\n") else text + "\n", encoding="utf-8",
                       newline="\n")
        print(out.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
