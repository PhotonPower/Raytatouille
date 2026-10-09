"""Systems for the optimization tests (#169): the M5 reference systems with the merit functions
of the M5 acceptance (#170) in their section "optimization", and variants of them."""

from __future__ import annotations

import json
from typing import Any

from conftest import REFERENCE_DIR

import raytatouille as rt


def singlet_merit() -> rt.System:
    """m5/singlet_optim: R1, R2 and the image distance variable; EFL 100 mm (weight 1e4) and
    the RMS spot generator (defaults: 3 rings, 6 arms, reference centroid)."""
    return rt.load(REFERENCE_DIR / "m5" / "singlet_optim.rtt.json")


def gap_merit() -> rt.System:
    """m5/two_lens_gap: the row D and the image distance variable; EFL 55 mm and the marginal
    ray height 0 on IMG (two equations, two unknowns, zero residual)."""
    return rt.load(REFERENCE_DIR / "m5" / "two_lens_gap.rtt.json")


def with_merit(relative: str, optimization: dict[str, Any] | None) -> rt.System:
    """The reference system `relative` with its section "optimization" replaced (None: none)."""
    data = json.loads((REFERENCE_DIR / relative).read_text(encoding="utf-8"))
    data.pop("optimization", None)
    if optimization is not None:
        data["optimization"] = optimization
    return rt.System.from_json(json.dumps(data))
