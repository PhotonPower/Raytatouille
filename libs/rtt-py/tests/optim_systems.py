"""Systems with a merit function for the optimization tests (#169): the M5 reference systems
with an optimization section added as JSON, the same merit functions as optim_cases.cpp builds
from the model types."""

from __future__ import annotations

import json
from typing import Any

from conftest import REFERENCE_DIR

import raytatouille as rt


def with_merit(relative: str, optimization: dict[str, Any]) -> rt.System:
    data = json.loads((REFERENCE_DIR / relative).read_text(encoding="utf-8"))
    data["optimization"] = optimization
    return rt.System.from_json(json.dumps(data))


def singlet_merit() -> rt.System:
    """m5/singlet_optim (R1, R2 and the image distance variable) with EFL 100 mm and the RMS
    spot generator (defaults: 3 rings, 6 arms, reference centroid)."""
    return with_merit("m5/singlet_optim.rtt.json", {
        "operands": [{"type": "efl", "path": "main", "target": 100.0}],
        "generators": [{"type": "rms_spot", "path": "main"}],
    })


def gap_merit() -> rt.System:
    """m5/two_lens_gap (the row D and the image distance variable) with EFL 40 mm and the
    marginal ray height 0 on IMG."""
    return with_merit("m5/two_lens_gap.rtt.json", {
        "operands": [
            {"type": "efl", "path": "main", "target": 40.0},
            {"type": "ray_y", "path": "main", "surface": "IMG", "py": 1.0, "target": 0.0},
        ],
    })
