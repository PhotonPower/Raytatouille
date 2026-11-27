"""The aiming of the optimization and of the tolerance trials in Python (#291; ADR 0030,
addendum #291): OptimizeOptions.aiming, optimize(aiming=), MeritFunction(aiming=) and the
aiming of sensitivity/monte_carlo, real by default (the numbers are checked in C++,
test_optimize_aiming.cpp and test_trial_aiming.cpp)."""
from __future__ import annotations

import json

import pytest

from conftest import REFERENCE_DIR

import raytatouille as rt


def _two_lenses() -> rt.System:
    """m1/two_lenses_stop_between (aimings differ) with a spot operand at the field edge and the
    image distance variable."""
    d = json.loads((REFERENCE_DIR / "m1" / "two_lenses_stop_between.rtt.json").read_text())
    d["root"]["children"][3]["pose"]["position"][2] = {"value": 80.0, "variable": True}
    d["optimization"] = {"operands": [{"type": "spot_rms", "path": "main", "field": 2,
                                       "target": 0.0}],
                         "generators": [{"type": "rms_spot", "path": "main"}]}
    return rt.System.from_json(json.dumps(d))


def test_merit_function_aims_like_the_spot() -> None:
    system = _two_lenses()
    cs = rt.compile(system)
    values = {}
    for aiming in (rt.trace.Aiming.REAL, rt.trace.Aiming.PARAXIAL):
        merit = rt.optim.MeritFunction(system, aiming=aiming)
        assert merit.aiming == aiming
        value = merit.evaluate(merit.start()).values[0]
        spot = rt.analysis.spot(cs, field=2, wavelength=cs.reference_wavelength, aiming=aiming)
        assert value == spot.stats.rms_centroid
        values[aiming] = value
    assert values[rt.trace.Aiming.REAL] != values[rt.trace.Aiming.PARAXIAL]  # content guard
    assert rt.optim.MeritFunction(system).aiming == rt.trace.Aiming.REAL


def test_optimize_takes_the_aiming() -> None:
    system = _two_lenses()
    a = rt.optim.optimize(system, aiming=rt.trace.Aiming.PARAXIAL)
    b = rt.optim.optimize(system,
                          options=rt.optim.OptimizeOptions(aiming=rt.trace.Aiming.PARAXIAL))
    c = rt.optim.optimize(system)
    assert a.patch == b.patch
    assert a.patch != c.patch  # content guard: another bundle, another image distance
    assert rt.optim.OptimizeOptions().aiming == rt.trace.Aiming.REAL


def test_tolerance_runs_take_the_aiming() -> None:
    d = json.loads((REFERENCE_DIR / "m1" / "two_lenses_stop_between.rtt.json").read_text())
    d["tolerances"] = {
        "items": [{"id": "l2-dz", "node": "L2", "kind": "shift_z", "limit": 0.05}],
        "criteria": [{"id": "spot", "operand": {"type": "spot_rms", "path": "main", "field": 2},
                      "max": 1.0}]}
    system = rt.System.from_json(json.dumps(d))
    real = rt.tolerance.sensitivity(system, confirm=False)
    paraxial = rt.tolerance.sensitivity(system, confirm=False, aiming=rt.trace.Aiming.PARAXIAL)
    assert real.nominal.criteria[0].value != pytest.approx(paraxial.nominal.criteria[0].value,
                                                            rel=1e-12)
    # None takes compensators.aiming, an explicit aiming replaces it
    from_options = rt.tolerance.sensitivity(
        system, confirm=False,
        compensators=rt.optim.OptimizeOptions(aiming=rt.trace.Aiming.PARAXIAL))
    assert from_options.nominal.criteria[0].value == paraxial.nominal.criteria[0].value
    replaced = rt.tolerance.sensitivity(
        system, confirm=False, aiming=rt.trace.Aiming.REAL,
        compensators=rt.optim.OptimizeOptions(aiming=rt.trace.Aiming.PARAXIAL))
    assert replaced.nominal.criteria[0].value == real.nominal.criteria[0].value
