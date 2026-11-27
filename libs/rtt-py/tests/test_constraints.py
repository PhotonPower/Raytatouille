"""Exact constraints of the merit function in Python (#196, ADR 0056): the operand attribute
constraint (file key "constraint": "equal", no weight), OptimizeOptions.constraint_tolerance,
OperandValue.constraint in the result and in the result format 0.1.25 (the numbers of the solver
are checked in C++, test_lm_constraints.cpp and the M5 acceptance)."""
from __future__ import annotations

import json
from typing import Any

import pytest

from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import model


def _constrained() -> rt.System:
    return rt.System.from_json(
        (REFERENCE_DIR / "m5" / "singlet_constraint.rtt.json").read_text(encoding="utf-8"))


def test_the_operand_carries_the_constraint() -> None:
    system = _constrained()
    op = system.optimization.operands[0]
    assert isinstance(op, model.FirstOrderOperand)
    assert op.constraint == model.ConstraintKind.EQUAL
    assert op.weight == 1.0  # a constraint has no weight; the model keeps 1
    assert op.target == 100.0
    d = json.loads(system.to_json())
    assert d["optimization"]["operands"][0] == {"type": "efl", "path": "main", "target": 100.0,
                                                "constraint": "equal"}
    plain = json.loads(system.to_json())
    del plain["optimization"]["operands"][0]["constraint"]
    assert rt.System.from_json(json.dumps(plain)).optimization.operands[0].constraint is None
    with pytest.raises(rt.ParseError):
        bad = json.loads(system.to_json())
        bad["optimization"]["operands"][0]["weight"] = 1.0
        rt.System.from_json(json.dumps(bad))


def test_optimize_holds_the_constraint() -> None:
    system = _constrained()
    options = rt.optim.OptimizeOptions(max_iterations=5)
    assert options.constraint_tolerance == 1e-12
    result = rt.optim.optimize(system, options=options)
    [efl] = result.operands
    assert efl.constraint
    assert efl.weight == 0.0 and efl.contribution == 0.0
    assert abs(efl.value - 100.0) <= 1e-12 * 100.0  # projected at the start, held in every step
    # the option reaches the run, also through optimize(aiming=), which copies the options
    loose = rt.optim.OptimizeOptions(max_iterations=0, constraint_tolerance=1e-3)
    copied = rt.optim.optimize(system, options=loose, aiming=rt.trace.Aiming.REAL)
    assert abs(copied.operands[0].value - 100.0) <= 1e-3 * 100.0
    assert abs(copied.operands[0].value - 100.0) > 1e-12 * 100.0  # content guard: loose enough
    # the result format: constraint from 0.1.25, False for an older file
    doc = json.loads(rt.results.to_json(result))
    assert doc["schema_version"] == "0.1.25"
    assert doc["data"]["operands"][0]["constraint"] is True
    assert rt.results.load_json(json.dumps(doc)).data["operands"][0]["constraint"] is True
    old = json.loads(rt.results.to_json(result))
    old["schema_version"] = "0.1.24"
    del old["data"]["operands"][0]["constraint"]
    assert rt.results.load_json(json.dumps(old)).data["operands"][0]["constraint"] is False


def test_errors_of_constraints() -> None:
    system = _constrained()
    d = json.loads(system.to_json())
    # only the constraint, no generator: nothing to minimize
    d["optimization"]["generators"] = []
    with pytest.raises(rt.OptimError) as e:
        rt.optim.optimize(rt.System.from_json(json.dumps(d)))
    assert [x.code for x in e.value.diagnostics] == ["optim.no_operands"]


def _compensated(constraint: bool) -> dict[str, Any]:
    """m1/two_lenses_stop_between with L1 shifted by 0.2 mm, L2 and the image as compensators
    and the EFL in the compensator merit, as a constraint or as a weighted operand."""
    d = json.loads((REFERENCE_DIR / "m1" / "two_lenses_stop_between.rtt.json").read_text())
    efl0 = rt.paraxial.first_order(rt.compile(rt.System.from_json(json.dumps(d)))).efl
    op: dict[str, Any] = {"type": "efl", "path": "main", "target": efl0}
    if constraint:
        op["constraint"] = "equal"
    d["tolerances"] = {
        "items": [{"id": "l1-dz", "node": "L1", "kind": "shift_z", "limit": 0.2,
                   "distribution": "uniform"}],
        "compensators": [{"node": "L2", "param": "/pose/position/2", "travel_min": -2.0,
                          "travel_max": 2.0},
                         {"node": "image", "param": "/pose/position/2", "travel_min": -5.0,
                          "travel_max": 5.0}],
        "compensator_merit": {"operands": [op], "generators": [{"type": "rms_spot",
                                                                "path": "main"}]},
        "criteria": [{"id": "efl", "operand": {"type": "efl", "path": "main"},
                      "max_change": 1e-6}]}
    result = rt.tolerance.sensitivity(rt.System.from_json(json.dumps(d)), confirm=False)
    return {"efl0": efl0, "plus": result.to_dict()["items"][0]["plus"]}


def test_a_constraint_in_the_compensator_merit() -> None:
    # ADR 0056, point 6: allowed there; the compensators hold the EFL of the trial exactly
    # (to eps_c = 1e-12 relative), so the criterion |EFL - nominal| <= 1e-6 holds.
    held = _compensated(True)
    value = held["plus"]["criteria"][0]["value"]
    assert abs(value - held["efl0"]) <= 1e-12 * held["efl0"]
    assert held["plus"]["criteria"][0]["satisfied"]
    # content guard: as a weighted operand beside the spot the EFL drifts and the criterion fails
    weighted = _compensated(False)
    assert abs(weighted["plus"]["criteria"][0]["value"] - weighted["efl0"]) > 1e-3
    assert not weighted["plus"]["criteria"][0]["satisfied"]


def test_constraints_of_the_compensator_merit_are_checked_once() -> None:
    # more constraints than compensators: an input error of the run, not of every trial
    d = json.loads((REFERENCE_DIR / "m1" / "two_lenses_stop_between.rtt.json").read_text())
    d["tolerances"] = {
        "items": [{"id": "l1-dz", "node": "L1", "kind": "shift_z", "limit": 0.2}],
        "compensators": [{"node": "image", "param": "/pose/position/2", "travel_min": -5.0,
                          "travel_max": 5.0}],
        "compensator_merit": {
            "operands": [{"type": "efl", "path": "main", "target": 38.0, "constraint": "equal"},
                         {"type": "bfl", "path": "main", "target": 30.0, "constraint": "equal"}],
            "generators": [{"type": "rms_spot", "path": "main"}]},
        "criteria": [{"id": "efl", "operand": {"type": "efl", "path": "main"}, "max": 40.0}]}
    with pytest.raises(rt.OptimError) as e:
        rt.tolerance.sensitivity(rt.System.from_json(json.dumps(d)), confirm=False)
    assert [x.code for x in e.value.diagnostics] == ["optim.constraints_too_many"]
