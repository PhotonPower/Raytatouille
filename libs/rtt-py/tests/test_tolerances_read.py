"""The tolerance set from Python (#253, ADR 0035): read-only types of the section "tolerances",
their defaults from the model, and the edit form through a patch."""

from __future__ import annotations

import json

import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import model

FILE = REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json"


def system() -> rt.System:
    return rt.load(FILE)


def test_set_of_the_reference_file() -> None:
    t = system().tolerances
    assert isinstance(t, model.ToleranceSet)
    assert t.tolerance_class == model.ToleranceClass.NORMAL
    assert len(t.overrides) == 5
    assert len(t.items) == 3
    assert len(t.compensators) == 1
    assert len(t.criteria) == 3
    assert isinstance(t.compensator_merit, model.Optimization)
    assert isinstance(t.compensator_merit.operands[0], model.SpotRmsOperand)
    assert isinstance(t.monte_carlo, model.MonteCarloDefaults)
    assert (t.monte_carlo.seed, t.monte_carlo.trials) == (20261012, 200)


def test_overrides_in_their_three_forms() -> None:
    cls, value, segment, pivot, none = system().tolerances.overrides
    assert isinstance(cls, model.ToleranceClassOverride)
    assert (cls.node, cls.tolerance_class) == ("objective", model.ToleranceClass.PRECISE)
    assert isinstance(value, model.ToleranceValueOverride)
    assert (value.surface, value.node, value.segment, value.kind, value.limit,
            value.distribution) == ("L1.S2", None, None, model.ToleranceKind.CURVATURE_RELATIVE,
                                    0.005, model.ToleranceDistribution.NORMAL)
    assert isinstance(segment, model.ToleranceValueOverride)
    assert (segment.node, segment.segment, segment.kind, segment.distribution) == (
        "L1", 1, model.ToleranceKind.INDEX, model.ToleranceDistribution.UNIFORM)
    assert isinstance(pivot, model.TolerancePivotOverride)
    assert (pivot.node, pivot.surface, pivot.pivot) == ("objective", None, [0.0, 0.0, 9.5])
    assert isinstance(none, model.ToleranceClassOverride)
    assert none.tolerance_class == model.ToleranceClass.NONE


def test_items_compensators_and_criteria() -> None:
    t = system().tolerances
    wedge, conic, row = t.items
    assert isinstance(wedge, model.ToleranceItem)
    assert (wedge.id, wedge.surface, wedge.node, wedge.parameter, wedge.param, wedge.segment,
            wedge.kind, wedge.limit, wedge.distribution) == (
        "S3-wedge", "L1.S3", None, None, None, None, model.ToleranceKind.TILT_X_DEG, 0.02,
        model.ToleranceDistribution.NORMAL)
    assert (conic.surface, conic.param, conic.kind, conic.distribution) == (
        "L1.S1", "/shape/base/conic", model.ToleranceKind.PARAM,
        model.ToleranceDistribution.UNIFORM)
    assert (row.parameter, row.surface, row.kind, row.distribution) == (
        "R1", None, model.ToleranceKind.PARAMETER, model.ToleranceDistribution.ENDPOINT)

    [focus] = t.compensators
    assert isinstance(focus, model.Compensator)
    assert (focus.node, focus.surface, focus.parameter, focus.param, focus.travel_min,
            focus.travel_max) == ("image", None, None, "/pose/position/2", -1.0, 1.0)

    axis, edge, efl = t.criteria
    assert isinstance(axis, model.Criterion)
    assert (axis.id, axis.min, axis.max, axis.max_change) == ("spot-axis", None, 0.02, None)
    assert isinstance(edge.operand, model.SpotRmsOperand)
    assert edge.operand.field == 2
    # ADR 0035, point 7: a criterion's operand has the model's target 0 and weight 1.
    assert (edge.operand.target, edge.operand.weight) == (0.0, 1.0)
    assert isinstance(efl.operand, model.FirstOrderOperand)
    assert efl.max_change == 0.5


def test_without_tolerances_the_set_is_empty() -> None:
    t = rt.load(REFERENCE_DIR / "m0" / "singlet.rtt.json").tolerances
    assert t.tolerance_class == model.ToleranceClass.NONE
    assert (t.overrides, t.items, t.compensators, t.criteria) == ([], [], [], [])
    assert (t.compensator_merit, t.monte_carlo) == (None, None)


def test_enum_members() -> None:
    assert [k.name for k in model.ToleranceKind] == [
        "RADIUS", "CURVATURE_RELATIVE", "CONIC", "THICKNESS", "PARAM", "PARAMETER", "DECENTER_X",
        "DECENTER_Y", "SHIFT_Z", "TILT_X_DEG", "TILT_Y_DEG", "TILT_Z_DEG", "INDEX", "ABBE"]
    assert [d.name for d in model.ToleranceDistribution] == ["NORMAL", "UNIFORM", "ENDPOINT"]
    assert [c.name for c in model.ToleranceClass] == ["NONE", "COARSE", "NORMAL", "PRECISE"]


@pytest.mark.parametrize("attribute", ["items", "tolerance_class", "monte_carlo"])
def test_the_set_is_read_only(attribute: str) -> None:
    t = system().tolerances
    with pytest.raises(AttributeError):
        setattr(t, attribute, None)


def test_a_patch_adds_an_item() -> None:
    # ADR 0035, point 10: the section is part of the edit form, so a patch reaches it.
    s = system()
    edited = rt.apply_patch(s, [{"op": "add", "path": "/tolerances/items/-", "value": {
        "id": "gap", "node": "image", "kind": "thickness", "limit": 0.1}}])
    item = edited.tolerances.items[-1]
    assert (item.id, item.node, item.kind) == ("gap", "image", model.ToleranceKind.THICKNESS)
    assert len(s.tolerances.items) == 3  # the input is not changed
    data = edited.to_dict()
    assert data["tolerances"]["items"][-1]["distribution"] == "normal"
    assert json.loads(edited.to_json())["tolerances"]["items"][-1] == {
        "id": "gap", "node": "image", "kind": "thickness", "limit": 0.1}
