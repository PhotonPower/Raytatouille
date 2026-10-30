"""Sensitivity from Python (#261, part 2a; ADR 0036, points 6 and 10): rt.tolerance.sensitivity,
its result types and the result format ToleranceSensitivity.

The system is the M7 reference achromat with the class switched off, so only its three explicit
items take part (a tilt, a conic and the table parameter R1), with the image position as
compensator and three criteria (spot-axis and spot-zone quadratic, efl linear). The numbers are
checked against C++ bit for bit in test_bitwise_tolerance.py; the tests here check the structure,
the relations between the values (ADR 0036, point 6) and the three places of a cancellation."""

from __future__ import annotations

import json
import math
import threading
from functools import cache
from typing import Any

import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR

import raytatouille as rt
from raytatouille import model, tolerance

FILE = REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json"
ITEMS = ["S3-wedge", "S1-conic", "R1"]
CRITERIA = ["spot-axis", "spot-zone", "efl"]
TS = tolerance.TrialStatus
KIND = tolerance.CriterionKind


def schott() -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    return lib


def explicit_only() -> rt.System:
    """The achromat without class items: the three explicit items (as in
    test_tolerance_expand.py, test_without_a_class_only_the_explicit_items)."""
    return rt.apply_patch(rt.load(FILE), [
        {"op": "replace", "path": "/tolerances/class", "value": "none"},
        {"op": "remove", "path": "/tolerances/overrides/4"},  # image none
        {"op": "remove", "path": "/tolerances/overrides/0"},  # objective precise
    ])


@cache
def full() -> tolerance.ToleranceSensitivity:
    return tolerance.sensitivity(explicit_only(), materials=schott())


def test_the_result_has_one_entry_per_item_and_criterion() -> None:
    s = full()
    assert s.status == tolerance.SensitivityStatus.COMPLETED
    assert s.nominal.status == TS.VALID
    assert [t.item.id for t in s.items] == ITEMS
    assert s.expansion is not None
    assert [i.id for i in s.expansion.items] == ITEMS
    assert [c.id for c in s.criteria] == CRITERIA
    assert [c.kind for c in s.criteria] == [KIND.QUADRATIC, KIND.QUADRATIC, KIND.LINEAR]
    for t in s.items:
        assert t.sigma_ratio == tolerance.sigma_ratio(t.item.distribution)
        for trial in (t.plus, t.minus):
            assert trial is not None
            assert trial.status == TS.VALID
            assert len(trial.criteria) == 3
            assert len(trial.compensators) == 1  # the image position
            assert trial.lm_status is not None
    for c in s.criteria:
        assert len(c.items) == 3
        assert c.complete and c.excluded == []
        assert c.nominal is not None and c.nominal_reason == ""
    assert len(s.limits) == 3


def test_the_contributions_follow_from_the_two_sides() -> None:
    """ADR 0036, point 6: linear delta = (Phi+ - Phi-)/2, shift = (Phi+ + Phi-)/2 - Phi_0;
    quadratic q = (Phi+^2 + Phi-^2)/2 - Phi_0^2; RSS and prediction from them."""
    s = full()
    for k, c in enumerate(s.criteria):
        phi0 = c.nominal
        assert phi0 is not None
        assert phi0 == s.nominal.criteria[k].value
        for i, contribution in enumerate(c.items):
            plus, minus = contribution.plus, contribution.minus
            assert plus is not None and minus is not None
            assert plus == s.items[i].plus.criteria[k].value  # type: ignore[union-attr]
            assert minus == s.items[i].minus.criteria[k].value  # type: ignore[union-attr]
            assert contribution.excluded == ""
            if c.kind == KIND.LINEAR:
                assert contribution.delta == pytest.approx((plus - minus) / 2, rel=1e-14)
                assert contribution.shift == pytest.approx((plus + minus) / 2 - phi0,
                                                           rel=1e-12, abs=1e-15)
                assert contribution.q is None
            else:
                assert contribution.q == pytest.approx(
                    (plus * plus + minus * minus) / 2 - phi0 * phi0, rel=1e-12, abs=1e-18)
                assert contribution.delta is None and contribution.shift is None
        assert c.rss is not None
        if c.kind == KIND.LINEAR:
            deltas = [x.delta for x in c.items if x.delta is not None]
            assert len(deltas) == 3
            assert c.rss == pytest.approx(math.sqrt(sum(d * d for d in deltas)), rel=1e-14)
            assert c.predicted_sigma == pytest.approx(math.sqrt(sum(
                (t.sigma_ratio * d) ** 2 for t, d in zip(s.items, deltas))), rel=1e-14)
            assert c.predicted_mean_square is None
        else:
            qs = [x.q for x in c.items if x.q is not None]
            assert len(qs) == 3
            total = phi0 * phi0 + sum(qs)
            assert c.rss == pytest.approx(math.sqrt(max(0.0, total)), rel=1e-12)
            assert c.predicted_sigma is None
            assert c.predicted_mean_square is not None


def test_the_inverse_and_the_confirmation_use_the_smallest_limit() -> None:
    s = full()
    for c in s.criteria:
        # Only a positive allowance is distributed: then one limit per item, else none and a note
        # (spot-zone: Phi_0 = 0.041 mm lies above its max 0.03, so A = 0).
        distributed = c.inverse.allowance is not None and c.inverse.allowance > 0.0
        assert len(c.inverse.limits) == (3 if distributed else 0)
        assert (c.inverse.note == "") == distributed
        assert all(0 <= i < 3 for i in c.inverse.unlimited)
    assert [c.inverse.note == "" for c in s.criteria] == [True, False, True]  # content guard
    for i in range(3):
        own = [x for c in s.criteria if c.inverse.limits and (x := c.inverse.limits[i]) is not None]
        assert s.limits[i] == (min(own) if own else None)
    assert any(x is not None for x in s.limits)  # content guard: some criterion limits an item
    confirmation = s.confirmation
    assert confirmation is not None
    assert confirmation.limits == [
        t.item.limit if limit is None else limit for t, limit in zip(s.items, s.limits)]
    assert len(confirmation.plus) == 3 and len(confirmation.minus) == 3
    assert [c.id for c in confirmation.criteria] == CRITERIA
    for c in confirmation.criteria:
        assert c.rss is not None
        assert c.inverse.limits == [] and c.inverse.allowance is None


def test_without_confirmation_the_first_part_is_the_same() -> None:
    s = tolerance.sensitivity(explicit_only(), materials=schott(), confirm=False)
    assert s.confirmation is None
    assert s.status == tolerance.SensitivityStatus.COMPLETED
    a = s.to_dict()
    b = full().to_dict()
    del a["confirmation"], b["confirmation"]
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_sigma_ratio_of_the_three_distributions() -> None:
    d = model.ToleranceDistribution
    assert tolerance.sigma_ratio(d.NORMAL) == pytest.approx(0.43981, abs=5e-6)
    assert tolerance.sigma_ratio(d.UNIFORM) == pytest.approx(1 / math.sqrt(3), rel=1e-15)
    assert tolerance.sigma_ratio(d.ENDPOINT) == 1.0


def test_bitwise_the_same_for_one_and_all_threads() -> None:
    one = tolerance.sensitivity(explicit_only(), materials=schott(), threads=1)
    assert one.to_json() == full().to_json()  # full(): all threads


def test_progress_reports_both_stages() -> None:
    seen: dict[str, set[int]] = {}
    lock = threading.Lock()

    def progress(done: int, total: int, stage: str) -> None:
        with lock:
            seen.setdefault(stage, set()).add(total)

    tolerance.sensitivity(explicit_only(), materials=schott(), progress=progress)
    assert set(seen) == {"sensitivity", "confirmation"}
    assert seen["sensitivity"] == {6} and seen["confirmation"] == {6}  # 2 trials per item


def test_cancelled_before_the_start() -> None:
    # Cancelled during Phi_0 (prepare): nominal is CANCELLED, no trial runs, no exception.
    token = rt.CancelToken()
    token.cancel()
    s = tolerance.sensitivity(explicit_only(), materials=schott(), cancel=token)
    assert s.status == tolerance.SensitivityStatus.CANCELLED
    assert s.nominal.status == TS.CANCELLED
    assert all(t.plus is None and t.minus is None for t in s.items)
    assert all(c.rss is None for c in s.criteria)
    assert s.limits == [None, None, None]
    assert s.confirmation is None


def cancel_in(stage: str) -> tolerance.ToleranceSensitivity:
    token = rt.CancelToken()

    def progress(done: int, total: int, at: str) -> None:
        if at == stage and done >= 1:
            token.cancel()

    return tolerance.sensitivity(explicit_only(), materials=schott(), threads=1, cancel=token,
                                 progress=progress)


def test_cancelled_during_the_trials() -> None:
    s = cancel_in("sensitivity")
    assert s.status == tolerance.SensitivityStatus.CANCELLED
    assert s.nominal.status == TS.VALID
    assert s.confirmation is None
    assert all(c.rss is None for c in s.criteria)
    assert s.limits == [None, None, None]
    # Only items with both sides are reduced; the others are excluded with a reason.
    assert any(t.plus is None or t.minus is None for t in s.items)
    for c in s.criteria:
        for t, contribution in zip(s.items, c.items):
            done = t.plus is not None and t.minus is not None
            assert (contribution.excluded == "") == done


def test_cancelled_during_the_confirmation_keeps_the_first_part() -> None:
    s = cancel_in("confirmation")
    assert s.status == tolerance.SensitivityStatus.CANCELLED
    assert s.confirmation is None
    a = s.to_dict()
    b = full().to_dict()
    for key in ("nominal", "items", "criteria", "limits", "expansion"):
        assert json.dumps(a[key], sort_keys=True) == json.dumps(b[key], sort_keys=True), key


@pytest.mark.parametrize("make", [full, lambda: cancel_in("sensitivity")],
                         ids=["completed", "cancelled"])
def test_round_trip_through_the_result_format(make: Any) -> None:
    s = make()
    loaded = rt.results.load_json(s.to_json())
    assert loaded.type == "ToleranceSensitivity"
    assert loaded.schema_version == rt.results.SCHEMA_VERSION
    assert json.dumps(loaded.data, sort_keys=True, default=str) == json.dumps(
        s.to_dict(), sort_keys=True, default=str)


@pytest.mark.parametrize("make", [full, lambda: cancel_in("sensitivity")],
                         ids=["completed", "cancelled"])
def test_the_file_meets_the_schema(make: Any) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((REFERENCE_DIR.parent.parent / "schema"
                         / "raytatouille-result.schema.json").read_text(encoding="utf-8"))
    jsonschema.validate(json.loads(make().to_json()), schema)


def test_an_invalid_set_raises_before_any_trial() -> None:
    bad = rt.apply_patch(explicit_only(), [
        {"op": "replace", "path": "/tolerances/items/0/limit", "value": -1.0}],
        check="structure_only")
    with pytest.raises(rt.CompileError):
        tolerance.sensitivity(bad, materials=schott())
