"""Monte Carlo from Python (#261, part 2b; ADR 0036, points 7 and 10): rt.tolerance.monte_carlo,
its result types, the TrialTable columns and the result format MonteCarloResult.

The system is the M7 reference achromat with the class switched off, as in
test_tolerance_sensitivity.py: three explicit items, the image position as compensator and three
criteria. N = 12 trials keep the run short; the seed is that of the file. The numbers are checked
against C++ bit for bit in test_bitwise_tolerance.py; here the statistics are recomputed from the
columns (ADR 0036, point 7) and the reproducibility, the cancellation and the options are
checked."""

from __future__ import annotations

import json
import math
import threading
from functools import cache
from typing import Any

import numpy as np
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR

import raytatouille as rt
from raytatouille import tolerance

FILE = REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json"
N = 12
VALID = tolerance.TrialStatus.VALID.value


def schott() -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    return lib


def explicit_only() -> rt.System:
    return rt.apply_patch(rt.load(FILE), [
        {"op": "replace", "path": "/tolerances/class", "value": "none"},
        {"op": "remove", "path": "/tolerances/overrides/4"},  # image none
        {"op": "remove", "path": "/tolerances/overrides/0"},  # objective precise
    ])


def file_seed() -> int:
    seed = json.loads(FILE.read_text(encoding="utf-8"))["tolerances"]["monte_carlo"]["seed"]
    assert isinstance(seed, int)
    return seed


def run(**options: Any) -> tolerance.MonteCarloResult:
    return tolerance.monte_carlo(explicit_only(), materials=schott(), trials=N, **options)


@cache
def full() -> tolerance.MonteCarloResult:
    return run()


def test_the_result_has_one_row_per_trial() -> None:
    r = full()
    assert r.status == tolerance.MonteCarloStatus.COMPLETED
    assert r.seed == file_seed() and r.requested == N
    assert r.nominal.status == tolerance.TrialStatus.VALID
    assert [i.id for i in r.items] == ["S3-wedge", "S1-conic", "R1"]
    assert r.expansion is not None and [i.id for i in r.expansion.items] == [i.id for i in r.items]
    assert [c.id for c in r.criteria] == ["spot-axis", "spot-zone", "efl"]
    t = r.trials
    assert len(t) == N
    assert t.trial.dtype == np.uint32 and t.trial.tolist() == list(range(N))
    assert t.values.shape == (N, 3) and t.values.dtype == np.float64
    assert t.criteria.shape == (N, 3) and t.satisfied.shape == (N, 3)
    assert t.satisfied.dtype == np.int8 and set(t.satisfied.ravel().tolist()) <= {-1, 0, 1}
    for name in ("compensator_start", "compensator_end", "at_bound", "infeasible"):
        assert getattr(t, name).shape == (N, 1), name
    assert t.status.shape == t.good.shape == t.lm_status.shape == t.mu.shape == (N,)
    assert t.good.dtype == np.bool_ and t.lm_iterations.dtype == np.int32
    assert len(t.reason) == N
    assert all((reason == "") == (status == VALID) for reason, status in zip(t.reason, t.status))
    # Content guards: the perturbations differ between the trials and lie within the limits
    # (ENDPOINT: exactly at +-T).
    limits = np.array([i.limit for i in r.items])
    assert len({tuple(row) for row in t.values.tolist()}) == N
    assert bool(np.all(np.abs(t.values) <= limits))
    assert bool(np.all(np.abs(t.values[:, 2]) == limits[2]))  # R1 is ENDPOINT


def wilson(good: int, total: int, z: float) -> tuple[float, float]:
    """Wilson score interval (NIST/SEMATECH e-Handbook, Sec. 7.2.4.1)."""
    p = good / total
    z2 = z * z
    centre = (p + z2 / (2 * total)) / (1 + z2 / total)
    half = z / (1 + z2 / total) * math.sqrt(p * (1 - p) / total + z2 / (4 * total * total))
    return centre - half, centre + half


def check_yield(y: tolerance.Yield, good: int, total: int, z: float) -> None:
    assert (y.good, y.total) == (good, total)
    assert y.fraction == pytest.approx(good / total, rel=1e-15)
    lower, upper = wilson(good, total, z)
    assert y.lower == (0.0 if good == 0 else pytest.approx(lower, rel=1e-12, abs=1e-15))
    assert y.upper == (1.0 if good == total else pytest.approx(upper, rel=1e-12))


def test_the_statistics_follow_from_the_columns() -> None:
    """ADR 0036, point 7: n = the values of the valid trials; mean, standard deviation (n - 1),
    min, max; percentiles by nearest rank r = ceil(p n); yields over all N with the Wilson
    interval."""
    r = full()
    t = r.trials
    valid = t.status == VALID
    check_yield(r.yield_, int(np.sum(t.good)), N, 1.959963984540054)
    assert t.good.tolist() == [bool(v) and all(s == 1 for s in row)
                               for v, row in zip(valid, t.satisfied.tolist())]
    for k, c in enumerate(r.criteria):
        column = t.criteria[:, k]
        assert c.with_value == int(np.sum(~np.isnan(column)))
        assert c.without_value == N - c.with_value
        values = column[valid & ~np.isnan(column)]
        assert c.count == len(values)
        assert c.count >= 2  # content guard: there are statistics
        assert c.mean == pytest.approx(float(np.mean(values)), rel=1e-12)
        assert c.standard_deviation == pytest.approx(float(np.std(values, ddof=1)), rel=1e-10)
        assert (c.min, c.max) == (float(np.min(values)), float(np.max(values)))
        ordered = np.sort(values)
        assert [x.p for x in c.percentiles] == [0.5, 0.9, 0.95, 0.99]
        for x in c.percentiles:
            assert x.rank == math.ceil(x.p * c.count)
            assert x.value == float(ordered[x.rank - 1])
        check_yield(c.yield_, int(np.sum(t.satisfied[:, k] == 1)), N, 1.959963984540054)


def test_good_trials_are_valid_with_every_criterion_satisfied() -> None:
    """In the reference case spot-zone (max 0.03 mm) already fails at the nominal (0.041 mm), so
    no trial is good. With max 0.0411 mm some are (review by the physik-reviewer, P3-3); the
    values of the seed lie at 0.04085 to 0.04141 mm, none closer than 2.6e-5 mm to the bound."""
    system = rt.apply_patch(explicit_only(), [
        {"op": "replace", "path": "/tolerances/criteria/1/max", "value": 0.0411}])
    r = tolerance.monte_carlo(system, materials=schott(), trials=N)
    t = r.trials
    valid = t.status == VALID
    expected = valid & np.all(t.satisfied == 1, axis=1)
    assert 0 < int(np.sum(expected)) < N  # content guard: good and bad trials
    assert t.good.tolist() == expected.tolist()
    assert bool(np.all(t.satisfied[:, 1] == (t.criteria[:, 1] <= 0.0411)))
    check_yield(r.yield_, int(np.sum(expected)), N, 1.959963984540054)
    check_yield(r.criteria[1].yield_, int(np.sum(t.satisfied[:, 1] == 1)), N, 1.959963984540054)


def test_numpy_integers_are_a_seed_and_a_number_of_trials() -> None:
    # Trial n depends only on the seed and n (ADR 0036, point 7, D2), so these are the first rows
    # of the full run.
    trials: Any = np.int64(2)
    seed: Any = np.uint64(file_seed())
    r = tolerance.monte_carlo(explicit_only(), materials=schott(), trials=trials, seed=seed)
    assert (r.seed, r.requested) == (file_seed(), 2)
    assert r.trials.values.tolist() == full().trials.values[:2].tolist()


def test_the_same_seed_gives_the_same_run_and_another_seed_another() -> None:
    assert run().to_json() == full().to_json()
    assert run(threads=1).to_json() == full().to_json()
    other = run(seed=file_seed() + 1)
    assert other.seed == file_seed() + 1
    assert not np.array_equal(other.trials.values, full().trials.values)


def test_without_keep_trials_only_the_invalid_trials_stay() -> None:
    """The kept rows are exactly the rows of the full run that are not valid. Without an invalid
    trial among the N the table is empty; that case is checked explicitly against the full run
    (review of the sketch, note a)."""
    lean = run(keep_trials=False)
    rows = full().trials
    invalid = rows.trial[rows.status != VALID]
    assert lean.trials.trial.tolist() == invalid.tolist()
    if len(invalid) == 0:
        assert len(lean.trials) == 0 and lean.trials.values.shape == (0, 3)
    else:
        mask = rows.status != VALID
        assert np.array_equal(lean.trials.values, rows.values[mask])
    # The statistics do not depend on keeping the trials.
    a, b = lean.to_dict(), full().to_dict()
    del a["trials"], b["trials"]
    assert json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True,
                                                                     default=str)


def test_invalid_trials_are_nan_and_minus_one_in_the_columns() -> None:
    """The reference case has no invalid trial, so the empty values of the columns (NaN, -1) are
    checked here: a wedge of up to 60 degrees on L1.S3 makes the chief ray miss the image in some
    trials, so that the compensator cannot run (review of #302, note c)."""
    system = rt.apply_patch(explicit_only(), [
        {"op": "replace", "path": "/tolerances/items/0/limit", "value": 60.0},
        {"op": "add", "path": "/tolerances/items/0/distribution", "value": "uniform"},
    ])
    with pytest.warns(rt.RaytatouilleWarning):
        r = tolerance.monte_carlo(system, materials=schott(), trials=6)
    t = r.trials
    invalid = t.status != VALID
    assert 0 < int(np.sum(invalid)) < len(t)  # content guard: both kinds of rows
    assert all(reason != "" for reason in np.array(t.reason)[invalid])
    assert bool(np.all(np.isnan(t.criteria[invalid])))
    assert bool(np.all(t.satisfied[invalid] == -1))
    assert bool(np.all(t.lm_status[invalid] == -1)) and bool(np.all(np.isnan(t.mu[invalid])))
    assert bool(np.all(t.lm_iterations[invalid] == 0))
    # No trial is good only because spot-zone (max 0.03 mm) already fails at the nominal
    # (0.041 mm), not because of the invalid rows; good trials are tested separately.
    assert not bool(np.any(t.good))
    assert not bool(np.any(np.isnan(t.criteria[~invalid])))
    assert bool(np.all(t.lm_status[~invalid] >= 0)) and not bool(np.any(np.isnan(t.mu[~invalid])))
    for c in r.criteria:
        assert (c.with_value, c.without_value) == (int(np.sum(~invalid)), int(np.sum(invalid)))
    with pytest.warns(rt.RaytatouilleWarning):
        lean = tolerance.monte_carlo(system, materials=schott(), trials=6, keep_trials=False)
    assert lean.trials.trial.tolist() == t.trial[invalid].tolist()
    loaded = rt.results.load_json(r.to_json())
    assert np.array_equal(loaded.data["trials"]["criteria"], t.criteria, equal_nan=True)
    assert np.array_equal(loaded.data["trials"]["mu"], t.mu, equal_nan=True)
    assert np.array_equal(loaded.data["trials"]["satisfied"], t.satisfied)
    assert np.array_equal(loaded.data["trials"]["lm_status"], t.lm_status)
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((REFERENCE_DIR.parent.parent / "schema"
                         / "raytatouille-result.schema.json").read_text(encoding="utf-8"))
    jsonschema.validate(json.loads(r.to_json()), schema)


def test_percentiles_and_z_are_options() -> None:
    r = run(percentiles=[0.25, 1.0], z=1.0)
    for c in r.criteria:
        assert [x.p for x in c.percentiles] == [0.25, 1.0]
        assert c.percentiles[1].rank == c.count  # p = 1 is the maximum
        assert c.percentiles[1].value == c.max
    check_yield(r.yield_, r.yield_.good, N, 1.0)


def test_progress_reports_one_stage() -> None:
    seen: dict[str, set[int]] = {}
    lock = threading.Lock()

    def progress(done: int, total: int, stage: str) -> None:
        with lock:
            seen.setdefault(stage, set()).add(total)

    run(progress=progress)
    assert seen == {"monte_carlo": {N}}


def test_cancelled_before_the_start_runs_no_trial() -> None:
    token = rt.CancelToken()
    token.cancel()
    r = run(cancel=token)
    assert r.status == tolerance.MonteCarloStatus.CANCELLED
    assert r.nominal.status == tolerance.TrialStatus.CANCELLED
    assert len(r.trials) == 0
    assert r.yield_.total == 0 and r.yield_.fraction is None


def test_a_cancelled_run_keeps_finished_trials_as_in_the_full_run() -> None:
    # Which trials finish depends on the timing; each finished one equals its row of the full
    # run (review of the sketch, note b).
    token = rt.CancelToken()

    def progress(done: int, total: int, stage: str) -> None:
        if done >= 1:
            token.cancel()

    r = run(threads=1, cancel=token, progress=progress)
    assert r.status == tolerance.MonteCarloStatus.CANCELLED
    rows = full().trials
    kept = r.trials.trial.tolist()
    assert 1 <= len(kept) < N
    assert kept == sorted(set(kept))
    assert np.array_equal(r.trials.values, rows.values[kept])
    assert np.array_equal(r.trials.criteria, rows.criteria[kept], equal_nan=True)
    assert r.yield_.total == len(kept)


@pytest.mark.parametrize("options,message", [
    ({"trials": 0}, "trials"),
    ({"trials": -1}, "trials"),
    ({"trials": 2**64}, "trials"),
    ({"seed": -1}, "seed"),
    ({"seed": 2**64}, "seed"),
    ({"percentiles": [0.0]}, "percentile"),
    ({"percentiles": [1.5]}, "percentile"),
    ({"z": 0.0}, "z"),
    ({"z": math.nan}, "z"),
])
def test_invalid_options_raise_value_error(options: dict[str, Any], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        tolerance.monte_carlo(explicit_only(), materials=schott(), **{"trials": N, **options})


def test_a_seed_or_number_of_trials_must_be_an_int() -> None:
    with pytest.raises(ValueError, match="seed must be an int"):
        tolerance.monte_carlo(explicit_only(), materials=schott(), trials=N,
                              seed=1.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="trials must be an int"):
        tolerance.monte_carlo(explicit_only(), materials=schott(),
                              trials=True)  # bool is an int to mypy, not to monte_carlo
    numpy_true: Any = np.True_
    for name in ("seed", "trials"):
        with pytest.raises(ValueError, match=f"{name} must be an int"):
            tolerance.monte_carlo(explicit_only(), materials=schott(),
                                  **{"trials": N, name: numpy_true})


def test_without_a_number_of_trials_it_is_an_error() -> None:
    system = rt.apply_patch(explicit_only(), [
        {"op": "remove", "path": "/tolerances/monte_carlo/trials"}], check="structure_only")
    with pytest.raises(ValueError, match="trials"):
        tolerance.monte_carlo(system, materials=schott())


@pytest.mark.parametrize("keep", [True, False], ids=["all_trials", "invalid_only"])
def test_round_trip_through_the_result_format(keep: bool) -> None:
    r = full() if keep else run(keep_trials=False)
    data = r.to_dict()
    assert "yield" in data and "yield_" not in data
    assert all("yield" in c and "yield_" not in c for c in data["criteria"])
    loaded = rt.results.load_json(r.to_json())
    assert loaded.type == "MonteCarloResult"
    assert loaded.schema_version == rt.results.SCHEMA_VERSION
    assert loaded.data["yield"]["good"] == r.yield_.good
    assert np.array_equal(loaded.data["trials"]["values"], r.trials.values)
    assert np.array_equal(loaded.data["trials"]["criteria"], r.trials.criteria, equal_nan=True)


@pytest.mark.parametrize("keep", [True, False], ids=["all_trials", "invalid_only"])
def test_the_file_meets_the_schema(keep: bool) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((REFERENCE_DIR.parent.parent / "schema"
                         / "raytatouille-result.schema.json").read_text(encoding="utf-8"))
    r = full() if keep else run(keep_trials=False)
    jsonschema.validate(json.loads(r.to_json()), schema)
