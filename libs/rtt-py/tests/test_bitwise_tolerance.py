"""Tolerancing from Python is bitwise equal to C++ (#261).

rtt_py_reference (tolerance_cases.cpp) expands the tolerance set of the M7 reference achromat
with the SCHOTT catalogue and runs the sensitivity of its three explicit items (class none);
flatten_items and flatten_sensitivity below build the same arrays from the Python results.
The abbe limits (a share of nu_d of the glass) are only checked exactly here, because nu_d is
not available in Python. expand runs on one thread, sensitivity on the given threads; the cases
are compared for both runs of the reference program (1 and 4 threads).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import model, tolerance

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def text(lines: list[str]) -> npt.NDArray[np.uint8]:
    return np.frombuffer("\n".join(lines).encode("utf-8"), dtype=np.uint8)


def optional_index(i: int | None) -> int:
    return 0 if i is None else i + 1


# Flattening; names and order as write_items() in tolerance_cases.cpp.


def flatten_items(r: tolerance.ToleranceItems) -> Arrays:
    return {
        "ids": text([i.id for i in r.items]),
        "limit": f8([i.limit for i in r.items]),
        "kind": u8([i.kind.value for i in r.items]),
        "distribution": u8([i.distribution.value for i in r.items]),
        "origin": u8([n for o in r.origins for n in (o.source.value, o.index)]),
        "class_ids": text([c.id for c in r.classes.items]),
        "class_ints": u8([n for c in r.classes.items for n in (
            c.kind.value, c.level.value, c.tolerance_class.value, c.skip.value,
            optional_index(c.value_override), optional_index(c.replaced_by))]),
        "overrides": u8([n for v in r.classes.value_overrides
                         for n in (v.index, v.effect.value, v.skip.value)]),
    }


def achromat(threads: int) -> Arrays:
    del threads  # expand is serial
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    system = rt.load(REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json")
    return flatten_items(tolerance.expand(system, materials=lib))


def value(x: float | None) -> float:
    return math.nan if x is None else x


def add_trial(v: list[float], t: tolerance.TrialResult | None) -> None:
    if t is None:
        return
    v.extend(value(c.value) for c in t.criteria)
    v.extend(c.end for c in t.compensators)
    v.append(value(t.mu))
    v.append(float(t.lm_iterations))


def criterion_values(cs: list[tolerance.CriterionSensitivity]) -> list[float]:
    v: list[float] = []
    for c in cs:
        v.extend(value(x) for x in (c.nominal, c.rss, c.predicted_sigma,
                                    c.predicted_mean_square, c.inverse.allowance,
                                    c.inverse.share))
        for i in c.items:
            v.extend(value(x) for x in (i.plus, i.minus, i.worst, i.delta, i.shift, i.q))
        v.extend(value(x) for x in c.inverse.limits)
    return v


def criterion_ints(cs: list[tolerance.CriterionSensitivity]) -> list[int]:
    v: list[int] = []
    for c in cs:
        v += [c.kind.value, int(c.complete), len(c.excluded), len(c.inverse.unlimited)]
        v += c.excluded + c.inverse.unlimited
    return v


def flatten_sensitivity(r: tolerance.ToleranceSensitivity) -> Arrays:
    trials: list[float] = []
    add_trial(trials, r.nominal)
    for i in r.items:
        add_trial(trials, i.plus)
        add_trial(trials, i.minus)
    status = [r.status.value, r.nominal.status.value, int(r.confirmation is not None)]
    status += [0 if t is None else t.status.value + 1 for i in r.items for t in (i.plus, i.minus)]
    confirmation: list[float] = []
    if r.confirmation is not None:
        confirmation += r.confirmation.limits
        for plus, minus in zip(r.confirmation.plus, r.confirmation.minus):
            add_trial(confirmation, plus)
            add_trial(confirmation, minus)
        confirmation += criterion_values(r.confirmation.criteria)
    return {
        "ids": text([i.item.id for i in r.items]),
        "sigma_ratio": f8([i.sigma_ratio for i in r.items]),
        "trials": f8(trials),
        "criteria": f8(criterion_values(r.criteria)),
        "criterion_ints": u8(criterion_ints(r.criteria)),
        "limits": f8([value(x) for x in r.limits]),
        "status": u8(status),
        "confirmation": f8(confirmation),
    }


def add_yield(v: list[float], y: tolerance.Yield) -> None:
    v += [float(y.good), float(y.total), value(y.fraction), value(y.lower), value(y.upper)]


def flatten_monte_carlo(r: tolerance.MonteCarloResult) -> Arrays:
    t = r.trials
    ints = [r.seed, r.requested, r.status.value, len(t)]
    statistics: list[float] = []
    add_yield(statistics, r.yield_)
    for c in r.criteria:
        ints += [c.with_value, c.without_value, c.count, len(c.percentiles)]
        statistics += [value(x) for x in (c.mean, c.standard_deviation, c.min, c.max)]
        for p in c.percentiles:
            ints.append(p.rank)
            statistics += [p.p, p.value]
        add_yield(statistics, c.yield_)
    rows: list[float] = []
    flags: list[int] = []
    for i in range(len(t)):
        lm = int(t.lm_status[i])
        flags += [int(t.trial[i]), int(t.status[i]), int(t.good[i]), 0 if lm < 0 else lm + 1]
        rows += t.values[i].tolist()
        for k in range(t.criteria.shape[1]):
            rows.append(float(t.criteria[i, k]))
            flags.append(int(t.satisfied[i, k]) + 1)
        for c in range(t.compensator_start.shape[1]):
            rows += [float(t.compensator_start[i, c]), float(t.compensator_end[i, c])]
            flags += [int(t.at_bound[i, c]), int(t.infeasible[i, c])]
        rows += [float(t.mu[i]), float(t.lm_iterations[i])]
    return {"ints": u8(ints), "statistics": f8(statistics), "rows": f8(rows), "flags": u8(flags)}


@cache
def monte_carlo_run() -> tolerance.MonteCarloResult:
    system = explicit_only(rt.load(REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json"))
    return tolerance.monte_carlo(system, materials=schott_library(), trials=12)


def achromat_monte_carlo(threads: int) -> Arrays:
    # Bitwise the same for every number of threads; C++ runs with 1 and 4.
    return flatten_monte_carlo(monte_carlo_run())


def explicit_only(system: rt.System) -> rt.System:
    # As explicit_only() in tolerance_cases.cpp and test_tolerance_sensitivity.py.
    return rt.apply_patch(system, [
        {"op": "replace", "path": "/tolerances/class", "value": "none"},
        {"op": "remove", "path": "/tolerances/overrides/4"},
        {"op": "remove", "path": "/tolerances/overrides/0"},
    ])


@cache
def sensitivity_run(threads: int) -> tolerance.ToleranceSensitivity:
    # Cached: each run takes seconds; the arrays are rebuilt from the same result.
    system = explicit_only(rt.load(REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json"))
    return tolerance.sensitivity(system, materials=schott_library(), threads=threads)


def achromat_sensitivity(threads: int) -> Arrays:
    return flatten_sensitivity(sensitivity_run(threads))


Case = Callable[[int], Arrays]

# Same calls as run_tolerance_cases() in tolerance_cases.cpp.
CASES: dict[str, Case] = {
    "achromat_expand": achromat,
    "achromat_sensitivity": achromat_sensitivity,
    "achromat_monte_carlo": achromat_monte_carlo,
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_tolerance_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = CASES[name](threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_have_content() -> None:
    """26 items (23 generated, 3 explicit); the two abbe limits differ (a crown and a flint);
    at least three different skip reasons; both single values applied."""
    a = CASES["achromat_expand"](1)
    assert len(a["limit"]) == 26
    abbe = a["limit"][a["kind"] == model.ToleranceKind.ABBE.value]
    assert len(abbe) == 2 and abbe[0] != abbe[1] and bool(np.all(abbe > 0.0))
    skips = a["class_ints"].reshape(-1, 6)[:, 3]
    assert len(set(skips.tolist()) - {tolerance.ClassSkip.NONE.value}) >= 3
    applied = a["overrides"].reshape(-1, 3)[:, 1]
    assert applied.tolist() == [tolerance.OverrideEffect.APPLIED.value] * 2


def test_the_sensitivity_case_has_content() -> None:
    """Completed with confirmation; every criterion complete with an RSS > 0; some item limited;
    the two sides differ; and a sign: efl grows with R1 (delta > 0). R1 = 61.5 mm is the convex
    front surface of the crown, phi_1 = (n - 1) / R1, so d phi_1 / d R1 = -(n - 1) / R1^2 < 0:
    a larger R1 lowers the power and raises efl = 1 / Phi; the image position (the compensator)
    does not change efl. An equal sign error on both sides would pass the other checks."""
    s = sensitivity_run(1)
    assert s.status == tolerance.SensitivityStatus.COMPLETED
    assert s.confirmation is not None
    assert all(c.complete and c.rss is not None and c.rss > 0.0 for c in s.criteria)
    assert any(x is not None for x in s.limits)
    assert any(c.items[i].plus != c.items[i].minus for c in s.criteria for i in range(3))
    efl = next(c for c in s.criteria if c.id == "efl")
    r1 = [t.item.id for t in s.items].index("R1")
    delta = efl.items[r1].delta
    assert delta is not None and delta > 0.0


def schott_library() -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    return lib


def test_the_monte_carlo_case_has_content() -> None:
    """Completed, all 12 trials kept, each criterion with statistics, and the perturbations of
    every trial different."""
    r = monte_carlo_run()
    assert r.status == tolerance.MonteCarloStatus.COMPLETED
    assert len(r.trials) == 12 and r.yield_.total == 12
    assert all(c.count >= 2 and c.standard_deviation is not None for c in r.criteria)
    assert len({tuple(row) for row in r.trials.values.tolist()}) == 12
