"""M7 acceptance from Python is bitwise equal to C++ (#262).

rtt_py_reference (m7_cases.cpp) runs the sensitivity and the Monte Carlo of the M7 acceptance
(test_m7_acceptance.cpp, cases 4, 6, 7, 8 and 13) on the reference systems under
tests/reference/m7; flatten_sensitivity and flatten_monte_carlo below build the same arrays from
the Python results. Python binds only the overloads with a System, so both sides restrict the
items in the system (keep_items). The cases are compared for both runs of the reference program
(1 and 4 threads), and each has a content guard: a comparison of two equal wrong results would
pass otherwise (AGENTS.md).
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
from raytatouille import tolerance

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def value(x: float | None) -> float:
    return math.nan if x is None else x


@cache
def schott() -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    return lib


def load(name: str) -> rt.System:
    return rt.load(REFERENCE_DIR / "m7" / name)


def keep_items(system: rt.System, ids: list[str]) -> rt.System:
    """The system with only the explicit items `ids` (as keep_items() in m7_cases.cpp)."""
    items = system.to_dict()["tolerances"]["items"]
    remove = [i for i, item in enumerate(items) if item["id"] not in ids]
    return rt.apply_patch(system, [{"op": "remove", "path": f"/tolerances/items/{i}"}
                                   for i in reversed(remove)])


# Flattening; names and order as write_sensitivity() and write_monte_carlo() in m7_cases.cpp.


def flatten_sensitivity(r: tolerance.ToleranceSensitivity) -> Arrays:
    v: list[float] = []
    for c in r.criteria:
        v += [value(x) for x in (c.nominal, c.rss, c.predicted_sigma, c.predicted_mean_square)]
        for i in c.items:
            v += [value(x) for x in (i.delta, i.q, i.plus, i.minus)]
    status = [r.status.value, r.nominal.status.value]
    status += [0 if t is None else t.status.value + 1 for i in r.items for t in (i.plus, i.minus)]
    return {"criteria": f8(v), "status": u8(status)}


def add_yield(f: list[float], u: list[int], y: tolerance.Yield) -> None:
    f += [value(y.fraction), value(y.lower), value(y.upper)]
    u += [y.good, y.total]


def flatten_monte_carlo(r: tolerance.MonteCarloResult) -> Arrays:
    t = r.trials
    n = len(t.trial)
    trial_ints = np.stack([t.trial.astype(np.uint64), t.status.astype(np.uint64),
                           t.good.astype(np.uint64)], axis=1).reshape(-1) if n else u8([])
    statistics: list[float] = []
    counts: list[int] = []
    for c in r.criteria:
        statistics += [value(x) for x in (c.mean, c.standard_deviation, c.min, c.max)]
        for p in c.percentiles:
            statistics.append(p.value)
            counts.append(p.rank)
        counts += [c.with_value, c.without_value, c.count]
        add_yield(statistics, counts, c.yield_)
    add_yield(statistics, counts, r.yield_)
    counts += [r.status.value, r.seed, r.requested]
    return {
        "trial_ints": u8(trial_ints),
        "values": f8(t.values).reshape(-1),
        "criteria": f8(t.criteria).reshape(-1),
        "compensators": f8(t.compensator_end).reshape(-1),
        "statistics": f8(statistics),
        "counts": u8(counts),
    }


# The runs, cached: each takes seconds, and the content guards reuse them.


@cache
def sensitivity_lens(threads: int) -> tolerance.ToleranceSensitivity:
    system = keep_items(load("decentred_lens_finite.rtt.json"), ["lens-dy", "image-dy"])
    return tolerance.sensitivity(system, materials=rt.MaterialLibrary(), confirm=False,
                                 threads=threads)


@cache
def sensitivity_bk7(threads: int) -> tolerance.ToleranceSensitivity:
    return tolerance.sensitivity(load("planoconvex_bk7.rtt.json"), materials=schott(),
                                 confirm=False, threads=threads)


@cache
def mc_lens(threads: int) -> tolerance.MonteCarloResult:
    return tolerance.monte_carlo(load("decentred_lens.rtt.json"), materials=rt.MaterialLibrary(),
                                 threads=threads)


@cache
def mc_compensator(threads: int, seed: int | None = None) -> tolerance.MonteCarloResult:
    return tolerance.monte_carlo(load("paraboloid_compensator.rtt.json"),
                                 materials=rt.MaterialLibrary(), trials=20, seed=seed,
                                 threads=threads)


@cache
def mc_classes(threads: int) -> tolerance.MonteCarloResult:
    return tolerance.monte_carlo(load("achromat_tolerances.rtt.json"), materials=schott(),
                                 trials=10, threads=threads)


Case = Callable[[int], Arrays]

# Same calls as run_m7_cases() in m7_cases.cpp.
CASES: dict[str, Case] = {
    "m7_sensitivity_lens": lambda t: flatten_sensitivity(sensitivity_lens(t)),
    "m7_sensitivity_bk7": lambda t: flatten_sensitivity(sensitivity_bk7(t)),
    "m7_mc_lens": lambda t: flatten_monte_carlo(mc_lens(t)),
    "m7_mc_compensator": lambda t: flatten_monte_carlo(mc_compensator(t)),
    "m7_mc_classes": lambda t: flatten_monte_carlo(mc_classes(t)),
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_m7_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = CASES[name](threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_sensitivity_cases_have_content() -> None:
    """Case 6: delta = (1 - m) T = 0.15 for the lens and -T = -0.05 for the detector (sign);
    case 4: efl-F < efl-d < efl-C (normal dispersion), all trials valid."""
    lens = sensitivity_lens(1)
    assert lens.status == tolerance.SensitivityStatus.COMPLETED
    image_y = next(c for c in lens.criteria if c.id == "image-y")
    deltas = [i.delta for i in image_y.items]
    assert deltas[0] == pytest.approx(0.15, abs=1e-12)
    assert deltas[1] == pytest.approx(-0.05, abs=1e-12)
    bk7 = sensitivity_bk7(1)
    assert bk7.status == tolerance.SensitivityStatus.COMPLETED
    nominal = [c.nominal for c in bk7.criteria]
    assert all(x is not None for x in nominal)
    assert nominal[0] < nominal[1] < nominal[2]  # type: ignore[operator]
    assert all(t is not None and t.status == tolerance.TrialStatus.VALID
               for i in bk7.items for t in (i.plus, i.minus))


def test_the_monte_carlo_cases_have_content() -> None:
    """Case 8: the yield of |y| <= T/2 lies strictly between 0 and N. Case 7: the compensator ran
    in every trial, and another seed draws other values. Case 13: every item of the expansion
    takes part."""
    lens = mc_lens(1)
    assert lens.status == tolerance.MonteCarloStatus.COMPLETED
    image_y = next(c for c in lens.criteria if c.id == "image-y")
    assert 0 < image_y.yield_.good < image_y.yield_.total == 200
    comp = mc_compensator(1)
    assert len(comp.trials.trial) == 20
    assert bool(np.all(comp.trials.lm_status >= 0))  # -1: no run
    other = mc_compensator(1, seed=26206)
    assert not np.array_equal(other.trials.values, comp.trials.values)
    classes = mc_classes(1)
    assert classes.expansion is not None
    assert len(classes.items) == len(classes.expansion.items) > 3
    assert classes.trials.values.shape == (10, len(classes.items))
