"""Optimization from Python is bitwise equal to C++ (#169).

rtt_py_reference (optim_cases.cpp) optimizes the M5 reference systems with the merit functions
of the M5 acceptance (#170) and evaluates the singlet's merit at its start; the flatten_*
functions below build the same arrays from the Python results. Both sides run with 1 and with 4
threads.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_EXE
from optim_systems import gap_merit, singlet_merit, solve_merit

import raytatouille as rt
from raytatouille.optim import OptimStatus

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def b1(values: Any) -> npt.NDArray[np.bool_]:
    return np.array(values, dtype=np.bool_)


# Flattening; names and order as write_result() and write_evaluation() in optim_cases.cpp.


def flatten_result(r: rt.optim.OptimResult) -> Arrays:
    h = r.history
    return {
        "k": np.array(h.k), "phi": np.array(h.phi), "mu": np.array(h.mu),
        "rho": np.array(h.rho), "step_norm": np.array(h.step_norm),
        "accepted": np.array(h.accepted), "evaluations": np.array(h.evaluations),
        "var_start": f8([v.start for v in r.variables]),
        "var_end": f8([v.end for v in r.variables]),
        "var_changed": b1([v.changed for v in r.variables]),
        "var_at_bound": b1([v.at_bound for v in r.variables]),
        "op_value": f8([o.value for o in r.operands]),
        "op_contribution": f8([o.contribution for o in r.operands]),
        "gen_rms": f8([g.rms for g in r.generators]),
        "gen_contribution": f8([g.contribution for g in r.generators]),
        "gen_rays": u8([n for g in r.generators for n in (g.rays_launched, g.rays_lost)]),
        "ints": u8([r.status.value, r.iterations, r.evaluations, r.failed_evaluations,
                    len(r.diagnostics)]),
        "patch": np.frombuffer(r.patch.encode("utf-8"), dtype=np.uint8),
    }


def flatten_evaluation(e: rt.optim.MeritEvaluation) -> Arrays:
    return {
        "values": np.array(e.values), "residuals": np.array(e.residuals),
        "gen_mean_square": f8([g.mean_square for g in e.generators]),
        "gen_rays": u8([n for g in e.generators for n in (g.rays_launched, g.rays_lost)]),
        "ints": u8([int(e.valid()), len(e.warnings)]),
    }


def merit_start(threads: int) -> Arrays:
    merit = rt.optim.MeritFunction(singlet_merit())
    return flatten_evaluation(merit.evaluate(merit.start(), threads=threads))


Case = Callable[[int], Arrays]

# Same calls as run_optim_cases() in optim_cases.cpp.
CASES: dict[str, Case] = {
    "singlet_optimize": lambda t: flatten_result(rt.optim.optimize(singlet_merit(), threads=t)),
    "gap_optimize": lambda t: flatten_result(rt.optim.optimize(gap_merit(), threads=t)),
    "solve_optimize": lambda t: flatten_result(rt.optim.optimize(solve_merit(), threads=t)),
    "singlet_merit_start": merit_start,
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_optim_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = CASES[name](threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_have_content() -> None:
    """The compared arrays carry content: both runs take accepted steps and change every
    variable, the merit falls, the patch is not empty, the start evaluation is valid."""
    for name in ("singlet_optimize", "gap_optimize", "solve_optimize"):
        a = CASES[name](1)
        assert len(a["k"]) >= 2 and bool(np.any(a["accepted"])), name
        assert bool(np.all(a["var_changed"])), name
        # phi of the accepted state never rises and falls overall (review of #195: strict).
        assert bool(np.all(np.diff(a["phi"]) <= 0.0)) and a["phi"][-1] < a["phi"][0], name
        assert bytes(a["patch"]) != b"[]", name
    start = CASES["singlet_merit_start"](1)
    assert int(start["ints"][0]) == 1 and len(start["residuals"]) == 1 + 2 * 3 * 6
    assert int(start["gen_rays"][0]) == 18 and float(start["gen_mean_square"][0]) > 0.0
    # Case 1 of the M5 acceptance (#170): the run converges (default options, so also by the
    # merit test), and the efl observer (weight 0) shows the EFL held by the parameter table.
    solve = CASES["solve_optimize"](1)
    converged = {OptimStatus.CONVERGED_GRADIENT.value, OptimStatus.CONVERGED_STEP.value,
                 OptimStatus.CONVERGED_MERIT.value}
    assert int(solve["ints"][0]) in converged
    assert abs(float(solve["op_value"][0]) - 100.0) <= 1e-8
    assert float(solve["op_contribution"][0]) == 0.0
