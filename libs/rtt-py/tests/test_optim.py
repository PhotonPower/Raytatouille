"""Optimization from Python (#169, ADR 0030): variables, merit function, optimize, the result
patch with an Editor, cancellation and progress, start errors."""

from __future__ import annotations

import json
import math
import threading

import numpy as np
import pytest
from conftest import REFERENCE_DIR
from optim_systems import gap_merit, singlet_merit, with_merit

import raytatouille as rt
from raytatouille.optim import OptimStatus

CONVERGED = {OptimStatus.CONVERGED_GRADIENT, OptimStatus.CONVERGED_STEP,
             OptimStatus.CONVERGED_MERIT}


def test_variables_in_the_order_of_adr_0030() -> None:
    # The row D of the table first, then the variable model Param (the image distance).
    gap = rt.optim.variables(rt.load(REFERENCE_DIR / "m5" / "two_lens_gap.rtt.json"))
    assert [v.pointer for v in gap][0] == "/parameters/0/value"
    assert (gap[0].row, gap[0].configuration, gap[0].param_index) == ("D", None, None)
    assert (gap[0].start, gap[0].min, gap[0].max) == (10.0, 1.0, None)
    assert len(gap) == 2 and gap[1].row == "" and gap[1].param_index is not None
    assert gap[1].pointer.endswith("/value") and gap[1].start == 45.0
    singlet = rt.optim.variables(singlet_merit())
    assert [v.start for v in singlet] == [103.36, -103.36, 107.3]


def test_merit_function_evaluates_the_operands() -> None:
    merit = rt.optim.MeritFunction(singlet_merit())
    assert merit.size == 1 + 2 * 3 * 6  # EFL, then (x, y) of 18 rays
    assert merit.generator_sizes == [36]
    start = merit.start()
    assert list(start) == [v.start for v in merit.variables]
    e = merit.evaluate(start)
    assert e.valid() and e.undefined == [""]
    efl = rt.paraxial.first_order(rt.compile(singlet_merit())).efl
    assert efl is not None and e.values[0] == efl  # exactly the analysis
    assert e.residuals[0] == e.values[0] - 100.0  # weight 1
    [g] = e.generators
    assert (g.rays_launched, g.rays_lost, g.undefined) == (18, 0, "")
    # The sum of squares of the generator residuals is its mean square (weight 1).
    assert math.isclose(float(np.sum(e.residuals[1:] ** 2)), g.mean_square, rel_tol=1e-12)
    with pytest.raises(ValueError, match="2 values for 3 variables"):
        merit.evaluate([1.0, 2.0])


def test_optimize_returns_a_new_system_and_its_patch() -> None:
    system = singlet_merit()
    before = system.to_json()
    result = rt.optim.optimize(system)
    assert result.status in CONVERGED
    assert system.to_json() == before  # the input is not changed
    assert result.system != system
    assert all(v.changed for v in result.variables)
    assert result.iterations == len(result.history) >= 1
    assert result.history.phi[-1] < result.history.phi[0] or result.iterations == 1
    # EFL 100 mm reached within a small fraction of the start deviation.
    [efl] = result.operands
    assert efl.pointer == "/optimization/operands/0" and abs(efl.value - 100.0) < 1e-3
    [spot] = result.generators
    assert spot.pointer == "/optimization/generators/0" and spot.rays_lost == 0
    # The patch turns the input into the result; an Editor applies it as one step.
    assert rt.apply_patch(system, result.patch).to_json() == result.system.to_json()
    editor = rt.Editor(system)
    editor.apply(result.patch)
    assert editor.system.to_json() == result.system.to_json()
    editor.undo()
    assert editor.system.to_json() == before
    ops = json.loads(result.patch)
    assert [op["path"] for op in ops] == [v.pointer for v in result.variables]
    assert [op["value"] for op in ops] == [v.end for v in result.variables]


def test_options_and_threads() -> None:
    none = rt.optim.optimize(singlet_merit(), options=rt.optim.OptimizeOptions(max_iterations=0))
    assert (none.status, none.patch, none.iterations) == (OptimStatus.MAX_ITERATIONS, "[]", 0)
    assert none.system == singlet_merit()
    defaults = rt.optim.OptimizeOptions()
    assert defaults.function_precision == rt.optim.MERIT_PRECISION == 1e-9
    assert defaults.max_iterations == 100 and defaults.gtol == 0.0
    one = rt.optim.optimize(gap_merit(), threads=1)
    four = rt.optim.optimize(gap_merit(), threads=4)
    assert one.patch == four.patch and one.evaluations == four.evaluations


def test_cancelled_before_the_start_returns_the_input() -> None:
    token = rt.CancelToken()
    token.cancel()
    system = singlet_merit()
    result = rt.optim.optimize(system, cancel=token)
    assert result.status == OptimStatus.CANCELLED  # a status, not an exception
    assert (result.patch, result.evaluations, result.iterations) == ("[]", 0, 0)
    assert result.system == system
    assert result.operands == [] and result.generators == []


def test_cancel_from_progress_keeps_the_last_accepted_state() -> None:
    token = rt.CancelToken()
    stages: list[str] = []
    lock = threading.Lock()

    def progress(done: int, total: int, stage: str) -> None:
        with lock:
            stages.append(stage)
        if stage == "optimize" and done >= 1:
            token.cancel()

    result = rt.optim.optimize(singlet_merit(), cancel=token, progress=progress)
    assert result.status == OptimStatus.CANCELLED
    assert {"jacobian", "optimize"} <= set(stages)
    # The state is a valid system: the patch still turns the input into it.
    assert rt.apply_patch(singlet_merit(), result.patch).to_json() == result.system.to_json()


def test_start_errors_carry_codes() -> None:
    with pytest.raises(rt.OptimError) as error:
        rt.optim.optimize(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json"))
    assert error.value.codes == ["optim.no_variables", "optim.no_operands"]
    assert isinstance(error.value, ValueError)
    with pytest.raises(rt.OptimError) as error:
        rt.optim.optimize(rt.load(REFERENCE_DIR / "m5" / "singlet_optim.rtt.json"))
    assert error.value.codes == ["optim.no_operands"]
    assert error.value.diagnostics[0].location == "/optimization"


def test_a_variable_at_its_bound_warns() -> None:
    # The image distance may not exceed 100 mm; the focus of EFL 100 mm lies behind it.
    data = json.loads(singlet_merit().to_json())
    image = data["root"]["children"][2]["pose"]["position"]
    image[2] = {"value": 99.0, "variable": True, "max": 100.0}
    system = rt.System.from_json(json.dumps(data))
    with pytest.warns(rt.RaytatouilleWarning) as caught:
        result = rt.optim.optimize(system)
    codes = [w.message.code for w in caught if isinstance(w.message, rt.RaytatouilleWarning)]
    assert "optim.parameter_at_bound" in codes
    [image_var] = [v for v in result.variables if v.pointer.startswith("/root/children/2")]
    assert image_var.at_bound and image_var.end == 100.0


def test_to_dict_leaves_out_the_system() -> None:
    result = rt.optim.optimize(gap_merit())
    data = result.to_dict()
    assert "system" not in data and data["patch"] == result.patch
    assert data["status"] == result.status.name
    assert list(data["history"]) == ["k", "phi", "mu", "rho", "step_norm", "accepted",
                                     "evaluations"]
    loaded = rt.results.load_json(result.to_json())
    assert loaded.type == "OptimResult" and loaded.data["patch"] == result.patch


def test_with_merit_helper_keeps_the_reference_untouched() -> None:
    # optim_systems builds the merit functions from the reference files; the files have none.
    assert rt.load(REFERENCE_DIR / "m5" / "singlet_optim.rtt.json").optimization.operands == []
    assert with_merit("m5/singlet_optim.rtt.json", {"operands": []}).optimization.operands == []
