"""Optimization (rtt-optim, ADR 0030): the variables of a system, its merit function and the
Levenberg-Marquardt run.

The merit function is part of the system (``System.optimization``: operands and generators;
variables are the Params and table rows marked variable). ``optimize`` returns an OptimResult
with a new System; the input is not changed. Its ``patch`` (RFC 6902) turns the input into the
result, so an Editor can apply it as one step and take it back with undo::

    import raytatouille as rt

    system = rt.load("singlet_optim.rtt.json")
    result = rt.optim.optimize(system)
    print(result.status, result.history.phi[-1])
    editor = rt.Editor(system)
    editor.apply(result.patch)   # the optimized system, one undo step
    editor.undo()                # back to the input

Results are bitwise the same for every number of ``threads``. ``cancel`` (a CancelToken) ends
the run with status CANCELLED and the last accepted state, not with an exception;
``progress(done, total, stage)`` is called with the stages "jacobian" and "optimize" (as in
raytatouille.trace). Units of the values: those of the field or row (mm, um, degree, ...).
"""

from __future__ import annotations

import warnings

from . import _core
from ._core import (
    MERIT_PRECISION,
    CoatingLibrary,
    GeneratorStats,
    GeneratorValue,
    MaterialLibrary,
    MeritEvaluation,
    MeritFunction,
    OperandValue,
    OptimIterations,
    OptimizeOptions,
    OptimResult,
    OptimStatus,
    Severity,
    System,
    Variable,
    VariableValue,
)
from .errors import OptimError, RaytatouilleWarning
from .trace import CancelToken, ProgressCallback

__all__ = [
    "MERIT_PRECISION",
    "GeneratorStats",
    "GeneratorValue",
    "MeritEvaluation",
    "MeritFunction",
    "OperandValue",
    "OptimError",
    "OptimIterations",
    "OptimResult",
    "OptimStatus",
    "OptimizeOptions",
    "Variable",
    "VariableValue",
    "optimize",
    "variables",
]


def variables(system: System) -> list[Variable]:
    """The unknowns of ``system`` in the order of ADR 0030, point 5: the variable rows of the
    parameter table in table order (a ``values`` row once per configuration), then the
    variable model Params. Each has the JSON pointer its value has in the edit form, the start
    value and the bounds (None if not set)."""
    return _core.collect_variables(system)


def optimize(
    system: System,
    *,
    materials: MaterialLibrary | None = None,
    coatings: CoatingLibrary | None = None,
    options: OptimizeOptions | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> OptimResult:
    """Optimizes the variables of ``system`` against its merit function (ADR 0030) with
    Levenberg-Marquardt. ``materials`` (None: only VACUUM, AIR and CONST:) and ``coatings``
    as in rt.compile; ``options`` the solver settings (None: the defaults of OptimizeOptions).
    The system is copied before the run; the libraries are used by reference, so do not add
    catalogues to them from another thread during the run (as for every analysis).

    The result has the new ``system``, the ``patch`` from the input to it, the ``history`` of
    the solves (normalised merit ``phi`` per solve), the final values of the ``operands``,
    ``generators`` and ``variables`` and the ``diagnostics``. Its warnings (e.g.
    optim.parameter_at_bound, optim.rays_lost, rays.lost) are also issued as
    RaytatouilleWarning.

    Raises CompileError if the system is invalid or does not compile at the start; OptimError
    with optim.no_variables and/or optim.no_operands, or merit.operand_unsupported; the
    exception of an analysis at the start (ParaxialError, AnalysisError, NoStopError), or
    ValueError naming the operand or generator whose value is not defined at the start.
    """
    result = _core.optimize(system, materials, coatings,
                            OptimizeOptions() if options is None else options, threads, cancel,
                            progress)
    for d in result.diagnostics:
        if d.severity == Severity.WARNING:
            warnings.warn(RaytatouilleWarning(d.message, d.code, d.location), stacklevel=2)
    return result
