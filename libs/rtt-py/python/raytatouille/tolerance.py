"""Tolerancing (rtt-tolerance, ADR 0035 and 0036): the items of a tolerance set, the class
table, the sensitivity run and the ratio sigma / T of a distribution.

The tolerance set is part of the system (``System.tolerances``: class, overrides, explicit
items, compensators, criteria; edited by patch like the rest of the model). ``expand`` turns it
into the complete item list: the items the classes generate, in tree order, with the limit of a
single-value override or of the class table, then the explicit items. It also tells where each
item comes from and why a class skips a kind at a target::

    import raytatouille as rt

    system = rt.load("achromat_tolerances.rtt.json")
    materials = rt.MaterialLibrary()
    materials.add_catalog("schott.agf")
    r = rt.tolerance.expand(system, materials=materials)
    for item, origin in zip(r.items, r.origins):
        print(item.id, item.kind, item.limit, origin.source)

Units of the limits: those of the kind (mm for lengths, degree for the tilts, dimensionless for
curvature_relative, conic, index and abbe). The class limit of abbe is relative: a share of
nu_d of the glass (SCHOTT TIE-29, table 5-1), so ``class_limit`` rejects ABBE and
``abbe_share`` gives the share; the limit of a generated abbe item is abbe_share(class) * nu_d.
"""

from __future__ import annotations

import operator
import warnings
from collections.abc import Sequence

from . import _core
from ._core import (
    ClassSkip,
    ClassTolerance,
    ClassTolerances,
    CoatingLibrary,
    CompensatorValue,
    Confirmation,
    CriterionKind,
    CriterionSensitivity,
    CriterionStatistics,
    CriterionValue,
    InverseSensitivity,
    ItemContribution,
    ItemOrigin,
    ItemOriginSource,
    ItemTrials,
    MaterialLibrary,
    MonteCarloResult,
    MonteCarloStatus,
    OptimizeOptions,
    OverrideEffect,
    OverrideStatus,
    Percentile,
    SensitivityStatus,
    Severity,
    System,
    ToleranceClass,
    ToleranceDistribution,
    ToleranceItems,
    ToleranceKind,
    ToleranceLevel,
    ToleranceSensitivity,
    TrialResult,
    TrialStatus,
    TrialTable,
    Yield,
)
from .errors import RaytatouilleWarning
from .trace import CancelToken, ProgressCallback

__all__ = [
    "ClassSkip",
    "ClassTolerance",
    "ClassTolerances",
    "CompensatorValue",
    "Confirmation",
    "CriterionKind",
    "CriterionSensitivity",
    "CriterionStatistics",
    "CriterionValue",
    "InverseSensitivity",
    "ItemContribution",
    "ItemOrigin",
    "ItemOriginSource",
    "ItemTrials",
    "MonteCarloResult",
    "MonteCarloStatus",
    "OverrideEffect",
    "OverrideStatus",
    "Percentile",
    "SensitivityStatus",
    "ToleranceItems",
    "ToleranceLevel",
    "ToleranceSensitivity",
    "TrialResult",
    "TrialStatus",
    "TrialTable",
    "Yield",
    "abbe_share",
    "class_limit",
    "expand",
    "monte_carlo",
    "sensitivity",
    "sigma_ratio",
]


def expand(system: System, *, materials: MaterialLibrary | None = None) -> ToleranceItems:
    """The items of ``system.tolerances`` (ADR 0036, point 9).

    ``materials`` resolves the glasses of the generated abbe items, whose limit is a share of
    their nu_d; None means the default library (VACUUM, AIR, CONST only). The system is
    copied; ``materials`` must not change during the call.

    Raises CompileError with all errors if the tolerance set is invalid (validate), or if the
    glass of a generated abbe item is unknown (material.unknown) or gives no nu_d
    (tolerances.material_unsupported).
    """
    return _core.tolerance_expand(system, materials)


def class_limit(tolerance_class: ToleranceClass, kind: ToleranceKind,
                level: ToleranceLevel) -> float:
    """The limit T of a generated item of ``kind`` at ``level`` in ``tolerance_class``, in the
    unit of the kind (ADR 0036, addendum #260).

    Raises ValueError for the class NONE, for a kind no class generates at that level, and for
    ABBE, whose class limit is relative (use ``abbe_share``)."""
    return _core.tolerance_class_limit(tolerance_class, kind, level)


def abbe_share(tolerance_class: ToleranceClass) -> float:
    """The share p of nu_d that limits abbe in ``tolerance_class`` (SCHOTT TIE-29, table 5-1:
    COARSE step 4, NORMAL step 3, PRECISE step 1, i.e. 0.008, 0.005, 0.002; dimensionless). The
    limit of a generated abbe item is p * nu_d of its glass. Raises ValueError for NONE."""
    return _core.tolerance_abbe_share(tolerance_class)


def sigma_ratio(distribution: ToleranceDistribution) -> float:
    """c = sigma / T of a distribution (ADR 0036, points 5 and 6), dimensionless: NORMAL
    (truncated at +-2 sigma = T) 0.43981, UNIFORM 1/sqrt(3), ENDPOINT 1."""
    return _core.tolerance_sigma_ratio(distribution)


def sensitivity(
    system: System,
    *,
    materials: MaterialLibrary | None = None,
    coatings: CoatingLibrary | None = None,
    confirm: bool = True,
    compensators: OptimizeOptions | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> ToleranceSensitivity:
    """Sensitivity, RSS and inverse sensitivity of ``system.tolerances`` (ADR 0036, point 6).

    The items are those of ``expand(system, materials=materials)`` (kept as ``expansion``).
    Each item runs at +T_i and -T_i, all others at 0, after the compensator run (Levenberg-
    Marquardt with ``compensators``, None: the defaults of OptimizeOptions), with the nominal
    reference. Per criterion the result has Phi(+-T_i), the contributions (linear: ``delta``
    and ``shift``; spot_rms and opd_rms quadratic: ``q``), the RSS, the predicted spread and
    the inverse sensitivity; ``limits`` has the smallest T'_i over the criteria per item. With
    ``confirm`` every item runs again at +-T'_i (``confirmation``). Values are in the unit of
    the criterion's operand, limits in the unit of the item's kind; indices (``excluded``,
    ``unlimited``) are 0-based. The names are those of rtt/tolerance/sensitivity.hpp.

    The run is bitwise the same for every number of ``threads``. ``cancel`` (a CancelToken)
    does not raise: the result has status CANCELLED. Cancelled during Phi_0, ``nominal`` is
    CANCELLED and nothing else is computed; during the trials, only items with both sides are
    reduced, without RSS, inverse and confirmation; during the confirmation, the first part is
    complete and ``confirmation`` is None. ``progress`` gets the stages "sensitivity" and
    "confirmation", one work item per trial.

    The warnings of the run (tolerances.rss_negative) are also issued as RaytatouilleWarning;
    those of the single trials stay in their TrialResult. The system is copied; the libraries
    are used by reference and must not change during the run.

    Raises CompileError for an invalid tolerance set or system (as expand), ValueError for an
    item whose target does not exist or does not fit its kind, OptimError for an operand of a
    criterion or of the compensator merit (compensator_merit, else optimization) that the merit
    function does not support, and whatever the progress callback raises.
    """
    result = _core.tolerance_sensitivity(
        system, materials, coatings, confirm,
        OptimizeOptions() if compensators is None else compensators, threads, cancel, progress)
    for d in result.warnings:
        if d.severity == Severity.WARNING:
            warnings.warn(RaytatouilleWarning(d.message, d.code, d.location), stacklevel=2)
    return result


def _uint64(name: str, value: int | None) -> int | None:
    """`value` as an int; ValueError for a seed or a number of trials that is not an integer in
    [0, 2^64) (review of #261, F2), where nanobind would raise OverflowError or TypeError. NumPy
    integers are accepted (operator.index), bool is not. The range of the number of trials,
    [1, 2^32], is checked in C++."""
    if value is None:
        return None
    try:
        if isinstance(value, bool):
            raise TypeError
        number = operator.index(value)
    except TypeError:
        raise ValueError(f"{name} must be an int, not {type(value).__name__}") from None
    if not 0 <= number < 2**64:
        raise ValueError(f"{name} must lie in [0, 2^64), got {number}")
    return number


def monte_carlo(
    system: System,
    *,
    materials: MaterialLibrary | None = None,
    coatings: CoatingLibrary | None = None,
    trials: int | None = None,
    seed: int | None = None,
    percentiles: Sequence[float] = (0.5, 0.9, 0.95, 0.99),
    z: float = 1.959963984540054,
    keep_trials: bool = True,
    compensators: OptimizeOptions | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> MonteCarloResult:
    """Monte Carlo of ``system.tolerances`` (ADR 0036, point 7).

    The items are those of ``expand(system, materials=materials)`` (kept as ``expansion``).
    Trial n (0-based) perturbs every item with the draws of ``seed`` and n (ADR 0036, point 5)
    and evaluates the criteria after the compensator run (``compensators``, None: the defaults
    of OptimizeOptions). ``trials`` and ``seed`` default to tolerances.monte_carlo of the file
    (the seed else 0; without a number of trials it is an error). ``percentiles`` are p in
    (0, 1] (nearest rank, r = ceil(p n)); ``z`` sets the Wilson interval of the yields (default
    95 %). With ``keep_trials`` every trial is kept in ``trials`` (a TrialTable of NumPy
    columns, one row per trial in ascending number), otherwise only those that are not valid.

    Values are in the unit of the criterion's operand, perturbations in the unit of the item's
    kind. The yield is ``yield_`` (``yield`` is a Python keyword; "yield" in the result format).
    The names are those of rtt/tolerance/monte_carlo.hpp.

    The run is bitwise the same for every number of ``threads`` and every repetition with the
    same seed. ``cancel`` does not raise: the result has status CANCELLED and holds the trials
    that finished, which depends on the timing (not reproducible). Cancelled during Phi_0, no
    trial runs. ``progress`` gets the stage "monte_carlo", one work item per trial.

    The warnings of the trials, merged to one per code and location with the number of trials
    and the first trial (``warnings``), are also issued as RaytatouilleWarning. The system is
    copied; the libraries are used by reference and must not change during the run.

    Raises ValueError for a seed or number of trials that is not an int in [0, 2^64), for a
    number of trials neither given nor in the file or outside [1, 2^32], a percentile outside
    (0, 1], a z that is not finite and positive, and as ``sensitivity`` (CompileError,
    ValueError, OptimError) for the tolerance set; and whatever the progress callback raises.
    """
    result = _core.tolerance_monte_carlo(
        system, materials, coatings, _uint64("trials", trials), _uint64("seed", seed), [float(p) for p in percentiles], z,
        keep_trials, OptimizeOptions() if compensators is None else compensators, threads,
        cancel, progress)
    for d in result.warnings:
        if d.severity == Severity.WARNING:
            warnings.warn(RaytatouilleWarning(d.message, d.code, d.location), stacklevel=2)
    return result
