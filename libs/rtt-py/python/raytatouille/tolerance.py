"""Tolerancing (rtt-tolerance, ADR 0035 and 0036): the items of a tolerance set and the class
table.

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

import warnings

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
    CriterionValue,
    InverseSensitivity,
    ItemContribution,
    ItemOrigin,
    ItemOriginSource,
    ItemTrials,
    MaterialLibrary,
    OptimizeOptions,
    OverrideEffect,
    OverrideStatus,
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
    "CriterionValue",
    "InverseSensitivity",
    "ItemContribution",
    "ItemOrigin",
    "ItemOriginSource",
    "ItemTrials",
    "OverrideEffect",
    "OverrideStatus",
    "SensitivityStatus",
    "ToleranceItems",
    "ToleranceLevel",
    "ToleranceSensitivity",
    "TrialResult",
    "TrialStatus",
    "abbe_share",
    "class_limit",
    "expand",
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
    item whose target does not exist or does not fit its kind, OptimError for a criterion or
    compensator operand the merit function does not support.
    """
    result = _core.tolerance_sensitivity(
        system, materials, coatings, confirm,
        OptimizeOptions() if compensators is None else compensators, threads, cancel, progress)
    for d in result.warnings:
        if d.severity == Severity.WARNING:
            warnings.warn(RaytatouilleWarning(d.message, d.code, d.location), stacklevel=2)
    return result
