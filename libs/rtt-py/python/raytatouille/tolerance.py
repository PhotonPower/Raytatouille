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

from . import _core
from ._core import (
    ClassSkip,
    ClassTolerance,
    ClassTolerances,
    ItemOrigin,
    ItemOriginSource,
    MaterialLibrary,
    OverrideEffect,
    OverrideStatus,
    System,
    ToleranceClass,
    ToleranceItems,
    ToleranceKind,
    ToleranceLevel,
)

__all__ = [
    "ClassSkip",
    "ClassTolerance",
    "ClassTolerances",
    "ItemOrigin",
    "ItemOriginSource",
    "OverrideEffect",
    "OverrideStatus",
    "ToleranceItems",
    "ToleranceLevel",
    "abbe_share",
    "class_limit",
    "expand",
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
