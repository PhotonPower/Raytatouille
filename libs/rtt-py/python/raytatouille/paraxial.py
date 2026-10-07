"""First-order (paraxial) data and Seidel sums of a compiled system (rtt-paraxial).

Indices are absolute: in AIR (Ciddor air, n about 1.00027) EFL = 1/power is smaller by the
factor n_air than in programs that compute relative to air.
"""

from __future__ import annotations

from . import _core
from ._core import (
    ChromaticPair,
    FirstOrder,
    MaterialLibrary,
    Pupil,
    RayStart,
    Seidel,
    SeidelSurfaces,
    SeidelTerms,
    first_order,
)
from ._util import SystemLike, chromatic_pair, compiled

__all__ = [
    "ChromaticPair",
    "FirstOrder",
    "Pupil",
    "RayStart",
    "Seidel",
    "SeidelSurfaces",
    "SeidelTerms",
    "first_order",
    "seidel",
]


def seidel(
    system: SystemLike,
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    pair: ChromaticPair | tuple[int, int] | None = None,
    materials: MaterialLibrary | None = None,
) -> Seidel:
    """Seidel sums S_I ... S_V per event and for the path, in mm (W040 = S_I / 8, ...).

    ``system`` is a System (compiled with ``materials``) or a CompiledSystem; for several
    analyses compile once with rt.compile(). ``path`` is an index or a name, ``wavelength`` an
    index (None: reference). ``pair`` = (first, second) wavelength indices gives the colour
    terms C_L and C_T (n(first) - n(second), e.g. F and C); without it they are 0. Conventions
    (marginal ray at the paraxial entrance pupil, chief ray of the largest field, signs) in
    rtt/paraxial/seidel.hpp.

    Raises ParaxialError if the path is not rotationally symmetric (``surface`` and
    ``location`` name the place where known), NoStopError (also a ParaxialError) if it has no
    stop.
    """
    return _core.seidel(
        compiled(system, materials), path, wavelength, chromatic_pair(pair)
    )
