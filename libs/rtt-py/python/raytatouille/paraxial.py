"""First-order (paraxial) data, prescription data and Seidel sums (rtt-paraxial).

Indices are absolute: in AIR (Ciddor air, n about 1.00027) EFL = 1/power is smaller by the
factor n_air than in programs that compute relative to air.
"""

from __future__ import annotations

from . import _core
from ._core import (
    ChromaticPair,
    FirstOrder,
    MaterialLibrary,
    Prescription,
    PrescriptionSurfaces,
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
    "Prescription",
    "PrescriptionSurfaces",
    "Pupil",
    "RayStart",
    "Seidel",
    "SeidelSurfaces",
    "SeidelTerms",
    "first_order",
    "prescription",
    "seidel",
]


def prescription(
    system: SystemLike,
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    materials: MaterialLibrary | None = None,
) -> Prescription:
    """Paraxial prescription data of a path, as in a prescription report.

    Per event of the path (``surfaces``, NumPy copies): vertex z, signed index n after the
    event, marginal ray y, u, i and chief ray y_bar, u_bar, i_bar (i = u + y c with u before
    the event, rad), and the Lagrange invariant after the event. System data: total_track
    (unfolded with mirrors), object_distance, paraxial_working_f_number = 1 / (2 |n' u'|),
    paraxial_image_na = |n' u'|, lagrange_invariant = n (u_bar y - u y_bar) and first_order
    (lateral and angular magnification). The rays are those of seidel().

    Without a stop or a usable entrance pupil the affected values are None, and NaN in the
    arrays, without an error (as first_order()). ``system`` is a System (compiled with
    ``materials``) or a CompiledSystem; ``path`` an index or a name, ``wavelength`` an index
    (None: reference). Conventions in rtt/paraxial/prescription.hpp.

    Raises ParaxialError if the path is not rotationally symmetric, the path or wavelength
    does not exist, or the field definition is invalid for the chief ray.
    """
    return _core.prescription(compiled(system, materials), path, wavelength)


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
