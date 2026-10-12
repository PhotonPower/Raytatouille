"""Analyses as data objects (rtt-analysis): spot diagram, ray fans, OPD, longitudinal and
lateral colour, distortion, field curvature and Seidel sums; path transmission and OPL
difference of two paths (#122) and the ghost ranking (#124); reports as data (#177): the
raytrace report, the system data report and the dimension report (CSV with to_csv(), see
raytatouille.reports); the OPD at arbitrary pupil points, the pupil function with refocus and
the Zernike fit with the Strehl approximation (ADR 0032; #202, #203, #209); the diffraction
PSF, the Strehl ratio through focus, the MTF and the encircled energy (ADR 0033; #205, #206, #209).

Every function takes a System (compiled on each call with ``materials``) or a CompiledSystem;
for several analyses compile once with rt.compile(). ``path`` is an index or a path name,
``field`` a field index, ``wavelength`` a wavelength index where None means the reference
wavelength (in spot(): polychromatic). ``threads`` limits the worker threads; results are
bitwise the same for every number of threads. The GIL is released during a computation; do not
change the columns of a RayBatch passed as ``start`` from another thread meanwhile (as for
rt.trace.trace).

The bundle and sweep analyses (spot, ray_fan, opd_map, opd_fan, opd_points, pupil_function,
longitudinal_colour,
lateral_colour, distortion, field_curvature) take ``cancel`` (a CancelToken; cancel() from
another thread raises raytatouille.errors.Cancelled after at most one block per worker: rays,
or one field point or wavelength for the sweeps) and ``progress(done, total, stage)``
(called with the GIL from any thread, at most every 50 ms; stages "aim" and "trace" per ray
bundle, "field" for the sweeps, "wavelength" for the colour analyses). Neither changes a
result; an exception of ``progress`` ends the analysis and is raised. ``progress`` must not
start a raytatouille computation: through work stealing it may be called again in the same
thread. The single-point functions (``*_at``) and seidel() take neither: they are fast.

Units and conventions as in C++ (rtt/analysis/*.hpp): coordinates on the image surface in mm
in its local x, y (the surface of the last path event), relative to the chief ray of the
reference wavelength; OPD in waves at the reference wavelength, W > 0 leading. Lists in the
results are read-only NumPy copies.

Pupil samplings for spot() are objects of raytatouille.trace or a shorthand string:
``"hexapolar:N"`` (rings), ``"grid:N"``, ``"fan_x:N"``, ``"fan_y:N"``, ``"random:N"`` or
``"random:N:SEED"`` (count, seed; seed 0 if omitted) and ``"single:PX,PY"``.
"""

from __future__ import annotations

import math
import re
import warnings
from collections.abc import Sequence
from typing import Any, TypeVar

import numpy as np
import numpy.typing as npt

from . import _core
from .errors import RaytatouilleWarning
from ._core import (
    ApertureKind,
    DimensionReport,
    DimensionSegments,
    DistortionPoint,
    DistortionSweep,
    EeCentre,
    EncircledEnergy,
    FanPoints,
    Field,
    FieldCurvaturePoint,
    FieldCurvatureSweep,
    Foci,
    GhostEntries,
    GhostRanking,
    GhostSystem,
    Interferogram,
    InterferogramArrival,
    InterferogramPoints,
    JonesBasis,
    JonesPupil,
    LateralColour,
    LongitudinalColour,
    MaterialLibrary,
    Mtf,
    MtfCurve,
    NollIndex,
    OpdFan,
    OpdMap,
    OpdPoints,
    OpdPupilPoints,
    OplDifferencePoints,
    PathOplDifference,
    PathRays,
    PathTransmission,
    Point2,
    Points2,
    Psf,
    PsfPeak,
    PsfWavelengths,
    PupilCells,
    PupilFunction,
    PupilMapping,
    RayFan,
    RayLosses,
    RaytraceReport,
    RaytraceRows,
    ReferenceSphere,
    SpotDiagram,
    SpotStatistics,
    StokesPupil,
    SystemReport,
    ThroughFocus,
    ZernikeFit,
)
from ._util import SystemLike, chromatic_pair, compiled
from .paraxial import ChromaticPair, seidel
from .trace import (
    Aiming,
    CancelToken,
    FanXPupil,
    FanYPupil,
    GridPupil,
    HexapolarPupil,
    ProgressCallback,
    PupilSampling,
    RandomPupil,
    RayBatch,
    SinglePupilPoint,
)

__all__ = [
    "ApertureKind",
    "DimensionReport",
    "DimensionSegments",
    "DistortionPoint",
    "DistortionSweep",
    "EeCentre",
    "EncircledEnergy",
    "FanPoints",
    "FieldCurvaturePoint",
    "FieldCurvatureSweep",
    "Foci",
    "GhostEntries",
    "GhostRanking",
    "Interferogram",
    "InterferogramArrival",
    "InterferogramPoints",
    "JonesBasis",
    "JonesPupil",
    "LateralColour",
    "LongitudinalColour",
    "MAX_PUPIL_GRID",
    "MAX_ZERNIKE_TERMS",
    "Mtf",
    "MtfCurve",
    "NollIndex",
    "OpdFan",
    "OpdMap",
    "OpdPoints",
    "OpdPupilPoints",
    "OplDifferencePoints",
    "PathOplDifference",
    "PathRays",
    "PathTransmission",
    "Point2",
    "Points2",
    "Psf",
    "PsfPeak",
    "PsfWavelengths",
    "PupilCells",
    "PupilFunction",
    "PupilMapping",
    "RayFan",
    "RayLosses",
    "RaytraceReport",
    "RaytraceRows",
    "ReferenceSphere",
    "SpotDiagram",
    "SpotStatistics",
    "StokesPupil",
    "SystemReport",
    "ThroughFocus",
    "ZernikeFit",
    "dimension_report",
    "distortion",
    "distortion_at",
    "encircled_energy",
    "field_curvature",
    "field_curvature_at",
    "ghost_ranking",
    "interferogram",
    "jones_pupil",
    "lateral_colour",
    "longitudinal_colour",
    "mtf",
    "noll_from_nm",
    "noll_to_nm",
    "opd_fan",
    "opd_map",
    "opd_points",
    "opl_difference",
    "path_transmission",
    "psf",
    "psf_at",
    "psf_positions",
    "pupil_cells",
    "pupil_function",
    "ray_fan",
    "raytrace_report",
    "refocus",
    "sampling",
    "seidel",
    "spot",
    "stokes_pupil",
    "system_report",
    "zernike",
    "zernike_fit",
    "zernike_radial",
]

#: Largest grid N of pupil_function and pupil_cells (memory bound, rtt/analysis/pupil.hpp).
MAX_PUPIL_GRID: int = _core.MAX_PUPIL_GRID
#: Largest number of Zernike terms (n <= 20 completely).
MAX_ZERNIKE_TERMS: int = _core.MAX_ZERNIKE_TERMS

_SHORTHAND = (
    "expected 'hexapolar:N', 'grid:N', 'fan_x:N', 'fan_y:N', 'random:N', 'random:N:SEED' or "
    "'single:PX,PY'"
)
_DIGITS = re.compile(r"[0-9]+")
_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
_INT32_MAX = 2**31 - 1
_UINT64_MAX = 2**64 - 1


def _integer(text: str, low: int, high: int) -> int:
    """Decimal digits only, low <= value <= high; ValueError otherwise."""
    if not _DIGITS.fullmatch(text):
        raise ValueError(text)
    value = int(text)
    if not low <= value <= high:
        raise ValueError(text)
    return value


def _number(text: str) -> float:
    """Finite decimal number; ValueError otherwise."""
    if not _NUMBER.fullmatch(text):
        raise ValueError(text)
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(text)
    return value


def sampling(rays: str | PupilSampling) -> PupilSampling:
    """Pupil sampling from a shorthand string (see module docstring); sampling objects are
    returned unchanged. Accepted values: rings >= 0, n >= 1, count and seed >= 0 (unsigned 64
    bit), finite PX, PY. Raises ValueError naming the shorthand for anything else."""
    if not isinstance(rays, str):
        return rays
    kind, _, arguments = rays.partition(":")
    parts = arguments.split(":") if arguments else []
    try:
        if kind == "hexapolar" and len(parts) == 1:
            return HexapolarPupil(rings=_integer(parts[0], 0, _INT32_MAX))
        if kind in ("grid", "fan_x", "fan_y") and len(parts) == 1:
            n = _integer(parts[0], 1, _INT32_MAX)
            if kind == "grid":
                return GridPupil(n=n)
            return FanXPupil(n=n) if kind == "fan_x" else FanYPupil(n=n)
        if kind == "random" and len(parts) in (1, 2):
            count = _integer(parts[0], 0, _UINT64_MAX)
            seed = _integer(parts[1], 0, _UINT64_MAX) if len(parts) == 2 else 0
            return RandomPupil(count=count, seed=seed)
        if kind == "single" and len(parts) == 1:
            values = parts[0].split(",")
            if len(values) == 2:
                return SinglePupilPoint(px=_number(values[0]), py=_number(values[1]))
    except ValueError as error:
        raise ValueError(f"invalid ray sampling {rays!r}: {_SHORTHAND}") from error
    raise ValueError(f"invalid ray sampling {rays!r}: {_SHORTHAND}")


_R = TypeVar("_R", SpotDiagram, RayFan, OpdMap, OpdFan, OpdPupilPoints, PupilFunction,
              ZernikeFit, PathTransmission, PathOplDifference, Interferogram, GhostRanking,
              SystemReport, JonesPupil, Psf, ThroughFocus, Mtf,
              EncircledEnergy)


def _warn(result: _R) -> _R:
    """Issues the warnings of `result` as RaytatouilleWarning (ADR 0022, 0023): rays.lost when
    more rays are lost than lost_warning_fraction (default one half; vignetting at the field
    edge is intended), stop.clips_beam when rays end vignetted at the stop. They stay in
    result.warnings, and result.losses has the counts in any case."""
    for d in result.warnings:
        warnings.warn(RaytatouilleWarning(d.message, d.code, d.location), stacklevel=3)
    return result


def spot(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    rays: str | PupilSampling = "hexapolar:6",
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> SpotDiagram:
    """Spot diagram of ``field``. ``wavelength`` None gives a POLYCHROMATIC spot over all
    system wavelengths with the model's wavelength weights (unlike the other analyses, where
    None means the reference wavelength). ``rays`` is the pupil sampling per wavelength.

    Raises AnalysisError if the chief ray does not reach the image surface or no ray arrives.

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.spot(
        compiled(system, materials), path, field, wavelength, sampling(rays), aiming,
        lost_warning_fraction, threads, cancel, progress
    ))


def ray_fan(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    points: int = 21,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> RayFan:
    """Tangential and sagittal ray fans of ``field`` with ``points`` points on [-1, 1]:
    transverse aberration relative to the chief ray of the reference wavelength, mm.

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.ray_fan(
        compiled(system, materials), path, field, wavelength, points, aiming,
        lost_warning_fraction, threads, cancel, progress
    ))


def opd_map(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    grid: int = 33,
    aiming: Aiming = Aiming.REAL,
    focus_shift: float = 0.0,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> OpdMap:
    """OPD map of ``field`` on a grid x grid pupil grid (points inside the unit circle), in
    waves at the reference wavelength, against the reference sphere centred on the chief ray.

    An exit pupil at infinity uses the limit of the reference sphere (radius inf).
    ``focus_shift`` (mm, finite) moves the centre of the reference sphere along the chief ray of
    the reference wavelength (through focus, ADR 0032 point 5; 0 gives the sphere of before).

    Raises AnalysisError e.g. for a path without stop.

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.opd_map(
        compiled(system, materials), path, field, wavelength, grid, aiming, focus_shift,
        lost_warning_fraction, threads, cancel, progress
    ))


def opd_fan(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    points: int = 21,
    aiming: Aiming = Aiming.REAL,
    focus_shift: float = 0.0,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> OpdFan:
    """Tangential (px = 0) and sagittal (py = 0) OPD fans of ``field``; as opd_map().

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.opd_fan(
        compiled(system, materials), path, field, wavelength, points, aiming, focus_shift,
        lost_warning_fraction, threads, cancel, progress
    ))


def opd_points(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    rays: str | PupilSampling = "hexapolar:6",
    aiming: Aiming = Aiming.REAL,
    focus_shift: float = 0.0,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> OpdPupilPoints:
    """OPD of ``field`` at the points of the pupil sampling ``rays`` (any sampling, also
    GaussPupil and PointListPupil, or a shorthand string), conventions as opd_map(); no
    statistics. The points come in the order of the sampling; check status before using w."""
    return _warn(_core.opd_points(
        compiled(system, materials), path, field, wavelength, sampling(rays), aiming, focus_shift,
        lost_warning_fraction, threads, cancel, progress
    ))


def pupil_cells(grid: int, obscuration: float = 0.0) -> PupilCells:
    """The N x N cells of the pupil function (ADR 0032, point 1): area fraction inside the unit
    circle and sample point (centre, or centroid of the part inside); rows from py = -1, px
    fastest. With ``obscuration`` = epsilon > 0 the cells cover the annulus epsilon <= |p| <= 1
    (addendum #229): the exact annular area, the centroid of the annular part, f = 0 inside the
    inner circle, and ``inner_rim`` = 1 for the cells on the inner rim (the pupil function traces
    them as virtual rays). ValueError for N outside [2, MAX_PUPIL_GRID] or obscuration outside
    [0, 1)."""
    return _core.pupil_cells(grid, obscuration)


def pupil_function(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    grid: int = 64,
    apodization: float = 0.0,
    focus_shift: float = 0.0,
    aiming: Aiming = Aiming.REAL,
    with_prt: bool = False,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    mapping_warning: float = 0.01,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> PupilFunction:
    """Pupil function of ``field`` at ``wavelength`` on grid x grid cells (ADR 0032, point 3):
    per cell the area fraction, sample point, amplitude sqrt(weight) exp(-apodization rho^2), W
    in waves at the reference wavelength and status; ``with_prt`` adds the accumulated PRT
    matrix per cell. Check ``fraction > 0`` first, then ``status``. The phase of a cell at
    wavelength lambda relative to the chief ray is -2 pi W lambda_ref / lambda.

    ``mapping`` is the affine mapping of the pupil onto the direction sines (None for an afocal
    image space); above ``mapping_warning`` its delta_max warns with pupil.mapping_nonlinear.

    Raises AnalysisError (location /paths/i) for a path with crystal modes, ValueError for
    options out of range. Warnings stay in ``warnings`` and are also issued as
    RaytatouilleWarning.
    """
    return _warn(_core.pupil_function(
        compiled(system, materials), path, field, wavelength, grid, apodization, focus_shift,
        aiming, with_prt, lost_warning_fraction, mapping_warning, threads, cancel, progress
    ))


def refocus(system: SystemLike, pupil: PupilFunction, focus_shift: float, *,
            materials: MaterialLibrary | None = None) -> PupilFunction:
    """``pupil`` for another ``focus_shift`` (mm) from its end states, without tracing again;
    bitwise the same as pupil_function() with that focus_shift. ``system`` must be the system
    that produced ``pupil`` (ValueError if it obviously is not)."""
    return _warn(_core.refocus(compiled(system, materials), pupil, focus_shift))


def zernike_fit(pupil: PupilFunction, *, terms: int = 37,
                incomplete_warning_fraction: float = 0.01) -> ZernikeFit:
    """Zernike fit of ``pupil`` with ``terms`` Noll terms (ADR 0032, point 6): coefficients in
    waves at the reference wavelength (index j - 1), residual RMS, sigma (piston removed),
    sigma_bf (about the diffraction focus) and the Strehl approximation after Marechal
    exp(-(2 pi sigma_bf lambda_ref / lambda)^2) (good down to about 0.1).

    Raises AnalysisError (location /paths/i) for too few arrived cells or singular normal
    equations, ValueError for options out of range. An incomplete pupil warns with
    zernike.pupil_incomplete.
    """
    return _warn(_core.zernike_fit(pupil, terms, incomplete_warning_fraction))


def zernike(j: int, rho: float, theta: float) -> float:
    """Orthonormal Zernike polynomial Z_j(rho, theta) after Noll (Niu/Tian Eq. (3)); theta in
    rad from +x counter-clockwise. ValueError for j outside [1, MAX_ZERNIKE_TERMS]."""
    return _core.zernike(j, rho, theta)


def zernike_radial(n: int, m: int, rho: float) -> float:
    """Radial polynomial R_n^m(rho) (Niu/Tian Eq. (4)) for 0 <= m <= n <= 20, n - m even."""
    return _core.zernike_radial(n, m, rho)


def noll_to_nm(j: int) -> NollIndex:
    """(n, m) of the Noll index j (m >= 0; sine for odd j)."""
    return _core.noll_to_nm(j)


def noll_from_nm(n: int, m: int, sine: bool = False) -> int:
    """Noll index of (n, m) with the sine (odd j) or the cosine (even j) for m != 0."""
    return _core.noll_from_nm(n, m, sine)



# ------------------------------------------------- diffraction (ADR 0033; #205, #206, P1b) ---


def psf(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    *,
    wavelengths: Sequence[int] | None = None,
    grid: int = 64,
    padding: int = 4,
    apodization: float = 0.0,
    focus_shift: float = 0.0,
    aiming: Aiming = Aiming.REAL,
    aliasing_factor: float = 3.0,
    undersampled_waves: float = 0.5,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    mapping_warning: float = 0.01,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> Psf:
    """Diffraction PSF of ``field`` by FFT of the pupil function (ADR 0033, points 2-4, 7).

    ``intensity`` is an (L, L) array with L = ``padding`` x ``grid`` (larger only for a shorter
    wavelength), row k2 + L // 2 and column k1 + L // 2, the fraction of the pupil energy per
    cell on the lattice x(k) = ``lattice`` @ k in mm along ``e_x``, ``e_y`` of the plane through
    C' normal to the chief ray (``psf_positions(result)`` gives all positions). ``strehl`` is
    the Strehl ratio at C'. ``wavelengths``: indices, or None for all system wavelengths with
    their model weights (each wavelength on its own lattice, weighted with w_i E_i).

    Raises AnalysisError (location /paths/i) for an afocal image space or a path with crystal
    modes, ValueError for options out of range. Warnings (rays.lost, stop.clips_beam,
    pupil.mapping_nonlinear, psf.aliasing, psf.pupil_undersampled) are issued as
    RaytatouilleWarning and stay in ``warnings``.
    """
    return _warn(_core.psf(
        compiled(system, materials), path, field,
        None if wavelengths is None else list(wavelengths), grid, padding, apodization,
        focus_shift, aiming, aliasing_factor, undersampled_waves, lost_warning_fraction,
        mapping_warning, threads, cancel, progress,
    ))


def psf_at(pupil: PupilFunction, points: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Intensity of the monochromatic PSF of ``pupil`` at ``points`` (shape (n, 2), mm along
    e_x, e_y of ``pupil.mapping``) by the direct sum over the cells, in Strehl units: 1 is the
    peak of the same pupil without aberration, the value at (0, 0) the Strehl ratio.

    Raises AnalysisError for an afocal pupil function or one without arrived cells, ValueError
    for points of another shape."""
    xy = np.ascontiguousarray(points, dtype=np.float64)
    if xy.ndim == 1 and xy.shape == (2,):
        xy = xy.reshape(1, 2)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"points must have the shape (n, 2), not {xy.shape}")
    return np.asarray(_core.psf_at(pupil, xy))


def psf_positions(result: Psf) -> npt.NDArray[np.float64]:
    """Positions of the lattice cells of ``result``, shape (L, L, 2) in mm along e_x, e_y
    relative to C', in the order of ``intensity``: [k2 + L // 2, k1 + L // 2] = lattice @ (k1,
    k2). Also available as ``Psf.positions()``."""
    k = np.arange(result.length, dtype=np.float64) - result.length // 2
    k1, k2 = np.meshgrid(k, k)  # k1 varies along a row, k2 along a column
    b = result.lattice
    out = np.empty((result.length, result.length, 2), dtype=np.float64)
    out[..., 0] = b[0, 0] * k1 + b[0, 1] * k2
    out[..., 1] = b[1, 0] * k1 + b[1, 1] * k2
    return out


def through_focus(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    wavelength: int | None = None,
    *,
    shifts: Sequence[float],
    grid: int = 64,
    apodization: float = 0.0,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    mapping_warning: float = 0.01,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> ThroughFocus:
    """Strehl ratio |sum u|^2 / |sum f A|^2 at ``wavelength`` for every focus shift in
    ``shifts`` (mm along the chief ray, ADR 0032, point 5): one trace, W recomputed per shift
    without tracing again (ADR 0033, points 8 and 11). ``arrived_fraction`` is the area of the
    cells that cross each shifted reference sphere.

    Raises ValueError for a non-finite shift, AnalysisError as psf()."""
    return _warn(_core.through_focus(
        compiled(system, materials), path, field, wavelength, [float(x) for x in shifts], grid,
        apodization, aiming, lost_warning_fraction, mapping_warning, threads, cancel, progress,
    ))


def _pupil_form(name: str, given: dict[str, Any]) -> None:
    """The form on a PupilFunction excludes the options of the system form."""
    named = [key for key, value in given.items() if value is not None]
    if named:
        raise ValueError(f"{', '.join(named)} only apply to a system, not to a PupilFunction "
                         f"({name} of a pupil function uses its own grid, field and wavelength)")


def mtf(
    system: SystemLike | PupilFunction,
    path: int | str | None = None,
    field: int | None = None,
    *,
    wavelengths: Sequence[int] | None = None,
    grid: int | None = None,
    map: bool = False,
    complex_otf: bool = False,
    apodization: float | None = None,
    focus_shift: float | None = None,
    aiming: Aiming | None = None,
    undersampled_waves: float | None = None,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float | None = None,
    mapping_warning: float | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> Mtf:
    """Diffraction MTF (ADR 0033, point 5) from the autocorrelation of the pupil function:
    ``tangential`` (shifts in py, meridional for a field on +-y) and ``sagittal`` (px), each
    with ``frequency`` (cycles per mm, up to the cut-off where the value is exactly 0) and
    ``modulation`` = |OTF|; ``complex_otf`` adds the OTF, ``map`` the 2D map over all shifts
    (one wavelength only). Polychromatic: the complex OTF of each wavelength interpolated onto
    the frequencies of the reference and summed with w_i E_i.

    ``system`` is a System or CompiledSystem (path 0, field 0, grid 64 and the defaults of
    pupil_function() unless given), or a PupilFunction (e.g. from refocus()): then only ``map``
    and ``complex_otf`` apply, the others raise ValueError.

    Raises AnalysisError (location /paths/i) for an afocal image space or a path with crystal
    modes, ValueError for options out of range. Warnings are issued as RaytatouilleWarning.
    """
    if isinstance(system, PupilFunction):
        _pupil_form("mtf", {"path": path, "field": field, "wavelengths": wavelengths,
                            "grid": grid, "apodization": apodization,
                            "focus_shift": focus_shift, "aiming": aiming,
                            "undersampled_waves": undersampled_waves, "materials": materials,
                            "threads": threads, "lost_warning_fraction": lost_warning_fraction,
                            "mapping_warning": mapping_warning, "cancel": cancel,
                            "progress": progress})
        return _warn(_core.mtf_of_pupil(system, map, complex_otf))
    return _warn(_core.mtf(
        compiled(system, materials), 0 if path is None else path, 0 if field is None else field,
        None if wavelengths is None else list(wavelengths), 64 if grid is None else grid,
        0.0 if apodization is None else apodization, 0.0 if focus_shift is None else focus_shift,
        Aiming.REAL if aiming is None else aiming, map, complex_otf,
        0.5 if undersampled_waves is None else undersampled_waves,
        0.5 if lost_warning_fraction is None else lost_warning_fraction,
        0.01 if mapping_warning is None else mapping_warning, threads, cancel, progress,
    ))


def encircled_energy(
    system: SystemLike | PupilFunction,
    radii: Sequence[float],
    path: int | str | None = None,
    field: int | None = None,
    *,
    centre: EeCentre = EeCentre.REFERENCE,
    wavelengths: Sequence[int] | None = None,
    grid: int | None = None,
    apodization: float | None = None,
    focus_shift: float | None = None,
    aiming: Aiming | None = None,
    undersampled_waves: float | None = None,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float | None = None,
    mapping_warning: float | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> EncircledEnergy:
    """Encircled energy at ``radii`` (mm in the PSF plane; ADR 0033, point 6): the fraction of
    the pupil energy inside circles about ``centre`` (EeCentre.REFERENCE: C';
    EeCentre.LINEAR_PHASE: the point given by the phase of the first harmonic, exact for a tilt,
    not the centroid). Radii up to ``radius_max`` (a quarter of the shortest period of the
    periodic PSF); larger, negative or non-finite radii raise ValueError.

    ``system`` as in mtf(): a System or CompiledSystem, or a PupilFunction (then only
    ``centre`` applies besides ``radii``). Unlike psf() and mtf(), ``radii`` is the second
    positional argument in both forms: encircled_energy(system, radii, path, field) and
    encircled_energy(pupil, radii)."""
    values = [float(r) for r in radii]
    if isinstance(system, PupilFunction):
        _pupil_form("encircled_energy", {
            "path": path, "field": field, "wavelengths": wavelengths, "grid": grid,
            "apodization": apodization, "focus_shift": focus_shift, "aiming": aiming,
            "undersampled_waves": undersampled_waves, "materials": materials,
            "threads": threads, "lost_warning_fraction": lost_warning_fraction,
            "mapping_warning": mapping_warning, "cancel": cancel, "progress": progress})
        return _warn(_core.encircled_energy_of_pupil(system, values, centre))
    return _warn(_core.encircled_energy(
        compiled(system, materials), 0 if path is None else path, 0 if field is None else field,
        values, centre, None if wavelengths is None else list(wavelengths),
        64 if grid is None else grid, 0.0 if apodization is None else apodization,
        0.0 if focus_shift is None else focus_shift, Aiming.REAL if aiming is None else aiming,
        0.5 if undersampled_waves is None else undersampled_waves,
        0.5 if lost_warning_fraction is None else lost_warning_fraction,
        0.01 if mapping_warning is None else mapping_warning, threads, cancel, progress,
    ))


def longitudinal_colour(
    system: SystemLike,
    path: int | str = 0,
    *,
    pair: ChromaticPair | tuple[int, int] | None = None,
    zone: float = 1.0,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> LongitudinalColour:
    """Longitudinal colour: focus(first) - focus(second) along the image-space propagation,
    mm, paraxially and with the real ray at pupil height ``zone`` in (0, 1]. ``pair`` =
    (first, second) wavelength indices; None: first and last system wavelength."""
    return _core.longitudinal_colour(
        compiled(system, materials), path, chromatic_pair(pair), zone, aiming, threads,
        cancel, progress,
    )


def lateral_colour(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    *,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> LateralColour:
    """Lateral colour of ``field``: chief ray per wavelength on the image surface and its
    offset from the chief ray of the reference wavelength, mm."""
    return _core.lateral_colour(
        compiled(system, materials), path, field, aiming, threads, cancel, progress
    )


def _field(field: Field | tuple[float, float]) -> Field:
    if isinstance(field, Field):
        return field
    x, y = field
    return Field(x=x, y=y)


def distortion(
    system: SystemLike,
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    samples: int = 11,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> DistortionSweep:
    """Distortion D = (h_real - h_par) / h_par in percent over ``samples`` relative fields
    0 ... 1 along +y of the largest field point; real chief-ray height on the image surface,
    paraxial height in its vertex plane (rtt/analysis/field.hpp)."""
    return _core.distortion(
        compiled(system, materials), path, wavelength, samples, aiming, threads, cancel,
        progress,
    )


def distortion_at(
    system: SystemLike,
    field: Field | tuple[float, float],
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
) -> DistortionPoint:
    """Distortion at one field value (x, y) in the units of the system's field type."""
    return _core.distortion_at(
        compiled(system, materials), path, _field(field), wavelength, aiming, threads
    )


def field_curvature(
    system: SystemLike,
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    samples: int = 11,
    delta: float = 1e-3,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> FieldCurvatureSweep:
    """Tangential and sagittal focus from the image-surface vertex along the image-space
    propagation, mm, over the field sweep; neighbour rays at +-delta (normalised pupil)."""
    return _core.field_curvature(
        compiled(system, materials), path, wavelength, samples, delta, aiming, threads,
        cancel, progress,
    )


def field_curvature_at(
    system: SystemLike,
    field: Field | tuple[float, float],
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    delta: float = 1e-3,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
) -> FieldCurvaturePoint:
    """Field curvature at one field value (x, y) in the units of the system's field type."""
    return _core.field_curvature_at(
        compiled(system, materials), path, _field(field), wavelength, delta, aiming, threads
    )


# ---------------------------------------------------------------- Jones pupil (#208, ADR 0034) ---

def jones_pupil(
    system: SystemLike,
    path: int | str = 0,
    *,
    start: RayBatch | None = None,
    chief: int | None = None,
    field: int | None = None,
    wavelength: int | None = None,
    rays: str | PupilSampling | None = None,
    aiming: Aiming | None = None,
    basis: JonesBasis = JonesBasis.DOUBLE_POLE,
    remove_global_phase: bool = False,
    field_amplitude: bool = False,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> JonesPupil:
    """Jones pupil of ``path`` (ADR 0034, points 5 and 6; #208): per pupil point the 2 x 2
    Jones matrix ``j`` (shape (n, 2, 2), rows the output, columns the input basis) and the maps
    transmission, diattenuation and physical retardance with orientation and ellipticity.

    Two forms, as path_transmission():

    - make_rays (convenience form; rotationally symmetric paths only): ``field`` (default 0)
      at ``wavelength`` (None: reference) with the pupil sampling ``rays`` (default "grid:33")
      and ``aiming`` (default REAL); the phase w is W of opd_points (reference sphere)::

          jp = rt.analysis.jones_pupil(compiled, field=1, rays="grid:17")

    - Start rays (for folded and tilted paths): ``start`` is a RayBatch of one wavelength from a
      source (P = identity), ``chief`` the index of the chief ray in it, which gives both
      frames; the phase is the OPL difference up to the image surface without a reference
      sphere (sphere.radius NaN). ``field``, ``wavelength``, ``rays`` and ``aiming`` are not
      allowed then (ValueError)::

          jp = rt.analysis.jones_pupil(compiled, "periscope", start=rays, chief=12)

    ``basis``: output basis, DOUBLE_POLE (default) or DIPOLE. ``remove_global_phase``: multiply
    every J by exp(-i psi), psi = arg det J(0, 0) / 2. ``field_amplitude``: J without the power
    factors sqrt(c) of the transmissions (the field in the image-space medium); points across
    an absorbing interface get j = NaN.

    Check ``status`` first: a point that did not arrive has w = 0 and NaN in j and the maps.
    Raises ValueError for an invalid path, chief or wavelength, mixed forms or empty start rays;
    ParaxialError (convenience form) for a path that is not rotationally symmetric;
    AnalysisError for a path with crystal modes (location /paths/i), a lost or grazing chief
    ray; Cancelled after a cancellation. Warnings (rays.lost, stop.clips_beam,
    pupil.basis_pole, pupil.field_amplitude_absorbing, pupil.chief_phase_unavailable) stay in
    ``warnings`` and are also issued as RaytatouilleWarning.
    """
    _start_form(start, field, rays, aiming, wavelength)
    if (start is None) != (chief is None):
        raise ValueError("start and chief go together: the chief ray is an index into start")
    cs = compiled(system, materials)
    if start is not None:
        assert chief is not None
        if chief < 0:
            raise ValueError(f"chief ray index {chief} is negative")
        return _warn(_core.jones_pupil_rays(
            cs, path, start, chief, basis, remove_global_phase, field_amplitude,
            lost_warning_fraction, threads, cancel, progress,
        ))
    return _warn(_core.jones_pupil(
        cs, path, 0 if field is None else field, wavelength,
        sampling("grid:33" if rays is None else rays),
        Aiming.REAL if aiming is None else aiming, basis, remove_global_phase, field_amplitude,
        lost_warning_fraction, threads, cancel, progress,
    ))



def stokes_pupil(jones: JonesPupil, e1: complex, e2: complex) -> StokesPupil:
    """Stokes map of a Jones pupil (ADR 0034, point 6, addendum #211) for the input state
    E_in = e1 x_in + e2 y_in, normalised to |E_in| = 1: per point s0 ... s3 of E_out = J E_in in
    the output basis and the ellipse (orientation, ellipticity). Scale: the power of the polarized
    input with every factor of the ray, S = weight / (||J||_F^2 / 2) Stokes(J E_in); s3 > 0 is
    right circular ((i x_in + y_in) / sqrt(2) gives s3 = s0 in free space). NaN where the
    Jones pupil has no J. The output is fully polarized; the unpolarized case is
    ``jones.transmission``.

    Raises ValueError for a zero or non-finite input state and for a Jones pupil with
    ``field_amplitude`` (the power scale needs the power-normalised J). The warnings of the
    Jones pupil are copied to ``warnings`` but not issued again.
    """
    return _core.stokes_pupil(jones, complex(e1), complex(e2))

# ----------------------------------------------- path evaluation and ghosts (#122 to #124) ---

def _start_form(
    start: RayBatch | None,
    field: int | None,
    rays: str | PupilSampling | None,
    aiming: Aiming | None,
    wavelength: int | None,
) -> None:
    """The main form (start rays) excludes the parameters of the convenience form."""
    if start is None:
        return
    given = [name for name, value in (("field", field), ("wavelength", wavelength),
                                      ("rays", rays), ("aiming", aiming)) if value is not None]
    if given:
        raise ValueError(f"{', '.join(given)} only apply without start rays (each start ray "
                         "carries its field, pupil point and wavelength)")


def path_transmission(
    system: SystemLike,
    path: int | str = 0,
    *,
    start: RayBatch | None = None,
    field: int | None = None,
    wavelength: int | None = None,
    rays: str | PupilSampling | None = None,
    aiming: Aiming | None = None,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> PathTransmission:
    """Transmission of ``path`` (#122): the final weight of every launched ray (lost rays 0)
    with mean, min and max. weight is the power for an unpolarized source (ADR 0021), so
    ``mean`` is the transmitted power fraction of a uniformly illuminated pupil (with start
    weights other than 1 the apodized transmitted power).

    Two forms:

    - Start rays (main form; also for folded and tilted paths such as interferometers, which
      the paraxial aiming of make_rays does not accept): ``start`` is a RayBatch in global
      coordinates; it is copied, not changed. Every ray counts as launched; a ray that is not
      ALIVE at the start counts as lost with its status. Start weight and start OPL are
      carried along, and each ray has its own wavelength index (``wl``). ``field``,
      ``wavelength``, ``rays`` and ``aiming`` are not allowed then (ValueError)::

          start = rt.trace.RayBatch(1)  # one ray at the origin along +z
          t = rt.analysis.path_transmission(compiled, "test arm", start=start)

    - make_rays (convenience form; rotationally symmetric paths only): ``field`` (default 0)
      at ``wavelength`` (None: reference) with the pupil sampling ``rays`` (default
      "hexapolar:6", shorthand as in spot()) and ``aiming`` (default REAL)::

          t = rt.analysis.path_transmission(compiled, "main", field=1, rays="hexapolar:12")

    Raises ValueError for an invalid path, empty start rays, a wavelength index that is not a
    system wavelength or an invalid status in the start rays, and for mixing the two forms;
    ParaxialError (convenience form) for a path that is not rotationally symmetric; NoStopError
    (convenience form) without a stop; Cancelled after a cancellation.

    ``lost_warning_fraction`` in [0, 1]: above this fraction of lost rays the result warns with
    rays.lost (ADR 0023). Warnings stay in ``warnings`` and are also issued as
    RaytatouilleWarning. ``cancel`` and ``progress`` as in the module docstring (stages "aim"
    in the convenience form, then "trace").
    """
    _start_form(start, field, rays, aiming, wavelength)
    cs = compiled(system, materials)
    if start is not None:
        return _warn(_core.path_transmission_rays(
            cs, path, start, lost_warning_fraction, threads, cancel, progress
        ))
    return _warn(_core.path_transmission(
        cs, path, 0 if field is None else field, wavelength,
        sampling("hexapolar:6" if rays is None else rays),
        Aiming.REAL if aiming is None else aiming, lost_warning_fraction, threads, cancel,
        progress,
    ))


def opl_difference(
    system: SystemLike,
    path_a: int | str,
    path_b: int | str,
    *,
    start: RayBatch | None = None,
    field: int | None = None,
    wavelength: int | None = None,
    rays: str | PupilSampling | None = None,
    aiming: Aiming | None = None,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> PathOplDifference:
    """Optical path difference OPL_b - OPL_a of two paths (#122) per start ray, in mm: a path
    length on the image surface, not a wavefront (waves: interferogram(), #210). Both
    paths must end on the same surface (ValueError otherwise); the two rays of one start ray
    may land at different points. ``points.delta`` is 0 and ``points.status`` not ALIVE for a
    ray lost on either path; ``chief`` is the delta of the first start ray at pupil (0, 0)
    that arrived on both paths, None if there is none.

    The two forms of path_transmission(): start rays (main form), e.g. for a Michelson
    interferometer::

        d = rt.analysis.opl_difference(compiled, "reference arm", "test arm", start=start)

    or ``field``, ``wavelength``, ``rays`` and ``aiming`` for make_rays on ``path_a``
    (convenience form; the same rays are traced on ``path_b``)::

        d = rt.analysis.opl_difference(compiled, 0, 1, field=0, rays="grid:9")

    Raises as path_transmission(). The warnings of both paths (path a first) are issued as
    RaytatouilleWarning and stay in ``warnings``.
    """
    _start_form(start, field, rays, aiming, wavelength)
    cs = compiled(system, materials)
    if start is not None:
        return _warn(_core.opl_difference_rays(
            cs, path_a, path_b, start, lost_warning_fraction, threads, cancel, progress
        ))
    return _warn(_core.opl_difference(
        cs, path_a, path_b, 0 if field is None else field, wavelength,
        sampling("hexapolar:6" if rays is None else rays),
        Aiming.REAL if aiming is None else aiming, lost_warning_fraction, threads, cancel,
        progress,
    ))


def interferogram(
    system: SystemLike,
    path_a: int | str,
    path_b: int | str,
    *,
    start: RayBatch | None = None,
    field: int | None = None,
    wavelength: int | None = None,
    rays: str | PupilSampling | None = None,
    aiming: Aiming | None = None,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> Interferogram:
    """Interferogram of two paths with the same image surface and image medium (#210).

    Evaluated at the hit point x_a of each start ray on ``path_a``; the wave of ``path_b`` is
    carried from its own hit x_b along its ray, OPL_b(x_a) = OPL_b(x_b) + n' d_b . (x_a - x_b):
    exact for plane waves, the first order otherwise (error about n' s^2 / (2 R_b) for the
    shear s across the ray, reported per point as ``points.shear`` and as ``max_shear``). Per
    start ray: ``delta`` = OPL_b(x_a) - OPL_a(x_a) in mm, ``waves``, the unwrapped ``phase``
    2 pi delta / lambda in rad, ``fringe`` = 1 + cos(phase) and ``intensity`` = w_a + w_b +
    2 sqrt(w_a w_b) cos(phase) with the final weights: scalar two-beam interference with equal
    polarisation; the phase is that of the optical path only (no phases of metal reflections
    or coatings). ``points.arrival`` (InterferogramArrival) says which paths a ray arrived on;
    for ONLY_A, ONLY_B and NONE the intensity is w_a, w_b or 0 and the other values are 0.

    The two forms of path_transmission(): start rays (main form), e.g. for a Michelson::

        g = rt.analysis.interferogram(compiled, "reference arm", "test arm", start=start)

    or ``field``, ``wavelength``, ``rays`` and ``aiming`` for make_rays on ``path_a``
    (convenience form; the same rays are traced on ``path_b``).

    Raises as opl_difference(), and ValueError if the image medium of either path is a crystal
    or the two image media differ. The warnings of both paths (path a first) are issued as
    RaytatouilleWarning and stay in ``warnings``.
    """
    _start_form(start, field, rays, aiming, wavelength)
    cs = compiled(system, materials)
    if start is not None:
        return _warn(_core.interferogram_rays(
            cs, path_a, path_b, start, lost_warning_fraction, threads, cancel, progress
        ))
    return _warn(_core.interferogram(
        cs, path_a, path_b, 0 if field is None else field, wavelength,
        sampling("hexapolar:6" if rays is None else rays),
        Aiming.REAL if aiming is None else aiming, lost_warning_fraction, threads, cancel,
        progress,
    ))


def ghost_ranking(
    ghosts: GhostSystem,
    field: int = 0,
    wavelength: int | None = None,
    *,
    rays: str | PupilSampling = "hexapolar:6",
    aiming: Aiming = Aiming.REAL,
    resolution_radius: float = 0.005,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> GhostRanking:
    """Ranks the ghosts of a GhostSystem (raytatouille.compile_with_ghosts) by their irradiance
    at the image relative to the useful image (#124, ADR 0027). Every ghost is traced with the
    start rays of the base path (make_rays for ``field`` at ``wavelength``, None: reference,
    with ``rays`` and ``aiming``), and the rank value is

        rho = (P_g / P_b) (r_b^2 + r0^2) / (r_g^2 + r0^2)

    with P the transmitted power (mean final weight over the launched rays), r the weighted
    RMS radius of the arrived rays about their centroid in mm and r0 = ``resolution_radius``
    in mm (> 0), the resolution radius of the detector: a model choice, not a physical
    constant. The RMS radius is geometric (no diffraction). A ghost with very few arrived rays
    gets nearly the full factor (r_b^2 + r0^2) / r0^2; check ``entries.rays_arrived``.
    ``entries.focus_offset`` and ``entries.paraxial_blur_radius`` are paraxial diagnostics in
    mm, NaN where they have no value (collimated ghost, no finite entrance pupil).

    Example::

        g = rt.compile_with_ghosts(rt.load("cooke_triplet.rtt.json"), "main", materials=lib)
        ranking = rt.analysis.ghost_ranking(g, field=0)
        strongest = ranking.entries.relative_irradiance[0]

    Raises ValueError for a GhostSystem without ghosts, an invalid field or wavelength,
    resolution_radius <= 0 or lost_warning_fraction outside [0, 1]; ParaxialError if the base
    path is not rotationally symmetric; NoStopError without a stop; AnalysisError if no ray of
    the base path arrives; Cancelled after a cancellation. Ghosts do not warn about lost rays:
    the warnings of the base path are issued as RaytatouilleWarning and stay in ``warnings``.
    """
    return _warn(_core.ghost_ranking(
        ghosts, field, wavelength, sampling(rays), aiming, resolution_radius,
        lost_warning_fraction, threads, cancel, progress,
    ))


def raytrace_report(
    system: SystemLike,
    path: int | str = 0,
    *,
    start: RayBatch,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
) -> RaytraceReport:
    """Raytrace report of ``path`` (#177): the start rays ``start`` (a RayBatch in global
    coordinates, copied, not changed) are traced with path recording, and ``rows`` has one
    entry per ray and slot (ray-major) up to the slot where each ray stopped: slot 0 is the
    start, slot s the state after event s - 1 at surface ``rows.surface`` (NO_SURFACE in slot
    0). Position (mm) and direction global and in the frame of the event's surface (``local_*``,
    NaN in slot 0), OPL from the start (mm), weight (power, unpolarized source = 1) and status
    as uint8 RayStatus values. The global columns are bitwise those of rt.trace.trace(...,
    record_path=True) for the same rays (#80). Meant for a few rays (chief and marginal rays):
    about 65 bytes per ray and slot.

    Example::

        rays = rt.trace.make_rays(compiled, rt.trace.SinglePupilPoint(0.0, 1.0), fields=[0])
        report = rt.analysis.raytrace_report(compiled, "main", start=rays)
        text = report.to_csv()

    Raises ValueError for an invalid path or empty start rays.
    """
    return _core.raytrace_report(compiled(system, materials), path, start, threads)


def system_report(
    system: SystemLike,
    path: int | str = 0,
    wavelength: int | None = None,
    *,
    materials: MaterialLibrary | None = None,
) -> SystemReport:
    """System data report of ``path`` at ``wavelength`` (None: reference) (#177): the system
    settings (wavelengths in um, reference wavelength, numbers of fields, surfaces and path
    events, stop surface or None) and the paraxial prescription (as
    raytatouille.paraxial.prescription, with the first-order data). If the prescription
    rejects the path (not rotationally symmetric, a crystal, a diffraction order != 0, a stop
    aperture that is not circular, ...), ``prescription`` is None and ``warnings`` has
    report.paraxial_unavailable with the reason; the warning is also issued as
    RaytatouilleWarning.

    Raises ValueError for an invalid path or wavelength.
    """
    return _warn(_core.system_report(compiled(system, materials), path, wavelength))


def dimension_report(
    system: SystemLike,
    *,
    materials: MaterialLibrary | None = None,
) -> DimensionReport:
    """Dimension report (#177): one entry of ``segments`` per segment (the glass between
    surfaces j and j + 1) of every lens and plate, in element order, from the compiled
    geometry (relative placement and configurations included). In mm: centre thickness (vertex
    distance along the z axis of surface j), the semi-diameters of both surfaces (the
    circumscribed radius of the aperture, its kind in ``aperture_first``/``aperture_second``
    as uint8 ApertureKind values), the diameter 2 h and the edge thickness at h, the larger of
    the two semi-diameters (where the body ends; the model has no mechanical diameter). NaN
    where a value is undefined: no aperture, a surface that does not reach h, a segment that is
    not ``coaxial``.
    """
    return _core.dimension_report(compiled(system, materials))
