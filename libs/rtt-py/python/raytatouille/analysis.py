"""Analyses as data objects (rtt-analysis): spot diagram, ray fans, OPD, longitudinal and
lateral colour, distortion, field curvature and Seidel sums; path transmission and OPL
difference of two paths (#122) and the ghost ranking (#124).

Every function takes a System (compiled on each call with ``materials``) or a CompiledSystem;
for several analyses compile once with rt.compile(). ``path`` is an index or a path name,
``field`` a field index, ``wavelength`` a wavelength index where None means the reference
wavelength (in spot(): polychromatic). ``threads`` limits the worker threads; results are
bitwise the same for every number of threads.

The bundle and sweep analyses (spot, ray_fan, opd_map, opd_fan, longitudinal_colour,
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
from typing import TypeVar

from . import _core
from .errors import RaytatouilleWarning
from ._core import (
    DistortionPoint,
    DistortionSweep,
    FanPoints,
    Field,
    FieldCurvaturePoint,
    FieldCurvatureSweep,
    Foci,
    GhostEntries,
    GhostRanking,
    GhostSystem,
    LateralColour,
    LongitudinalColour,
    MaterialLibrary,
    OpdFan,
    OpdMap,
    OpdPoints,
    OplDifferencePoints,
    PathOplDifference,
    PathRays,
    PathTransmission,
    Point2,
    Points2,
    RayFan,
    RayLosses,
    ReferenceSphere,
    SpotDiagram,
    SpotStatistics,
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
    "DistortionPoint",
    "DistortionSweep",
    "FanPoints",
    "FieldCurvaturePoint",
    "FieldCurvatureSweep",
    "Foci",
    "GhostEntries",
    "GhostRanking",
    "LateralColour",
    "LongitudinalColour",
    "OpdFan",
    "OpdMap",
    "OpdPoints",
    "OplDifferencePoints",
    "PathOplDifference",
    "PathRays",
    "PathTransmission",
    "Point2",
    "Points2",
    "RayFan",
    "RayLosses",
    "ReferenceSphere",
    "SpotDiagram",
    "SpotStatistics",
    "distortion",
    "distortion_at",
    "field_curvature",
    "field_curvature_at",
    "ghost_ranking",
    "lateral_colour",
    "longitudinal_colour",
    "opd_fan",
    "opd_map",
    "opl_difference",
    "path_transmission",
    "ray_fan",
    "sampling",
    "seidel",
    "spot",
]

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


_R = TypeVar("_R", SpotDiagram, RayFan, OpdMap, OpdFan, PathTransmission, PathOplDifference,
              GhostRanking)


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
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
    lost_warning_fraction: float = 0.5,
    cancel: CancelToken | None = None,
    progress: ProgressCallback | None = None,
) -> OpdMap:
    """OPD map of ``field`` on a grid x grid pupil grid (points inside the unit circle), in
    waves at the reference wavelength, against the reference sphere centred on the chief ray.

    An exit pupil at infinity uses the limit of the reference sphere (radius inf).

    Raises AnalysisError e.g. for a path without stop.

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.opd_map(
        compiled(system, materials), path, field, wavelength, grid, aiming,
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
        compiled(system, materials), path, field, wavelength, points, aiming,
        lost_warning_fraction, threads, cancel, progress
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
    length on the image surface, not a wavefront (waves belong to an interferogram). Both
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
