"""Analyses as data objects (rtt-analysis): spot diagram, ray fans, OPD, longitudinal and
lateral colour, distortion, field curvature and Seidel sums.

Every function takes a System (compiled on each call with ``materials``) or a CompiledSystem;
for several analyses compile once with rt.compile(). ``path`` is an index or a path name,
``field`` a field index, ``wavelength`` a wavelength index where None means the reference
wavelength (in spot(): polychromatic). ``threads`` limits the worker threads; results are
bitwise the same for every number of threads.

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
    LateralColour,
    LongitudinalColour,
    MaterialLibrary,
    OpdFan,
    OpdMap,
    OpdPoints,
    Point2,
    Points2,
    RayFan,
    ReferenceSphere,
    SpotDiagram,
    SpotStatistics,
)
from ._util import SystemLike, chromatic_pair, compiled
from .paraxial import ChromaticPair, seidel
from .trace import (
    Aiming,
    FanXPupil,
    FanYPupil,
    GridPupil,
    HexapolarPupil,
    PupilSampling,
    RandomPupil,
    SinglePupilPoint,
)

__all__ = [
    "DistortionPoint",
    "DistortionSweep",
    "FanPoints",
    "FieldCurvaturePoint",
    "FieldCurvatureSweep",
    "Foci",
    "LateralColour",
    "LongitudinalColour",
    "OpdFan",
    "OpdMap",
    "OpdPoints",
    "Point2",
    "Points2",
    "RayFan",
    "ReferenceSphere",
    "SpotDiagram",
    "SpotStatistics",
    "distortion",
    "distortion_at",
    "field_curvature",
    "field_curvature_at",
    "lateral_colour",
    "longitudinal_colour",
    "opd_fan",
    "opd_map",
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


_R = TypeVar("_R", SpotDiagram, RayFan, OpdMap, OpdFan)


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
        lost_warning_fraction, threads
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
) -> RayFan:
    """Tangential and sagittal ray fans of ``field`` with ``points`` points on [-1, 1]:
    transverse aberration relative to the chief ray of the reference wavelength, mm.

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.ray_fan(
        compiled(system, materials), path, field, wavelength, points, aiming,
        lost_warning_fraction, threads
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
) -> OpdMap:
    """OPD map of ``field`` on a grid x grid pupil grid (points inside the unit circle), in
    waves at the reference wavelength, against the reference sphere centred on the chief ray.

    Raises AnalysisError e.g. for a path without stop or an exit pupil at infinity.

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.opd_map(
        compiled(system, materials), path, field, wavelength, grid, aiming,
        lost_warning_fraction, threads
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
) -> OpdFan:
    """Tangential (px = 0) and sagittal (py = 0) OPD fans of ``field``; as opd_map().

    ``lost_warning_fraction`` in [0, 1] (ValueError otherwise): above this fraction of lost
    rays the result warns with rays.lost (ADR 0023); losses has the counts in any case.
    Warnings stay in ``warnings`` and are also issued as RaytatouilleWarning.
    """
    return _warn(_core.opd_fan(
        compiled(system, materials), path, field, wavelength, points, aiming,
        lost_warning_fraction, threads
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
) -> LongitudinalColour:
    """Longitudinal colour: focus(first) - focus(second) along the image-space propagation,
    mm, paraxially and with the real ray at pupil height ``zone`` in (0, 1]. ``pair`` =
    (first, second) wavelength indices; None: first and last system wavelength."""
    return _core.longitudinal_colour(
        compiled(system, materials), path, chromatic_pair(pair), zone, aiming, threads
    )


def lateral_colour(
    system: SystemLike,
    path: int | str = 0,
    field: int = 0,
    *,
    aiming: Aiming = Aiming.REAL,
    materials: MaterialLibrary | None = None,
    threads: int | None = None,
) -> LateralColour:
    """Lateral colour of ``field``: chief ray per wavelength on the image surface and its
    offset from the chief ray of the reference wavelength, mm."""
    return _core.lateral_colour(compiled(system, materials), path, field, aiming, threads)


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
) -> DistortionSweep:
    """Distortion D = (h_real - h_par) / h_par in percent over ``samples`` relative fields
    0 ... 1 along +y of the largest field point; real chief-ray height on the image surface,
    paraxial height in its vertex plane (rtt/analysis/field.hpp)."""
    return _core.distortion(
        compiled(system, materials), path, wavelength, samples, aiming, threads
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
) -> FieldCurvatureSweep:
    """Tangential and sagittal focus from the image-surface vertex along the image-space
    propagation, mm, over the field sweep; neighbour rays at +-delta (normalised pupil)."""
    return _core.field_curvature(
        compiled(system, materials), path, wavelength, samples, delta, aiming, threads
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
