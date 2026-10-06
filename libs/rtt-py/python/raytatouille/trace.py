"""Ray generation and sequential tracing (rtt-trace).

``make_rays`` creates a RayBatch for a pupil sampling, ``trace`` traces it in place. The
columns of a RayBatch are NumPy views without a copy. ``trace(..., record_path=True)`` also
records the path of every ray (RayPaths) for drawing rays in a layout.
"""

from __future__ import annotations

from typing import Literal, Union, overload

import numpy as np
import numpy.typing as npt

from . import _core
from ._core import (
    DEFAULT_MAX_RECORDED_RAYS,
    NO_SURFACE,
    Aiming,
    CompiledSystem,
    FanXPupil,
    FanYPupil,
    GridPupil,
    HexapolarPupil,
    RandomPupil,
    RayBatch,
    RayPaths,
    RayStatus,
    SinglePupilPoint,
    TraceStats,
    make_rays,
)

#: Any of the pupil samplings accepted by make_rays.
PupilSampling = Union[
    SinglePupilPoint, HexapolarPupil, GridPupil, FanXPupil, FanYPupil, RandomPupil
]

__all__ = [
    "DEFAULT_MAX_RECORDED_RAYS",
    "NO_SURFACE",
    "Aiming",
    "FanXPupil",
    "FanYPupil",
    "GridPupil",
    "HexapolarPupil",
    "PupilSampling",
    "RandomPupil",
    "RayBatch",
    "RayPaths",
    "RayStatus",
    "SinglePupilPoint",
    "TraceStats",
    "make_rays",
    "trace",
]


@overload
def trace(system: CompiledSystem, rays: RayBatch, *, path: int | str = 0,
          threads: int | None = None, record_path: Literal[False] = False,
          record_rays: npt.ArrayLike | None = None,
          max_recorded_rays: int = DEFAULT_MAX_RECORDED_RAYS) -> TraceStats: ...


@overload
def trace(system: CompiledSystem, rays: RayBatch, *, path: int | str = 0,
          threads: int | None = None, record_path: Literal[True],
          record_rays: npt.ArrayLike | None = None,
          max_recorded_rays: int = DEFAULT_MAX_RECORDED_RAYS) -> tuple[TraceStats, RayPaths]: ...


@overload
def trace(system: CompiledSystem, rays: RayBatch, *, path: int | str = 0,
          threads: int | None = None, record_path: bool = False,
          record_rays: npt.ArrayLike | None = None,
          max_recorded_rays: int = DEFAULT_MAX_RECORDED_RAYS,
          ) -> TraceStats | tuple[TraceStats, RayPaths]: ...


def trace(system: CompiledSystem, rays: RayBatch, *, path: int | str = 0,
          threads: int | None = None, record_path: bool = False,
          record_rays: npt.ArrayLike | None = None,
          max_recorded_rays: int = DEFAULT_MAX_RECORDED_RAYS,
          ) -> TraceStats | tuple[TraceStats, RayPaths]:
    """Traces ``rays`` in place along ``path`` (index or name) with the sequential tracer and
    returns the counts per status. ``threads`` limits the worker threads (None: all); the result
    is bitwise the same for every number of threads. The GIL is released; do not read or change
    the columns of ``rays`` from another thread meanwhile.

    With ``record_path=True`` the path of the rays is recorded as well and the result is
    ``(TraceStats, RayPaths)``; the rays end bitwise as without recording. ``record_rays``
    selects the rays to record (1-D integer indices into ``rays``, in this order, no
    duplicates); without it all rays are recorded, at most ``max_recorded_rays`` (default
    10000, about 65 bytes per ray and slot), otherwise ValueError. Select a subset explicitly
    (e.g. every k-th ray) rather than the first rays: those are only the inner rings of a
    hexapolar bundle.

    Raises ValueError for an unknown path name, a wavelength index that is not a system
    wavelength, an invalid status or an invalid selection, IndexError for an unknown path
    index."""
    if not record_path:
        if record_rays is not None:
            raise ValueError("record_rays needs record_path=True")
        stats: TraceStats = _core.trace(system, rays, path=path, threads=threads)
        return stats
    selection = None
    if record_rays is not None:
        indices = np.asarray(record_rays)
        if indices.ndim == 1 and indices.size == 0:
            raise ValueError("record_rays is empty")
        if indices.ndim != 1 or not np.issubdtype(indices.dtype, np.integer):
            raise ValueError("record_rays must be a 1-D array of integer ray indices")
        if indices.dtype.kind == "u" and int(indices.max()) > np.iinfo(np.int64).max:
            raise ValueError("record_rays has an index of 2**63 or more, larger than any batch")
        selection = np.ascontiguousarray(indices, dtype=np.int64)
    result: tuple[TraceStats, RayPaths] = _core.trace_recorded(
        system, rays, path=path, threads=threads, record_rays=selection,
        max_recorded_rays=max_recorded_rays)
    return result
