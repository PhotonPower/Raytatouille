"""Ray generation and sequential tracing (rtt-trace).

``make_rays`` creates a RayBatch for a pupil sampling, ``trace`` traces it in place. The
columns of a RayBatch are NumPy views without a copy.
"""

from typing import Union

from ._core import (
    NO_SURFACE,
    Aiming,
    FanXPupil,
    FanYPupil,
    GridPupil,
    HexapolarPupil,
    RandomPupil,
    RayBatch,
    RayStatus,
    SinglePupilPoint,
    TraceStats,
    make_rays,
    trace,
)

#: Any of the pupil samplings accepted by make_rays.
PupilSampling = Union[
    SinglePupilPoint, HexapolarPupil, GridPupil, FanXPupil, FanYPupil, RandomPupil
]

__all__ = [
    "NO_SURFACE",
    "Aiming",
    "FanXPupil",
    "FanYPupil",
    "GridPupil",
    "HexapolarPupil",
    "PupilSampling",
    "RandomPupil",
    "RayBatch",
    "RayStatus",
    "SinglePupilPoint",
    "TraceStats",
    "make_rays",
    "trace",
]
