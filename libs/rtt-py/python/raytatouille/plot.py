"""Simple plots of analysis results with matplotlib, without analysis logic.

matplotlib is optional (``pip install raytatouille[plot]``) and imported only when a function
is called. Every function draws into ``ax`` (None: a new figure) and returns the Axes. Rays
that did not arrive (status not ALIVE) are left out.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from ._core import (
    DistortionSweep,
    FieldCurvatureSweep,
    LongitudinalColour,
    OpdFan,
    OpdMap,
    RayFan,
    RayStatus,
    SpotDiagram,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = [
    "distortion",
    "field_curvature",
    "longitudinal_colour",
    "opd_fan",
    "opd_map",
    "ray_fan",
    "spot",
]


def _axes(ax: Axes | None) -> Axes:
    if ax is not None:
        return ax
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "raytatouille.plot needs matplotlib: pip install raytatouille[plot]"
        ) from error
    _, new_ax = plt.subplots()
    return new_ax


def _arrived(status: Any) -> Any:
    return np.asarray(status) == int(RayStatus.ALIVE)


def spot(diagram: SpotDiagram, ax: Axes | None = None) -> Axes:
    """Spot points (one colour per wavelength) and the chief ray (+), local x, y in mm."""
    ax = _axes(ax)
    for wl in np.unique(diagram.wavelengths):
        sel = diagram.wavelengths == wl
        ax.scatter(diagram.x[sel], diagram.y[sel], s=4, label=f"wavelength {int(wl)}")
    ax.plot([diagram.chief.x], [diagram.chief.y], "k+", label="chief ray")
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
    ax.set_title(f"Spot, field {diagram.field}")
    ax.legend()
    return ax


def ray_fan(fan: RayFan, ax: Axes | None = None) -> Axes:
    """Tangential ey(py) and sagittal ex(px) aberration in mm over the normalised pupil."""
    ax = _axes(ax)
    t, s = fan.tangential, fan.sagittal
    at, as_ = _arrived(t.status), _arrived(s.status)
    ax.plot(t.p[at], t.ey[at], label="tangential ey(py)")
    ax.plot(s.p[as_], s.ex[as_], label="sagittal ex(px)")
    ax.set_xlabel("normalised pupil coordinate")
    ax.set_ylabel("transverse aberration / mm")
    ax.set_title(f"Ray fan, field {fan.field}, wavelength {fan.wavelength}")
    ax.legend()
    return ax


def opd_map(opd: OpdMap, ax: Axes | None = None) -> Axes:
    """OPD in waves over the normalised pupil as coloured points."""
    ax = _axes(ax)
    p = opd.points
    arrived = _arrived(p.status)
    points = ax.scatter(p.px[arrived], p.py[arrived], c=p.w[arrived], s=12)
    ax.figure.colorbar(points, ax=ax, label="W / waves")
    ax.set_aspect("equal")
    ax.set_xlabel("px")
    ax.set_ylabel("py")
    ax.set_title(f"OPD, field {opd.field}, RMS {opd.rms:.4g}, PV {opd.pv:.4g} waves")
    return ax


def opd_fan(fan: OpdFan, ax: Axes | None = None) -> Axes:
    """Tangential W(py) and sagittal W(px) in waves."""
    ax = _axes(ax)
    t, s = fan.tangential, fan.sagittal
    at, as_ = _arrived(t.status), _arrived(s.status)
    ax.plot(t.py[at], t.w[at], label="tangential W(py)")
    ax.plot(s.px[as_], s.w[as_], label="sagittal W(px)")
    ax.set_xlabel("normalised pupil coordinate")
    ax.set_ylabel("W / waves")
    ax.set_title(f"OPD fan, field {fan.field}, wavelength {fan.wavelength}")
    ax.legend()
    return ax


def distortion(sweep: DistortionSweep, ax: Axes | None = None) -> Axes:
    """Distortion in percent (x) over the relative field (y)."""
    ax = _axes(ax)
    ax.plot(sweep.percent, sweep.fraction)
    ax.set_xlabel("distortion / %")
    ax.set_ylabel("relative field")
    ax.set_title("Distortion")
    return ax


def field_curvature(sweep: FieldCurvatureSweep, ax: Axes | None = None) -> Axes:
    """Tangential and sagittal focus from the image-surface vertex (x, mm) over the relative
    field (y)."""
    ax = _axes(ax)
    ax.plot(sweep.tangential, sweep.fraction, label="tangential")
    ax.plot(sweep.sagittal, sweep.fraction, label="sagittal")
    ax.set_xlabel("focus from the image surface / mm")
    ax.set_ylabel("relative field")
    ax.set_title("Field curvature")
    ax.legend()
    return ax


def longitudinal_colour(colour: LongitudinalColour, ax: Axes | None = None) -> Axes:
    """Paraxial and real focus (global z, mm) per wavelength index."""
    ax = _axes(ax)
    foci = colour.foci
    ax.plot(foci.paraxial_z, foci.wavelength, "o-", label="paraxial")
    ax.plot(foci.real_z, foci.wavelength, "s-", label="real zone ray")
    ax.set_xlabel("focus z / mm")
    ax.set_ylabel("wavelength index")
    ax.set_title("Longitudinal colour")
    ax.legend()
    return ax
