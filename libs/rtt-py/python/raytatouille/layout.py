"""Geometry export for drawing a system (rtt-compile, #81): surfaces with global pose, shape,
aperture and media, vectorised sag and normals, profile polylines in a section plane and
closed outlines of the glass segments of lenses.

Units and frames: lengths in mm. "Local" is the coordinate system of a surface (its vertex at
the origin, +z along its axis); "global" are the system coordinates (right-handed, optical axis
+z). ``surface`` is an index into ``CompiledSystem.surface_ids`` or a surface id, ``element`` an
index into ``elements(system)``.

Section planes are given as ``"yz"`` (x = 0; profiles run along +y), ``"xz"`` (y = 0; along
+x) or a tuple ``(point, normal)`` in global coordinates (along z_local x normal). Profiles
exist for planes parallel to the local z axis of the surface (meridional sections, tilts in the
plane, fold mirrors, planes offset from the axis); other planes raise ValueError (general
sections come with the 3D meshes, M10). A surface without aperture is cut at the domain of its
shape (e.g. the hemisphere of a sphere): give lens surfaces apertures for drawings.
"""

from __future__ import annotations

from typing import NamedTuple, Union

import numpy as np
import numpy.typing as npt

from . import _core
from ._core import CompiledElement, CompiledSystem

__all__ = [
    "CompiledElement",
    "SectionPlane",
    "SurfaceLayout",
    "elements",
    "normal",
    "outlines",
    "profile",
    "sag",
    "surfaces",
]

FloatArray = npt.NDArray[np.float64]

#: "yz", "xz" or (point, normal) in global coordinates, mm.
SectionPlane = Union[str, tuple[npt.ArrayLike, npt.ArrayLike]]


class SurfaceLayout(NamedTuple):
    """One surface of the layout.

    ``rotation`` (3, 3) and ``translation`` (3,) map local to global coordinates:
    p_global = rotation @ p_local + translation; ``translation`` is the global vertex.
    ``shape`` is "plane", "conic" or "even_asphere" with ``curvature`` c = 1/R in 1/mm,
    ``conic_constant`` k and ``coefficients`` A4, A6, ... in mm^(1 - n); ``max_radius`` is the
    domain of the shape in mm (None: unbounded). ``aperture`` is None or (type, parameters in
    mm) with type "circular" (radius, inner_radius), "rectangular" (half_width_x, half_width_y)
    or "elliptical" (semi_axis_x, semi_axis_y), centred on the local axis.

    ``medium_front`` and ``medium_back`` index CompiledSystem.media: the media before and after
    the surface in the surface order of its element, as a refraction through the element's
    surfaces in order from the environment would see them (ADR 0017), independent of the paths.
    The media of a path are in its events; a reflection or a path running backwards may see the
    surface differently. For plates and prisms of one material with more than two surfaces this
    is only the convention of the surface order (each refraction toggles inside/environment);
    which side is glass follows from CompiledElement.media and CompiledElement.segmented."""

    id: str
    element: int
    kind: str
    rotation: FloatArray
    translation: FloatArray
    shape: str
    curvature: float
    conic_constant: float
    coefficients: list[float]
    max_radius: float | None
    aperture: tuple[str, dict[str, float]] | None
    medium_front: int
    medium_back: int

    @property
    def to_global(self) -> FloatArray:
        """Homogeneous 4 x 4 matrix local -> global."""
        m = np.eye(4)
        m[:3, :3] = self.rotation
        m[:3, 3] = self.translation
        return m


def surfaces(system: CompiledSystem) -> list[SurfaceLayout]:
    """All surfaces in tree order (index = surface index)."""
    out = []
    for entry in _core.layout_surfaces(system):
        (sid, element, kind, rotation, translation, (shape, c, k, coefficients), max_radius,
         aperture, front, back) = entry
        out.append(SurfaceLayout(sid, int(element), kind, rotation, translation, shape, float(c),
                                 float(k), list(coefficients), max_radius, aperture, int(front),
                                 int(back)))
    return out


def elements(system: CompiledSystem) -> list[CompiledElement]:
    """All elements in tree order; SurfaceLayout.element indexes this list."""
    result: list[CompiledElement] = _core.layout_elements(system)
    return result


def _surface_index(system: CompiledSystem, surface: int | str) -> int:
    ids = system.surface_ids
    if isinstance(surface, str):
        if surface not in ids:
            raise ValueError(f"unknown surface id {surface!r}")
        return ids.index(surface)
    index = int(surface)
    if not 0 <= index < len(ids):
        raise IndexError(f"surface index {index} out of range ({len(ids)} surfaces)")
    return index


def _samples(samples: int) -> int:
    if isinstance(samples, bool) or not isinstance(samples, (int, np.integer)):
        raise TypeError(f"samples must be an integer, got {type(samples).__name__}")
    if samples < 2:
        raise ValueError("a profile needs at least 2 samples")
    return int(samples)


def _plane(plane: SectionPlane) -> tuple[FloatArray, FloatArray]:
    if isinstance(plane, str):
        # Normals chosen so that the profiles run along +y ("yz") and +x ("xz"): t = z x n.
        normals = {"yz": [1.0, 0.0, 0.0], "xz": [0.0, -1.0, 0.0]}
        if plane not in normals:
            raise ValueError(f"unknown section plane {plane!r}; use 'yz', 'xz' or (point, normal)")
        return np.zeros(3), np.array(normals[plane])
    point, normal = plane
    return (np.ascontiguousarray(point, dtype=np.float64),
            np.ascontiguousarray(normal, dtype=np.float64))


def _points(x: npt.ArrayLike, y: npt.ArrayLike) -> tuple[FloatArray, FloatArray, tuple[int, ...]]:
    bx, by = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
    return np.ascontiguousarray(bx).ravel(), np.ascontiguousarray(by).ravel(), bx.shape


def sag(system: CompiledSystem, surface: int | str, x: npt.ArrayLike,
        y: npt.ArrayLike) -> FloatArray:
    """Sag z(x, y) of the surface at local points, mm; x and y broadcast, the result has their
    shape. NaN outside the domain of the shape; the aperture is not applied."""
    fx, fy, shape = _points(x, y)
    z: FloatArray = _core.layout_sag(system, _surface_index(system, surface), fx, fy)
    return z.reshape(shape)


def normal(system: CompiledSystem, surface: int | str, x: npt.ArrayLike, y: npt.ArrayLike,
           frame: str = "local") -> FloatArray:
    """Unit normal (-dz/dx, -dz/dy, 1)/|...| at local points (x, y), +z at the vertex; in
    ``frame`` "local" or "global" (rotated by the surface pose). Shape (..., 3)."""
    if frame not in ("local", "global"):
        raise ValueError("frame must be 'local' or 'global'")
    fx, fy, shape = _points(x, y)
    n: FloatArray = _core.layout_normal(system, _surface_index(system, surface), fx, fy,
                                        frame == "global")
    return n.reshape(*shape, 3)


def profile(system: CompiledSystem, surface: int | str, plane: SectionPlane = "yz",
            samples: int = 101) -> list[FloatArray]:
    """Profile polylines of the surface in the section plane, global coordinates, each (M, 3)
    with M = ``samples`` points from aperture edge to aperture edge; an annulus gives two, a
    plane that misses the surface none."""
    point, normal_ = _plane(plane)
    result: list[FloatArray] = _core.layout_profile(system, _surface_index(system, surface),
                                                    point, normal_, _samples(samples))
    return result


def outlines(system: CompiledSystem, element: int, plane: SectionPlane = "yz",
             samples: int = 101) -> list[FloatArray]:
    """Closed outlines (last point = first) of the glass segments of a lens or segmented plate
    in the section plane, global coordinates, (M, 3) each: profile of the front surface, edge,
    profile of the back surface backwards, edge. Rims ending at different places along the
    section line are joined by a step (parallel to the axis at the outermost rim as seen from the
    glass, then along the section line): the cylindrical edge at the larger radius, a bore at the
    smaller one. Only lenses and segmented plates have outlines; all plates of one material
    (windows, prisms), mirrors, thin elements, stops and detectors give an empty list; draw
    their profiles."""
    count = len(_core.layout_elements(system))
    if not 0 <= int(element) < count:
        raise IndexError(f"element index {element} out of range ({count} elements)")
    point, normal_ = _plane(plane)
    result: list[FloatArray] = _core.layout_outlines(system, int(element), point, normal_,
                                                     _samples(samples))
    return result
