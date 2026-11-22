"""Design helpers and lens catalogues (ADR 0053): starting systems from a few numbers.

``singlet``, ``doublet`` (cemented) and ``triplet`` (three singlets with air gaps) return a
complete, validated System: a stop (its own node in front of the first lens, or a lens surface),
the lenses, a detector, the entrance pupil diameter, fields and wavelengths. Lengths in mm,
wavelengths in um, field angles in degree; radii with the sign convention of the system format
(centre of curvature in +z: positive), ``None`` is plano. Materials are references
``CATALOG:NAME``, resolved when the system is compiled::

    import raytatouille as rt

    system = rt.design.doublet(
        (61.5, -35.2, -136.9), (6.0, 3.0), ("SCHOTT:N-BK7", "SCHOTT:F2"),
        epd=20.0, stop_distance=5.0, semi_diameter=12.7, image_distance=96.614)
    catalog = rt.design.load_lens_catalog("lenses.json")
    lens = rt.design.from_catalog(catalog["ACH-25-100"], epd=20.0, stop_distance=5.0,
                                  image_distance=96.614)

Names in the built system: lenses ``L1``, ``L2``, ``L3`` with surfaces ``L1.S1``, ...; the stop
node ``stop`` with surface ``STO``; the detector ``image`` with surface ``IMG``; one path
``main`` with automatic events. Layout (ADR 0053, point 1): the first lens node at z =
``stop_distance`` (or 0 with ``stop_surface``), the surfaces relative to their node at the
cumulative thicknesses, the next node at z + (last surface) + gap, the detector at z_last +
(last surface + ``image_distance``). Aspheres, tilts, coatings and air-spaced doublets are not
part of the helpers; change the result with ``rt.apply_patch``.
"""

from __future__ import annotations

import json
import math
import warnings
from collections.abc import Sequence
from typing import Any

from ._core import (
    SCHEMA_VERSION,
    LensCatalog,
    LensEntry,
    Severity,
    System,
    load_lens_catalog,
    parse_lens_catalog,
    validate,
)
from .errors import DesignError, RaytatouilleWarning

#: Default wavelengths in um (F, d, C lines); the reference is index 1 (d).
DEFAULT_WAVELENGTHS: tuple[float, ...] = (0.4861, 0.5876, 0.6563)

Radius = float | None


def _finite(value: float, what: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite, got {value!r}")
    return number


def _positive(value: float, what: str) -> float:
    number = _finite(value, what)
    if number <= 0.0:
        raise ValueError(f"{what} must be > 0, got {value!r}")
    return number


def _count(values: Sequence[Any], count: int, what: str) -> list[Any]:
    items = list(values)
    if len(items) != count:
        raise ValueError(f"{what} needs {count} values, got {len(items)}")
    return items


def _radius(value: Radius, what: str) -> float | None:
    if value is None or (isinstance(value, float) and math.isinf(value)):
        return None
    number = _finite(value, what)
    if number == 0.0:
        raise ValueError(f"{what} must not be 0 (None is plano)")
    return number


def _material(value: str, what: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{what} must be a non-empty reference CATALOG:NAME, got {value!r}")
    return value


class _Lens:
    """One lens node: radii, thicknesses and materials of its glass segments."""

    def __init__(self, radii: list[float | None], thicknesses: list[float],
                 materials: list[str]) -> None:
        self.radii = radii
        self.thicknesses = thicknesses
        self.materials = materials


def _build(lenses: list[_Lens], gaps: list[float], *, epd: float, image_distance: float,
           stop_distance: float | None, stop_surface: str | None, stop_radius: float | None,
           semi_diameter: float | None, fields_y: Sequence[float], wavelengths: Sequence[float],
           reference: int, name: str) -> System:
    epd = _positive(epd, "epd")
    image_distance = _finite(image_distance, "image_distance")
    if (stop_distance is None) == (stop_surface is None):
        raise ValueError("give exactly one of stop_distance and stop_surface")
    if semi_diameter is not None:
        semi_diameter = _positive(semi_diameter, "semi_diameter")
    waves = [_positive(w, "wavelength") for w in wavelengths]
    if not waves:
        raise ValueError("wavelengths must not be empty")
    if isinstance(reference, bool) or not isinstance(reference, int) or not (
            0 <= reference < len(waves)):
        raise ValueError(f"reference must be an index 0..{len(waves) - 1}, got {reference!r}")
    fields = [_finite(f, "field angle") for f in fields_y]
    if not fields:
        raise ValueError("fields_y must not be empty")

    ids = [f"L{i + 1}.S{j + 1}" for i, lens in enumerate(lenses) for j in range(len(lens.radii))]
    if stop_surface is not None:
        if stop_surface not in ids:
            raise ValueError(f"stop_surface {stop_surface!r} is none of {', '.join(ids)}")
        if stop_radius is None:
            raise ValueError("stop_surface needs stop_radius (no default for a surface stop)")
        stop_radius = _positive(stop_radius, "stop_radius")
        z = 0.0
    else:
        assert stop_distance is not None
        # A ray from the stop plane does not hit a surface at distance 0 (#151, #177).
        z = _positive(stop_distance, "stop_distance")
        stop_radius = epd / 2.0 if stop_radius is None else _positive(stop_radius, "stop_radius")

    children: list[dict[str, Any]] = []
    if stop_surface is None:
        children.append({"type": "stop", "name": "stop", "surfaces": [
            {"id": "STO", "aperture": {"type": "circular", "radius": stop_radius}}]})
    last = 0.0
    for i, lens in enumerate(lenses):
        if i > 0:
            z = z + last + gaps[i - 1]  # left to right (ADR 0053, point 1)
        surfaces: list[dict[str, Any]] = []
        position = 0.0
        for j, radius in enumerate(lens.radii):
            if j > 0:
                position = position + lens.thicknesses[j - 1]
            surface: dict[str, Any] = {"id": f"L{i + 1}.S{j + 1}"}
            if j > 0:
                surface["pose"] = {"position": [0.0, 0.0, position]}
            if radius is not None:
                surface["shape"] = {"base": {"type": "conic", "radius": radius}}
            if surface["id"] == stop_surface:
                surface["aperture"] = {"type": "circular", "radius": stop_radius}
            elif semi_diameter is not None:
                surface["aperture"] = {"type": "circular", "radius": semi_diameter}
            surfaces.append(surface)
        last = position
        node: dict[str, Any] = {"type": "lens", "name": f"L{i + 1}"}
        if z != 0.0:
            node["pose"] = {"position": [0.0, 0.0, z]}
        node["material"] = lens.materials[0] if len(lens.materials) == 1 else lens.materials
        node["surfaces"] = surfaces
        children.append(node)
    image_z = z + (last + image_distance)  # relative to the last node (ADR 0053, point 1)
    children.append({"type": "detector", "name": "image",
                     "pose": {"position": [0.0, 0.0, image_z]}, "surfaces": [{"id": "IMG"}]})

    aperture: dict[str, Any] = {"type": "epd", "value": epd}
    if stop_surface is not None:
        aperture["stop"] = {"surface": stop_surface}
    data: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "name": name,
        "units": {"length": "mm", "wavelength": "um"},
        "wavelengths": [{"um": w, "reference": True} if k == reference else {"um": w}
                        for k, w in enumerate(waves)],
        "aperture": aperture,
        "fields": {"points": [{"y": f} for f in fields]},
        "root": {"type": "assembly", "name": "system", "children": children},
        "paths": [{"name": "main", "events": "auto"}],
    }
    system = System.from_json(json.dumps(data))
    diagnostics = validate(system)
    errors = [d for d in diagnostics if d.severity == Severity.ERROR]
    if errors:
        raise DesignError("the built system is invalid: " + "; ".join(str(d) for d in errors),
                          errors)
    for d in diagnostics:
        warnings.warn(RaytatouilleWarning(d.message, d.code, d.location), stacklevel=3)
    return system


def singlet(radii: Sequence[Radius], thickness: float, material: str, *, epd: float,
            image_distance: float, stop_distance: float | None = None,
            stop_surface: str | None = None, stop_radius: float | None = None,
            semi_diameter: float | None = None, fields_y: Sequence[float] = (0.0,),
            wavelengths: Sequence[float] = DEFAULT_WAVELENGTHS, reference: int = 1,
            name: str = "Singlet") -> System:
    """A single lens L1 with a stop and a detector (ADR 0053, point 1).

    ``radii`` = (R1, R2) in mm (None: plano); ``thickness`` > 0 in mm; ``material`` a
    reference CATALOG:NAME. ``epd``: entrance pupil diameter in mm. ``image_distance``: mm from
    the vertex of L1.S2 to the image plane. The stop: exactly one of ``stop_distance`` (mm > 0, a
    stop node at z = 0 in front of the lens) and ``stop_surface`` (a surface id such as
    "L1.S1"); ``stop_radius`` in mm, default epd/2 with ``stop_distance``, required with
    ``stop_surface``. ``semi_diameter``: aperture radius in mm of every other lens surface.
    ``fields_y``: field angles in degree; ``wavelengths`` in um with the 0-based index
    ``reference``; ``name`` of the system. Raises ValueError for wrong arguments and DesignError
    if the built system is invalid; warnings of validate arrive as RaytatouilleWarning.
    """
    r = [_radius(v, f"radii[{k}]") for k, v in enumerate(_count(radii, 2, "radii"))]
    lens = _Lens(r, [_positive(thickness, "thickness")], [_material(material, "material")])
    return _build([lens], [], epd=epd, image_distance=image_distance, stop_distance=stop_distance,
                  stop_surface=stop_surface, stop_radius=stop_radius, semi_diameter=semi_diameter,
                  fields_y=fields_y, wavelengths=wavelengths, reference=reference, name=name)


def doublet(radii: Sequence[Radius], thicknesses: Sequence[float], materials: Sequence[str], *,
            epd: float, image_distance: float, stop_distance: float | None = None,
            stop_surface: str | None = None, stop_radius: float | None = None,
            semi_diameter: float | None = None, fields_y: Sequence[float] = (0.0,),
            wavelengths: Sequence[float] = DEFAULT_WAVELENGTHS, reference: int = 1,
            name: str = "Doublet") -> System:
    """A cemented doublet L1 (three surfaces, a material per segment, ADR 0017) with a stop and
    a detector (ADR 0053, point 1).

    ``radii`` = (R1, R2, R3) in mm (None: plano); ``thicknesses`` = (t1, t2) > 0 in mm;
    ``materials`` = (glass 1, glass 2) as CATALOG:NAME. ``image_distance``: mm from the vertex
    of L1.S3 to the image plane. The other arguments as in ``singlet``.
    """
    r = [_radius(v, f"radii[{k}]") for k, v in enumerate(_count(radii, 3, "radii"))]
    t = [_positive(v, f"thicknesses[{k}]")
         for k, v in enumerate(_count(thicknesses, 2, "thicknesses"))]
    m = [_material(v, f"materials[{k}]") for k, v in enumerate(_count(materials, 2, "materials"))]
    return _build([_Lens(r, t, m)], [], epd=epd, image_distance=image_distance,
                  stop_distance=stop_distance, stop_surface=stop_surface, stop_radius=stop_radius,
                  semi_diameter=semi_diameter, fields_y=fields_y, wavelengths=wavelengths,
                  reference=reference, name=name)


def triplet(radii: Sequence[Radius], thicknesses: Sequence[float], gaps: Sequence[float],
            materials: Sequence[str], *, epd: float, image_distance: float,
            stop_distance: float | None = None, stop_surface: str | None = None,
            stop_radius: float | None = None, semi_diameter: float | None = None,
            fields_y: Sequence[float] = (0.0,), wavelengths: Sequence[float] = DEFAULT_WAVELENGTHS,
            reference: int = 1, name: str = "Triplet") -> System:
    """Three singlets L1, L2, L3 with air gaps, a stop and a detector (ADR 0053, point 1).

    ``radii``: 6 values in mm (None: plano), two per lens; ``thicknesses``: 3 centre thicknesses
    > 0 in mm; ``gaps``: 2 air gaps >= 0 in mm from the rear vertex of a lens to the front
    vertex of the next; ``materials``: 3 references CATALOG:NAME. ``image_distance``: mm from
    the vertex of L3.S2 to the image plane. The other arguments as in ``singlet``; a stop between
    the lenses is a surface stop, e.g. ``stop_surface="L2.S1"``.
    """
    r = [_radius(v, f"radii[{k}]") for k, v in enumerate(_count(radii, 6, "radii"))]
    t = [_positive(v, f"thicknesses[{k}]")
         for k, v in enumerate(_count(thicknesses, 3, "thicknesses"))]
    g = []
    for k, v in enumerate(_count(gaps, 2, "gaps")):
        gap = _finite(v, f"gaps[{k}]")
        if gap < 0.0:
            raise ValueError(f"gaps[{k}] must be >= 0, got {v!r}")
        g.append(gap)
    m = [_material(v, f"materials[{k}]") for k, v in enumerate(_count(materials, 3, "materials"))]
    lenses = [_Lens(r[2 * k:2 * k + 2], [t[k]], [m[k]]) for k in range(3)]
    return _build(lenses, g, epd=epd, image_distance=image_distance, stop_distance=stop_distance,
                  stop_surface=stop_surface, stop_radius=stop_radius, semi_diameter=semi_diameter,
                  fields_y=fields_y, wavelengths=wavelengths, reference=reference, name=name)


def from_catalog(entry: LensEntry, *, epd: float, image_distance: float,
                 stop_distance: float | None = None, stop_surface: str | None = None,
                 stop_radius: float | None = None, semi_diameter: float | None = None,
                 fields_y: Sequence[float] = (0.0,),
                 wavelengths: Sequence[float] = DEFAULT_WAVELENGTHS, reference: int = 1,
                 name: str | None = None) -> System:
    """The system of a catalogue lens (ADR 0053, point 3): the helper of ``entry.kind`` with
    the radii, thicknesses, gaps and materials of the entry. ``semi_diameter`` defaults to
    ``entry.diameter / 2``, ``name`` to ``entry.name``. The stop is required as with the
    helpers. The other arguments as in ``singlet``.
    """
    if semi_diameter is None and entry.diameter is not None:
        semi_diameter = entry.diameter / 2.0
    common: dict[str, Any] = {
        "epd": epd, "image_distance": image_distance, "stop_distance": stop_distance,
        "stop_surface": stop_surface, "stop_radius": stop_radius, "semi_diameter": semi_diameter,
        "fields_y": fields_y, "wavelengths": wavelengths, "reference": reference,
        "name": entry.name if name is None else name,
    }
    if entry.kind == "singlet":
        return singlet(entry.radii, entry.thicknesses[0], entry.materials[0], **common)
    if entry.kind == "doublet":
        return doublet(entry.radii, entry.thicknesses, entry.materials, **common)
    return triplet(entry.radii, entry.thicknesses, entry.gaps, entry.materials, **common)


__all__ = [
    "DEFAULT_WAVELENGTHS",
    "LensCatalog",
    "LensEntry",
    "doublet",
    "from_catalog",
    "load_lens_catalog",
    "parse_lens_catalog",
    "singlet",
    "triplet",
]
