"""Reports as CSV (#177). The reports themselves come from raytatouille.analysis
(raytrace_report, system_report, dimension_report); as data they also have to_dict() and
to_json() (result format, ADR 0023).

``to_csv(report)`` (also ``report.to_csv()``) gives CSV text: comma-separated, a header line,
lines ending in "\\n", fields with a comma, quote or line break quoted as in RFC 4180. Values:

- floats in the shortest representation that reads back to the same double (``repr``), so
  ``float(field)`` gives the value bit for bit; +-infinity as "Infinity" and "-Infinity";
- NaN and None as an empty field;
- status columns as RayStatus names (e.g. ALIVE), aperture kinds as ApertureKind names (e.g.
  CIRCULAR), booleans as "true" and "false";
- a surface index NO_SURFACE (slot 0 of the raytrace report) as an empty field.

Layout per report:

- RaytraceReport: one line per ray and slot (ray-major), the columns of ``rows`` in their
  order: ray, slot, surface, x, y, z, dx, dy, dz, local_x, local_y, local_z, local_dx,
  local_dy, local_dz, opl, weight, status.
- DimensionReport: one line per segment, the columns of ``segments`` in their order: element,
  first_surface, coaxial, centre_thickness, semi_diameter_first, semi_diameter_second,
  aperture_first, aperture_second, edge_thickness, diameter.
- SystemReport: two columns ``key,value``, one line per value of ``to_dict()`` in its order.
  The key is the path to the value: object keys joined with ".", list and array indices as
  "[i]" (e.g. ``wavelengths_um[0]``, ``prescription.first_order.efl``,
  ``prescription.surfaces.z[3]``, ``warnings[0].code``). A value that is None (e.g. ``stop``,
  ``prescription``, ``prescription.first_order.entrance_pupil``) is one line with an empty
  value; an empty list has no line. Keys are stable within the result format 0.1.x: new
  values only add lines.

Example::

    import raytatouille as rt

    report = rt.analysis.dimension_report(compiled)
    with open("dimensions.csv", "w", encoding="utf-8", newline="") as f:
        f.write(report.to_csv())
"""

from __future__ import annotations

import csv
import io
import math
from typing import Any, Union

import numpy as np

from ._core import (
    NO_SURFACE,
    ApertureKind,
    DimensionReport,
    RayStatus,
    RaytraceReport,
    SystemReport,
)
from .results import to_dict

__all__ = ["RAYTRACE_COLUMNS", "DIMENSION_COLUMNS", "to_csv"]

Report = Union[RaytraceReport, SystemReport, DimensionReport]

#: Columns of the raytrace report CSV, in order (the attributes of RaytraceRows).
RAYTRACE_COLUMNS: tuple[str, ...] = (
    "ray", "slot", "surface", "x", "y", "z", "dx", "dy", "dz", "local_x", "local_y", "local_z",
    "local_dx", "local_dy", "local_dz", "opl", "weight", "status")
#: Columns of the dimension report CSV, in order (the attributes of DimensionSegments).
DIMENSION_COLUMNS: tuple[str, ...] = (
    "element", "first_surface", "coaxial", "centre_thickness", "semi_diameter_first",
    "semi_diameter_second", "aperture_first", "aperture_second", "edge_thickness", "diameter")


def _field(v: Any) -> str:
    """One CSV field of a plain value (as to_dict gives it)."""
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        x = float(v)
        if math.isnan(x):
            return ""
        if math.isinf(x):
            return "Infinity" if x > 0 else "-Infinity"
        return repr(x)
    if isinstance(v, str):
        return v
    raise TypeError(f"reports: cannot write {type(v).__name__} as a CSV field")


def _table(columns: Any, names: tuple[str, ...]) -> list[list[str]]:
    """Header and one line per entry of a column class."""
    formats = {"status": lambda v: RayStatus(int(v)).name,
               "aperture_first": lambda v: ApertureKind(int(v)).name,
               "aperture_second": lambda v: ApertureKind(int(v)).name,
               "surface": lambda v: "" if int(v) == NO_SURFACE else str(int(v))}
    arrays = [getattr(columns, name) for name in names]
    lines = [list(names)]
    for i in range(len(columns)):
        lines.append([formats.get(name, _field)(a[i]) for name, a in zip(names, arrays)])
    return lines


def _flatten(key: str, v: Any, out: list[list[str]]) -> None:
    """Lines key,value of the system report (see the module docstring)."""
    if isinstance(v, dict):
        for k, x in v.items():
            _flatten(f"{key}.{k}" if key else k, x, out)
    elif isinstance(v, np.ndarray):
        for index in np.ndindex(v.shape):
            _flatten(key + "".join(f"[{i}]" for i in index), v[index], out)
    elif isinstance(v, list):
        for i, x in enumerate(v):
            _flatten(f"{key}[{i}]", x, out)
    else:
        out.append([key, _field(v)])


def to_csv(report: Report) -> str:
    """``report`` (RaytraceReport, SystemReport or DimensionReport) as CSV text; layout and
    values as in the module docstring.

    Raises TypeError for any other object.
    """
    if isinstance(report, RaytraceReport):
        lines = _table(report.rows, RAYTRACE_COLUMNS)
    elif isinstance(report, DimensionReport):
        lines = _table(report.segments, DIMENSION_COLUMNS)
    elif isinstance(report, SystemReport):
        lines = [["key", "value"]]
        _flatten("", to_dict(report), lines)
    else:
        raise TypeError(f"reports: {type(report).__name__} is not a report (RaytraceReport, "
                        "SystemReport, DimensionReport)")
    text = io.StringIO()
    csv.writer(text, lineterminator="\n").writerows(lines)
    return text.getvalue()
