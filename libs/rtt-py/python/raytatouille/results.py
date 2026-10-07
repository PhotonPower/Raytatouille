"""Result format raytatouille-result (ADR 0023): results as data and as JSON.

Every result object has ``to_dict()`` and ``to_json()``; ``load_json`` reads the JSON back.
The round trip gives the data, not the objects: ``load_json(r.to_json()).data`` equals
``r.to_dict()`` bit for bit (arrays with dtype and shape, floats including NaN and -0.0).

JSON text::

    {"format": "raytatouille-result", "schema_version": "0.1.0", "type": "SpotDiagram",
     "data": {...}}

Values in ``data``:

- numbers, strings, booleans and null as in JSON; floats in the shortest representation that
  reads back to the same double; enums as their names (e.g. "MISSED");
- a NumPy array as {"dtype": "float64", "shape": [n, ...], "values": [...]} with the values
  in C order; NaN and +-infinity in float arrays as the strings "NaN", "Infinity",
  "-Infinity"; complex values as [re, im];
- a float scalar NaN or +-infinity as {"float": "NaN"} (a plain string stays a string);
- a complex scalar as {"complex": [re, im]};
- nested results (e.g. the chief ray Point2) as objects with their attributes.

The JSON schema is schema/raytatouille-result.schema.json; this module is the reference.

Example::

    import raytatouille as rt

    spot = rt.analysis.spot(system, field=1)
    text = spot.to_json()
    data = rt.results.load_json(text).data   # dict with NumPy arrays
"""

from __future__ import annotations

import dataclasses
import enum
import json
import math
import re
from typing import Any, NamedTuple

import numpy as np

__all__ = [
    "FORMAT",
    "SCHEMA_VERSION",
    "TYPES",
    "Result",
    "encode_json",
    "load_json",
    "to_dict",
    "to_json",
]

FORMAT = "raytatouille-result"
SCHEMA_VERSION = "0.1.0"
_VERSION = re.compile(r"0\.1\.[0-9]+")

# Attributes per bound class (raytatouille._core), in output order. "name=method()" calls a
# method. Top-level types are listed in TYPES; the others only appear nested.
_RAY_BATCH = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "wl", "opl", "weight",
              "field", "pupil_x", "pupil_y", "last_surface", "status", "prt=prt_matrices()")
_FIELDS: dict[str, tuple[str, ...]] = {
    # analysis
    "Point2": ("x", "y"),
    "Points2": ("x", "y"),
    "RayLosses": ("launched", "by_status", "worst_surface", "worst_surface_count"),
    "SpotStatistics": ("centroid", "rms_centroid", "rms_chief", "geo_centroid", "geo_chief"),
    "SpotDiagram": ("field", "wavelength", "image_surface", "chief", "stats", "rays_launched",
                    "rays_arrived", "vignetted_fraction", "losses", "warnings", "x", "y",
                    "wavelengths", "weight"),
    "FanPoints": ("p", "ex", "ey", "status"),
    "RayFan": ("field", "wavelength", "image_surface", "chief", "tangential", "sagittal",
               "losses", "warnings"),
    "ReferenceSphere": ("centre", "radius"),
    "OpdPoints": ("px", "py", "w", "status"),
    "OpdMap": ("field", "wavelength", "sphere", "points", "rms", "pv", "arrived", "vignetted",
               "losses", "warnings"),
    "OpdFan": ("field", "wavelength", "sphere", "tangential", "sagittal", "losses", "warnings"),
    "Foci": ("wavelength", "paraxial_z", "real_z"),
    "LongitudinalColour": ("pair", "foci", "paraxial", "real"),
    "LateralColour": ("field", "chief", "offset"),
    "DistortionPoint": ("fraction", "field", "real_height", "paraxial_height", "percent"),
    "DistortionSweep": ("fraction", "field_x", "field_y", "real_height", "paraxial_height",
                        "percent"),
    "FieldCurvaturePoint": ("fraction", "field", "tangential", "sagittal", "astigmatism"),
    "FieldCurvatureSweep": ("fraction", "field_x", "field_y", "tangential", "sagittal",
                            "astigmatism"),
    # paraxial
    "ChromaticPair": ("first", "second"),
    "Pupil": ("z", "diameter"),
    "FirstOrder": ("object_index", "image_index", "image_direction", "power", "efl",
                   "front_focal_length", "rear_focal_length", "ffl", "bfl", "front_focal_z",
                   "rear_focal_z", "front_principal_z", "rear_principal_z", "image_z",
                   "lateral_magnification", "angular_magnification", "entrance_pupil",
                   "exit_pupil"),
    "SeidelTerms": ("s1", "s2", "s3", "s4", "s5", "c_l", "c_t"),
    "RayStart": ("z", "y", "u"),
    "SeidelSurfaces": ("surface", "y", "y_bar", "a", "a_bar", "lagrange", "s1", "s2", "s3",
                       "s4", "s5", "c_l", "c_t"),
    "Seidel": ("surfaces", "sum", "lagrange", "marginal", "chief", "chromatic"),
    # trace
    "TraceStats": ("rays",),
    "RayBatch": _RAY_BATCH,
    "RayPaths": ("slots", "n_events", "ray_count", "ray_indices", "event_surfaces", "position",
                 "direction", "opl", "weight", "status", "count", "lost_at"),
    # layout (#81)
    "CompiledElement": ("name", "first_surface", "surface_count", "media", "segmented"),
    # model
    "Field": ("x", "y", "weight"),
    "Diagnostic": ("severity", "code", "location", "message"),
}

#: Types that to_json writes and load_json accepts as "type".
TYPES: tuple[str, ...] = (
    "SpotDiagram", "RayFan", "OpdMap", "OpdFan", "LongitudinalColour", "LateralColour",
    "DistortionSweep", "DistortionPoint", "FieldCurvatureSweep", "FieldCurvaturePoint",
    "FirstOrder", "Seidel", "TraceStats", "RayBatch", "RayPaths",
    "Diattenuation", "Diattenuations", "Retardance", "Retardances",
    "GlassInfo", "GlassMap", "SurfaceLayout", "CompiledElement", "Diagnostic",
)

_DTYPES = frozenset({"bool", "int8", "int16", "int32", "int64", "uint8", "uint16", "uint32",
                     "uint64", "float32", "float64", "complex64", "complex128"})
_SPECIAL = {"NaN": math.nan, "Infinity": math.inf, "-Infinity": -math.inf}


class Result(NamedTuple):
    """A result read by load_json."""

    type: str
    schema_version: str
    data: dict[str, Any]


# ------------------------------------------------------------------------------ to_dict ---

def _value(v: Any) -> Any:
    if v is None or isinstance(v, str):
        return v
    if isinstance(v, enum.Enum):  # before int: nanobind enums may be int subclasses
        return v.name
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return float(v)
    if isinstance(v, (complex, np.complexfloating)):
        return complex(v)
    if isinstance(v, np.ndarray):
        return np.array(v, copy=True, order="C")
    if isinstance(v, tuple) and hasattr(v, "_fields"):  # NamedTuple
        return {f: _value(getattr(v, f)) for f in v._fields}
    if dataclasses.is_dataclass(v) and not isinstance(v, type):
        return {f.name: _value(getattr(v, f.name)) for f in dataclasses.fields(v)}
    fields = _FIELDS.get(type(v).__name__)
    if fields is not None:
        out: dict[str, Any] = {}
        for spec in fields:
            key, _, method = spec.partition("=")
            out[key] = _value(getattr(v, method[:-2])() if method else getattr(v, key))
        return out
    if isinstance(v, (list, tuple)):
        return [_value(x) for x in v]
    if isinstance(v, dict) and all(isinstance(k, str) for k in v):  # e.g. aperture parameters
        return {k: _value(x) for k, x in v.items()}
    raise TypeError(f"results: cannot convert {type(v).__name__}")


def to_dict(result: Any) -> dict[str, Any]:
    """The data of `result` (one of TYPES) as plain Python values; see the module docstring.

    Raises TypeError for an object that is not a result.
    """
    name = type(result).__name__
    if name not in TYPES:
        raise TypeError(f"results: {name} is not a result type ({', '.join(TYPES)})")
    data = _value(result)
    assert isinstance(data, dict)
    return data


# -------------------------------------------------------------------------------- JSON ---

def _float_item(x: float) -> float | str:
    if math.isfinite(x):
        return x
    return "NaN" if math.isnan(x) else ("Infinity" if x > 0 else "-Infinity")


def _encode(v: Any) -> Any:
    if isinstance(v, dict):
        return {k: _encode(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_encode(x) for x in v]
    if isinstance(v, np.ndarray):
        dtype = str(v.dtype)
        if dtype not in _DTYPES:
            raise ValueError(f"results: arrays of dtype {dtype} cannot be written "
                             f"({', '.join(sorted(_DTYPES))})")
        flat = v.ravel(order="C").tolist()
        if v.dtype.kind == "f":
            values: list[Any] = [_float_item(x) for x in flat]
        elif v.dtype.kind == "c":
            values = [[_float_item(z.real), _float_item(z.imag)] for z in flat]
        else:
            values = flat
        return {"dtype": dtype, "shape": list(v.shape), "values": values}
    if isinstance(v, float):
        return v if math.isfinite(v) else {"float": _float_item(v)}
    if isinstance(v, complex):
        return {"complex": [_float_item(v.real), _float_item(v.imag)]}
    return v


def encode_json(type_name: str, data: dict[str, Any], indent: int | None = None) -> str:
    """JSON text of `data` (as to_dict gives it) under the result type `type_name`.

    Raises ValueError for an unknown type or an array of a dtype outside the format.
    """
    if type_name not in TYPES:
        raise ValueError(f"results: unknown result type '{type_name}'")
    envelope = {"format": FORMAT, "schema_version": SCHEMA_VERSION, "type": type_name,
                "data": _encode(data)}
    return json.dumps(envelope, allow_nan=False, indent=indent, ensure_ascii=False)


def to_json(result: Any, indent: int | None = None) -> str:
    """`result` (one of TYPES) as JSON text in the format raytatouille-result."""
    return encode_json(type(result).__name__, to_dict(result), indent)


def _required(type_name: str) -> list[str]:
    """Keys that `data` of `type_name` must have (as written by to_dict)."""
    fields = _FIELDS.get(type_name)
    if fields is not None:
        return [spec.partition("=")[0] for spec in fields]
    from . import layout, materials, polar  # here: they import this module

    cls: Any = {"Diattenuation": polar.Diattenuation, "Diattenuations": polar.Diattenuations,
                "Retardance": polar.Retardance, "Retardances": polar.Retardances,
                "GlassMap": materials.GlassMap, "GlassInfo": materials.GlassInfo,
                "SurfaceLayout": layout.SurfaceLayout}[type_name]
    if dataclasses.is_dataclass(cls):
        return [f.name for f in dataclasses.fields(cls)]
    return list(cls._fields)


def _number(x: Any, where: str) -> float:
    if isinstance(x, str):
        if x not in _SPECIAL:
            raise ValueError(f"results: {where}: '{x}' is not a number, NaN or (-)Infinity")
        return _SPECIAL[x]
    if type(x) not in (int, float):  # bool is no number here
        raise ValueError(f"results: {where}: {x!r} is not a number")
    return float(x)


def _integer(x: Any, dtype: str, where: str) -> int:
    info = np.iinfo(dtype)
    if type(x) is not int or not info.min <= x <= info.max:
        raise ValueError(f"results: {where}: {x!r} is not a {dtype} value")
    return x


def _array(v: dict[str, Any], where: str) -> np.ndarray[Any, Any]:
    dtype, shape, values = v["dtype"], v["shape"], v["values"]
    if not isinstance(dtype, str) or dtype not in _DTYPES:
        raise ValueError(f"results: {where}: unsupported dtype {dtype!r}")
    if (not isinstance(shape, list) or not all(type(n) is int and n >= 0 for n in shape)
            or not isinstance(values, list)):
        raise ValueError(f"results: {where}: shape must be a list of sizes, values a list")
    if len(values) != math.prod(shape):
        raise ValueError(f"results: {where}: {len(values)} values for shape {shape}")
    kind = np.dtype(dtype).kind
    if kind == "f":
        items: list[Any] = [_number(x, where) for x in values]
    elif kind == "c":
        if not all(isinstance(z, list) and len(z) == 2 for z in values):
            raise ValueError(f"results: {where}: complex values must be [re, im]")
        items = [complex(_number(z[0], where), _number(z[1], where)) for z in values]
    elif kind == "b":
        if not all(type(x) is bool for x in values):
            raise ValueError(f"results: {where}: bool values must be true or false")
        items = values
    else:
        items = [_integer(x, dtype, where) for x in values]
    return np.array(items, dtype=dtype).reshape(shape)


def _decode(v: Any, where: str) -> Any:
    if isinstance(v, dict):
        keys = set(v)
        if keys & {"dtype", "float", "complex"}:  # reserved keys: only the exact encodings
            if keys == {"dtype", "shape", "values"}:
                return _array(v, where)
            if keys == {"float"}:
                if not isinstance(v["float"], str) or v["float"] not in _SPECIAL:
                    raise ValueError(f"results: {where}: float must be NaN, Infinity or "
                                     "-Infinity")
                return _SPECIAL[v["float"]]
            if keys == {"complex"}:
                z = v["complex"]
                if not (isinstance(z, list) and len(z) == 2):
                    raise ValueError(f"results: {where}: complex must be [re, im]")
                return complex(_number(z[0], where), _number(z[1], where))
            raise ValueError(f"results: {where}: the keys dtype, float and complex are reserved "
                             "for arrays {dtype, shape, values}, {float} and {complex}")
        return {k: _decode(x, f"{where}/{k}") for k, x in v.items()}
    if isinstance(v, list):
        return [_decode(x, f"{where}/{i}") for i, x in enumerate(v)]
    return v


def _no_constant(name: str) -> Any:
    raise ValueError(f"results: {name} is not standard JSON; write it as a string (ADR 0023)")


def _unique_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    keys = [k for k, _ in pairs]
    if len(set(keys)) != len(keys):
        raise ValueError("results: duplicate key in a JSON object (ADR 0020)")
    return dict(pairs)


def load_json(text: str) -> Result:
    """Reads JSON text written by to_json. The data come back as to_dict gave them; NaN is
    canonical (sign and payload are not kept).

    Raises ValueError for text that is not a raytatouille-result of schema version 0.1.x, an
    unknown type, missing keys of the type, or malformed values (a NaN literal, a duplicate
    key, an array whose values do not match its dtype or shape, a misused reserved key).
    """
    envelope = json.loads(text, parse_constant=_no_constant, object_pairs_hook=_unique_keys)
    if not isinstance(envelope, dict):
        raise ValueError("results: the JSON text must be an object")
    if set(envelope) != {"format", "schema_version", "type", "data"}:
        raise ValueError("results: expected exactly format, schema_version, type and data")
    if envelope["format"] != FORMAT:
        raise ValueError(f"results: format must be '{FORMAT}', not {envelope['format']!r}")
    version = envelope["schema_version"]
    if not isinstance(version, str) or not _VERSION.fullmatch(version):
        raise ValueError(f"results: schema_version {version!r} is not supported (0.1.x)")
    type_name = envelope["type"]
    if not isinstance(type_name, str) or type_name not in TYPES:
        raise ValueError(f"results: unknown result type {type_name!r}")
    if not isinstance(envelope["data"], dict):
        raise ValueError("results: data must be an object")
    missing = [k for k in _required(type_name) if k not in envelope["data"]]
    if missing:
        raise ValueError(f"results: data of {type_name} misses {', '.join(missing)}")
    data = _decode(envelope["data"], "/data")
    return Result(type_name, version, data)
