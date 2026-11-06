"""Checks schema/raytatouille-result.schema.json against the example results (ADR 0023).

raytatouille.results.load_json is the reference. These tests make sure that the schema
accepts every example under tests/reference/results (one per result type, written by
to_json) and rejects the structural mistakes that load_json rejects.
"""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = json.loads(
    (ROOT / "schema" / "raytatouille-result.schema.json").read_text(encoding="utf-8")
)
EXAMPLES = sorted((ROOT / "tests" / "reference" / "results").glob("*.result.json"))
VALIDATOR = jsonschema.Draft202012Validator(SCHEMA)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_schema_is_valid():
    jsonschema.Draft202012Validator.check_schema(SCHEMA)


def test_one_example_per_type():
    types = sorted(load(p)["type"] for p in EXAMPLES)
    assert types == sorted(SCHEMA["properties"]["type"]["enum"])


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.name)
def test_example_matches_schema(path: Path):
    VALIDATOR.validate(load(path))


def _broken(mutate, name="OpdFan.result.json"):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "results" / name))
    mutate(doc)
    return doc


BROKEN = {
    "format": lambda d: d.update(format="raytatouille-system"),
    "version_minor": lambda d: d.update(schema_version="0.2.0"),
    "version_major": lambda d: d.update(schema_version="1.0.0"),
    "unknown_type": lambda d: d.update(type="Spot"),
    "extra_key": lambda d: d.update(extra=1),
    "data_not_object": lambda d: d.update(data=[]),
    "missing_field": lambda d: d["data"].pop("sphere"),
    "dtype": lambda d: d["data"]["tangential"]["w"].update(dtype="object"),
    "float_text": lambda d: d["data"]["tangential"]["w"]["values"].__setitem__(0, "nan"),
    "int_as_float": lambda d: d["data"]["tangential"]["status"]["values"].__setitem__(0, 1.5),
    "array_without_shape": lambda d: d["data"]["tangential"]["w"].pop("shape"),
    "special_float": lambda d: d["data"]["sphere"].update(radius={"float": "nan"}),
}


@pytest.mark.parametrize("name", sorted(BROKEN))
def test_schema_rejects(name: str):
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(_broken(BROKEN[name]))


def test_complex_values_are_pairs():
    doc = _broken(lambda d: None, "Diattenuations.result.json")
    VALIDATOR.validate(doc)
    doc["data"]["axis"]["values"][0] = 1.0
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


WAVE_KEYS = ["wave_x", "wave_y", "wave_z", "mode_index"]


def test_older_ray_batch_without_wave_is_valid():
    # RayBatch gained wave_x/y/z and mode_index in 0.1.5 (#134, ADR 0023 addendum): files of
    # an older version stay valid without them.
    for version in ("0.1.0", "0.1.4"):
        doc = _broken(lambda d: [d["data"].pop(k) for k in WAVE_KEYS], "RayBatch.result.json")
        doc["schema_version"] = version
        VALIDATOR.validate(doc)


@pytest.mark.parametrize("key", WAVE_KEYS)
@pytest.mark.parametrize("version", ["0.1.5", "0.1.12"])
def test_ray_batch_from_0_1_5_needs_wave(key: str, version: str):
    doc = _broken(lambda d: d["data"].pop(key), "RayBatch.result.json")
    doc["schema_version"] = version
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


# RayBatch gained q, RayPaths polarization, prt and q in 0.1.9 (#207, ADR 0034, ADR 0023
# addendum): files of an older version stay valid without them, from 0.1.9 on they are required.
ADDED_0_1_9 = [("RayBatch.result.json", "q"), ("RayPaths.result.json", "polarization"),
               ("RayPaths.result.json", "prt"), ("RayPaths.result.json", "q")]


@pytest.mark.parametrize("name,key", ADDED_0_1_9)
def test_older_files_without_the_0_1_9_keys_are_valid(name: str, key: str):
    for version in ("0.1.5", "0.1.8"):
        doc = _broken(lambda d: d["data"].pop(key), name)
        doc["schema_version"] = version
        VALIDATOR.validate(doc)


@pytest.mark.parametrize("name,key", ADDED_0_1_9)
@pytest.mark.parametrize("version", ["0.1.9", "0.1.12"])
def test_files_from_0_1_9_need_the_new_keys(name: str, key: str, version: str):
    doc = _broken(lambda d: d["data"].pop(key), name)
    doc["schema_version"] = version
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


# PupilCells gained obscuration and inner_rim, PupilFunction obscuration in 0.1.12 (#229 C,
# Python #209 P1b, ADR 0023 addendum): older files stay valid without them.
ADDED_0_1_12 = [("PupilCells.result.json", "obscuration"), ("PupilCells.result.json", "inner_rim"),
                ("PupilFunction.result.json", "obscuration")]


@pytest.mark.parametrize("name,key", ADDED_0_1_12)
def test_older_files_without_the_0_1_12_keys_are_valid(name: str, key: str):
    for version in ("0.1.8", "0.1.11"):
        doc = _broken(lambda d: d["data"].pop(key), name)
        doc["schema_version"] = version
        VALIDATOR.validate(doc)


@pytest.mark.parametrize("name,key", ADDED_0_1_12)
@pytest.mark.parametrize("version", ["0.1.12", "0.1.20"])
def test_files_from_0_1_12_need_the_new_keys(name: str, key: str, version: str):
    doc = _broken(lambda d: d["data"].pop(key), name)
    doc["schema_version"] = version
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


# ToleranceItems (#261): an item names exactly one anchor (surface, node or parameter), a class
# entry exactly one of surface and node; the others are null.


def _set_anchors(doc: dict, where: str, values: dict) -> None:
    entry = doc["data"]["items"][0] if where == "item" else doc["data"]["classes"]["items"][0]
    entry.update(values)


@pytest.mark.parametrize("where,values", [
    ("item", {"surface": None, "node": None, "parameter": None}),
    ("item", {"surface": "L1.S1", "node": "L1"}),
    ("item", {"node": "L1", "parameter": "R1"}),
    ("class", {"surface": None, "node": None}),
    ("class", {"surface": "L1.S1", "node": "L1"}),
])
def test_tolerance_items_need_exactly_one_anchor(where: str, values: dict):
    doc = _broken(lambda d: _set_anchors(d, where, values), "ToleranceItems.result.json")
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


@pytest.mark.parametrize("value", [[], 5, {}], ids=["list", "number", "empty_object"])
def test_the_sensitivity_expansion_is_null_or_tolerance_items(value: object):
    # expansion is anyOf [null, $defs/tolerance_items]; the definition is an object with the
    # required keys, so neither another type nor an empty object passes (review of #302).
    doc = _broken(lambda d: d["data"].update(expansion=value), "ToleranceSensitivity.result.json")
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


# SurfaceLayout gained curvature_x and conic_constant_x in 0.1.17 (#315, ADR 0042, point 5,
# ADR 0023 addendum): older files stay valid without them, from 0.1.17 on they are required.
ADDED_0_1_17 = ["curvature_x", "conic_constant_x"]


def _surface_layout_without(key: str, version: str) -> dict:
    doc = _broken(lambda d: None, "SurfaceLayout.result.json")
    doc["data"].update(curvature_x=doc["data"]["curvature"],
                       conic_constant_x=doc["data"]["conic_constant"])
    doc["data"].pop(key)
    doc["schema_version"] = version
    return doc


@pytest.mark.parametrize("key", ADDED_0_1_17)
def test_older_surface_layouts_without_the_x_section_are_valid(key: str):
    for version in ("0.1.0", "0.1.16"):
        VALIDATOR.validate(_surface_layout_without(key, version))


@pytest.mark.parametrize("key", ADDED_0_1_17)
@pytest.mark.parametrize("version", ["0.1.17", "0.1.20"])
def test_surface_layouts_from_0_1_17_need_the_x_section(key: str, version: str):
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(_surface_layout_without(key, version))


# SurfaceLayout gained normalization_radius and axicon_angle_deg in 0.1.18 (#314, ADR 0043,
# point 4, ADR 0023 addendum): older files stay valid without them, from 0.1.18 on they are
# required; both are a number or null.
ADDED_0_1_18 = ["normalization_radius", "axicon_angle_deg"]


def _surface_layout_0_1_18(version: str) -> dict:
    doc = _broken(lambda d: None, "SurfaceLayout.result.json")
    doc["data"].update(curvature_x=doc["data"]["curvature"],
                       conic_constant_x=doc["data"]["conic_constant"],
                       normalization_radius=None, axicon_angle_deg=None)
    doc["schema_version"] = version
    return doc


@pytest.mark.parametrize("key", ADDED_0_1_18)
def test_older_surface_layouts_without_the_rotational_values_are_valid(key: str):
    for version in ("0.1.0", "0.1.17"):
        doc = _surface_layout_0_1_18(version)
        doc["data"].pop(key)
        VALIDATOR.validate(doc)


@pytest.mark.parametrize("key", ADDED_0_1_18)
@pytest.mark.parametrize("version", ["0.1.18", "0.1.20"])
def test_surface_layouts_from_0_1_18_need_the_rotational_values(key: str, version: str):
    doc = _surface_layout_0_1_18(version)
    VALIDATOR.validate(doc)
    doc["data"].pop(key)
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


def test_the_rotational_values_are_numbers_or_null():
    doc = _surface_layout_0_1_18("0.1.18")
    doc["data"].update(shape="axicon", axicon_angle_deg=-5.0)
    VALIDATOR.validate(doc)
    doc["data"]["axicon_angle_deg"] = "-5"
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)
