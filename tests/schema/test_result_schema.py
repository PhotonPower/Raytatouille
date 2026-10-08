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
