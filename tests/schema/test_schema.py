"""Checks that schema/raytatouille.schema.json agrees with the reference files.

The C++ parser in rtt-io is the authority. These tests make sure that the schema
accepts every canonical reference file and rejects typical structural mistakes.
"""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = json.loads((ROOT / "schema" / "raytatouille.schema.json").read_text(encoding="utf-8"))
REFERENCE_FILES = sorted((ROOT / "tests" / "reference").rglob("*.rtt.json"))
VALIDATOR = jsonschema.Draft202012Validator(SCHEMA)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_schema_is_valid():
    jsonschema.Draft202012Validator.check_schema(SCHEMA)


def test_reference_files_exist():
    assert len(REFERENCE_FILES) >= 3


@pytest.mark.parametrize("path", REFERENCE_FILES, ids=lambda p: p.name)
def test_reference_file_matches_schema(path: Path):
    VALIDATOR.validate(load(path))


def _broken(mutate):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    mutate(doc)
    return doc


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.update(colour="red"),
        lambda d: d["units"].update(length="inch"),
        lambda d: d.update(schema_version="1.0.0"),
        lambda d: d["root"]["children"][1]["surfaces"][0]["shape"]["base"].pop("radius"),
        lambda d: d["root"]["children"][1].update(type="prism"),
        lambda d: d["root"]["children"][1]["pose"].update(position=[0.0, 5.0]),
        lambda d: d["paths"][0].update(events="manual"),
        lambda d: d["aperture"].pop("value"),
    ],
    ids=["unknown-key", "units", "version", "missing-radius", "element-type",
         "position-length", "events", "aperture-value"],
)
def test_structural_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_broken(mutate))
