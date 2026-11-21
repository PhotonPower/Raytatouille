"""Checks that schema/raytatouille-lenses.schema.json agrees with the lens catalogues.

The C++ parser in rtt-io (parse_lens_catalog) is the authority (ADR 0053). These tests make sure
that the schema accepts every catalogue under tests/catalogs/lenses and rejects the structural
mistakes that the parser rejects (libs/rtt-io/tests/test_lens_catalog.cpp).
"""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = json.loads(
    (ROOT / "schema" / "raytatouille-lenses.schema.json").read_text(encoding="utf-8")
)
CATALOGS = sorted((ROOT / "tests" / "catalogs" / "lenses").glob("*.json"))
VALIDATOR = jsonschema.Draft202012Validator(SCHEMA)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_schema_is_valid():
    jsonschema.Draft202012Validator.check_schema(SCHEMA)


def test_catalogs_exist():
    assert len(CATALOGS) >= 1


@pytest.mark.parametrize("path", CATALOGS, ids=lambda p: p.name)
def test_catalog_matches_schema(path: Path):
    VALIDATOR.validate(load(path))


def _broken(mutate):
    doc = copy.deepcopy(load(ROOT / "tests" / "catalogs" / "lenses" / "demo.json"))
    mutate(doc)
    return doc


def _singlet(doc):
    return doc["lenses"][0]


def _doublet(doc):
    return doc["lenses"][1]


def _triplet(doc):
    return doc["lenses"][2]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.update(colour="red"),
        lambda d: d.update(format="other"),
        lambda d: d.update(schema_version="0.2.0"),
        lambda d: d.update(catalog="A:B"),
        lambda d: d.update(catalog=""),
        lambda d: d.pop("catalog"),
        lambda d: d.update(lenses={}),
        lambda d: _singlet(d).pop("name"),
        lambda d: _singlet(d).update(name="A:B"),
        lambda d: _singlet(d).update(kind="quartet"),
        lambda d: _singlet(d).update(colour="red"),
        lambda d: _singlet(d).update(radii=[51.68]),
        lambda d: _singlet(d).update(radii=[0.0, None]),
        lambda d: _singlet(d).update(radii=[51.68, "x"]),
        lambda d: _singlet(d).update(thicknesses=[0.0]),
        lambda d: _singlet(d).update(gaps=[]),
        lambda d: _singlet(d).update(materials=[""]),
        lambda d: _singlet(d).update(diameter=0),
        lambda d: _doublet(d).update(thicknesses=[6.0]),
        lambda d: _doublet(d).update(materials=["SCHOTT:N-BK7"]),
        lambda d: _triplet(d).pop("gaps"),
        lambda d: _triplet(d).update(gaps=[5.86]),
        lambda d: _triplet(d).update(gaps=[5.86, -1.0]),
    ],
)
def test_schema_rejects_structural_errors(mutate):
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(_broken(mutate))


def test_schema_accepts_a_plano_radius_and_no_optional_keys():
    doc = _broken(lambda d: [_singlet(d).pop(k) for k in ("description", "diameter")])
    _singlet(doc)["radii"] = [None, None]
    VALIDATOR.validate(doc)
