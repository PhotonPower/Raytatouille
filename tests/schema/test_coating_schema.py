"""Checks that schema/raytatouille-coatings.schema.json agrees with the coating catalogues.

The C++ parser in rtt-coating (CoatingLibrary::add_catalog) is the authority (ADR 0019).
These tests make sure that the schema accepts every catalogue under tests/catalogs/coatings
and rejects the structural mistakes that the parser rejects.
"""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = json.loads(
    (ROOT / "schema" / "raytatouille-coatings.schema.json").read_text(encoding="utf-8")
)
CATALOGS = sorted((ROOT / "tests" / "catalogs" / "coatings").glob("*.json"))
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
    doc = copy.deepcopy(load(ROOT / "tests" / "catalogs" / "coatings" / "demo.json"))
    mutate(doc)
    return doc


def _layer(doc):
    return doc["coatings"][0]["layers"][0]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.update(colour="red"),
        lambda d: d.update(format="other"),
        lambda d: d.update(schema_version="0.2.0"),
        lambda d: d.update(catalog="A:B"),
        lambda d: d.pop("catalog"),
        lambda d: d["coatings"][0].pop("name"),
        lambda d: d["coatings"][0].update(layers=[]),
        lambda d: d["coatings"][0].pop("design_wavelength_um"),
        lambda d: d["coatings"][0].update(design_wavelength_um=0),
        lambda d: _layer(d).update(thickness_um=0.1),
        lambda d: _layer(d).pop("qwot"),
        lambda d: _layer(d).update(qwot=-1),
        lambda d: _layer(d).update(material=""),
        lambda d: _layer(d).update(n=1.5),
        lambda d: d["coatings"][0].update(name="A:B"),
        lambda d: d.update(catalog=""),
    ],
    ids=["unknown-key", "format", "version", "catalog-colon", "missing-catalog",
         "missing-name", "no-layers", "qwot-without-design-wavelength",
         "design-wavelength-zero", "both-thicknesses", "no-thickness", "negative-qwot",
         "empty-material", "unknown-layer-key", "coating-name-colon", "empty-catalog"],
)
def test_structural_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_broken(mutate))
