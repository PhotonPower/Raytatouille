"""aperture.reading in the JSON schema (#296, ADR 0040, point 1)."""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = json.loads((ROOT / "schema" / "raytatouille.schema.json").read_text(encoding="utf-8"))
VALIDATOR = jsonschema.Draft202012Validator(SCHEMA)
REAL = ROOT / "tests" / "reference" / "r3" / "immersion_real_na.rtt.json"


def real(mutate=None) -> dict:
    data = copy.deepcopy(json.loads(REAL.read_text(encoding="utf-8")))
    if mutate is not None:
        mutate(data)
    return data


def test_both_readings_are_accepted():
    assert VALIDATOR.is_valid(real())
    assert VALIDATOR.is_valid(real(lambda d: d["aperture"].update(reading="paraxial")))


@pytest.mark.parametrize(
    "mutate",
    [lambda d: d["aperture"].update(reading="exact"),
     lambda d: d["aperture"].update(reading=1),
     lambda d: d["aperture"].update(type="epd", value=10.0)],
    ids=["unknown", "number", "other-type"],
)
def test_wrong_readings_are_rejected(mutate):
    assert not VALIDATOR.is_valid(real(mutate))
