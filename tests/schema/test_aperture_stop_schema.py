"""aperture.stop in the JSON schema (#53, ADR 0039, point 9): the reference form of ADR 0038."""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = json.loads((ROOT / "schema" / "raytatouille.schema.json").read_text(encoding="utf-8"))
VALIDATOR = jsonschema.Draft202012Validator(SCHEMA)
SURFACE_STOP = ROOT / "tests" / "reference" / "m2" / "cooke_triplet_surface_stop.rtt.json"


def surface_stop(mutate=None) -> dict:
    data = copy.deepcopy(json.loads(SURFACE_STOP.read_text(encoding="utf-8")))
    if mutate is not None:
        mutate(data)
    return data


def test_the_surface_stop_reference_is_accepted():
    data = surface_stop()
    assert data["aperture"]["stop"] == {"surface": "L2.S1"}
    assert VALIDATOR.is_valid(data)


@pytest.mark.parametrize(
    "mutate",
    [lambda d: d["aperture"].update(stop="L2.S1"),
     lambda d: d["aperture"].update(stop="origin"),
     lambda d: d["aperture"].update(stop=1.0),
     lambda d: d["aperture"].update(stop={}),
     lambda d: d["aperture"].update(stop={"node": "L2"}),
     lambda d: d["aperture"].update(stop={"surface": "L2.S1", "node": "L2"}),
     lambda d: d["aperture"].update(stop={"surface": 3})],
    ids=["text", "origin", "number", "empty", "node", "extra-key", "surface-number"],
)
def test_other_forms_are_rejected(mutate):
    assert not VALIDATOR.is_valid(surface_stop(mutate))
