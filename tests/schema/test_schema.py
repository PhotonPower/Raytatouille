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


def _lens(doc):
    return doc["root"]["children"][1]


def _with_order(order):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    doc["paths"] = [{"name": "main", "events": [{"surface": "L1.S1", "kind": "transmit",
                                                  "order": order}]}]
    return doc


def test_event_order_is_an_int():
    # Event::order is an int; the C++ reader rejects values outside it (#35 E).
    for order in (2147483647, -2147483648, 0, 1):
        assert VALIDATOR.is_valid(_with_order(order)), order
    for order in (2147483648, -2147483649, 4294967296, 1.5):
        assert not VALIDATOR.is_valid(_with_order(order)), order


def test_weights_are_not_negative():
    # model::validate rejects negative weights (wavelengths.weight_invalid,
    # fields.weight_invalid); the schema says the same (#35). Zero stays allowed.
    def with_weights(wavelength, field):
        doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
        doc["wavelengths"][0]["weight"] = wavelength
        doc["fields"]["points"][0]["weight"] = field
        return doc

    for wavelength, field in ((0.0, 0.0), (1.0, 2.5)):
        assert VALIDATOR.is_valid(with_weights(wavelength, field)), (wavelength, field)
    for wavelength, field in ((-1.0, 1.0), (1.0, -0.5)):
        assert not VALIDATOR.is_valid(with_weights(wavelength, field)), (wavelength, field)


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _lens(d).update(material="SCHOTT:N-BK7"),
        lambda d: _lens(d).update(material=["SCHOTT:N-BK7"]),
        lambda d: _lens(d).update(material=["SCHOTT:N-BK7", "SCHOTT:F2"]),
    ],
    ids=["material-shorthand", "material-list-one", "material-list-two"],
)
def test_segment_materials_are_accepted(mutate):
    VALIDATOR.validate(_broken(mutate))


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _lens(d).update(material=[]),
        lambda d: _lens(d).update(material=["SCHOTT:N-BK7", 2]),
        lambda d: _lens(d).update(material=[["SCHOTT:N-BK7"]]),
        lambda d: _lens(d).update(material={"name": "SCHOTT:N-BK7"}),
        lambda d: _lens(d).update(material=True),
        lambda d: d.update(schema_version="0.1.0"),
    ],
    ids=["material-list-empty", "material-list-number", "material-list-nested",
         "material-object", "material-bool", "old-version"],
)
def test_segment_material_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_broken(mutate))


def test_achromat_reference_uses_a_material_list():
    doc = load(ROOT / "tests" / "reference" / "m2" / "achromat.rtt.json")
    assert _lens(doc)["material"] == ["SCHOTT:N-BK7", "SCHOTT:F2"]
    VALIDATOR.validate(doc)


CALCITE = {"ordinary": "BIREFRINGENT:CALCITE", "extraordinary": "BIREFRINGENT:CALCITE-E"}


def _first_surface(doc):
    return _lens(doc)["surfaces"][0]


def _event(doc, **event):
    doc["paths"] = [{"name": "main", "events": [{"surface": "L1.S1", **event}]}]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _lens(d).update(material=CALCITE, optic_axis=[0.0, 1.0, 1.0]),
        lambda d: _first_surface(d).update(
            phases=[{"type": "linear_grating", "lines_per_mm": 300.0}],
            diffraction_efficiency=[{"order": 1, "efficiency": 0.8},
                                    {"order": -1, "efficiency": 0.1}]),
        # Values and duplicates are semantics, checked by rtt validate like in the parser.
        lambda d: _first_surface(d).update(diffraction_efficiency=[]),
        lambda d: _first_surface(d).update(
            diffraction_efficiency=[{"order": 1, "efficiency": 1.5}]),
        lambda d: _event(d, kind="reflect", order=1),
        lambda d: _event(d, kind="ordinary", order=-2),
    ],
    ids=["crystal", "efficiency", "efficiency-empty", "efficiency-above-1", "reflect-order",
         "ordinary-order"],
)
def test_schema_0_3_keys_are_accepted(mutate):
    # ADR 0025, ADR 0026.
    VALIDATOR.validate(_broken(mutate))


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _event(d, kind="diffract", order=1),
        lambda d: _lens(d).update(material={"ordinary": "BIREFRINGENT:CALCITE"}),
        lambda d: _lens(d).update(material={**CALCITE, "axis": [0.0, 0.0, 1.0]}),
        lambda d: _lens(d).update(material={**CALCITE, "extraordinary": 1}),
        lambda d: _lens(d).update(optic_axis=[0.0, 1.0]),
        lambda d: _lens(d).update(optic_axis="z"),
        lambda d: _first_surface(d).update(diffraction_efficiency=[{"order": 1}]),
        lambda d: _first_surface(d).update(
            diffraction_efficiency=[{"order": 1.5, "efficiency": 0.5}]),
        lambda d: _first_surface(d).update(
            diffraction_efficiency=[{"order": 1, "efficiency": 0.5, "phase": 0.0}]),
        lambda d: _first_surface(d).update(diffraction_efficiency={"order": 1, "efficiency": 0.5}),
        lambda d: d.update(schema_version="0.2.0"),
    ],
    ids=["diffract", "crystal-missing-extraordinary", "crystal-extra-key",
         "crystal-not-a-string", "axis-length", "axis-type", "efficiency-missing",
         "efficiency-order-float", "efficiency-extra-key", "efficiency-not-a-list",
         "version-0.2"],
)
def test_schema_0_3_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_broken(mutate))
