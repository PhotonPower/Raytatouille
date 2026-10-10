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


# ---- schema 0.4 (#162): relative poses (ADR 0028), parameter table (ADR 0029) ----

def _tour():
    return copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "feature_tour.rtt.json"))


def _with(mutate):
    doc = _tour()
    mutate(doc)
    return doc


def _a_s2_z(d):
    return d["root"]["children"][2]["surfaces"][1]["pose"]["position"]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _a_s2_z(d).__setitem__(2, {"param": "A_S2_Z"}),
        lambda d: _a_s2_z(d).__setitem__(2, {"value": 6.0, "min": 0.0, "max": 10.0}),
        lambda d: d["root"]["children"][2].update(
            pose={"reference": "relative_to_preceding", "order": "rotate_first"}),
        lambda d: d["root"]["children"][2].update(pose={"reference": "relative_to_sibling"}),
        lambda d: d["parameters"].append({"name": "_x9", "values": [1.0, 2.0]}),
    ],
    ids=["bound-param", "bounds", "relative-rotate-first", "sibling", "values-row"],
)
def test_0_4_forms_are_accepted(mutate):
    VALIDATOR.validate(_with(mutate))


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _a_s2_z(d).__setitem__(2, {"value": 6.0, "pickup": "2 * 3"}),
        lambda d: _a_s2_z(d).__setitem__(2, {"param": "A_S2_Z", "value": 6.0}),
        lambda d: _a_s2_z(d).__setitem__(2, {"param": 3}),
        lambda d: d["root"]["children"][2].update(pose={"reference": "relative"}),
        lambda d: d["root"]["children"][2].update(pose={"order": "zyx"}),
        lambda d: d.update(configurations=[]),
        lambda d: d.update(configurations=["near"]),
        lambda d: d["parameters"].append({"name": "X"}),
        lambda d: d["parameters"].append({"name": "X", "value": 1.0, "expression": "1"}),
        lambda d: d["parameters"].append({"name": "X", "expression": "1", "min": 0.0}),
        lambda d: d["parameters"].append({"name": "X", "values": 1.0}),
        lambda d: d["parameters"].append({"name": "X", "value": 1.0, "weight": 1.0}),
    ],
    ids=["pickup", "param-with-value", "param-not-string", "reference-value", "order-value",
         "empty-configurations", "configuration-string", "row-without-form", "row-two-forms",
         "bounds-at-expression", "values-not-array", "row-unknown-key"],
)
def test_0_4_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_with(mutate))


# ---- schema 0.4 (#162 part B): the merit function (ADR 0030) ----

def _ops(d):
    return d["optimization"]["operands"]


def _gens(d):
    return d["optimization"]["generators"]


def test_tour_has_every_operand_and_generator_type():
    doc = _tour()
    assert {o["type"] for o in _ops(doc)} == {
        "efl", "bfl", "image_fnumber", "magnification", "ray_x", "ray_y", "spot_rms", "opd_rms",
        "param_value"}
    assert {g["type"] for g in _gens(doc)} == {"rms_spot", "rms_wavefront"}


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.pop("optimization"),
        lambda d: d["optimization"].pop("generators"),
        lambda d: d["optimization"].pop("operands"),
        lambda d: _ops(d).append({"type": "ray_y", "path": "p", "surface": "S", "target": 1.0,
                                  "weight": 0.0, "px": -2.5, "occurrence": 0}),
        lambda d: _ops(d).append({"type": "spot_rms", "path": "p", "wavelength": 1,
                                  "reference": "centroid", "target": 0.0}),
        lambda d: _gens(d).append({"type": "rms_wavefront", "path": "p", "wavelengths": [1, 0]}),
    ],
    ids=["no-section", "no-generators", "no-operands", "ray-all-keys", "spot-wavelength",
         "wavefront-selection"],
)
def test_merit_function_forms_are_accepted(mutate):
    VALIDATOR.validate(_with(mutate))


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _ops(d).append({"type": "focus", "path": "p", "target": 0.0}),
        lambda d: _ops(d).append({"type": "efl", "path": "p"}),
        lambda d: _ops(d).append({"type": "efl", "target": 0.0}),
        lambda d: _ops(d).append({"type": "efl", "path": "p", "target": 0.0, "grid": 3}),
        lambda d: _ops(d).append({"type": "param_value", "parameter": "X", "path": "p",
                                  "target": 0.0}),
        lambda d: _ops(d).append({"type": "spot_rms", "path": "p", "polychromatic": True,
                                  "wavelength": 0, "target": 0.0}),
        lambda d: _ops(d).append({"type": "ray_x", "path": "p", "target": 0.0}),
        lambda d: _ops(d).append({"type": "opd_rms", "path": "p", "field": -1, "target": 0.0}),
        lambda d: _ops(d).append({"type": "opd_rms", "path": "p", "field": 65536,
                                  "target": 0.0}),
        lambda d: _ops(d).append({"type": "opd_rms", "path": "p", "grid": 1.5, "target": 0.0}),
        lambda d: _ops(d).append({"type": "spot_rms", "path": "p", "reference": "best",
                                  "target": 0.0}),
        lambda d: _gens(d).append({"type": "rms_spot", "path": "p", "target": 0.0}),
        lambda d: _gens(d).append({"type": "rms_wavefront", "path": "p", "reference": "chief"}),
        lambda d: _gens(d).append({"type": "rms_spot", "path": "p", "fields": []}),
        lambda d: _gens(d).append({"type": "rms_spot"}),
        lambda d: d["optimization"].update(extra=[]),
    ],
    ids=["unknown-type", "operand-without-target", "operand-without-path", "foreign-key",
         "param-value-with-path", "polychromatic-with-wavelength", "ray-without-surface",
         "negative-index", "index-too-large", "grid-not-integer", "reference-value",
         "generator-with-target", "wavefront-with-reference", "empty-selection",
         "generator-without-path", "section-unknown-key"],
)
def test_merit_function_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_with(mutate))

# ---- schema 0.5 (#178, ADR 0031): ideal lens and ideal cylinder lens ----

def _ideal(mutate):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "r2" / "ideal_lens.rtt.json"))
    mutate(doc)
    return doc


def _il(d):
    return d["root"]["children"][1]["surfaces"][0]["interaction"]


def _cl(d):
    return d["root"]["children"][2]["surfaces"][0]["interaction"]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _il(d).pop("object_distance"),
        lambda d: _il(d).update(focal_length={"value": 50.0, "variable": True, "min": 20.0}),
        lambda d: _il(d).update(focal_length={"param": "F"}, object_distance={"param": "OBJ"}),
        lambda d: _il(d).update(focal_length=-50.0, object_distance=-30.0),
        lambda d: _cl(d).update(axis_deg=30.0),
    ],
    ids=["infinity", "variable-f", "bound", "negative", "cylinder-axis"],
)
def test_0_5_ideal_lens_forms_are_accepted(mutate):
    VALIDATOR.validate(_ideal(mutate))


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _il(d).pop("focal_length"),
        lambda d: _il(d).update(focal_length="50"),
        lambda d: _il(d).update(axis_deg=30.0),
        lambda d: _cl(d).update(axis_deg={"value": 30.0}),
        lambda d: _il(d).update(object_distance=True),
        lambda d: _il(d).update(weight=1.0),
    ],
    ids=["missing-f", "f-string", "axis-at-ideal-lens", "axis-not-number", "distance-bool",
         "unknown-key"],
)
def test_0_5_ideal_lens_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_ideal(mutate))
