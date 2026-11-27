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


# ---- schema 0.6 (#253, ADR 0035): the tolerance set ----

def _tol(mutate):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m7" / "achromat_tolerances.rtt.json"))
    mutate(doc)
    return doc


def _t(d):
    return d["tolerances"]


def test_tolerance_reference_has_every_form():
    t = _t(_tol(lambda d: None))
    assert set(t) == {"class", "overrides", "items", "compensators", "compensator_merit",
                      "criteria", "monte_carlo"}
    forms = [("class" in o, "kind" in o, "pivot" in o) for o in t["overrides"]]
    assert (True, False, False) in forms
    assert (False, True, False) in forms
    assert (False, False, True) in forms
    assert {a for i in t["items"] for a in ("surface", "node", "parameter") if a in i} == {
        "surface", "parameter"}


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d.pop("tolerances"),
        lambda d: d.update(tolerances={}),
        lambda d: _t(d).update({"class": "none"}),
        lambda d: _t(d)["items"].append({"id": "x", "node": "L1", "kind": "index", "segment": 0,
                                         "limit": 0.001}),
        lambda d: _t(d)["items"].append({"id": "x", "node": "L1", "kind": "shift_z",
                                         "limit": 0.01, "distribution": "normal"}),
        lambda d: _t(d)["overrides"].append({"surface": "L1.S1", "pivot": [0.0, 0.0, 0.0]}),
        lambda d: _t(d)["compensators"].append({"node": "image", "param": "/pose/position/2",
                                                "travel_min": -0.5, "travel_max": 0.5}),
        lambda d: _t(d).pop("compensator_merit"),
        lambda d: _t(d)["compensator_merit"].update(generators=[{"type": "rms_spot",
                                                                 "path": "main"}]),
        lambda d: _t(d)["criteria"].append({"id": "y", "operand": {
            "type": "ray_y", "path": "main", "surface": "IMG", "configuration": "c"},
            "min": -0.1, "max": 0.1}),
        lambda d: _t(d).update(monte_carlo={}),
        lambda d: _t(d).update(monte_carlo={"seed": 18446744073709551615, "trials": 0}),
    ],
    ids=["no-section", "empty-section", "class-none", "segment", "default-distribution",
         "surface-pivot", "node-param-compensator", "no-compensator-merit", "merit-generator",
         "ray-criterion", "empty-monte-carlo", "uint64-range"],
)
def test_0_6_tolerance_forms_are_accepted(mutate):
    VALIDATOR.validate(_tol(mutate))


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: _t(d).update(extra=1),
        lambda d: _t(d).update({"class": "fine"}),
        lambda d: _t(d).update(items={}),
        lambda d: _t(d)["items"][0].pop("surface"),
        lambda d: _t(d)["items"][0].update(node="L1"),
        lambda d: _t(d)["items"][0].pop("id"),
        lambda d: _t(d)["items"][0].pop("limit"),
        lambda d: _t(d)["items"][0].update(kind="wedge"),
        lambda d: _t(d)["items"][0].update(distribution="gauss"),
        lambda d: _t(d)["items"][0].update(segment=-1),
        lambda d: _t(d)["items"][0].update(pivot=[0.0, 0.0, 1.0]),
        lambda d: _t(d)["overrides"][0].update(kind="index"),
        lambda d: _t(d)["overrides"][3].update(kind="tilt_x_deg", limit=0.1),
        lambda d: _t(d)["overrides"].append({"node": "L1"}),
        lambda d: _t(d)["overrides"].append({"surface": "L1.S1", "class": "coarse"}),
        lambda d: _t(d)["overrides"][1].update(parameter="BFD"),
        lambda d: _t(d)["overrides"][1].update(node="L1"),
        lambda d: _t(d)["overrides"][3].update(pivot=[0.0, 0.0]),
        lambda d: _t(d)["compensators"][0].pop("travel_min"),
        lambda d: _t(d)["compensators"][0].update(segment=0),
        lambda d: _t(d)["compensator_merit"]["operands"][0].pop("target"),
        lambda d: _t(d)["criteria"][0]["operand"].update(target=0.0),
        lambda d: _t(d)["criteria"][0]["operand"].update(weight=1.0),
        lambda d: _t(d)["criteria"][0]["operand"].update(type="rms_spot"),
        lambda d: _t(d)["criteria"][0].pop("operand"),
        lambda d: _t(d)["criteria"][0].update(target=0.0),
        lambda d: _t(d)["monte_carlo"].update(seed=-1),
        lambda d: _t(d)["monte_carlo"].update(seed=18446744073709551616),
        lambda d: _t(d)["monte_carlo"].update(trials=1.5),
        lambda d: _t(d)["monte_carlo"].update(threads=4),
        lambda d: d.update(schema_version="0.5.0"),
    ],
    ids=["section-unknown-key", "class-name", "items-not-array", "item-without-anchor",
         "item-two-anchors", "item-without-id", "item-without-limit", "kind-name",
         "distribution-name", "negative-segment", "item-pivot", "class-and-kind",
         "pivot-and-kind", "override-without-form", "class-at-surface", "override-parameter",
         "override-two-anchors", "pivot-length", "compensator-without-travel",
         "compensator-segment", "merit-operand-without-target", "criterion-target",
         "criterion-weight", "criterion-generator", "criterion-without-operand",
         "criterion-unknown-key", "negative-seed", "seed-too-large", "trials-not-integer",
         "monte-carlo-unknown-key", "old-version"],
)
def test_0_6_tolerance_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_tol(mutate))


# ---- schema 0.7 (#282, ADR 0037): the medium of the object space ----

def _obj(mutate):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m1" / "singlet_const.rtt.json"))
    doc["object"] = {"at_infinity": False, "distance": 50.0}
    mutate(doc)
    return doc


def test_0_7_object_medium_is_accepted():
    assert VALIDATOR.is_valid(_obj(lambda d: d["object"].update(medium="CONST:1.53")))
    assert VALIDATOR.is_valid(_obj(lambda d: None))


@pytest.mark.parametrize(
    "mutate",
    [lambda d: d["object"].update(medium=1.53),
     lambda d: d["object"].update(medium=None),
     lambda d: d["object"].update(medium=["CONST:1.5"]),
     lambda d: d["object"].update(unknown_key=1)],
    ids=["number", "null", "list", "unknown-key"],
)
def test_0_7_object_medium_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_obj(mutate))


# ---- schema 0.8 (#288, ADR 0038): the reference of the object distance ----

def test_0_8_object_reference_is_accepted():
    assert VALIDATOR.is_valid(_obj(lambda d: d["object"].update(reference="origin")))
    assert VALIDATOR.is_valid(_obj(lambda d: d["object"].update(reference={"surface": "L1.S1"})))


@pytest.mark.parametrize(
    "mutate",
    [lambda d: d["object"].update(reference="vertex"),
     lambda d: d["object"].update(reference=1.0),
     lambda d: d["object"].update(reference={}),
     lambda d: d["object"].update(reference={"node": "L1"}),
     lambda d: d["object"].update(reference={"surface": "L1.S1", "node": "L1"}),
     lambda d: d["object"].update(reference={"surface": 3})],
    ids=["other-text", "number", "empty", "node", "extra-key", "surface-number"],
)
def test_0_8_object_reference_errors_are_rejected(mutate):
    assert not VALIDATOR.is_valid(_obj(mutate))


# ---- schema 0.9 (#315, ADR 0042): the anamorphic base shapes ----

def _base(base):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    _lens(doc)["surfaces"][0]["shape"]["base"] = base
    return doc


@pytest.mark.parametrize(
    "base",
    [{"type": "biconic", "radius_x": 40.0, "radius_y": -60.0},
     {"type": "biconic", "radius_x": 40.0, "radius_y": {"param": "RX"}, "conic_x": -0.5,
      "conic_y": {"value": 0.2, "variable": True}},
     {"type": "toroid", "radius": 48.0, "radius_of_rotation": 55.0},
     {"type": "toroid", "radius": 48.0, "conic": -1.2, "coefficients": [2e-6, -3e-9],
      "radius_of_rotation": -30.0},
     {"type": "cylinder", "radius": -75.0},
     {"type": "cylinder", "radius": -75.0, "conic": 0.3, "coefficients": [1e-6]}],
    ids=["biconic", "biconic-all", "toroid", "toroid-all", "cylinder", "cylinder-all"],
)
def test_0_9_anamorphic_shapes_are_accepted(base):
    assert VALIDATOR.is_valid(_base(base))


@pytest.mark.parametrize(
    "base",
    [{"type": "biconic", "radius_x": 40.0},
     {"type": "biconic", "radius_x": 40.0, "radius_y": 40.0, "radius": 40.0},
     {"type": "biconic", "radius_x": 40.0, "radius_y": 40.0, "coefficients": [1e-6]},
     {"type": "toroid", "radius": 48.0},
     {"type": "toroid", "radius_of_rotation": 48.0},
     {"type": "toroid", "radius": 48.0, "radius_of_rotation": 55.0, "conic_x": 0.0},
     {"type": "cylinder"},
     {"type": "cylinder", "radius": 10.0, "radius_of_rotation": 55.0}],
    ids=["biconic-no-y", "biconic-radius", "biconic-coefficients", "toroid-no-rotation",
         "toroid-no-radius", "toroid-conic-x", "cylinder-no-radius", "cylinder-rotation"],
)
def test_0_9_anamorphic_shape_errors_are_rejected(base):
    assert not VALIDATOR.is_valid(_base(base))


# ---- schema 0.9 (#314, ADR 0043): odd asphere and axicon ----

@pytest.mark.parametrize(
    "base",
    [{"type": "odd_asphere", "curvature": 0.02},
     {"type": "odd_asphere", "curvature": 0.0, "conic": -0.5, "coefficients": [1e-3, 0.0, 2e-6]},
     {"type": "odd_asphere", "curvature": {"param": "C"},
      "conic": {"value": 0.2, "variable": True}, "coefficients": [{"value": 1e-3}]},
     {"type": "axicon", "angle_deg": -5.0},
     {"type": "axicon", "angle_deg": {"value": 10.0, "variable": True}}],
    ids=["odd", "odd-all", "odd-params", "axicon", "axicon-param"],
)
def test_0_9_rotational_shapes_are_accepted(base):
    assert VALIDATOR.is_valid(_base(base))


@pytest.mark.parametrize(
    "base",
    [{"type": "odd_asphere"},
     {"type": "odd_asphere", "radius": 50.0},
     {"type": "odd_asphere", "curvature": 0.02, "radius": 50.0},
     {"type": "odd_asphere", "curvature": "0.02"},
     {"type": "odd_asphere", "curvature": 0.02, "coefficients": 1e-3},
     {"type": "axicon"},
     {"type": "axicon", "angle": 5.0},
     {"type": "axicon", "angle_deg": 5.0, "conic": 0.0}],
    ids=["odd-empty", "odd-radius", "odd-both", "odd-string", "odd-scalar-coefficients",
         "axicon-empty", "axicon-angle", "axicon-conic"],
)
def test_0_9_rotational_shape_errors_are_rejected(base):
    assert not VALIDATOR.is_valid(_base(base))


# ---- schema 0.9 (#314, ADR 0043): Q-con and Q-bfs ----

@pytest.mark.parametrize(
    "base",
    [{"type": "qcon", "curvature": 0.02, "normalization_radius": 12.7},
     {"type": "qcon", "curvature": 0.0, "conic": -0.5, "normalization_radius": 5.0,
      "coefficients": [1e-3, -2e-4]},
     {"type": "qbfs", "curvature": {"param": "C"}, "conic": {"value": 0.2, "variable": True},
      "normalization_radius": {"value": 12.7}, "coefficients": [{"value": 1e-3}]}],
    ids=["qcon", "qcon-all", "qbfs-params"],
)
def test_0_9_forbes_shapes_are_accepted(base):
    assert VALIDATOR.is_valid(_base(base))


@pytest.mark.parametrize(
    "base",
    [{"type": "qcon", "curvature": 0.02},
     {"type": "qbfs", "normalization_radius": 5.0},
     {"type": "qcon", "radius": 50.0, "normalization_radius": 5.0},
     {"type": "qbfs", "curvature": 0.02, "normalization_radius": 5.0, "angle_deg": 1.0},
     {"type": "qcon", "curvature": 0.02, "normalization_radius": "5"}],
    ids=["qcon-no-rho", "qbfs-no-curvature", "qcon-radius", "qbfs-angle", "qcon-string"],
)
def test_0_9_forbes_shape_errors_are_rejected(base):
    assert not VALIDATOR.is_valid(_base(base))


def test_active_flag_of_nodes():
    """Schema 0.10 (#380, ADR 0049): "active" is a boolean on assemblies and elements; the schema
    accepts false also at the root, which rtt validate rejects (node.root_inactive)."""
    doc = _broken(lambda d: _lens(d).update(active=False))
    assert VALIDATOR.is_valid(doc)
    doc["root"]["active"] = True
    assert VALIDATOR.is_valid(doc)
    assert not VALIDATOR.is_valid(_broken(lambda d: _lens(d).update(active="no")))
    assert not VALIDATOR.is_valid(_broken(lambda d: d.update(schema_version="0.9.0")))


# ---- schema 0.11 (#320, ADR 0050, point 1): aperture dimensions are Params ----

def _aperture(aperture):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    _lens(doc)["surfaces"][0]["aperture"] = aperture
    return doc


@pytest.mark.parametrize(
    "aperture",
    [{"type": "circular", "radius": {"value": 12.7, "variable": True, "min": 5.0},
      "inner_radius": {"param": "R_IN"}},
     {"type": "rectangular", "half_width_x": {"value": 5.0}, "half_width_y": 4.0},
     {"type": "elliptical", "semi_axis_x": 6.0, "semi_axis_y": {"param": "B"}}],
    ids=["circular", "rectangular", "elliptical"],
)
def test_0_11_aperture_params_are_accepted(aperture):
    assert VALIDATOR.is_valid(_aperture(aperture))


@pytest.mark.parametrize(
    "aperture",
    [{"type": "circular", "radius": "12.7"},
     {"type": "circular", "radius": {"variable": True}},
     {"type": "rectangular", "half_width_x": {"param": "A", "value": 1.0}, "half_width_y": 1.0},
     {"type": "elliptical", "semi_axis_x": [6.0], "semi_axis_y": 1.0}],
    ids=["string", "no-value", "param-and-value", "list"],
)
def test_0_11_aperture_param_errors_are_rejected(aperture):
    assert not VALIDATOR.is_valid(_aperture(aperture))


def test_0_11_version_pattern():
    assert VALIDATOR.is_valid(_broken(lambda d: None))
    assert not VALIDATOR.is_valid(_broken(lambda d: d.update(schema_version="0.10.0")))


# ---- schema 0.11 (#320, ADR 0050, point 2): regular_polygon and polygon apertures ----

@pytest.mark.parametrize(
    "aperture",
    [{"type": "regular_polygon", "sides": 6, "circumradius": 12.7},
     {"type": "regular_polygon", "sides": 5, "circumradius": {"value": 5.0, "variable": True},
      "rotation_deg": 10.0},
     {"type": "polygon", "vertices": [[0.0, 0.0], [4.0, 0.0], [4.0, 2.0], [2.0, 2.0],
                                      [2.0, 4.0], [0.0, 4.0]]}],
    ids=["hexagon", "pentagon-param", "l-shape"],
)
def test_0_11_polygon_apertures_are_accepted(aperture):
    assert VALIDATOR.is_valid(_aperture(aperture))


@pytest.mark.parametrize(
    "aperture",
    [{"type": "regular_polygon", "circumradius": 5.0},
     {"type": "regular_polygon", "sides": 6.5, "circumradius": 5.0},
     {"type": "regular_polygon", "sides": 6, "circumradius": 5.0, "rotation_deg": "30"},
     {"type": "polygon"},
     {"type": "polygon", "vertices": [[0.0, 0.0, 1.0]]},
     {"type": "polygon", "vertices": [[0.0, {"value": 1.0}]]},
     {"type": "hexagonal", "circumradius": 5.0}],
    ids=["no-sides", "fractional-sides", "string-rotation", "no-vertices", "three-coordinates",
         "param-vertex", "hexagonal"],
)
def test_0_11_polygon_aperture_errors_are_rejected(aperture):
    assert not VALIDATOR.is_valid(_aperture(aperture))


# ---- schema 0.11 (#320, ADR 0050, point 3): obscurations ----

def _obscurations(obscurations):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    _lens(doc)["surfaces"][0]["obscurations"] = obscurations
    return doc


@pytest.mark.parametrize(
    "obscurations",
    [[],
     [{"type": "circular", "radius": 2.0}],
     [{"type": "circular", "radius": {"value": 2.0, "variable": True},
       "decenter_x": {"param": "DX"}, "decenter_y": -1.0}],
     [{"type": "rectangular", "half_width_x": 12.0, "half_width_y": 0.5, "rotation_deg": 45.0},
      {"type": "elliptical", "semi_axis_x": 3.0, "semi_axis_y": 1.5, "decenter_x": 5.0}]],
    ids=["empty", "disc", "disc-params", "spider-and-clip"],
)
def test_0_11_obscurations_are_accepted(obscurations):
    assert VALIDATOR.is_valid(_obscurations(obscurations))


@pytest.mark.parametrize(
    "obscurations",
    [{"type": "circular", "radius": 2.0},
     [{"type": "circular"}],
     [{"type": "circular", "radius": 2.0, "rotation_deg": 10.0}],
     [{"type": "rectangular", "half_width_x": 1.0}],
     [{"type": "elliptical", "semi_axis_x": 1.0, "semi_axis_y": 1.0, "rotation_deg": "10"}],
     [{"type": "spider", "width": 1.0, "arms": 4}]],
    ids=["not-a-list", "no-radius", "rotated-circle", "one-half-width", "string-rotation",
         "spider"],
)
def test_0_11_obscuration_errors_are_rejected(obscurations):
    assert not VALIDATOR.is_valid(_obscurations(obscurations))


# ---- schema 0.11 (#316, ADR 0051, point 3): zernike_fringe_sag ----

def _terms(terms):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    _lens(doc)["surfaces"][0]["shape"] = {"base": {"type": "plane"}, "terms": terms}
    return doc


@pytest.mark.parametrize(
    "terms",
    [[{"type": "zernike_fringe_sag", "normalization_radius": 10.0}],
     [{"type": "zernike_fringe_sag", "normalization_radius": {"value": 10.0},
       "coefficients": [0.0, {"value": 1e-4, "variable": True}]},
      {"type": "zernike_sag", "normalization_radius": 10.0, "coefficients": [0.0, 1e-4]}]],
    ids=["fringe", "fringe-and-noll"],
)
def test_0_11_fringe_terms_are_accepted(terms):
    assert VALIDATOR.is_valid(_terms(terms))


@pytest.mark.parametrize(
    "terms",
    [[{"type": "zernike_fringe_sag"}],
     [{"type": "zernike_fringe_sag", "normalization_radius": 10.0, "order": "fringe"}],
     [{"type": "zernike_standard_sag", "normalization_radius": 10.0}]],
    ids=["no-radius", "extra-key", "unknown-type"],
)
def test_0_11_fringe_term_errors_are_rejected(terms):
    assert not VALIDATOR.is_valid(_terms(terms))


# ---- schema 0.11 (#316, ADR 0051, point 4): xy_polynomial ----

@pytest.mark.parametrize(
    "terms",
    [[{"type": "xy_polynomial", "normalization_radius": 10.0}],
     [{"type": "xy_polynomial", "normalization_radius": {"value": 10.0},
       "coefficients": [0.0, {"value": 1e-4, "variable": True}, 2e-4]},
      {"type": "zernike_sag", "normalization_radius": 10.0, "coefficients": [0.0, 1e-4]}]],
    ids=["xy", "xy-and-noll"],
)
def test_0_11_xy_terms_are_accepted(terms):
    assert VALIDATOR.is_valid(_terms(terms))


@pytest.mark.parametrize(
    "terms",
    [[{"type": "xy_polynomial"}],
     [{"type": "xy_polynomial", "normalization_radius": 10.0, "degree": 4}],
     [{"type": "xy_polynomial", "normalization_radius": 10.0, "coefficients": "a"}]],
    ids=["no-radius", "extra-key", "bad-coefficients"],
)
def test_0_11_xy_term_errors_are_rejected(terms):
    assert not VALIDATOR.is_valid(_terms(terms))


# ---- schema 0.11 (#316, ADR 0051, point 5): grid_sag ----

_GRID = {"type": "grid_sag", "nx": 3, "ny": 3, "spacing_x": 1.0, "spacing_y": 1.0,
         "sag": [0.0] * 9}


@pytest.mark.parametrize(
    "terms",
    [[_GRID],
     [dict(_GRID, decenter_x=0.5, decenter_y=-0.5, dz_dx=[0.0] * 9, dz_dy=[0.0] * 9,
           d2z_dxdy=[0.0] * 9)],
     # Structurally readable; sizes, spacings, lengths and the derivative rule are validate's.
     [dict(_GRID, nx=2)],
     [dict(_GRID, spacing_x=0.0)],
     [dict(_GRID, dz_dx=[0.0] * 9)]],
    ids=["grid", "grid-with-derivatives", "nx-too-small", "zero-spacing", "partial-derivatives"],
)
def test_0_11_grid_terms_are_accepted(terms):
    assert VALIDATOR.is_valid(_terms(terms))


@pytest.mark.parametrize(
    "terms",
    [[{k: v for k, v in _GRID.items() if k != "sag"}],
     [dict(_GRID, nx=3.5)],
     [dict(_GRID, sag=[{"value": 0.0}] * 9)],
     [dict(_GRID, normalization_radius=1.0)]],
    ids=["no-sag", "nx-not-integer", "param-values", "extra-key"],
)
def test_0_11_grid_term_errors_are_rejected(terms):
    assert not VALIDATOR.is_valid(_terms(terms))


# ---- schema 0.11 (#316, ADR 0051, point 7): tolerance kind zernike ----

def _with_items(items):
    doc = copy.deepcopy(load(ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json"))
    doc["tolerances"] = {"items": items}
    return doc


def test_0_11_zernike_tolerance_is_accepted():
    assert VALIDATOR.is_valid(_with_items(
        [{"id": "z", "surface": "L1.S1", "term": 4, "kind": "zernike", "limit": 2e-4}]))
    # The range of term and the kind it belongs to are validate's.
    assert VALIDATOR.is_valid(_with_items(
        [{"id": "z", "surface": "L1.S1", "term": 1, "kind": "conic", "limit": 2e-4}]))


@pytest.mark.parametrize("term", [-1, 2.5, "4"], ids=["negative", "fraction", "string"])
def test_0_11_zernike_term_errors_are_rejected(term):
    assert not VALIDATOR.is_valid(_with_items(
        [{"id": "z", "surface": "L1.S1", "term": term, "kind": "zernike", "limit": 2e-4}]))


def test_0_12_object_slope_is_accepted():
    doc = load(ROOT / "tests" / "reference" / "m1" / "tilted_object_ideal_lens.rtt.json")
    assert doc["object"]["slope"] == [0.0, -1.0]
    VALIDATOR.validate(doc)
    doc["object"]["slope"] = [{"value": 0.1, "variable": True}, {"param": "tilt"}]
    VALIDATOR.validate(doc)


@pytest.mark.parametrize("mutate", [
    lambda o: o.__setitem__("slope", [1.0]),
    lambda o: o.__setitem__("slope", [1.0, 2.0, 3.0]),
    lambda o: o.__setitem__("slope", 1.0),
    lambda o: o.__setitem__("slope", [1.0, "a"]),
])
def test_0_12_object_slope_errors_are_rejected(mutate):
    doc = load(ROOT / "tests" / "reference" / "m1" / "tilted_object_ideal_lens.rtt.json")
    mutate(doc["object"])
    with pytest.raises(jsonschema.ValidationError):
        VALIDATOR.validate(doc)


def test_0_13_constraint_is_accepted():
    # m5/singlet_constraint (#196, ADR 0056): the efl as an exact constraint, without a weight
    doc = load(ROOT / "tests" / "reference" / "m5" / "singlet_constraint.rtt.json")
    assert doc["optimization"]["operands"][0]["constraint"] == "equal"
    VALIDATOR.validate(doc)


@pytest.mark.parametrize("mutate", [
    lambda o: o.__setitem__("weight", 1.0),        # a constraint has no weight
    lambda o: o.__setitem__("constraint", "less"),  # only "equal" (ADR 0056, point 6)
])
def test_0_13_constraint_rejects(mutate):
    doc = load(ROOT / "tests" / "reference" / "m5" / "singlet_constraint.rtt.json")
    mutate(doc["optimization"]["operands"][0])
    assert not VALIDATOR.is_valid(doc)


def test_0_13_criterion_has_no_constraint():
    doc = load(ROOT / "tests" / "reference" / "m5" / "singlet_constraint.rtt.json")
    doc["tolerances"] = {"criteria": [{"id": "efl", "max": 101.0, "operand": {
        "type": "efl", "path": "main", "constraint": "equal"}}]}
    assert not VALIDATOR.is_valid(doc)


def test_0_13_constraint_in_the_compensator_merit_is_accepted():
    doc = load(ROOT / "tests" / "reference" / "m5" / "singlet_constraint.rtt.json")
    doc["tolerances"] = {"compensator_merit": {"operands": [
        {"type": "efl", "path": "main", "target": 100.0, "constraint": "equal"}]}}
    VALIDATOR.validate(doc)
