"""The merit function from Python (#162 part B, ADR 0030): read-only types of the section
"optimization", their defaults from the model, and validate's merit codes."""

from __future__ import annotations

import json

from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import model


def tour() -> rt.System:
    return rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")


def test_operands_of_the_tour() -> None:
    ops = tour().optimization.operands
    efl, bfl, fno, mag, ray_x, ray_y, spot, opd, value = ops
    assert isinstance(efl, model.FirstOrderOperand)
    assert (efl.quantity, efl.path, efl.target, efl.weight, efl.configuration, efl.wavelength) == (
        model.FirstOrderQuantity.EFL, "first order", 50.0, 1.0, None, None)
    assert isinstance(bfl, model.FirstOrderOperand)
    assert (bfl.quantity, bfl.configuration, bfl.wavelength, bfl.target, bfl.weight) == (
        model.FirstOrderQuantity.BFL, "far", 1, 30.0, 0.5)
    assert isinstance(fno, model.FirstOrderOperand)
    assert fno.quantity == model.FirstOrderQuantity.IMAGE_F_NUMBER
    assert isinstance(mag, model.FirstOrderOperand)
    assert mag.quantity == model.FirstOrderQuantity.MAGNIFICATION
    assert isinstance(ray_x, model.RayOperand)
    assert (ray_x.coordinate, ray_x.path, ray_x.surface, ray_x.occurrence, ray_x.field, ray_x.px,
            ray_x.py, ray_x.wavelength) == (model.RayCoordinate.X, "ghost G.S2-A.S2", "A.S2", 1, 1,
                                            0.5, 0.0, None)
    assert isinstance(ray_y, model.RayOperand)
    assert (ray_y.coordinate, ray_y.occurrence, ray_y.field, ray_y.py) == (
        model.RayCoordinate.Y, None, 0, 1.0)
    assert isinstance(spot, model.SpotRmsOperand)
    assert (spot.field, spot.wavelength, spot.polychromatic, spot.reference, spot.rings) == (
        1, None, True, model.SpotReference.CHIEF, 4)
    assert isinstance(opd, model.OpdRmsOperand)
    assert (opd.field, opd.wavelength, opd.grid) == (0, 0, 17)
    assert isinstance(value, model.ParamValueOperand)
    assert (value.parameter, value.configuration, value.target) == ("FOCUS", "far", 0.25)


def test_generators_of_the_tour() -> None:
    spots, waves = tour().optimization.generators
    assert isinstance(spots, model.SpotGenerator)
    assert (spots.path, spots.configuration, spots.fields, spots.wavelengths, spots.reference,
            spots.rings, spots.arms, spots.weight) == (
        "first order", None, [0, 1], [0], model.SpotReference.CHIEF, 4, 8, 2.0)
    assert isinstance(waves, model.WavefrontGenerator)
    assert (waves.configuration, waves.fields, waves.wavelengths, waves.rings, waves.arms,
            waves.weight) == ("near", None, None, 3, 6, 1.0)


def test_defaults_come_from_the_model() -> None:
    # ADR 0030 and #162 B: the defaults live in rtt/model/optimization.hpp only.
    data = json.loads((REFERENCE_DIR / "m0" / "feature_tour.rtt.json").read_text(encoding="utf-8"))
    data["optimization"] = {
        "operands": [{"type": "spot_rms", "path": "first order", "target": 0.0},
                     {"type": "opd_rms", "path": "first order", "target": 0.0}],
        "generators": [{"type": "rms_spot", "path": "first order"}]}
    s = rt.System.from_json(json.dumps(data))
    spot, opd = s.optimization.operands
    assert isinstance(spot, model.SpotRmsOperand)
    assert (spot.field, spot.polychromatic, spot.reference, spot.rings, spot.weight) == (
        0, False, model.SpotReference.CENTROID, 6, 1.0)
    assert isinstance(opd, model.OpdRmsOperand)
    assert opd.grid == 33
    (gen,) = s.optimization.generators
    assert isinstance(gen, model.SpotGenerator)
    assert (gen.rings, gen.arms, gen.weight, gen.reference) == (3, 6, 1.0,
                                                               model.SpotReference.CENTROID)


def test_no_merit_function() -> None:
    s = rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")
    assert s.optimization.operands == [] and s.optimization.generators == []
    assert s.to_dict()["optimization"] == {"operands": [], "generators": []}
    assert "optimization" not in json.loads(s.to_json())


def test_edit_form_and_validate() -> None:
    # The edit form holds every value; an operand added by a patch is checked by validate.
    s = rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")
    patched = rt.apply_patch(s, [{"op": "add", "path": "/optimization/operands/-",
                                  "value": {"type": "efl", "path": "main", "target": 100.0}}])
    efl = patched.optimization.operands[0]
    assert isinstance(efl, model.FirstOrderOperand) and efl.target == 100.0
    bad = rt.System.from_json(json.dumps(
        {**json.loads(patched.to_json()),
         "optimization": {"operands": [{"type": "efl", "path": "none", "target": 1.0}]}}))
    codes = [(d.code, d.location) for d in rt.validate(bad)]
    assert ("merit.unknown_path", "/optimization/operands/0/path") in codes
