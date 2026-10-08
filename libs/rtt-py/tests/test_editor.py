"""Changing a system with JSON Patch, undo and redo in Python (#82, ADR 0024 points 3-5)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from conftest import REFERENCE_DIR, reference_files

import raytatouille as rt


def canonical(system: rt.System) -> str:
    return system.to_json()


def load(relative: str) -> rt.System:
    return rt.load(REFERENCE_DIR / relative)


# ---------------------------------------------------------------- apply_patch -----


def test_apply_patch_returns_a_new_system_and_leaves_the_input() -> None:
    s = load("m0/singlet.rtt.json")
    before = canonical(s)
    changed = rt.apply_patch(s, [{"op": "replace", "path": "/name", "value": "renamed"}])
    assert changed.name == "renamed"
    assert canonical(s) == before
    # JSON text works as well.
    assert rt.apply_patch(s, '[{"op": "replace", "path": "/name", "value": "x"}]').name == "x"


def test_apply_patch_with_inverse() -> None:
    s = load("m0/singlet.rtt.json")
    r = rt.apply_patch_with_inverse(s, [{"op": "add", "path": "/wavelengths/-",
                                         "value": {"um": 0.55}}])
    assert isinstance(r, rt.PatchResult)
    assert len(r.system.wavelengths) == len(s.wavelengths) + 1
    assert r.inverse == [{"op": "remove", "path": f"/wavelengths/{len(s.wavelengths)}"}]
    assert all(isinstance(d, rt.Diagnostic) for d in r.diagnostics)
    back = rt.apply_patch(r.system, r.inverse, check="structure_only")
    assert canonical(back) == canonical(s)


def test_edit_error() -> None:
    s = load("m0/singlet.rtt.json")
    with pytest.raises(rt.EditError) as info:
        rt.apply_patch(s, [{"op": "replace", "path": "/name", "value": "x"},
                           {"op": "remove", "path": "/root/children/9"}])
    e = info.value
    assert isinstance(e, rt.RaytatouilleError) and isinstance(e, ValueError)
    assert (e.code, e.location, e.op_index) == ("edit.path_not_found", "/root/children/9", 1)
    assert e.diagnostics == []
    # A new error of validate(): its code and place, all its diagnostics, no operation index.
    with pytest.raises(rt.EditError) as dup:
        rt.apply_patch(s, [{"op": "replace", "path": "/root/children/1/surfaces/1/id",
                            "value": "L1.S1"}])
    assert dup.value.code == "surface.id_duplicate"
    assert dup.value.op_index is None
    assert dup.value.diagnostics and all(d.code == "surface.id_duplicate"
                                         for d in dup.value.diagnostics)
    assert dup.value.location == dup.value.diagnostics[0].location
    # NaN is invalid JSON: the C++ parser reports it as an invalid patch (one error path).
    with pytest.raises(rt.EditError) as nan:
        rt.apply_patch(s, [{"op": "replace", "path": "/environment/temperature_c",
                            "value": float("nan")}])
    assert nan.value.code == "edit.patch_invalid"
    with pytest.raises(ValueError):
        rt.apply_patch(s, [], check="sometimes")  # type: ignore[arg-type]
    # A single operation instead of a list is a clear TypeError.
    with pytest.raises(TypeError, match="list of operations"):
        rt.apply_patch(s, {"op": "remove", "path": "/name"})  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="list of operations"):
        rt.Editor(s).apply({"op": "remove", "path": "/name"})  # type: ignore[arg-type]


# --------------------------------------------------------------------- Editor -----


def test_editor_undo_and_redo_are_bit_identical() -> None:
    s = load("m2/cooke_triplet.rtt.json")
    states = [canonical(s)]
    ed = rt.Editor(s)
    assert not ed.can_undo and not ed.can_redo
    ed.set("/name", "edited")
    states.append(canonical(ed.system))
    ed.insert("/fields/points/0", {"y": 3.0})
    states.append(canonical(ed.system))
    ed.remove("/wavelengths/2")
    states.append(canonical(ed.system))
    ed.move("/root/children/0", "/root/children/-")
    states.append(canonical(ed.system))
    ed.apply([{"op": "copy", "from": "/fields/points/1", "path": "/fields/points/-"},
              {"op": "test", "path": "/name", "value": "edited"}])
    states.append(canonical(ed.system))
    assert len(ed.history) == 5
    for k in range(5, 0, -1):
        ed.undo()
        assert canonical(ed.system) == states[k - 1]
    assert not ed.can_undo and ed.can_redo
    for k in range(1, 6):
        ed.redo()
        assert canonical(ed.system) == states[k]
    assert not ed.can_redo
    with pytest.raises(IndexError, match="nothing to redo"):
        ed.redo()


def test_a_new_command_after_undo_drops_redo() -> None:
    ed = rt.Editor(load("m0/singlet.rtt.json"))
    ed.set("/name", "a")
    ed.undo()
    assert ed.can_redo
    ed.set("/name", "b")
    assert not ed.can_redo
    assert ed.history == [[{"op": "replace", "path": "/name", "value": "b"}]]
    ed.undo()
    with pytest.raises(IndexError, match="nothing to undo"):
        ed.undo()


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_editor_on_every_reference_system(file: Path) -> None:
    s = rt.load(file)
    ed = rt.Editor(s)
    first = ed.system.to_dict()["aperture"]["value"]["value"]
    ed.set("/aperture/value", first + 1.0)  # a Param object and a number: writes .../value
    ed.insert("/fields/points/-", {"y": 0.5})
    ed.remove("/root/pose")
    ed.move("/wavelengths/0", "/wavelengths/-")
    after = canonical(ed.system)
    for _ in range(4):
        ed.undo()
    assert canonical(ed.system) == canonical(s)
    for _ in range(4):
        ed.redo()
    assert canonical(ed.system) == after


def test_set_on_a_param_keeps_variable_and_pickup() -> None:
    tour = load("m0/feature_tour.rtt.json")
    ed = rt.Editor(tour)
    ed.set("/object/distance", 300.0)
    assert ed.history[-1] == [{"op": "replace", "path": "/object/distance/value", "value": 300.0}]
    distance = ed.system.object_space.distance
    assert (distance.value, distance.variable) == (300.0, True)
    z = "/root/children/2/surfaces/1/pose/position/2"
    ed.set(z, 7.0)
    assert ed.system.json_at(z) == {"value": 7.0, "variable": False, "pickup": "2 * 3"}
    # A whole Param object replaces it; a bool is no number, so it replaces the object too and
    # the strict reader rejects it.
    ed.set("/object/distance", {"value": 1.0, "variable": False})
    assert not ed.system.object_space.distance.variable
    with pytest.raises(rt.EditError) as info:
        ed.set("/object/distance", True)
    assert info.value.code == "edit.invalid_value"


def test_set_on_a_missing_optional_member_adds_it() -> None:
    ed = rt.Editor(load("m0/feature_tour.rtt.json"))
    image = ed.system.locate_node("image")
    assert image is not None
    ed.set(image + "/surfaces/0/aperture", {"type": "circular", "radius": 5.0})
    assert ed.history[-1][0]["op"] == "add"
    surface = ed.system.root.children[5]
    assert isinstance(surface, rt.model.Element)
    aperture = surface.surfaces[0].aperture
    assert isinstance(aperture, rt.model.CircularAperture) and aperture.radius == 5.0


def test_set_does_not_append_to_an_array_or_create_parents() -> None:
    # set is replace, or add only for a missing member of an object (ADR 0024 point 3).
    ed = rt.Editor(load("m2/cooke_triplet.rtt.json"))
    with pytest.raises(rt.EditError) as past_end:
        ed.set("/fields/points/3", {"y": 1.0})  # 3 points: index 3 does not exist
    assert past_end.value.code == "edit.path_not_found"
    with pytest.raises(rt.EditError) as no_parent:
        ed.set("/environment/nothing/here", 1.0)
    assert no_parent.value.code == "edit.path_not_found"
    assert ed.history == [] and len(ed.system.fields.points) == 3


def test_a_text_patch_goes_to_the_strict_parser() -> None:
    ed = rt.Editor(load("m0/singlet.rtt.json"))
    # RFC 6902 A.13: a duplicate key is an invalid patch, not "the last one wins".
    with pytest.raises(rt.EditError) as duplicate:
        ed.apply('[{"op": "replace", "path": "/name", "value": "x", "op": "remove"}]')
    assert duplicate.value.code == "edit.patch_invalid"
    with pytest.raises(rt.EditError) as broken:
        ed.apply('[{"op": ')
    assert broken.value.code == "edit.patch_invalid"
    ed.apply('[{"op": "replace", "path": "/name", "value": "x"}]')
    assert ed.history == [[{"op": "replace", "path": "/name", "value": "x"}]]


def test_a_system_that_cannot_be_written() -> None:
    s = load("m0/singlet.rtt.json")
    s.environment.temperature_c = float("nan")  # possible through the mutable Environment
    with pytest.raises(rt.EditError) as info:
        rt.apply_patch(s, [{"op": "test", "path": "/name", "value": s.name}])
    assert info.value.code == "edit.base_not_representable"
    with pytest.raises(ValueError):
        rt.Editor(s)


def test_an_error_leaves_the_editor_unchanged() -> None:
    ed = rt.Editor(load("m0/singlet.rtt.json"))
    ed.set("/name", "a")
    state, history = canonical(ed.system), ed.history
    with pytest.raises(rt.EditError):
        ed.apply([{"op": "replace", "path": "/name", "value": "b"},
                  {"op": "remove", "path": "/paths/7"}])
    with pytest.raises(rt.EditError):
        ed.set("/root/children/1/surfaces/1/id", "L1.S1")  # a new error
    assert canonical(ed.system) == state and ed.history == history
    assert ed.can_undo and not ed.can_redo


def test_a_failing_undo_or_redo_leaves_state_and_stacks() -> None:
    # H4 of the second review of ADR 0024. The inverse of an accepted patch does not fail, so a
    # stored entry is broken on purpose (white box: the stacks hold (patch, inverse) pairs).
    broken = [{"op": "remove", "path": "/root/children/99"}]
    ed = rt.Editor(load("m0/singlet.rtt.json"))
    ed.set("/name", "a")
    ed.set("/name", "b")
    patch, _ = ed._undo[-1]  # noqa: SLF001
    ed._undo[-1] = (patch, broken)  # noqa: SLF001
    diagnostics = [d.code for d in ed.diagnostics]
    with pytest.raises(rt.EditError):
        ed.undo()
    assert ed.system.name == "b" and len(ed.history) == 2 and not ed.can_redo
    assert ed._undo[-1] == (patch, broken) and len(ed._undo) == 2  # noqa: SLF001
    assert [d.code for d in ed.diagnostics] == diagnostics

    ed = rt.Editor(load("m0/singlet.rtt.json"))
    ed.set("/name", "a")
    ed.undo()
    _, inverse = ed._redo[-1]  # noqa: SLF001
    ed._redo[-1] = (broken, inverse)  # noqa: SLF001
    with pytest.raises(rt.EditError):
        ed.redo()
    assert ed.system.name != "a" and ed.history == [] and ed.can_redo
    assert ed._redo == [(broken, inverse)]  # noqa: SLF001


def test_an_empty_system_can_be_built() -> None:
    ed = rt.Editor(rt.System())
    assert any(d.severity == rt.Severity.ERROR for d in ed.diagnostics)
    ed.insert("/wavelengths/-", {"um": 0.5876, "reference": True})
    ed.set("/root/name", "system")
    assert ed.system.wavelengths[0].um == 0.5876
    assert ed.can_undo


def test_editor_works_on_its_own_copy() -> None:
    s = load("m0/singlet.rtt.json")
    ed = rt.Editor(s)
    s.name = "outside"
    assert ed.system.name != "outside"
    copy = ed.system
    copy.name = "copy"
    assert ed.system.name != "copy"


# ------------------------------------------------------------- history replay -----


def test_export_and_replay() -> None:
    s = load("m2/cooke_triplet.rtt.json")
    ed = rt.Editor(s)
    ed.set("/name", "edited")
    ed.remove("/fields/points/2")
    ed.undo()
    ed.insert("/fields/points/0", {"y": 1.0})
    data = ed.export_history()
    assert json.loads(json.dumps(data)) == data  # plain JSON
    assert data["base"] == canonical(s)
    assert data["patches"] == ed.history
    again = rt.Editor.from_history(data)
    assert canonical(again.system) == canonical(ed.system)
    assert again.history == ed.history
    # Replaying the patches on the base system as script commands gives the same model.
    system = rt.System.from_json(data["base"])
    for p in data["patches"]:
        system = rt.apply_patch(system, p)
    assert canonical(system) == canonical(ed.system)


def test_replay_uses_the_full_check() -> None:
    # H5: a changed history whose patch adds an error is rejected on replay.
    data: dict[str, Any] = rt.Editor(load("m0/singlet.rtt.json")).export_history()
    data["patches"] = [[{"op": "replace", "path": "/root/children/1/surfaces/1/id",
                         "value": "L1.S1"}]]
    with pytest.raises(rt.EditError) as info:
        rt.Editor.from_history(data)
    assert info.value.code == "surface.id_duplicate"


def test_replay_rejects_data_that_is_no_history() -> None:
    for data in ({}, {"base": "not json", "patches": []}, {"base": "{}", "patches": []},
                 {"base": rt.System().to_json()}):
        with pytest.raises(ValueError, match="not an exported Editor history"):
            rt.Editor.from_history(data)


def test_replay_rejects_another_format_version() -> None:
    data = rt.Editor(load("m0/singlet.rtt.json")).export_history()
    major, minor, patch = rt.System().schema_version.split(".")
    other = f"{major}.{int(minor) + 1}.0"
    data["base"] = data["base"].replace(f'"schema_version": "{rt.System().schema_version}"',
                                        f'"schema_version": "{other}"')
    with pytest.raises(ValueError, match="schema"):
        rt.Editor.from_history(data)
    # A new patch version is a compatible addition and is allowed.
    data2 = rt.Editor(load("m0/singlet.rtt.json")).export_history()
    data2["base"] = data2["base"].replace(f'"schema_version": "{rt.System().schema_version}"',
                                          f'"schema_version": "{major}.{minor}.{int(patch) + 1}"')
    assert rt.Editor.from_history(data2).history == []
