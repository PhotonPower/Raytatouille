"""The flag active of assemblies and elements (#380, ADR 0049, schema 0.10): read from the
file, shown on the typed tree, switched with the Editor and taken back with undo."""

from __future__ import annotations

import json

import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import model


def singlet() -> rt.System:
    return rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")


def test_active_is_true_by_default_and_read_from_the_file() -> None:
    system = singlet()
    children = system.root.children
    assert system.root.active and all(c.active for c in children)
    data = json.loads(system.to_json())
    assert "active" not in data["root"]["children"][1]  # canonical form: only if false
    data["root"]["children"][1]["active"] = False
    switched = rt.System.from_json(json.dumps(data))
    lens = switched.root.children[1]
    assert isinstance(lens, model.Element) and not lens.active
    assert json.loads(switched.to_json())["root"]["children"][1]["active"] is False


def test_editor_switches_an_element_off_and_back() -> None:
    system = singlet()
    before = system.to_json()
    path = system.locate_node("L1")
    assert path == "/root/children/1"
    assert system.to_dict()["root"]["children"][1]["active"] is True  # edit form: always
    editor = rt.Editor(system)
    editor.set(path + "/active", False)
    lens = editor.system.root.children[1]
    assert isinstance(lens, model.Element) and lens.name == "L1" and not lens.active
    editor.undo()
    assert editor.system.to_json() == before


def test_root_cannot_be_inactive_and_old_files_cannot_hold_the_key() -> None:
    data = json.loads(singlet().to_json())
    data["root"]["active"] = False
    codes = [d.code for d in rt.validate(rt.System.from_json(json.dumps(data)))]
    assert "node.root_inactive" in codes
    old = json.loads(singlet().to_json())
    old["schema_version"] = "0.9.0"
    old["root"]["children"][1]["active"] = False
    with pytest.raises(ValueError, match="0.10"):
        rt.System.from_json(json.dumps(old))
