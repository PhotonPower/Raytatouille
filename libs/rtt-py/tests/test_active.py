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


def test_inactive_nodes_in_tolerancing_and_variables() -> None:
    """PR 4 of ADR 0049: an explicit item at an inactive node is left out with
    tolerances.target_inactive (ToleranceItems.warnings), class items there are skipped with
    ClassSkip.INACTIVE, and a variable Param there is no variable."""
    from raytatouille import optim, tolerance

    data = json.loads(
        rt.load(REFERENCE_DIR / "m7" / "decentred_lens_finite.rtt.json").to_json()
    )
    data["root"]["children"][1]["active"] = False  # the lens
    data["tolerances"]["class"] = "normal"
    system = rt.System.from_json(json.dumps(data))
    items = tolerance.expand(system)
    assert [w.code for w in items.warnings] == ["tolerances.target_inactive"]
    assert items.warnings[0].location == "/tolerances/items/0"
    assert "lens-dy" not in [i.id for i in items.items]
    at_lens = [c for c in items.classes.items if c.node == "lens"]
    assert at_lens and all(c.skip == tolerance.ClassSkip.INACTIVE for c in at_lens)

    z = {"value": 0.0, "variable": True}
    data["root"]["children"][1]["pose"] = {"position": [0.0, 0.0, z]}
    with_variable = rt.System.from_json(json.dumps(data))
    assert optim.variables(with_variable) == []


def test_layout_reports_and_result_format_carry_active() -> None:
    """PR 5 of ADR 0049: SurfaceLayout.active, CompiledElement.active and the column active of
    the dimension report; the result format 0.1.20 writes them, older files read as true."""
    from raytatouille import analysis, layout, results

    data = json.loads(singlet().to_json())
    data["root"]["children"][1]["active"] = False  # the lens L1
    cs = rt.compile(rt.System.from_json(json.dumps(data)))
    surfaces = layout.surfaces(cs)
    assert [s.active for s in surfaces] == [True, False, False, True]
    elements = layout.elements(cs)
    assert [e.active for e in elements] == [True, False, True]
    segments = analysis.dimension_report(cs).segments
    assert list(segments.active) == [False]  # the lens is still listed
    on = analysis.dimension_report(rt.compile(singlet())).segments
    assert list(on.active) == [True]

    text = surfaces[1].to_json()
    assert json.loads(text)["schema_version"] == "0.1.20"
    assert json.loads(text)["data"]["active"] is False
    assert json.loads(results.to_json(elements[1]))["data"]["active"] is False
    # An older file has no key active: it reads as true.
    old = json.loads(text)
    old["schema_version"] = "0.1.19"
    del old["data"]["active"]
    assert results.load_json(json.dumps(old)).data["active"] is True
    report = json.loads(analysis.dimension_report(cs).to_json())
    report["schema_version"] = "0.1.19"
    del report["data"]["segments"]["active"]
    loaded = results.load_json(json.dumps(report)).data["segments"]
    assert list(loaded["active"]) == [True]
