"""The edit form in Python: System.to_dict, json_at and locate_* (ADR 0024 points 2 and 7)."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema  # type: ignore[import-untyped]  # no py.typed; test tool (ADR 0018)
import pytest
from conftest import REFERENCE_DIR, REPO_ROOT, reference_files

import raytatouille as rt

SCHEMA = json.loads((REPO_ROOT / "schema" / "raytatouille.schema.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_the_edit_form_is_valid_against_the_schema(file: Path) -> None:
    # ADR 0024 point 7: the edit form is a second spelling of the same schema version.
    data = rt.load(file).to_dict()
    jsonschema.Draft202012Validator(SCHEMA).validate(data)
    assert rt.System.from_json(json.dumps(data)) == rt.load(file)


def test_to_dict_keeps_the_writer_order() -> None:
    data = rt.load(REFERENCE_DIR / "m0" / "singlet.rtt.json").to_dict()
    assert list(data)[:3] == ["schema_version", "name", "units"]
    assert data["aperture"]["value"] == {"value": 20.0, "variable": False}


def test_json_at() -> None:
    s = rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")
    data = s.to_dict()
    assert s.json_at("") == data
    assert s.json_at("/root/children/2/surfaces/0/id") == "A.S1"
    assert s.json_at("/object/distance") == {"value": 250.0, "variable": True}
    assert s.json_at("/wavelengths/1/um") == 1.064
    # RFC 6901: "~1" is "/", "~0" is "~" (no such key here: KeyError).
    with pytest.raises(KeyError):
        s.json_at("/environment/a~1b")
    with pytest.raises(KeyError):
        s.json_at("/wavelengths/9")
    with pytest.raises(KeyError):
        s.json_at("/name/x")  # a string has no members
    for bad in ("name", "/wavelengths/01", "/wavelengths/-", "/a~2"):
        with pytest.raises(ValueError):
            s.json_at(bad)


def test_locate() -> None:
    tour = rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")
    assert tour.locate_surface("QWP") == "/root/children/0/children/1/surfaces/0"
    assert tour.locate_node("quarter wave plate") == "/root/children/0/children/1"
    assert tour.locate_node("system") == "/root"
    assert tour.locate_path("ghost G.S2-A.S2") == "/paths/1"
    assert tour.locate_surface("nope") is None
    assert tour.locate_node("nope") is None and tour.locate_path("nope") is None
    # Element and surface may share a name: separate name spaces (ADR 0024 point 2).
    mirror = rt.load(REFERENCE_DIR / "m1" / "paraboloid_mirror.rtt.json")
    node, surface = mirror.locate_node("M"), mirror.locate_surface("M")
    assert node is not None and surface is not None and node != surface
    assert mirror.json_at(surface + "/id") == "M"
    assert mirror.json_at(node + "/name") == "M"


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_locate_matches_every_surface_and_node(file: Path) -> None:
    s = rt.load(file)
    data = s.to_dict()

    def walk(node: dict[str, object], pointer: str) -> None:
        assert s.locate_node(str(node["name"])) == pointer
        if node["type"] == "assembly":
            children = node["children"]
            assert isinstance(children, list)
            for i, child in enumerate(children):
                walk(child, f"{pointer}/children/{i}")
        else:
            surfaces = node["surfaces"]
            assert isinstance(surfaces, list)
            for k, surface in enumerate(surfaces):
                assert s.locate_surface(surface["id"]) == f"{pointer}/surfaces/{k}"

    walk(data["root"], "/root")
    for i, path in enumerate(data["paths"]):
        assert s.locate_path(path["name"]) == f"/paths/{i}"
