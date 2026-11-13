"""The JSON schemas shipped with the package (#386): equal to the files under schema/."""

from __future__ import annotations

import json

import pytest
from conftest import REPO_ROOT

import raytatouille as rt


@pytest.mark.parametrize(
    ("name", "file"),
    [
        ("system", "raytatouille.schema.json"),
        ("result", "raytatouille-result.schema.json"),
        ("coatings", "raytatouille-coatings.schema.json"),
    ],
)
def test_schema_equals_the_file_in_the_repository(name: str, file: str) -> None:
    expected = json.loads((REPO_ROOT / "schema" / file).read_text(encoding="utf-8"))
    assert rt.schema(name) == expected
    assert expected["$schema"]  # content guard: a real schema, not an empty dict


def test_default_is_the_system_schema() -> None:
    assert rt.schema() == rt.schema("system")
    assert "root" in rt.schema()["properties"]


def test_every_call_gives_a_new_dict() -> None:
    first = rt.schema()
    first["title"] = "changed"
    assert rt.schema()["title"] != "changed"


def test_unknown_name() -> None:
    with pytest.raises(ValueError, match="system, result, coatings"):
        rt.schema("model")
