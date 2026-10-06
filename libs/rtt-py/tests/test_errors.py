"""C++ errors arrive as Python exceptions with the C++ message (raytatouille.errors)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import raytatouille as rt


def test_parse_error_has_the_json_pointer(reference_dir: Path) -> None:
    data = json.loads((reference_dir / "m1" / "singlet_const.rtt.json").read_text())
    data["wavelengths"] = 3
    with pytest.raises(rt.ParseError) as info:
        rt.System.from_json(json.dumps(data))
    error = info.value
    assert isinstance(error, rt.RaytatouilleError) and isinstance(error, ValueError)
    assert error.pointer == "/wavelengths"
    assert str(error).startswith("/wavelengths: ")


def test_duplicate_key_is_a_parse_error(reference_dir: Path) -> None:
    # json.dumps cannot write a duplicate key, so it is inserted into the text (#68).
    text = (reference_dir / "m1" / "singlet_const.rtt.json").read_text()
    head, sep, tail = text.partition('"schema_version"')
    assert sep
    with pytest.raises(rt.ParseError) as info:
        rt.System.from_json(head + '"schema_version": "9.9.9", ' + sep + tail)
    assert info.value.pointer == "/schema_version"
    assert str(info.value) == "/schema_version: duplicate key 'schema_version'"


def test_missing_file_is_an_os_error(tmp_path: Path) -> None:
    with pytest.raises(OSError):
        rt.load(tmp_path / "missing.rtt.json")
    with pytest.raises(OSError):
        rt.save(rt.System(), tmp_path / "no" / "such" / "dir" / "x.rtt.json")


def test_compile_error_has_the_diagnostics(reference_dir: Path) -> None:
    system = rt.load(reference_dir / "m1" / "singlet_const.rtt.json")
    system.environment.medium = "NOPE:GLASS"
    with pytest.raises(rt.CompileError) as info:
        rt.compile(system)
    error = info.value
    assert isinstance(error, rt.RaytatouilleError) and isinstance(error, ValueError)
    assert [d.location for d in error.diagnostics] == ["/environment/medium"]
    assert error.diagnostics[0].severity == rt.Severity.ERROR
    assert "NOPE:GLASS" in str(error)


def test_unknown_material() -> None:
    with pytest.raises(rt.UnknownMaterial) as info:
        rt.MaterialLibrary().index("NOPE:GLASS", 0.5876)
    error = info.value
    assert isinstance(error, rt.RaytatouilleError) and isinstance(error, KeyError)
    assert "NOPE:GLASS" in str(error)
    assert not str(error).startswith("'")


def test_agf_error_has_file_and_line(tmp_path: Path) -> None:
    catalog = tmp_path / "broken.agf"
    catalog.write_text("CC broken test catalogue\nNM GLASS 2 0 1.5 60\nCD 1.0 x\n")
    with pytest.raises(rt.AgfError) as info:
        rt.MaterialLibrary().add_catalog(catalog)
    error = info.value
    assert isinstance(error, rt.RaytatouilleError) and isinstance(error, ValueError)
    assert Path(error.file).name == "broken.agf"
    assert error.line == 3


def test_paraxial_error(reference_dir: Path) -> None:
    cs = rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))
    with pytest.raises(rt.ParaxialError):
        rt.paraxial.first_order(cs, path=5)
    with pytest.raises(rt.ParaxialError):
        rt.paraxial.first_order(cs, wavelength=9)
    with pytest.raises(ValueError, match="unknown path 'nope'"):
        rt.paraxial.first_order(cs, path="nope")


def test_all_errors_derive_from_the_base_class() -> None:
    for cls in (rt.ParseError, rt.CompileError, rt.ParaxialError, rt.UnknownMaterial,
                rt.AgfError):
        assert issubclass(cls, rt.RaytatouilleError)
