"""Loading, saving and validating systems (rt.load, rt.save, rt.validate)."""

from __future__ import annotations

from pathlib import Path

import pytest
from conftest import reference_files

import raytatouille as rt


def test_there_are_reference_files() -> None:
    assert len(reference_files()) >= 3


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_file_python_file_is_byte_identical(file: Path, tmp_path: Path) -> None:
    # Reference case of #32: file -> Python -> file is bitwise the same.
    system = rt.load(file)
    out = tmp_path / file.name
    rt.save(system, out)
    assert out.read_bytes() == file.read_bytes()
    system.save(tmp_path / "method.rtt.json")
    assert (tmp_path / "method.rtt.json").read_bytes() == file.read_bytes()
    assert system.to_json().encode("utf-8") == file.read_bytes()


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_json_text_round_trip(file: Path) -> None:
    system = rt.load(file)
    again = rt.System.from_json(system.to_json())
    assert again == system


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_reference_files_are_valid(file: Path) -> None:
    errors = [d for d in rt.validate(rt.load(file)) if d.severity == rt.Severity.ERROR]
    assert errors == []


def test_load_accepts_str_and_path(reference_dir: Path) -> None:
    file = reference_dir / "m1" / "singlet_const.rtt.json"
    assert rt.load(str(file)) == rt.load(file)


def test_environment_can_be_changed(reference_dir: Path) -> None:
    system = rt.load(reference_dir / "m1" / "singlet_const.rtt.json")
    assert system.environment.medium == "AIR"
    assert system.environment.temperature_c == 20.0
    assert system.environment.pressure_atm == 1.0
    system.environment.medium = "VACUUM"
    system.environment.temperature_c = 25.0
    assert system.environment.medium == "VACUUM"
    assert '"medium": "VACUUM"' in system.to_json()
    system.name = "renamed"
    assert rt.System.from_json(system.to_json()).name == "renamed"


def test_wavelengths_are_read_only_copies(reference_dir: Path) -> None:
    system = rt.load(reference_dir / "m1" / "singlet_const.rtt.json")
    wavelengths = system.wavelengths
    assert [w.um for w in wavelengths] == [0.4861, 0.5876, 0.6563]
    assert [w.reference for w in wavelengths] == [False, True, False]
    assert system.schema_version == "0.2.0"


def test_validate_reports_errors_without_raising() -> None:
    system = rt.System()
    diagnostics = rt.validate(system)
    assert any(d.severity == rt.Severity.ERROR for d in diagnostics)
    for d in diagnostics:
        assert d.location.startswith("/")
        assert d.message
        assert str(d).startswith("error" if d.severity == rt.Severity.ERROR else "warning")
