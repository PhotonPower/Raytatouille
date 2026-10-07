"""Stable diagnostic codes, error places and the warning channel (ADR 0022, #86)."""

from __future__ import annotations

import re
import warnings
from collections.abc import Callable
from pathlib import Path

import pytest

import raytatouille as rt


def variant(reference_dir: Path, tmp_path: Path, old: str, new: str) -> rt.System:
    """m1/singlet_const with one text replacement, loaded from a temporary file."""
    text = (reference_dir / "m1" / "singlet_const.rtt.json").read_text(encoding="utf-8")
    assert text.count(old) == 1, old
    path = tmp_path / "variant.rtt.json"
    path.write_text(text.replace(old, new), encoding="utf-8")
    return rt.load(path)


CONIC = '"base": {"type": "conic", "radius": 51.68}'
BARE_ASPHERE = '"base": {"type": "even_asphere", "radius": 51.68, "coefficients": []}'
FIELDS = '{"y": 5.0}'


def test_diagnostics_have_codes(reference_dir: Path, tmp_path: Path) -> None:
    system = variant(reference_dir, tmp_path, '{"um": 0.4861}',
                     '{"um": 0.4861, "reference": true}')
    [d] = rt.validate(system)
    assert d.code == "wavelengths.reference_count"
    assert d.severity == rt.Severity.ERROR and d.location == "/wavelengths"
    assert str(d).startswith("error [wavelengths.reference_count] /wavelengths: ")
    with pytest.raises(rt.CompileError) as info:
        rt.compile(system)
    assert info.value.codes == ["wavelengths.reference_count"]
    assert [x.code for x in info.value.diagnostics] == info.value.codes


def test_registry_lists_every_code_with_its_severity() -> None:
    codes = rt.diagnostics.CODES
    assert list(codes) == sorted(codes)
    info = codes["material.unknown"]
    assert info.code == "material.unknown" and info.severity == rt.Severity.ERROR
    assert codes["aperture.na_not_physical"].severity == rt.Severity.WARNING
    assert all(re.fullmatch(r"[a-z0-9_]+\.[a-z0-9_]+", c) for c in codes)


def test_registry_and_documentation_agree(reference_dir: Path) -> None:
    """docs/diagnostics.md lists exactly the registered codes, each with its severity and its
    producer."""
    doc = (reference_dir.parent.parent / "docs" / "diagnostics.md").read_text(encoding="utf-8")
    rows = re.findall(r"^\| `([a-z0-9_.]+)` \| (Fehler|Warnung) \| [^|]* \| ([a-z]+) \|", doc,
                      flags=re.MULTILINE)
    severity = {"Fehler": rt.Severity.ERROR, "Warnung": rt.Severity.WARNING}
    documented = {code: (severity[s], producer) for code, s, producer in rows}
    assert len(documented) == len(rows), "a code is documented twice"
    assert documented == {c: (i.severity, i.producer) for c, i in rt.diagnostics.CODES.items()}


def test_compile_reports_warnings_as_data_and_as_python_warnings(
        reference_dir: Path, tmp_path: Path) -> None:
    system = variant(reference_dir, tmp_path, CONIC, BARE_ASPHERE)
    with pytest.warns(rt.errors.RaytatouilleWarning) as caught:
        cs = rt.compile(system)
    [w] = caught
    assert w.filename == __file__  # points at the line that called compile()
    warning = w.message
    assert isinstance(warning, rt.errors.RaytatouilleWarning)
    assert isinstance(warning, UserWarning)
    assert warning.code == "shape.asphere_without_coefficients"
    assert warning.location == "/root/children/1/surfaces/0/shape/base/coefficients"
    assert str(warning).startswith("warning [shape.asphere_without_coefficients] ")
    [d] = cs.diagnostics
    assert (d.code, d.severity, d.location) == (warning.code, rt.Severity.WARNING,
                                                warning.location)
    # The unchanged system compiles without warnings.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json")).diagnostics == []


def test_surface_locations(reference_dir: Path) -> None:
    cs = rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))
    assert cs.surface_ids == ["STO", "L1.S1", "L1.S2", "IMG"]
    assert cs.surface_locations == ["/root/children/0/surfaces/0", "/root/children/1/surfaces/0",
                                    "/root/children/1/surfaces/1", "/root/children/2/surfaces/0"]


@pytest.mark.parametrize("file", ["m1/paraboloid_mirror.rtt.json", "m0/michelson.rtt.json"])
def test_no_stop_is_one_error_class(reference_dir: Path, file: str) -> None:
    cs = rt.compile(rt.load(reference_dir / file))
    calls: list[Callable[[], object]] = [
        lambda: rt.analysis.spot(cs, path=0, field=0),
        lambda: rt.analysis.ray_fan(cs, path=0, field=0),
        lambda: rt.analysis.opd_map(cs, path=0, field=0),
        lambda: rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=2), path=0),
    ]
    for call in calls:
        with pytest.raises(rt.NoStopError) as info:
            call()
        e = info.value
        # Existing handlers keep working: it is all of these at once.
        assert isinstance(e, (rt.ParaxialError))
        assert isinstance(e, rt.AnalysisError) and isinstance(e, ValueError)
        assert e.path_name == cs.path_names[0] and e.location == "/paths/0"


def test_paraxial_error_names_the_surface(reference_dir: Path) -> None:
    cs = rt.compile(rt.load(reference_dir / "m0" / "michelson.rtt.json"))
    with pytest.raises(rt.ParaxialError) as info:
        rt.paraxial.first_order(cs, path=0)
    e = info.value
    assert e.surface == "BS"
    assert e.location == cs.surface_locations[cs.surface_ids.index("BS")]


def test_paraxial_error_points_at_the_field(reference_dir: Path, tmp_path: Path) -> None:
    cs = rt.compile(variant(reference_dir, tmp_path, FIELDS, FIELDS + ', {"y": 95.0}'))
    with pytest.raises(rt.ParaxialError) as info:
        rt.paraxial.seidel(cs, path=0)
    assert info.value.surface is None and info.value.location == "/fields/points/3"


def test_analysis_error_describes_the_lost_ray(reference_dir: Path, tmp_path: Path) -> None:
    # 80 degree: the chief ray passes the stop and misses the lens (as in test_error_context.cpp).
    cs = rt.compile(variant(reference_dir, tmp_path, FIELDS, FIELDS + ', {"y": 80.0}'))
    with pytest.raises(rt.AnalysisError) as info:
        rt.analysis.spot(cs, path=0, field=3)
    e = info.value
    assert not isinstance(e, rt.NoStopError)
    assert e.surface == "STO" and e.location == "/root/children/0/surfaces/0"
    assert e.ray_status == rt.trace.RayStatus.MISSED
    assert e.field == 3 and e.wavelength == cs.reference_wavelength
