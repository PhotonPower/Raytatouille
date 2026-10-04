"""Compiling systems (rt.compile, MaterialLibrary) and first-order data (rt.paraxial)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

import raytatouille as rt


def singlet(reference_dir: Path, medium: str = "VACUUM") -> rt.CompiledSystem:
    system = rt.load(reference_dir / "m1" / "singlet_const.rtt.json")
    system.environment.medium = medium
    return rt.compile(system)


def test_compiled_system_properties(reference_dir: Path) -> None:
    cs = singlet(reference_dir, "AIR")
    assert cs.wavelengths_um == [0.4861, 0.5876, 0.6563]
    assert cs.reference_wavelength == 1
    assert cs.wavelength_weights == [1.0, 1.0, 1.0]
    assert cs.temperature_c == 20.0
    assert cs.field_count == 3
    assert cs.path_names == ["main"]
    assert cs.find_path("main") == 0
    assert cs.find_path("nope") is None
    assert cs.surface_ids[0] == "STO"
    assert cs.media[0].reference == "AIR"
    # Ciddor air (#25): n - 1 about 2.7e-4 in the visible.
    assert all(2.6e-4 < n.real - 1.0 < 2.9e-4 for n in cs.media[0].index)
    glass = [m for m in cs.media if m.reference.startswith("CONST:")]
    assert [n.real for n in glass[0].index] == [1.5168] * 3


def test_material_library_index(catalog_dir: Path) -> None:
    lib = rt.MaterialLibrary()
    assert lib.index("VACUUM", 0.5876) == 1.0
    assert lib.index("CONST:1.5,0.01", 0.5876) == complex(1.5, 0.01)
    lib.add_catalog(catalog_dir / "schott.agf")
    # N-BK7 relative to air at the d line (587.56 nm in air) is 1.5168 (SCHOTT data sheet).
    n_air = lib.index("AIR", 0.58756).real
    n_rel = lib.index("SCHOTT:N-BK7", 0.58756 * n_air).real / n_air
    assert abs(n_rel - 1.5168) < 5e-6


def test_compile_with_a_catalogue(reference_dir: Path, catalog_dir: Path) -> None:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    cs = rt.compile(rt.load(reference_dir / "m2" / "achromat.rtt.json"), lib)
    assert {m.reference for m in cs.media} == {"AIR", "SCHOTT:N-BK7", "SCHOTT:F2"}


def test_first_order_of_the_reference_singlet(reference_dir: Path) -> None:
    # Plano-convex lens, R1 = 51.68 mm, n = 1.5168, plane rear side, in VACUUM:
    # phi = (n - 1) / R1, EFL = R1 / (n - 1) = 100 mm (Greivenkamp, OPTI-201/202, Sec. 9,
    # p. 9-2; the plane surface has no power). BFL = EFL - d / n with d = 4 mm (lens from
    # z = 5 to z = 9), since H' lies at the curved vertex shifted by -d/n (p. 9-12).
    cs = singlet(reference_dir)
    fo = rt.paraxial.first_order(cs, path="main")
    assert fo.efl is not None and fo.bfl is not None
    assert fo.efl == pytest.approx(100.0, rel=1e-12)
    assert fo.bfl == pytest.approx(100.0 - 4.0 / 1.5168, rel=1e-12)
    assert fo.object_index == 1.0
    assert fo.image_direction == 1
    assert fo.entrance_pupil is not None
    assert fo.entrance_pupil.diameter == pytest.approx(20.0, rel=1e-12)
    assert fo.lateral_magnification is None  # object at infinity


def test_first_order_path_and_wavelength_arguments(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    by_name = rt.paraxial.first_order(cs, path="main", wavelength=None)
    by_index = rt.paraxial.first_order(cs, 0, cs.reference_wavelength)
    assert by_name.efl == by_index.efl
    assert rt.paraxial.first_order(cs).efl == by_name.efl


def test_first_order_in_air_is_absolute(reference_dir: Path) -> None:
    # Since #25 the indices are absolute: in AIR, EFL = R1 / (n - n_air) (docs/architecture.md,
    # Engine 1).
    cs = singlet(reference_dir, "AIR")
    n_air = cs.media[0].index[cs.reference_wavelength].real
    fo = rt.paraxial.first_order(cs)
    assert fo.efl == pytest.approx(51.68 / (1.5168 - n_air), rel=1e-12)
    assert fo.rear_focal_length == pytest.approx(n_air * 51.68 / (1.5168 - n_air), rel=1e-12)
    assert not math.isclose(fo.efl or 0.0, 100.0, rel_tol=1e-5)
