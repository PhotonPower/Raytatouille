"""Compiling systems (rt.compile, MaterialLibrary) and first-order data (rt.paraxial)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
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
    # z = 5 to z = 9): H' lies d/n before the plane rear vertex (p. 9-12, BFD = -y_k / u'_k
    # with y_k = 1 - (n - 1) d / (n R1) and u'_k = -1 / EFL).
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


def test_prescription_of_the_reference_singlet(reference_dir: Path) -> None:
    # Values by hand as in libs/rtt-paraxial/tests/test_prescription.cpp: stop (r = 10) at
    # z = 0, plano-convex lens n = 1.5168, R1 = 51.68 from z = 5 to 9, image at z = 106.363,
    # EPD 20, maximum field 5 deg. i = u + y c (Sasian, OPTI 517 L4, p. 24); F/# = EFL / EPD =
    # 1 / (2 n' |u'|) and NA = n' |u'| (Greivenkamp, OPTI-502 Sec. 9, p. 9-34, 9-35); Lagrange
    # invariant n (u_bar y - u y_bar) = 10 tan 5 deg in the stop (p. 9-41).
    cs = singlet(reference_dir)
    p = rt.paraxial.prescription(cs, path="main")
    q = p.surfaces
    t = math.tan(math.radians(5.0))
    c = 1.0 / 51.68
    assert len(q) == 4
    assert q.surface.tolist() == [0, 1, 2, 3]
    assert q.z == pytest.approx([0.0, 5.0, 9.0, 106.363], rel=1e-12, abs=1e-12)
    assert q.n.tolist() == [1.0, 1.5168, 1.0, 1.0]
    assert q.y[1] == 10.0 and q.u[0] == 0.0
    assert q.i[1] == pytest.approx(10.0 * c, rel=1e-12)
    assert q.i_bar[1] == pytest.approx(t + 5.0 * t * c, rel=1e-12)
    assert p.total_track == pytest.approx(106.363, rel=1e-12)
    assert p.object_distance is None
    assert p.paraxial_working_f_number == pytest.approx(5.0, rel=1e-12)
    assert p.paraxial_image_na == pytest.approx(0.1, rel=1e-12)
    assert p.lagrange_invariant == pytest.approx(10.0 * t, rel=1e-12)
    assert q.lagrange == pytest.approx(np.full(4, 10.0 * t), rel=1e-12)
    assert p.marginal_start is not None and p.marginal_start.y == 10.0
    assert p.chief_start is not None and p.chief_start.u == pytest.approx(t, rel=1e-12)
    assert p.first_order.efl == pytest.approx(100.0, rel=1e-12)
    # The same rays as seidel() (one internal construction, #84).
    s = rt.paraxial.seidel(cs)
    assert q.y.tobytes() == s.surfaces.y.tobytes()
    assert q.y_bar.tobytes() == s.surfaces.y_bar.tobytes()
    assert not q.y.flags.writeable


def test_prescription_without_a_stop_has_nan_and_none(reference_dir: Path) -> None:
    # Paraboloid R = -200 without a stop, EPD 60 (tests/reference/m1/paraboloid_mirror):
    # marginal ray only; F/# = EFL / EPD = 100 / 60, n' u' = -30 (n' - n) c with c = -1/200.
    system = rt.load(reference_dir / "m1" / "paraboloid_mirror.rtt.json")
    system.environment.medium = "VACUUM"
    p = rt.paraxial.prescription(system)  # a System is compiled first
    q = p.surfaces
    assert q.n.tolist() == [-1.0, -1.0]
    assert not np.isnan(q.y).any()
    assert np.isnan(q.y_bar).all() and np.isnan(q.u_bar).all() and np.isnan(q.i_bar).all()
    assert np.isnan(q.lagrange).all()
    assert p.chief_start is None and p.lagrange_invariant is None
    assert p.total_track == pytest.approx(100.0, rel=1e-12)
    assert p.paraxial_working_f_number == pytest.approx(100.0 / 60.0, rel=1e-12)
    assert p.paraxial_image_na == pytest.approx(0.3, rel=1e-12)
    with pytest.raises(rt.ParaxialError):
        rt.paraxial.prescription(system, wavelength=7)


def test_first_order_in_air_is_absolute(reference_dir: Path) -> None:
    # Since #25 the indices are absolute: in AIR, EFL = R1 / (n - n_air) (docs/architecture.md,
    # Engine 1).
    cs = singlet(reference_dir, "AIR")
    n_air = cs.media[0].index[cs.reference_wavelength].real
    fo = rt.paraxial.first_order(cs)
    assert fo.efl == pytest.approx(51.68 / (1.5168 - n_air), rel=1e-12)
    assert fo.rear_focal_length == pytest.approx(n_air * 51.68 / (1.5168 - n_air), rel=1e-12)
    assert not math.isclose(fo.efl or 0.0, 100.0, rel_tol=1e-5)
