"""Pupil function, refocus, OPD at pupil points and the Zernike fit from Python (#209, P1a;
ADR 0032). The C++ results themselves are tested in rtt-analysis; here the bindings: columns,
options, warnings, errors with location and the agreement of the Python calls with each other.
Python = C++ bit for bit is test_bitwise_pupil.py."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille import analysis as an

MAPPING_WARNING = r"ignore:warning \[pupil\.mapping_nonlinear\]:raytatouille.errors.RaytatouilleWarning"


def paraboloid(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m2" / "paraboloid_stop.rtt.json"))


def singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def test_pupil_cells_columns(reference_dir: Path) -> None:
    c = an.pupil_cells(8)
    assert c.grid == 8 and len(c) == 64
    for column in (c.fraction, c.px, c.py):
        assert column.shape == (64,) and column.dtype == np.float64
        assert not column.flags.writeable
    a = 2.0 / 8
    assert abs(float(np.sum(c.fraction)) * a * a - math.pi) <= 1e-13 * math.pi
    with pytest.raises(ValueError):
        an.pupil_cells(1)
    with pytest.raises(ValueError):
        an.pupil_cells(an.MAX_PUPIL_GRID + 1)


def test_opd_points_equals_opd_map_on_the_grid(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    m = an.opd_map(cs, field=1, grid=9)
    p = an.opd_points(cs, field=1, rays=rt.trace.GridPupil(9))
    assert p.points.w.tobytes() == m.points.w.tobytes()
    assert p.points.status.tobytes() == m.points.status.tobytes()
    assert np.any(p.points.w != 0.0)  # content guard
    # The shorthand string works too, and focus_shift moves the sphere along the chief ray.
    q = an.opd_points(cs, field=1, rays="grid:9", focus_shift=0.5)
    assert q.points.w.tobytes() != p.points.w.tobytes()
    assert abs(float(q.sphere.centre[2] - p.sphere.centre[2]) - 0.5) < 1e-3
    with pytest.raises(ValueError):
        an.opd_map(cs, focus_shift=math.nan)


@pytest.mark.filterwarnings(MAPPING_WARNING)
def test_pupil_function_columns_and_mapping(reference_dir: Path) -> None:
    pf = an.pupil_function(paraboloid(reference_dir), grid=16, with_prt=True)
    n2 = 16 * 16
    for name in ("fraction", "px", "py", "amplitude", "w"):
        column = getattr(pf, name)
        assert column.shape == (n2,) and column.dtype == np.float64, name
        assert not column.flags.writeable, name
    assert pf.status.shape == (n2,) and pf.status.dtype == np.uint8
    assert pf.prt.shape == (n2, 3, 3) and pf.prt.dtype == np.complex128
    assert pf.grid == 16 and pf.path == 0 and pf.field == 0 and pf.wavelength == 0
    assert pf.wavelength_um == pf.reference_wavelength_um
    # Cells with f = 0 are not traced: VIGNETTED with A = W = 0 (check fraction first).
    outside = pf.fraction == 0.0
    assert np.all(pf.status[outside] == int(rt.trace.RayStatus.VIGNETTED))
    assert np.all(pf.amplitude[outside] == 0.0)
    arrived = (pf.fraction > 0.0) & (pf.status == int(rt.trace.RayStatus.ALIVE))
    assert np.all(np.abs(pf.w[arrived]) < 1e-6)  # the paraboloid on the axis
    assert pf.arrived_fraction == 1.0
    mapping = pf.mapping
    assert mapping is not None
    assert mapping.m.shape == (2, 2)
    assert mapping.s0.tolist() == [0.0, 0.0]
    assert mapping.centre.shape == (3,) and mapping.d_c.tolist() == pytest.approx([0, 0, -1])
    assert mapping.delta_max > 0.01  # F/1.7: above the default threshold
    assert [w.code for w in pf.warnings] == ["pupil.mapping_nonlinear"]
    without = an.pupil_function(paraboloid(reference_dir), grid=16)
    assert without.prt.shape == (0, 3, 3)


def test_pupil_function_warns_and_refocus_agrees(reference_dir: Path) -> None:
    cs = paraboloid(reference_dir)
    with pytest.warns(rt.errors.RaytatouilleWarning, match="pupil.mapping_nonlinear"):
        base = an.pupil_function(cs, grid=16)
    with pytest.warns(rt.errors.RaytatouilleWarning):
        direct = an.pupil_function(cs, grid=16, focus_shift=0.3)
    with pytest.warns(rt.errors.RaytatouilleWarning):
        moved = an.refocus(cs, base, 0.3)
    for name in ("w", "amplitude", "status"):
        assert getattr(moved, name).tobytes() == getattr(direct, name).tobytes(), name
    assert moved.focus_shift == 0.3
    assert np.any(moved.w != base.w)  # content guard
    with pytest.raises(ValueError):
        an.refocus(cs, base, math.inf)
    with pytest.raises(ValueError):
        an.refocus(singlet(reference_dir), base, 0.3)  # another system


@pytest.mark.filterwarnings(MAPPING_WARNING)
def test_zernike_fit_and_polynomials(reference_dir: Path) -> None:
    pf = an.pupil_function(rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json")),
                           field=2, grid=24)
    fit = an.zernike_fit(pf, terms=15)
    assert fit.coefficients.shape == (15,) and not fit.coefficients.flags.writeable
    assert 0.0 < fit.strehl <= 1.0
    assert fit.sigma >= fit.sigma_bf >= 0.0
    assert fit.cells > 15 and fit.warnings == []
    assert an.noll_to_nm(11).n == 4 and an.noll_to_nm(11).m == 0
    assert an.noll_from_nm(3, 1, sine=True) == 7
    assert an.zernike(4, 1.0, 0.0) == pytest.approx(math.sqrt(3.0))
    assert an.zernike_radial(4, 0, 1.0) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        an.zernike_fit(pf, terms=0)
    with pytest.raises(ValueError):
        an.zernike(0, 0.5, 0.0)


def test_errors_with_location(reference_dir: Path) -> None:
    """AnalysisError with a location (ADR 0032 point 9, ADR 0034): a path with crystal modes
    gives /paths/i without a lost ray; the chief ray lost at a stop keeps its surface (#202)."""
    calcite = rt.compile(rt.load(reference_dir / "m4" / "calcite_walkoff.rtt.json"))
    with pytest.raises(rt.AnalysisError) as error:
        an.pupil_function(calcite, path="o", grid=4)
    assert error.value.location == "/paths/0"
    assert error.value.surface is None
    assert error.value.ray_status is None
    # The case with a lost ray, for comparison: the chief ray of field 2 (5 deg) crosses the
    # first lens surface 0.4 mm off the axis; with an aperture of 1e-6 mm it ends VIGNETTED there.
    s = rt.load(reference_dir / "m1" / "singlet_const.rtt.json")
    ed = rt.Editor(s)
    lens = ed.system.locate_surface("L1.S1")
    assert lens is not None
    ed.set(lens + "/aperture/radius", 1e-6)
    with pytest.raises(rt.AnalysisError) as lost:
        an.pupil_function(rt.compile(ed.system), field=2, grid=4)
    assert lost.value.surface == "L1.S1" and lost.value.location == lens


def test_options_out_of_range(reference_dir: Path) -> None:
    cs = paraboloid(reference_dir)
    for kwargs in ({"grid": 1}, {"apodization": -1.0}, {"focus_shift": math.nan},
                   {"lost_warning_fraction": 2.0}, {"mapping_warning": -0.1}):
        with pytest.raises(ValueError):
            an.pupil_function(cs, **kwargs)
