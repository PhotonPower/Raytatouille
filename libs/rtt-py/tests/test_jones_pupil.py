"""Jones pupil from Python (#208 B2; ADR 0034, points 5, 6 and 8). The C++ results themselves are
tested in rtt-analysis; here the bindings: columns, options, both forms, errors with location and
the agreement of the Python calls with each other. Python = C++ bit for bit is
test_bitwise_jones.py."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille import analysis as an

MAPS = ("px", "py", "w", "transmission", "diattenuation", "diattenuation_orientation",
        "diattenuation_ellipticity", "retardance", "retardance_orientation",
        "retardance_ellipticity")
FRAMES = ("x_in", "y_in", "k_in", "x_out", "y_out", "k_out")


def singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def periscope(reference_dir: Path) -> rt.CompiledSystem:
    """m3/metal_periscope: two metal mirrors (CONST:1.2,7.26) at 45 deg in the y-z plane, path
    "periscope" to the detector; folded, so only the form with start rays applies."""
    return rt.compile(rt.load(reference_dir / "m3" / "metal_periscope.rtt.json"))


def parallel_rays() -> rt.trace.RayBatch:
    """5 x 5 rays along +z from z = 0 at x, y in {-3, -1.5, 0, 1.5, 3}, labelled (x, y) / 3; the
    chief ray (0, 0) has index 12 (as jones_cases.cpp)."""
    rays = rt.trace.RayBatch(25)
    for i in range(25):
        x = -3.0 + 1.5 * (i % 5)
        y = -3.0 + 1.5 * (i // 5)
        rays.pos_x[i] = x
        rays.pos_y[i] = y
        rays.pupil_x[i] = x / 3.0
        rays.pupil_y[i] = y / 3.0
    return rays


def chief_index(jp: an.JonesPupil) -> int:
    hits = np.flatnonzero((jp.px == 0.0) & (jp.py == 0.0))
    assert hits.size == 1
    return int(hits[0])


def test_columns_of_the_convenience_form(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    jp = an.jones_pupil(cs, field=1, rays="grid:9")
    n = jp.px.shape[0]
    assert n > 40
    assert jp.j.shape == (n, 2, 2) and jp.j.dtype == np.complex128
    assert not jp.j.flags.writeable
    for name in MAPS:
        column = getattr(jp, name)
        assert column.shape == (n,) and column.dtype == np.float64, name
        assert not column.flags.writeable, name
    assert jp.status.shape == (n,) and jp.status.dtype == np.uint8
    for name in FRAMES:
        assert getattr(jp, name).shape == (3,), name
    assert jp.path == 0 and jp.field == 1 and jp.wavelength == cs.reference_wavelength
    assert jp.basis == an.JonesBasis.DOUBLE_POLE
    assert jp.remove_global_phase is False and jp.field_amplitude is False
    assert math.isfinite(jp.sphere.radius)
    assert jp.losses.launched == n and jp.warnings == []
    c = chief_index(jp)
    assert jp.w[c] == 0.0
    assert np.all(np.isfinite(jp.j[jp.status == int(rt.trace.RayStatus.ALIVE)]))
    # Content guard: the refractions of the singlet diattenuate off the axis (Fresnel).
    assert float(np.nanmax(jp.diattenuation)) > 1e-6
    # The frame is orthonormal and right-handed (x_out x y_out = k_out).
    assert np.allclose(np.cross(jp.x_out, jp.y_out), jp.k_out, atol=1e-15)


def test_options_of_the_convenience_form(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    base = an.jones_pupil(cs, field=1, rays="grid:9")
    dipole = an.jones_pupil(cs, field=1, rays="grid:9", basis=an.JonesBasis.DIPOLE)
    assert dipole.basis == an.JonesBasis.DIPOLE
    c = chief_index(base)
    # The dipole basis is the chief frame at the chief ray and differs elsewhere (ADR 0034, 5).
    assert np.allclose(dipole.j[c], base.j[c], atol=1e-15)
    assert float(np.nanmax(np.abs(dipole.j - base.j))) > 1e-9
    assert dipole.diattenuation.tobytes() == base.diattenuation.tobytes()
    # Field amplitudes: j divided by a positive factor per point, the maps unchanged (B1).
    field = an.jones_pupil(cs, field=1, rays="grid:9", field_amplitude=True)
    assert field.field_amplitude is True
    assert field.transmission.tobytes() == base.transmission.tobytes()
    assert field.diattenuation.tobytes() == base.diattenuation.tobytes()
    alive = base.status == int(rt.trace.RayStatus.ALIVE)
    ratio = np.abs(base.j[alive, 0, 0]) / np.abs(field.j[alive, 0, 0])
    assert np.all(ratio > 0.0)
    assert float(np.max(np.abs(ratio - 1.0))) > 1e-4  # content guard: the factors are not 1
    with pytest.raises(ValueError):
        an.jones_pupil(cs, lost_warning_fraction=2.0)


def test_form_with_given_start_rays(reference_dir: Path) -> None:
    cs = periscope(reference_dir)
    jp = an.jones_pupil(cs, "periscope", start=parallel_rays(), chief=12)
    assert math.isnan(jp.sphere.radius)  # no reference sphere in this form
    assert jp.px.shape == (25,) and jp.status.tolist() == [int(rt.trace.RayStatus.ALIVE)] * 25
    # Two reflections in the same plane of incidence: J = diag(r_s^2, r_p^2) up to the phase,
    # |r_s| != |r_p| for the metal, so D > 0 at every point (Lam, Eq. (3.4)).
    assert float(np.max(np.abs(jp.j[:, 0, 1]))) <= 1e-9
    assert float(np.max(np.abs(jp.j[:, 1, 0]))) <= 1e-9
    assert float(np.min(np.abs(np.abs(jp.j[:, 0, 0]) - np.abs(jp.j[:, 1, 1])))) > 1e-3
    assert float(np.min(jp.diattenuation)) > 1e-2
    assert np.allclose(jp.k_out, [0.0, 0.0, 1.0], atol=1e-12)  # out along +z again
    # The global phase: det J(0, 0) real and positive afterwards.
    turned = an.jones_pupil(cs, "periscope", start=parallel_rays(), chief=12,
                            remove_global_phase=True)
    det = np.linalg.det(turned.j[12])
    assert abs(math.atan2(det.imag, det.real)) <= 1e-12
    assert abs(turned.global_phase) > 1e-3  # content guard: there was a phase to remove


def test_forms_and_errors(reference_dir: Path) -> None:
    cs = periscope(reference_dir)
    start = parallel_rays()
    with pytest.raises(ValueError):
        an.jones_pupil(cs, "periscope", start=start, chief=12, field=0)  # mixed forms
    with pytest.raises(ValueError):
        an.jones_pupil(cs, "periscope", start=start)  # no chief ray
    with pytest.raises(ValueError):
        an.jones_pupil(cs, "periscope", chief=12)  # chief ray without start rays
    with pytest.raises(ValueError):
        an.jones_pupil(cs, "periscope", start=start, chief=25)  # out of range
    with pytest.raises(rt.ParaxialError):
        an.jones_pupil(cs, "periscope")  # the convenience form needs a symmetric path
    calcite = rt.compile(rt.load(reference_dir / "m4" / "calcite_walkoff.rtt.json"))
    with pytest.raises(rt.AnalysisError) as error:
        an.jones_pupil(calcite, path="o", rays="grid:4")
    assert error.value.location == "/paths/0"


def test_warnings_are_issued(reference_dir: Path) -> None:
    """The warnings stay in ``warnings`` and are issued as RaytatouilleWarning (ADR 0022, 0023):
    rays.lost with a lost ray and lost_warning_fraction 0. The other codes of the Jones pupil
    (pupil.basis_pole, pupil.field_amplitude_absorbing, pupil.chief_phase_unavailable) have
    their C++ cases in test_warnings.cpp; the path through raytatouille is the same."""
    cs = periscope(reference_dir)
    start = parallel_rays()
    start.dir_z[0] = -1.0  # runs away from the mirror: lost, not a pole
    with pytest.warns(rt.RaytatouilleWarning):
        jp = an.jones_pupil(cs, "periscope", start=start, chief=12, lost_warning_fraction=0.0)
    assert [w.code for w in jp.warnings] == ["rays.lost"]
    assert jp.status[0] != int(rt.trace.RayStatus.ALIVE) and jp.w[0] == 0.0
    assert math.isnan(jp.j[0, 0, 0].real)
