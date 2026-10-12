"""Stokes map of the pupil from Python (#211 B; ADR 0034, point 6). The C++ results themselves
are tested in rtt-analysis (test_stokes_pupil.cpp); here the bindings: columns, the normalised
input state, NaN without J and the errors. Python = C++ bit for bit is test_bitwise_jones.py."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille import analysis as an

COLUMNS = ("px", "py", "s0", "s1", "s2", "s3", "orientation", "ellipticity")


def periscope_jones(reference_dir: Path, *, field_amplitude: bool = False,
                    lost: bool = False) -> an.JonesPupil:
    """The Jones pupil of m3/metal_periscope (two metal mirrors at 45 deg) for 5 x 5 parallel
    rays along +z; chief ray index 12 (as test_jones_pupil.py). ``lost``: ray 0 runs away from
    the mirror (no J there)."""
    cs = rt.compile(rt.load(reference_dir / "m3" / "metal_periscope.rtt.json"))
    rays = rt.trace.RayBatch(25)
    for i in range(25):
        x = -3.0 + 1.5 * (i % 5)
        y = -3.0 + 1.5 * (i // 5)
        rays.pos_x[i] = x
        rays.pos_y[i] = y
        rays.pupil_x[i] = x / 3.0
        rays.pupil_y[i] = y / 3.0
    if lost:
        rays.dir_z[0] = -1.0
    return an.jones_pupil(cs, "periscope", start=rays, chief=12, field_amplitude=field_amplitude)


def test_columns_and_input_state(reference_dir: Path) -> None:
    jp = periscope_jones(reference_dir)
    st = an.stokes_pupil(jp, 3.0, 3.0j)  # normalised to (1, i) / sqrt(2)
    assert st.e1 == pytest.approx(math.sqrt(0.5)) and st.e2 == pytest.approx(1j * math.sqrt(0.5))
    for name in COLUMNS:
        column = getattr(st, name)
        assert column.shape == (25,) and column.dtype == np.float64, name
        assert not column.flags.writeable, name
    assert st.status.tolist() == jp.status.tolist()
    assert st.px.tobytes() == jp.px.tobytes()
    assert st.x_out.tolist() == jp.x_out.tolist()
    assert st.path == jp.path and st.wavelength_um == jp.wavelength_um and st.warnings == []
    # Fully polarized output: s0^2 = s1^2 + s2^2 + s3^2.
    assert np.allclose(st.s0**2, st.s1**2 + st.s2**2 + st.s3**2, rtol=1e-12, atol=0.0)
    # Power scale: s0(x) + s0(y) = 2 weight (ADR 0034, addendum #211).
    sx = an.stokes_pupil(jp, 1.0, 0.0)
    sy = an.stokes_pupil(jp, 0.0, 1.0)
    assert np.allclose(sx.s0 + sy.s0, 2.0 * jp.transmission, rtol=1e-12, atol=0.0)
    # Content guard: the metal mirrors turn the circular input elliptical.
    assert float(np.max(np.abs(np.abs(st.s3) - st.s0))) > 1e-3


def test_points_without_j_and_errors(reference_dir: Path) -> None:
    jp = periscope_jones(reference_dir, lost=True)
    st = an.stokes_pupil(jp, 1.0, 0.0)
    assert jp.status[0] != int(rt.trace.RayStatus.ALIVE)
    for name in ("s0", "s1", "s2", "s3", "orientation", "ellipticity"):
        column = getattr(st, name)
        assert math.isnan(column[0]) and np.all(np.isfinite(column[1:])), name
    with pytest.raises(ValueError):
        an.stokes_pupil(jp, 0.0, 0.0)
    with pytest.raises(ValueError):
        an.stokes_pupil(jp, math.nan, 1.0)
    with pytest.raises(ValueError):
        an.stokes_pupil(periscope_jones(reference_dir, field_amplitude=True), 1.0, 0.0)
