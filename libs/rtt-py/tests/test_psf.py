"""Diffraction from Python (#209, P1b): PSF, psf_at and Psf.positions(), the Strehl ratio
through focus, the MTF and the encircled energy of the aberration-free paraboloid against the
circular pupil (ADR 0033, Airy, Wyant & Creath Eq. (64)); the forms on a system and on a pupil
function, errors and run control. Bitwise equality with C++ is in test_bitwise_psf.py; the
tolerances of the C++ tests are in libs/rtt-analysis/tests/test_psf.cpp and test_mtf_ee.cpp.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import analysis as an

J11 = 3.8317059702075125  # first zero of J1 (DLMF Table 10.21.1)


def paraboloid(epd: float) -> rt.CompiledSystem:
    """m2/paraboloid_stop (f = 100 mm, AIR) with the entrance pupil diameter `epd`."""
    data = json.loads(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json").to_json())
    data["aperture"]["value"] = epd
    return rt.compile(rt.System.from_json(json.dumps(data)))


@pytest.fixture(scope="module")
def f10() -> rt.CompiledSystem:
    return paraboloid(10.0)


@pytest.fixture(scope="module")
def f50() -> rt.CompiledSystem:
    return paraboloid(2.0)


def unit(cs: rt.CompiledSystem) -> float:
    """lambda F#_w in mm from the pupil mapping, F#_w = 1 / (2 |M_11|)."""
    p = an.pupil_function(cs, grid=16)
    assert p.mapping is not None
    return p.wavelength_um * 1e-3 / (2.0 * abs(float(p.mapping.m[0, 0])))


def mtf_circle(x: float) -> float:
    return 2.0 / math.pi * (math.acos(x) - x * math.sqrt(1.0 - x * x))


def test_psf_lattice_energy_and_positions(f10: rt.CompiledSystem) -> None:
    psf = an.psf(f10, grid=16)
    n = psf.length
    assert (psf.grid, n) == (16, 64)  # L = Q N with the default padding 4
    assert psf.intensity.shape == (n, n)
    # W = 0: the peak is at C' and the Strehl ratio is 1; the cells hold the energy fraction.
    assert psf.intensity[n // 2, n // 2] == psf.peak.value
    assert (psf.peak.k1, psf.peak.k2) == (0, 0)
    assert psf.strehl == pytest.approx(1.0, abs=1e-12)
    assert psf.intensity.sum() == pytest.approx(1.0 - psf.edge_energy_deficit, abs=1e-12)
    assert len(psf.wavelengths) == 1 and psf.wavelengths.weight.tolist() == [1.0]
    # Psf.positions() against the lattice: [k2 + L // 2, k1 + L // 2] = lattice @ (k1, k2).
    xy = psf.positions()
    assert xy.shape == (n, n, 2)
    assert np.array_equal(xy, an.psf_positions(psf))
    b = psf.lattice
    for k1, k2 in ((0, 0), (1, 0), (0, 1), (-3, 5), (n // 2 - 1, -n // 2)):
        expected = b @ np.array([k1, k2], dtype=float)
        assert np.allclose(xy[k2 + n // 2, k1 + n // 2], expected, rtol=0.0, atol=1e-15)
    assert b[0, 0] != 0.0 and b[1, 1] != 0.0  # content guard: a real lattice


def test_airy_zero_through_psf_at(f10: rt.CompiledSystem) -> None:
    # First zero of the Airy pattern at j11 / pi lambda F#_w (#205); sampled on 1e-3 steps.
    s = unit(f10)
    pupil = an.pupil_function(f10, grid=128)
    r = np.linspace(1.15, 1.30, 151) * s
    values = an.psf_at(pupil, np.column_stack([r, np.zeros_like(r)]))
    assert values.shape == r.shape
    zero = float(r[int(np.argmin(values))]) / s
    assert abs(zero - J11 / math.pi) <= 1e-3
    assert float(values.min()) < 1e-4
    assert an.psf_at(pupil, [0.0, 0.0])[0] == pytest.approx(1.0, abs=1e-12)  # Strehl units
    with pytest.raises(ValueError, match="shape"):
        an.psf_at(pupil, np.zeros((3, 3)))


def test_through_focus_and_psf_agree(f50: rt.CompiledSystem) -> None:
    shifts = [0.0, 0.5, 1.0]
    tf = an.through_focus(f50, shifts=shifts, grid=64)
    assert tf.shifts.tolist() == shifts
    assert tf.strehl[0] == pytest.approx(1.0, abs=1e-12)
    assert tf.strehl[1] > tf.strehl[2] > 0.0  # content guard: defocus lowers S
    assert tf.arrived_fraction.tolist() == [1.0, 1.0, 1.0]
    for shift, strehl in zip(shifts, tf.strehl):
        assert an.psf(f50, grid=64, focus_shift=shift).strehl == strehl  # bitwise (#205)


def test_mtf_of_the_circular_pupil(f10: rt.CompiledSystem) -> None:
    # Relative 1e-3 on x in [0.1, 0.9] at N = 256 (#206: e_max 2.5e-4 there).
    mtf = an.mtf(f10, grid=256)
    nu_c = 1.0 / unit(f10)
    t = mtf.tangential
    assert t.frequency[-1] == pytest.approx(nu_c, rel=1e-12)
    assert (t.modulation[0], t.modulation[-1]) == (1.0, 0.0)
    for x in np.linspace(0.1, 0.9, 9):
        value = float(np.interp(x * nu_c, t.frequency, t.modulation))
        assert abs(value / mtf_circle(float(x)) - 1.0) <= 1e-3
    assert np.allclose(mtf.sagittal.modulation, t.modulation, rtol=0.0, atol=1e-12)
    assert t.otf.size == 0 and mtf.map.shape == (0, 0)
    full = an.mtf(f10, grid=16, map=True, complex_otf=True)
    assert full.map.shape == (31, 31) and full.otf_map.shape == (31, 31)
    assert full.tangential.otf.shape == (17,)
    assert full.map[15, 15] == 1.0  # l = 0 in the middle


def test_encircled_energy_and_the_pupil_form(f10: rt.CompiledSystem) -> None:
    s = unit(f10)
    r0 = J11 / math.pi * s
    ee = an.encircled_energy(f10, [0.5 * s, r0, 2.0 * s], grid=128)
    assert 0.83 < ee.energy[1] < 0.84  # 83.8 % inside the first zero
    assert ee.energy[0] < ee.energy[1] < ee.energy[2] < 1.0
    assert (ee.offset.x, ee.offset.y) == (0.0, 0.0)
    with pytest.raises(ValueError):
        an.encircled_energy(f10, [2.0 * ee.radius_max], grid=128)
    # On a pupil function: the same values as the system form with the same pupil options.
    pupil = an.pupil_function(f10, grid=32)
    assert np.array_equal(an.mtf(pupil).tangential.modulation,
                          an.mtf(f10, grid=32).tangential.modulation)
    on_pupil = an.encircled_energy(pupil, [r0], centre=an.EeCentre.LINEAR_PHASE)
    assert on_pupil.energy.tolist() == an.encircled_energy(
        f10, [r0], grid=32, centre=an.EeCentre.LINEAR_PHASE).energy.tolist()
    with pytest.raises(ValueError, match="only apply to a system"):
        an.mtf(pupil, grid=32)
    with pytest.raises(ValueError, match="only apply to a system"):
        an.encircled_energy(pupil, [r0], field=0)


def test_errors_and_run_control(f10: rt.CompiledSystem) -> None:
    with pytest.raises(ValueError):
        an.psf(f10, grid=16, padding=1)
    with pytest.raises(ValueError):
        an.through_focus(f10, shifts=[math.inf], grid=16)
    stages: list[str] = []
    an.psf(f10, grid=16, progress=lambda done, total, stage: stages.append(stage))
    assert stages
    token = rt.CancelToken()
    token.cancel()
    for call in (lambda: an.psf(f10, grid=16, cancel=token),
                 lambda: an.mtf(f10, grid=16, cancel=token),
                 lambda: an.encircled_energy(f10, [0.001], grid=16, cancel=token),
                 lambda: an.through_focus(f10, shifts=[0.0], grid=16, cancel=token)):
        with pytest.raises(rt.Cancelled):
            call()


def test_threads_give_the_same_result(f10: rt.CompiledSystem) -> None:
    one = an.psf(f10, grid=16, threads=1)
    four = an.psf(f10, grid=16, threads=4)
    assert one.intensity.tobytes() == four.intensity.tobytes()
    assert an.mtf(f10, grid=16, threads=1).tangential.modulation.tobytes() == \
        an.mtf(f10, grid=16, threads=4).tangential.modulation.tobytes()


def test_result_types_round_trip(f10: rt.CompiledSystem) -> None:
    for obj in (an.psf(f10, grid=8), an.through_focus(f10, shifts=[0.0, 0.1], grid=8),
                an.mtf(f10, grid=8, map=True, complex_otf=True),
                an.encircled_energy(f10, [0.001], grid=8)):
        loaded = rt.results.load_json(obj.to_json())
        assert loaded.type == type(obj).__name__
        assert loaded.schema_version == rt.results.SCHEMA_VERSION
    # Without the map and the complex OTF (review of #247): the empty arrays keep their shapes,
    # (0, 0) and (0,), through to_json and load_json.
    plain = rt.results.load_json(an.mtf(f10, grid=8).to_json())
    assert plain.data["map"].shape == (0, 0) and plain.data["otf_map"].shape == (0, 0)
    assert plain.data["tangential"]["otf"].shape == (0,)
    assert plain.data["tangential"]["modulation"].shape == (9,)


def annular(epd: float = 10.0) -> rt.CompiledSystem:
    """The paraboloid with a stop aperture of radius 5 mm and inner radius 2 mm at D = 10 mm:
    epsilon = 2 / 5 = 0.4 (as test_annular.cpp, #229 C)."""
    data = json.loads(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json").to_json())
    data["aperture"]["value"] = epd
    data["root"]["children"][0]["surfaces"][0]["aperture"] = {
        "type": "circular", "radius": 5.0, "inner_radius": 2.0}
    return rt.compile(rt.System.from_json(json.dumps(data)))


def test_annular_cells_and_pupil_function(f10: rt.CompiledSystem) -> None:
    # pupil_cells(N, eps): the exact area of the annulus, pi (1 - eps^2) in units of the cell
    # area a^2 = (2 / N)^2 (ADR 0032, addendum #229); eps = 0 is bitwise pupil_cells(N).
    n, eps = 32, 0.4
    ring = an.pupil_cells(n, eps)
    assert ring.obscuration == eps
    assert ring.fraction.sum() * (2.0 / n) ** 2 == pytest.approx(math.pi * (1 - eps**2),
                                                                rel=1e-12)
    assert 0 < int(ring.inner_rim.sum()) < n * n and ring.inner_rim.dtype == np.uint8
    full = an.pupil_cells(n)
    assert full.obscuration == 0.0 and not full.inner_rim.any()
    assert full.fraction.tobytes() == an.pupil_cells(n, 0.0).fraction.tobytes()
    with pytest.raises(ValueError):
        an.pupil_cells(n, 1.0)
    # The pupil function of the obscured paraboloid carries eps from the stop, and its MTF lies
    # clearly below that of the full pupil at x = 0.25 (the annulus transfers less there).
    pupil = an.pupil_function(annular(), grid=32)
    assert pupil.obscuration == pytest.approx(eps, abs=1e-12)
    assert an.pupil_function(f10, grid=32).obscuration == 0.0
    ring_mtf = an.mtf(pupil).tangential.modulation
    open_mtf = an.mtf(f10, grid=32).tangential.modulation
    assert ring_mtf[8] < open_mtf[8] - 0.05


def test_older_pupil_files_read_without_obscuration(f10: rt.CompiledSystem) -> None:
    # PupilCells and PupilFunction gained their obscuration keys in 0.1.12 (ADR 0023): an older
    # file reads as a pupil without obscuration, a newer one must have them.
    cases = ((an.pupil_cells(8, 0.4), ("obscuration", "inner_rim")),
             (an.pupil_function(f10, grid=8), ("obscuration",)))
    for obj, keys in cases:
        envelope = json.loads(obj.to_json())
        for key in keys:
            del envelope["data"][key]
        envelope["schema_version"] = "0.1.11"
        loaded = rt.results.load_json(json.dumps(envelope))
        assert loaded.data["obscuration"] == 0.0
        if "inner_rim" in keys:
            assert loaded.data["inner_rim"].dtype == np.uint8
            assert not loaded.data["inner_rim"].any()
            assert loaded.data["inner_rim"].shape == loaded.data["fraction"].shape
        assert list(loaded.data) == list(obj.to_dict())
        envelope["schema_version"] = "0.1.12"
        with pytest.raises(ValueError, match="misses"):
            rt.results.load_json(json.dumps(envelope))


@pytest.mark.filterwarnings(
    r"ignore:warning \[pupil\.mapping_nonlinear\]:raytatouille.errors.RaytatouilleWarning")
def test_mapping_limit_rejects_the_fft_analyses() -> None:
    """ADR 0033, addendum #286: the paraboloid at F/1.7 has delta_max about 0.0225; above
    mapping_limit the FFT analyses raise AnalysisError with code pupil.mapping_too_nonlinear at
    /paths/0; the default 0.05 and 1 (off) do not. The tolerances and the ideal lens at NA 0.5
    are in libs/rtt-analysis/tests/test_mapping_limit.cpp."""
    cs = rt.compile(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json"))
    pupil = an.pupil_function(cs, grid=16, mapping_limit=0.02)  # never rejects itself
    assert pupil.mapping_limit == 0.02
    assert pupil.mapping is not None
    assert 0.02 < pupil.mapping.delta_max < 0.03
    calls = (lambda: an.psf(cs, grid=16, mapping_limit=0.02),
             lambda: an.through_focus(cs, shifts=[0.0], grid=16, mapping_limit=0.02),
             lambda: an.mtf(cs, grid=16, mapping_limit=0.02),
             lambda: an.encircled_energy(cs, [0.001], grid=16, mapping_limit=0.02),
             lambda: an.psf_at(pupil, [[0.0, 0.0]]),
             lambda: an.mtf(pupil),
             lambda: an.encircled_energy(pupil, [0.001]))
    for call in calls:
        with pytest.raises(rt.AnalysisError) as info:
            call()
        assert info.value.code == "pupil.mapping_too_nonlinear"
        assert info.value.location == "/paths/0"
        assert "delta_max" in str(info.value)
    assert an.psf(cs, grid=16).strehl > 0.0  # default 0.05
    assert an.psf(cs, grid=16, mapping_limit=1.0).strehl > 0.0
    assert an.pupil_function(cs, grid=16).mapping_limit == 0.05
    with pytest.raises(ValueError, match="mapping_limit"):
        an.psf(cs, grid=16, mapping_limit=1.5)
    with pytest.raises(ValueError, match="mapping_limit"):
        an.mtf(pupil, mapping_limit=0.5)  # an option of the system form


def test_analysis_errors_without_a_code_carry_none() -> None:
    """ADR 0022, addendum #286: AnalysisError.code is None except for
    pupil.mapping_too_nonlinear, also for its subclass NoStopError (a path without a stop)."""
    cs = rt.compile(rt.load(REFERENCE_DIR / "m1" / "paraboloid_mirror.rtt.json"))
    with pytest.raises(rt.AnalysisError) as info:
        an.psf(cs, grid=16)
    assert isinstance(info.value, rt.NoStopError)
    assert info.value.code is None


def test_older_pupil_files_read_without_mapping_limit(f10: rt.CompiledSystem) -> None:
    # PupilFunction gained mapping_limit in 0.1.14 (ADR 0032, addendum #286): an older file was
    # computed without the limit and reads with 1 (off); a newer one must have it.
    obj = an.pupil_function(f10, grid=8)
    envelope = json.loads(obj.to_json())
    assert envelope["schema_version"] == rt.results.SCHEMA_VERSION
    del envelope["data"]["mapping_limit"]
    envelope["schema_version"] = "0.1.13"
    loaded = rt.results.load_json(json.dumps(envelope))
    assert loaded.data["mapping_limit"] == 1.0
    assert list(loaded.data) == list(obj.to_dict())
    envelope["schema_version"] = "0.1.14"
    with pytest.raises(ValueError, match="misses mapping_limit"):
        rt.results.load_json(json.dumps(envelope))


def test_encircled_energy_with_relative_radii(f10: rt.CompiledSystem) -> None:
    # #267: the whole range up to radius_max in one run. An int k gives the relative radii
    # i / (k - 1); both forms equal the run with the radii in mm (the two-call way) bit for bit.
    probe = an.encircled_energy(f10, [0.0], grid=32)
    r_max = probe.radius_max
    curve = an.encircled_energy(f10, 5, grid=32)
    t = [0.0, 0.25, 0.5, 0.75, 1.0]
    assert curve.radius_max == r_max
    assert curve.radii.tolist() == [x * r_max for x in t]
    assert curve.radii[-1] == r_max
    relative = an.encircled_energy(f10, t, grid=32, relative=True)
    assert relative.radii.tolist() == curve.radii.tolist()
    assert relative.energy.tolist() == curve.energy.tolist()
    in_mm = an.encircled_energy(f10, curve.radii.tolist(), grid=32)
    assert in_mm.energy.tolist() == curve.energy.tolist()
    assert curve.energy[0] == 0.0 and bool(np.all(np.diff(curve.energy) > 0.0))  # content guard
    assert curve.energy[-1] > 0.9
    assert an.encircled_energy(f10, np.int64(3), grid=32).radii[-1] == r_max
    # The pupil form.
    pupil = an.pupil_function(f10, grid=32)
    on_pupil = an.encircled_energy(pupil, 5)
    assert on_pupil.radii.tolist() == curve.radii.tolist()
    assert on_pupil.energy.tolist() == curve.energy.tolist()
    assert an.encircled_energy(pupil, t, relative=True).energy.tolist() == curve.energy.tolist()
    # Errors: a fraction above 1, fewer than 2 radii, bool.
    with pytest.raises(ValueError, match="above 1"):
        an.encircled_energy(f10, [0.5, math.nextafter(1.0, 2.0)], grid=32, relative=True)
    with pytest.raises(ValueError, match="above 1"):
        an.encircled_energy(pupil, [1.5], relative=True)
    for k in (0, 1, -3):
        with pytest.raises(ValueError, match=">= 2"):
            an.encircled_energy(f10, k, grid=32)
    for flag in (True, np.True_):
        with pytest.raises(TypeError, match="bool"):
            an.encircled_energy(f10, flag, grid=32)  # type: ignore[arg-type]
    # A float is neither a number of radii nor a sequence (review of #267, P3-2).
    with pytest.raises(TypeError, match="a sequence of radii or an int, not float"):
        an.encircled_energy(f10, 5.0, grid=32)  # type: ignore[arg-type]
