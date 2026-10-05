"""Polarization quantities and coatings from Python (#62): raytatouille.polar, CoatingLibrary.

Semantics of RayBatch.prt and RayBatch.weight: ADR 0021 (P power-normalised, weight = power for
an unpolarized source). The expected values are analytic; the bitwise comparison with C++ is in
test_bitwise.py.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pytest

import raytatouille as rt
from raytatouille.trace import RayStatus

EPSILON = 1e-4  # extinction ratio of the analyzer in m3/polarizer_qwp.rtt.json


def polarizer_system(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m3" / "polarizer_qwp.rtt.json"))


def traced(cs: rt.CompiledSystem, path: str,
           sampling: rt.trace.PupilSampling | None = None) -> rt.trace.RayBatch:
    rays = rt.trace.make_rays(cs, sampling or rt.trace.HexapolarPupil(3), path=path)
    rt.trace.trace(cs, rays, path=path)
    assert np.all(np.array(rays.status) == int(RayStatus.ALIVE))
    return rays


def on_axis(rays: rt.trace.RayBatch) -> npt.NDArray[np.bool_]:
    result: npt.NDArray[np.bool_] = np.array(rays.field) == 0
    return result


# --------------------------------------------------------------------------------- batch ---


def launch_directions(rays: rt.trace.RayBatch) -> npt.NDArray[np.float64]:
    return np.stack([np.array(rays.dir_x), np.array(rays.dir_y), np.array(rays.dir_z)], axis=1)


def test_initial_directions_are_the_launch_directions(reference_dir: Path,
                                                      catalog_dir: Path) -> None:
    # k0 = Re(P^T k) (ADR 0021: P = P_T + k k0^T with k^T P_T = 0) equals the launch direction up
    # to rounding, both where the path keeps the direction (polarizer train, 3 non-trivial
    # factors: POL1, QWP, POL2) and where it deflects (ar_singlet: refraction at two coated
    # surfaces, so k != k0 and an implementation returning k would fail). Rounding estimate,
    # not a strict bound: a few eps per factor for entries of order 1, and P^T k sums 3 terms,
    # about 3 * 3 * 2.2e-16 = 2e-15; tolerance 1e-14. (Measured on MSVC: 1.1e-15 for the
    # polarizer train; the estimate was written after that first measurement.)
    cs = polarizer_system(reference_dir)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(3), path="main")
    launch = launch_directions(rays)
    rt.trace.trace(cs, rays, path="main")
    k0 = rt.polar.initial_directions(rays)
    assert k0.shape == (len(rays), 3)
    assert k0.dtype == np.float64
    assert np.max(np.abs(k0 - launch)) <= 1e-14

    coatings = rt.CoatingLibrary()
    coatings.add_catalog(catalog_dir / "coatings")
    singlet = rt.compile(rt.load(reference_dir / "m3" / "ar_singlet.rtt.json"), coatings=coatings)
    rays = rt.trace.make_rays(singlet, rt.trace.HexapolarPupil(3))
    launch = launch_directions(rays)
    rt.trace.trace(singlet, rays)
    assert np.all(np.array(rays.status) == int(RayStatus.ALIVE))
    after = launch_directions(rays)
    assert np.max(np.abs(after - launch)) > 1e-3  # the lens deflects
    assert np.max(np.abs(rt.polar.initial_directions(rays) - launch)) <= 1e-14


def test_transverse_polarization(reference_dir: Path) -> None:
    # Projection of e perpendicular to k0 of every ray, normalised: (e - k0 (k0 . e)) / |...|.
    # The fields of m3/polarizer_qwp are collimated at t = 0, 3 and 6 deg in y, so
    # k0 = (0, sin t, cos t) and the projection of y is (0, cos^2 t, -sin t cos t) / cos t =
    # (0, cos t, -sin t). The result is a unit vector transverse to k0 (tolerance 1e-15: one
    # normalisation, a few roundings), x stays x, and e parallel to k0 raises ValueError. The
    # comparisons with the field angle inherit the rounding of k0 (estimate of
    # test_initial_directions): tolerance 1e-14 there.
    rays = traced(polarizer_system(reference_dir), "main")
    t = np.radians(np.array([0.0, 3.0, 6.0]))[np.array(rays.field)]
    k0 = rt.polar.initial_directions(rays)
    zeros = np.zeros(len(rays))
    assert np.max(np.abs(k0 - np.stack([zeros, np.sin(t), np.cos(t)], axis=1))) <= 1e-14
    y = rt.polar.transverse_polarization(rays, np.array([0.0, 1.0, 0.0]))
    assert y.shape == (len(rays), 3) and y.dtype == np.complex128
    assert np.max(np.abs(np.linalg.norm(y, axis=1) - 1.0)) <= 1e-15
    assert np.max(np.abs(np.einsum("ij,ij->i", k0, y))) <= 1e-15
    assert np.max(np.abs(y - np.stack([zeros, np.cos(t), -np.sin(t)], axis=1))) <= 1e-14
    x = rt.polar.transverse_polarization(rays, np.array([1.0, 0.0, 0.0]))
    assert np.max(np.abs(x - np.array([1.0, 0.0, 0.0]))) <= 1e-15
    with pytest.raises(ValueError, match="parallel"):
        rt.polar.transverse_polarization(rays, np.array([0.0, 0.0, 1.0]))


def test_crossed_polarizers_with_quarter_wave_plate(reference_dir: Path) -> None:
    # Rays k = (0, sin t, cos t) of the field angle t. In the transverse basis e1 = x,
    # e2 = k x x = (0, cos t, -sin t): the first polarizer passes e1, the analyzer axis y
    # projects onto e2 (crossed with e1 for every t), and the quarter-wave axis (1, 1, 0)/sqrt 2
    # projects to angle phi from e1 with tan(phi) = cos t (45 deg only on axis). A quarter-wave
    # plate at phi turns e1 into |e2 component|^2 = sin^2(2 phi) / 2, so the analyzer with
    # extinction ratio eps passes sin^2(2 phi)/2 + eps (1 - sin^2(2 phi)/2) of the polarized x
    # light and half of that of unpolarized light (ideal elements: rtt/polar/ideal.hpp, #60).
    # On axis: (1 + eps)/2 and (1 + eps)/4. Polarized y is blocked by the perfect first
    # polarizer: 0. Tolerance 1e-12.
    cs = polarizer_system(reference_dir)
    rays = traced(cs, "main")
    axis = on_axis(rays)
    weight = np.array(rays.weight)
    unpolarized = rt.polar.transmission(rays)
    assert unpolarized.tobytes() == weight.tobytes()
    cos_t = rt.polar.initial_directions(rays)[:, 2]
    sin2phi = 2.0 * cos_t / (1.0 + cos_t ** 2)
    through = 0.5 * sin2phi ** 2 + EPSILON * (1.0 - 0.5 * sin2phi ** 2)
    assert np.max(np.abs(weight - 0.5 * through)) <= 1e-12
    assert np.max(np.abs(weight[axis] - (1.0 + EPSILON) / 4.0)) <= 1e-12
    assert np.ptp(weight) > 1e-6  # the case of the bitwise test needs varying weights
    x = rt.polar.transverse_polarization(rays, np.array([1.0, 0.0, 0.0]))
    y = rt.polar.transverse_polarization(rays, np.array([0.0, 1.0, 0.0]))
    assert np.max(np.abs(rt.polar.transmission(rays, x) - through)) <= 1e-12
    assert np.max(np.abs(rt.polar.transmission(rays, y))) <= 1e-12
    # Collimated fields in y only: x is transverse for every ray up to rounding (|k0_x| ~ 1e-18,
    # far below the check of 1e-12), so a single state for all rays is accepted and gives the
    # same result up to rounding (1e-15).
    single = rt.polar.transmission(rays, np.array([1.0, 0.0, 0.0], dtype=np.complex128))
    assert np.max(np.abs(single - rt.polar.transmission(rays, x))) <= 1e-15
    with pytest.raises(ValueError, match="rays"):
        rt.polar.transmission(rays, np.zeros((0, 3), dtype=np.complex128))


def test_diattenuation_of_the_polarizer_train(reference_dir: Path) -> None:
    # The perfect first polarizer gives a rank-1 transverse part: D = 1, minimum = 0, and
    # maximum^2 = (1 + eps)/2 is the transmission of x (power, without s; ADR 0021). The axis of
    # maximum transmission is x up to a phase. Tolerance 1e-12.
    rays = traced(polarizer_system(reference_dir), "main")
    axis = on_axis(rays)
    d = rt.polar.diattenuation(rays)
    assert d.value.shape == (len(rays),)
    assert d.axis.shape == (len(rays), 3)
    assert np.max(np.abs(d.value[axis] - 1.0)) <= 1e-12
    assert np.max(np.abs(d.minimum[axis])) <= 1e-12
    assert np.max(np.abs(d.maximum[axis] ** 2 - (1.0 + EPSILON) / 2.0)) <= 1e-12
    assert np.max(np.abs(np.abs(d.axis[axis, 0]) - 1.0)) <= 1e-12


def test_stokes_after_the_quarter_wave_plate_is_circular(reference_dir: Path) -> None:
    # Path "circular" ends after the quarter-wave plate: x light becomes circular. On axis the
    # fast axis is (1, 1, 0)/sqrt 2 with phase e^{-i pi/4}, the slow axis (-1, 1, 0)/sqrt 2 with
    # e^{+i pi/4} (rtt/polar/ideal.hpp, #60), so x = (f - s)/sqrt 2 leaves as
    # (cos(pi/4), -i sin(pi/4), 0) = (1, -i)/sqrt 2 and S3 = 2 Im(E1 E2*) = +1 with e1 = x,
    # e2 = z x x = y: right circular, as Lam's (i, 1) up to the phase -i (docs/architecture.md,
    # Händigkeit und Stokes). S0 = 1, S1 = S2 = 0. Tolerance 1e-12.
    rays = traced(polarizer_system(reference_dir), "circular")
    axis = on_axis(rays)
    x = rt.polar.transverse_polarization(rays, np.array([1.0, 0.0, 0.0]))
    s = rt.polar.stokes(rays, x, np.array([1.0, 0.0, 0.0]))
    assert s.shape == (len(rays), 4)
    assert np.max(np.abs(s[axis, 0] - 1.0)) <= 1e-12
    assert np.max(np.abs(s[axis, 3] - 1.0)) <= 1e-12
    assert np.max(np.abs(s[axis, 1:3])) <= 1e-12


def test_retardance_of_the_quarter_wave_plate_and_nan_for_deflected_paths(
        reference_dir: Path) -> None:
    # Path "retarder only": P is the quarter-wave plate, k = k0 for every ray, and the ideal
    # retarder is diag(e^{-i pi/4}, e^{+i pi/4}) in a transverse basis for every direction
    # (rtt/polar/ideal.hpp): retardance pi/2 for all rays. Michelson arms leave along -y while
    # entering along +z: retardance undefined -> NaN. Tolerance 1e-12.
    rays = traced(polarizer_system(reference_dir), "retarder only")
    r = rt.polar.retardance(rays)
    assert np.max(np.abs(r.value - math.pi / 2.0)) <= 1e-12
    assert r.fast_axis.shape == (len(rays), 3)
    cs = rt.compile(rt.load(reference_dir / "m0" / "michelson.rtt.json"))
    arm = rt.trace.RayBatch(1)
    rt.trace.trace(cs, arm, path="reference arm")
    assert np.array(arm.status)[0] == int(RayStatus.ALIVE)
    deflected = rt.polar.retardance(arm)
    assert math.isnan(deflected.value[0])
    assert np.all(np.isnan(deflected.fast_axis[0].real))


def test_michelson_arms_split_the_power(reference_dir: Path) -> None:
    # Ideal 50/50 beam splitter: each arm transmits once and reflects once, weight 1/4.
    cs = rt.compile(rt.load(reference_dir / "m0" / "michelson.rtt.json"))
    for path in ("reference arm", "test arm"):
        rays = rt.trace.RayBatch(1)
        rt.trace.trace(cs, rays, path=path)
        assert abs(rays.weight[0] - 0.25) <= 1e-15


def test_batch_input_checks(reference_dir: Path) -> None:
    rays = traced(polarizer_system(reference_dir), "main")
    n = len(rays)
    with pytest.raises(ValueError, match="unit"):
        rt.polar.transmission(rays, np.array([2.0, 0.0, 0.0], dtype=np.complex128))
    with pytest.raises(ValueError, match="transverse"):
        rt.polar.transmission(rays, np.array([0.0, 0.0, 1.0], dtype=np.complex128))
    with pytest.raises(ValueError, match="rays"):
        rt.polar.transmission(rays, np.zeros((n + 1, 3), dtype=np.complex128))
    with pytest.raises(ValueError, match="finite"):
        rt.polar.transmission(rays, np.array([math.nan, 0.0, 0.0], dtype=np.complex128))
    with pytest.raises(ValueError, match="parallel"):
        rt.polar.transverse_polarization(rays, np.array([0.0, 0.0, 1.0]))
    with pytest.raises(ValueError, match="parallel"):
        x = rt.polar.transverse_polarization(rays, np.array([1.0, 0.0, 0.0]))
        rt.polar.stokes(rays, x, np.array([0.0, 0.0, 1.0]))


# ------------------------------------------------------------------------- single matrix ---


def test_single_matrix_functions() -> None:
    # Ideal mirror at 45 deg: P = prt_matrix with (-1, +1) is the field of an ideal conductor
    # -(I - 2 N N^T) on transverse fields; Q^-1 P = -I there, physical retardance 0 (#57).
    # Transmission (0.6, 0.8): D = (0.64 - 0.36) / 1.0 = 0.28 (Lam, Eq. (4.3)). A retarder
    # diag(e^{-i pi/4}, e^{+i pi/4}) at normal incidence has retardance pi/2 (Lam, Eq. (4.4)).
    # Stokes of (1, i, 0)/sqrt(2) with axis x: S0 = 1, S3 = 2 Im(E1 E2*) = -1 (#60).
    # Tolerances: 1e-14 for values built from a few tens of roundings of quantities of order 1
    # (sqrt(0.5), products, a 3x3 SVD that is backward stable), so that the result does not
    # depend on the compiler (MSVC, GCC, clang); 1e-12 for the retardances (phase of eigenvalues
    # after a polar decomposition, as in the C++ tests of #57).
    k_in = np.array([0.0, 0.0, 1.0])
    n = np.array([0.0, -math.sqrt(0.5), math.sqrt(0.5)])
    k_out = k_in - 2.0 * np.dot(k_in, n) * n
    p = rt.polar.prt_matrix(k_in, k_out, n, -1.0, 1.0)
    assert p.shape == (3, 3) and p.dtype == np.complex128
    conductor = -(np.eye(3) - 2.0 * np.outer(n, n))
    for e in (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])):
        assert np.max(np.abs(p @ e - conductor @ e)) <= 1e-14
    q = rt.polar.geometric_transform(k_in, k_out, n, True)
    assert q.shape == (3, 3) and q.dtype == np.float64
    assert abs(rt.polar.physical_retardance(p, q, k_in).value) <= 1e-12

    t = rt.polar.prt_matrix(k_in, k_in, k_in, 0.6, 0.8)
    d = rt.polar.diattenuation(t, k_in, k_in)
    assert abs(d.value - 0.28) <= 1e-14
    assert abs(d.maximum - 0.8) <= 1e-14 and abs(d.minimum - 0.6) <= 1e-14
    assert d.axis.shape == (3,)

    phase = complex(math.cos(math.pi / 4.0), math.sin(math.pi / 4.0))
    m = rt.polar.prt_matrix(k_in, k_in, k_in, phase.conjugate(), phase)
    assert abs(rt.polar.retardance(m, k_in).value - math.pi / 2.0) <= 1e-12

    e = np.array([1.0, 1.0j, 0.0]) / math.sqrt(2.0)
    s = rt.polar.stokes(e, np.array([1.0, 0.0, 0.0]), k_in)
    assert s.shape == (4,)
    assert abs(s[0] - 1.0) <= 1e-14 and abs(s[3] + 1.0) <= 1e-14


def test_single_matrix_input_checks() -> None:
    z = np.array([0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="unit"):
        rt.polar.prt_matrix(np.array([0.0, 0.0, 2.0]), z, z, 1.0, 1.0)
    with pytest.raises(ValueError, match="3"):
        rt.polar.prt_matrix(np.array([0.0, 1.0]), z, z, 1.0, 1.0)
    with pytest.raises(ValueError, match="finite"):
        rt.polar.diattenuation(np.full((3, 3), math.nan, dtype=np.complex128), z, z)
    with pytest.raises(ValueError, match="3x3"):
        rt.polar.retardance(np.eye(2, dtype=np.complex128), z)


# ------------------------------------------------------------------------------- coatings ---


def test_compile_with_coatings(reference_dir: Path, catalog_dir: Path) -> None:
    # ADR 0019: CoatingRef needs a CoatingLibrary; without one compile reports a CompileError.
    system = rt.load(reference_dir / "m3" / "ar_singlet.rtt.json")
    with pytest.raises(rt.CompileError):
        rt.compile(system)
    coatings = rt.CoatingLibrary()
    coatings.add_catalog(catalog_dir / "coatings" / "demo.json")
    assert "DEMO:AR_MGF2" in coatings
    assert "DEMO:NONE" not in coatings
    cs = rt.compile(system, coatings=coatings)
    assert cs.field_count == 3
    with pytest.raises(ValueError):
        coatings.add_catalog(catalog_dir / "coatings" / "demo.json")  # registered twice


def test_coating_catalog_error(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text('{"format": "raytatouille-coatings", "schema_version": "0.1.0", '
                   '"catalog": "BAD", "coatings": [{"name": "X", "layers": []}]}')
    with pytest.raises(rt.CoatingCatalogError) as info:
        rt.CoatingLibrary().add_catalog(bad)
    assert info.value.file.endswith("bad.json")
    assert info.value.pointer.startswith("/coatings/0")
    assert isinstance(info.value, rt.RaytatouilleError)


def test_absorbing_ar_plate_from_python(reference_dir: Path, catalog_dir: Path) -> None:
    # Same reference as the C++ test of #61: weight = (1 - R)^2 exp(-4 pi kappa d / lambda) with
    # r = (n0 ns - nc^2)/(n0 ns + nc^2) (quarter-wave layer, #58); kappa enters r, t and the
    # power factors only in O(kappa^2). Tolerance 1e-10 relative.
    coatings = rt.CoatingLibrary()
    coatings.add_catalog(catalog_dir / "coatings")
    cs = rt.compile(rt.load(reference_dir / "m3" / "absorbing_ar_plate.rtt.json"),
                    coatings=coatings)
    rays = rt.trace.RayBatch(1)
    rt.trace.trace(cs, rays)
    assert np.array(rays.status)[0] == int(RayStatus.ALIVE)
    assert np.array(rays.last_surface)[0] == 2  # IMG
    nc2 = 1.38 * 1.38
    r = (1.52 - nc2) / (1.52 + nc2)
    expected = (1.0 - r * r) ** 2 * math.exp(-4.0 * math.pi * 1e-6 * 10.0 / 0.55e-3)
    assert abs(rays.weight[0] - expected) <= 1e-10 * expected
