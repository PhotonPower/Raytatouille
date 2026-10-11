"""Accumulated geometric transformation Q and P/Q per slot from Python (#207, ADR 0034, points 2,
3 and 8 with the addendum "Q nur auf Anforderung").

The reference values of Q per event and the physics are tested in C++
(libs/rtt-trace/tests/test_accumulated_q.cpp); the bitwise comparison with rtt_py_reference is in
test_bitwise_q.py. Here the Python API: the request, shapes, dtypes, read-only views, the
relation between RayPaths and RayBatch, and one analytic check through the bindings.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille.trace import RayStatus


def periscope(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m3" / "metal_periscope.rtt.json"))


def fresnel_reflection(n_i: complex, n_t: complex, xi: float) -> tuple[complex, complex]:
    """r_s and r_p after Byrnes, arXiv:1603.02720v5, Eq. (6), Convention A (docs/quellen.md),
    with q_j = sqrt(n_j^2 - xi^2), Im q >= 0."""
    def q(n: complex) -> complex:
        root = complex(np.sqrt(complex(n * n - xi * xi)))
        return root if root.imag >= 0 else -root
    qi, qt = q(n_i), q(n_t)
    rs = (qi - qt) / (qi + qt)
    rp = (n_t * n_t * qi - n_i * n_i * qt) / (n_t * n_t * qi + n_i * n_i * qt)
    return rs, rp


def axial_rays(n: int = 2, *, with_q: bool = True) -> rt.trace.RayBatch:
    rays = rt.trace.RayBatch(n)
    for i in range(n):
        rays.pos_x[i] = 0.5 * i
    if with_q:
        rays.enable_q()
    return rays


def test_q_only_on_request() -> None:
    """A new batch has no Q: has_q is False, q_matrices() is None and q() raises; enable_q()
    adds Q = I for every ray and keeps values on a second call."""
    rays = rt.trace.RayBatch(2)
    assert not rays.has_q
    assert rays.q_matrices() is None
    with pytest.raises(ValueError, match="enable_q"):
        rays.q(0, 0)
    rays.enable_q()
    assert rays.has_q
    q = rays.q_matrices()
    assert q is not None
    assert q.shape == (2, 3, 3) and q.dtype == np.float64
    assert np.array_equal(q, np.broadcast_to(np.eye(3), (2, 3, 3)))
    rays.q(1, 0)[1] = 0.25
    rays.enable_q()
    q = rays.q_matrices()
    assert q is not None and q[1, 1, 0] == 0.25
    with pytest.raises(IndexError):
        rays.q(3, 0)


def test_q_of_one_reflection_is_householder(reference_dir: Path) -> None:
    """Path "single" of m3/metal_periscope: one reflection at the mirror tilted 45 deg about x,
    normal N = (0, -sin 45, cos 45). Q = I - 2 N N^T (ADR 0034, point 1), 1e-15; Q k0 = k."""
    cs = periscope(reference_dir)
    rays = axial_rays()
    stats = rt.trace.trace(cs, rays, path="single")
    assert stats.count(RayStatus.ALIVE) == 2
    q = rays.q_matrices()
    assert q is not None
    n = np.array([0.0, -np.sin(np.pi / 4), np.cos(np.pi / 4)])
    expected = np.eye(3) - 2.0 * np.outer(n, n)
    assert np.max(np.abs(q - expected)) <= 1e-15
    k = np.stack([rays.dir_x, rays.dir_y, rays.dir_z], axis=1)
    assert np.max(np.abs(q @ np.array([0.0, 0.0, 1.0]) - k)) <= 1e-15


def test_without_the_request_the_trace_is_the_same(reference_dir: Path) -> None:
    """Without enable_q() the batch stays without Q, and every other column is bit for bit that
    of the trace with Q (addendum #207 to ADR 0034)."""
    cs = periscope(reference_dir)
    with_q = axial_rays()
    without = axial_rays(with_q=False)
    rt.trace.trace(cs, with_q, path="periscope")
    rt.trace.trace(cs, without, path="periscope")
    assert not without.has_q and without.q_matrices() is None
    for name in ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "weight",
                 "status"):
        assert np.array_equal(getattr(without, name), getattr(with_q, name)), name
    assert np.array_equal(without.prt_matrices(), with_q.prt_matrices())


def test_ray_paths_with_polarization(reference_dir: Path) -> None:
    """record_path with polarization=True records P per slot, and Q if the batch carries Q;
    the last valid slot is bitwise the batch, slot 0 the input. The views are read-only."""
    cs = periscope(reference_dir)
    rays = axial_rays(3)
    # Ray 2 runs along x, parallel to the plane of the first mirror: Missed at event 0.
    rays.dir_x[2], rays.dir_z[2] = 1.0, 0.0
    stats, paths = rt.trace.trace(cs, rays, path="periscope", record_path=True,
                                  polarization=True)
    assert stats.count(RayStatus.ALIVE) == 2
    assert paths.polarization
    prt, q = paths.prt, paths.q
    assert prt is not None and q is not None
    assert prt.shape == (3, paths.slots, 3, 3) and prt.dtype == np.complex128
    assert q.shape == (3, paths.slots, 3, 3) and q.dtype == np.float64
    assert not prt.flags.writeable and not q.flags.writeable
    batch_p = rays.prt_matrices()
    batch_q = rays.q_matrices()
    assert batch_q is not None
    for r in range(3):
        last = int(paths.count[r]) - 1
        assert np.array_equal(prt[r, last], batch_p[r])
        assert np.array_equal(q[r, last], batch_q[r])
        assert np.array_equal(q[r, 0], np.eye(3))
    assert int(paths.lost_at[2]) >= 0  # content guard: a lost ray
    assert np.all(np.isnan(q[2, int(paths.count[2]):]))
    # Slot 1 holds Q of the first reflection, not I (content guard).
    assert np.max(np.abs(q[0, 1] - np.eye(3))) > 0.5


def test_ray_paths_polarization_without_q_and_without_request(reference_dir: Path) -> None:
    """polarization=True with a batch without Q: P per slot, q is None. Without polarization
    both are None and polarization is False."""
    cs = periscope(reference_dir)
    _, paths = rt.trace.trace(cs, axial_rays(with_q=False), path="periscope", record_path=True,
                              polarization=True)
    assert paths.polarization and paths.prt is not None and paths.q is None
    _, plain = rt.trace.trace(cs, axial_rays(), path="periscope", record_path=True)
    assert not plain.polarization and plain.prt is None and plain.q is None


def test_polarization_needs_record_path(reference_dir: Path) -> None:
    cs = periscope(reference_dir)
    with pytest.raises(ValueError, match="record_path"):
        rt.trace.trace(cs, axial_rays(), path="periscope", polarization=True)


def test_physical_retardance_of_a_path_from_the_batch(reference_dir: Path) -> None:
    """Path "periscope" (two 45 deg metal mirrors in the same plane): the physical retardance
    from P and Q of the batch equals |arg((-r_p / r_s)^2)| (ADR 0034, point 4) with the Fresnel
    amplitudes from rt.polar, 1e-10 as in the C++ test."""
    cs = periscope(reference_dir)
    rays = axial_rays(1)
    rt.trace.trace(cs, rays, path="periscope")
    q = rays.q_matrices()
    assert q is not None
    rs, rp = fresnel_reflection(1.0, complex(1.2, 7.26), float(np.sin(np.pi / 4)))
    single = -rp / rs
    expected = abs(np.angle(single * single))
    assert expected > 0.1  # content guard: a real retardance
    r = rt.polar.physical_retardance(rays.prt_matrices()[0], q[0], [0.0, 0.0, 1.0])
    assert abs(r.value - expected) <= 1e-10
