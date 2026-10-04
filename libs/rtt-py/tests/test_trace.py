"""Ray generation and tracing (rt.trace): RayBatch as NumPy views, lifetime, threads."""

from __future__ import annotations

import gc
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille.trace import RayStatus

COLUMNS = {
    "pos_x": np.float64,
    "pos_y": np.float64,
    "pos_z": np.float64,
    "dir_x": np.float64,
    "dir_y": np.float64,
    "dir_z": np.float64,
    "wl": np.uint16,
    "opl": np.float64,
    "weight": np.float64,
    "field": np.uint16,
    "pupil_x": np.float64,
    "pupil_y": np.float64,
    "last_surface": np.uint32,
    "status": np.uint8,
}


def compiled_singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def test_new_batch_has_the_documented_default_state() -> None:
    rays = rt.trace.RayBatch(4)
    assert len(rays) == 4 and rays.size == 4
    for name, dtype in COLUMNS.items():
        column = getattr(rays, name)
        assert column.dtype == dtype, name
        assert column.shape == (4,), name
    assert np.all(rays.pos_x == 0.0) and np.all(rays.dir_z == 1.0)
    assert np.all(rays.weight == 1.0) and np.all(rays.opl == 0.0)
    assert np.all(rays.status == RayStatus.ALIVE)
    assert np.all(rays.last_surface == rt.trace.NO_SURFACE)
    p = rays.prt_matrices()
    assert p.shape == (4, 3, 3) and p.dtype == np.complex128
    assert np.array_equal(p, np.broadcast_to(np.eye(3), (4, 3, 3)))


def test_empty_batch() -> None:
    rays = rt.trace.RayBatch(0)
    assert len(rays) == 0
    assert rays.pos_x.shape == (0,)
    assert rays.prt_matrices().shape == (0, 3, 3)


def test_columns_are_views_without_copy() -> None:
    rays = rt.trace.RayBatch(3)
    weight = rays.weight
    weight[1] = 0.25
    assert rays.weight[1] == 0.25
    assert np.shares_memory(rays.weight, weight)
    rays.prt(0, 0)[2] = 2.0 + 1.0j
    assert rays.prt(0, 0)[2] == 2.0 + 1.0j
    assert rays.prt_matrices()[2, 0, 0] == 2.0 + 1.0j
    # prt_matrices() is a documented copy.
    copy = rays.prt_matrices()
    copy[0, 0, 0] = 5.0
    assert rays.prt(0, 0)[0] == 1.0


def test_a_view_keeps_the_batch_alive() -> None:
    rays = rt.trace.RayBatch(1000)
    rays.pos_z[:] = np.arange(1000.0)
    view = rays.pos_z
    del rays
    gc.collect()
    _ = [rt.trace.RayBatch(1000) for _ in range(10)]  # reuse freed memory, if any
    assert np.array_equal(view, np.arange(1000.0))


def test_size_is_fixed_from_python() -> None:
    rays = rt.trace.RayBatch(2)
    assert not hasattr(rays, "resize")


def test_prt_element_out_of_range_is_an_index_error() -> None:
    rays = rt.trace.RayBatch(1)
    with pytest.raises(IndexError):
        rays.prt(3, 0)
    with pytest.raises(IndexError):
        rays.prt(0, -1)


def test_make_rays_and_trace(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=6))
    # 3 fields x (1 + 3 * 6 * 7) rays, fields as outer loop.
    assert len(rays) == 3 * 127
    assert np.array_equal(rays.field, np.repeat(np.arange(3, dtype=np.uint16), 127))
    assert np.all(rays.wl == cs.reference_wavelength)
    assert np.all(rays.status == RayStatus.ALIVE)
    stats = rt.trace.trace(cs, rays, path="main")
    assert sum(stats.rays) == len(rays)
    assert list(stats.rays) == list(np.bincount(rays.status, minlength=7))
    # Every ray reaches the image surface IMG at z = 106.363 mm, also the rim rays (|p| = 1)
    # aimed at the stop edge: apertures pass within kApertureTolerance = kAimTolerance (#50).
    assert stats.count(RayStatus.ALIVE) == len(rays)
    assert np.all(rays.status == RayStatus.ALIVE)
    assert np.all(rays.last_surface == cs.surface_ids.index("IMG"))
    assert np.allclose(rays.pos_z, 106.363, rtol=0.0, atol=1e-9)
    norm = np.hypot(np.hypot(rays.dir_x, rays.dir_y), rays.dir_z)
    assert np.allclose(norm, 1.0, rtol=0.0, atol=1e-14)


def test_make_rays_arguments(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    rays = rt.trace.make_rays(
        cs,
        rt.trace.FanYPupil(n=5),
        path=0,
        fields=[2],
        wavelength=0,
        aiming=rt.trace.Aiming.PARAXIAL,
    )
    assert len(rays) == 5
    assert np.all(rays.field == 2) and np.all(rays.wl == 0)
    assert np.array_equal(rays.pupil_y, np.linspace(-1.0, 1.0, 5))
    single = rt.trace.make_rays(cs, rt.trace.SinglePupilPoint(px=0.5, py=-0.25), fields=[0])
    assert (single.pupil_x[0], single.pupil_y[0]) == (0.5, -0.25)
    random = rt.trace.make_rays(cs, rt.trace.RandomPupil(count=50, seed=7), fields=[0])
    again = rt.trace.make_rays(cs, rt.trace.RandomPupil(count=50, seed=7), fields=[0])
    assert np.array_equal(random.pupil_x, again.pupil_x)
    assert len(rt.trace.make_rays(cs, rt.trace.GridPupil(n=3), fields=[0])) == 5


def test_result_does_not_depend_on_the_number_of_threads(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    results = []
    for threads in (1, 4, None):
        rays = rt.trace.make_rays(cs, rt.trace.RandomPupil(count=2000, seed=3))
        rt.trace.trace(cs, rays, threads=threads)
        results.append(np.stack([rays.pos_x, rays.pos_y, rays.pos_z, rays.opl]).tobytes())
    assert results[0] == results[1] == results[2]
    with pytest.raises(ValueError):
        rt.trace.trace(cs, rt.trace.RayBatch(1), threads=0)


def test_invalid_rays_raise_at_the_api_boundary(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    rays = rt.trace.RayBatch(2)
    rays.status[0] = 200
    with pytest.raises(ValueError):
        rt.trace.trace(cs, rays)
    rays = rt.trace.RayBatch(2)
    rays.wl[1] = 9
    with pytest.raises(ValueError):
        rt.trace.trace(cs, rays)
    with pytest.raises(ValueError, match="unknown path 'nope'"):
        rt.trace.trace(cs, rt.trace.RayBatch(1), path="nope")
    with pytest.raises(IndexError):
        rt.trace.trace(cs, rt.trace.RayBatch(1), path=3)
    with pytest.raises(ValueError):
        rt.trace.make_rays(cs, rt.trace.HexapolarPupil(), fields=[7])
