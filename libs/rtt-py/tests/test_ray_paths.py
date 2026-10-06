"""Ray path recording from Python (#80): rt.trace.trace(..., record_path=True) -> RayPaths.

The C++ reference values and the bitwise comparison with rtt_py_reference are in
libs/rtt-trace/tests/test_ray_paths.cpp and test_bitwise.py; here the Python API: shapes,
dtypes, read-only views, the selection and the input checks.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pytest

import raytatouille as rt
from raytatouille.trace import RayStatus


def singlet(reference_dir: Path) -> rt.CompiledSystem:
    system = rt.load(reference_dir / "m1" / "singlet_const.rtt.json")
    env = system.environment
    env.medium = "VACUUM"  # n = 1 exactly, so the OPL below is exact
    system.environment = env
    return rt.compile(system)


def axial_ray() -> rt.trace.RayBatch:
    rays = rt.trace.RayBatch(1)
    rays.pos_z[0] = -10.0
    return rays


def test_axial_ray_through_the_singlet(reference_dir: Path) -> None:
    # As the C++ reference test: the vertices z = 0, 5, 9, 106.363 mm, direction (0, 0, 1),
    # OPL 10, 15, 15 + 4 n, 15 + 4 n + 97.363 mm with n = 1.5168. Tolerance 1e-12 mm.
    cs = singlet(reference_dir)
    rays = axial_ray()
    stats, paths = rt.trace.trace(cs, rays, record_path=True)
    assert stats.count(RayStatus.ALIVE) == 1
    assert isinstance(paths, rt.trace.RayPaths)
    assert (paths.ray_count, paths.slots, paths.n_events, len(paths)) == (1, 5, 4, 1)
    assert paths.position.shape == (1, 5, 3) and paths.position.dtype == np.float64
    assert paths.direction.shape == (1, 5, 3)
    assert paths.opl.shape == (1, 5) and paths.weight.shape == (1, 5)
    assert paths.status.shape == (1, 5) and paths.status.dtype == np.uint8
    assert paths.count.dtype == np.uint32 and paths.lost_at.dtype == np.int32
    assert paths.ray_indices.dtype == np.uint64 and paths.event_surfaces.dtype == np.uint32
    z = np.array([-10.0, 0.0, 5.0, 9.0, 106.363])
    n = 1.5168
    opl = np.array([0.0, 10.0, 15.0, 15.0 + 4.0 * n, 15.0 + 4.0 * n + 97.363])
    expected = np.stack([np.zeros(5), np.zeros(5), z], axis=1)
    assert np.max(np.abs(paths.position[0] - expected)) <= 1e-12
    assert np.max(np.abs(paths.direction[0] - np.array([0.0, 0.0, 1.0]))) <= 1e-12
    assert np.max(np.abs(paths.opl[0] - opl)) <= 1e-12
    assert np.all(paths.status[0] == int(RayStatus.ALIVE))
    assert paths.count[0] == 5 and paths.lost_at[0] == -1
    ids = [cs.surface_ids[int(i)] for i in paths.event_surfaces]
    assert ids == ["STO", "L1.S1", "L1.S2", "IMG"]


def test_last_slot_is_the_final_state_and_lost_rays_are_marked(reference_dir: Path) -> None:
    # Rays at y = 0 ... 19.5 mm: ALIVE up to the stop radius of 10 mm, VIGNETTED above it at the
    # stop (event 0); slot count - 1 equals the final batch bitwise, later slots are NaN.
    cs = singlet(reference_dir)
    rays = rt.trace.RayBatch(40)
    rays.pos_z[:] = -10.0
    rays.pos_y[:] = 0.5 * np.arange(40)
    plain = rt.trace.RayBatch(40)
    plain.pos_z[:] = -10.0
    plain.pos_y[:] = 0.5 * np.arange(40)
    rt.trace.trace(cs, plain)
    _, paths = rt.trace.trace(cs, rays, record_path=True)
    for name in ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "weight",
                 "status", "last_surface"):
        assert np.array(getattr(rays, name)).tobytes() == np.array(getattr(plain, name)).tobytes()
    last = paths.count.astype(np.int64) - 1
    rows = np.arange(40)
    final = np.stack([np.array(plain.pos_x), np.array(plain.pos_y), np.array(plain.pos_z)],
                     axis=1)
    assert paths.position[rows, last].tobytes() == final.tobytes()
    assert paths.opl[rows, last].tobytes() == np.array(plain.opl).tobytes()
    assert paths.status[rows, last].tobytes() == np.array(plain.status).tobytes()
    vignetted = np.array(plain.status) == int(RayStatus.VIGNETTED)
    assert np.any(vignetted) and np.any(~vignetted)
    assert np.all(paths.lost_at[vignetted] == 0) and np.all(paths.count[vignetted] == 2)
    assert np.all(np.isnan(paths.position[vignetted, 2:]))
    assert np.all(paths.status[vignetted, 2:] == int(RayStatus.VIGNETTED))
    assert np.all(paths.lost_at[~vignetted] == -1) and np.all(paths.count[~vignetted] == 5)


def test_record_rays_selects_rays_of_a_bundle(reference_dir: Path) -> None:
    # The selection is recorded in the given order and equals the full record bitwise.
    cs = singlet(reference_dir)
    full_rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(4))
    selected_rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(4))
    n = len(full_rays)
    _, full = rt.trace.trace(cs, full_rays, record_path=True)
    selection = np.array([n - 1, 0, n // 2, 7])
    _, part = rt.trace.trace(cs, selected_rays, record_path=True, record_rays=selection, threads=1)
    assert part.ray_indices.tolist() == selection.tolist()
    for name in ("position", "direction", "opl", "weight", "status", "count", "lost_at"):
        selected: npt.NDArray[np.generic] = getattr(part, name)
        assert selected.tobytes() == getattr(full, name)[selection].tobytes(), name


def test_views_are_read_only_and_keep_the_paths_alive(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    _, paths = rt.trace.trace(cs, axial_ray(), record_path=True)
    before = sys.getrefcount(paths)
    position = paths.position
    assert sys.getrefcount(paths) == before + 1  # the view holds a reference to the RayPaths
    del position
    assert sys.getrefcount(paths) == before
    with pytest.raises(ValueError):
        paths.position[0, 0, 0] = 1.0


def test_empty_batch(reference_dir: Path) -> None:
    # No rays: arrays of shape (0, S, 3) and (0, S), event surfaces still listed.
    cs = singlet(reference_dir)
    stats, paths = rt.trace.trace(cs, rt.trace.RayBatch(0), record_path=True)
    assert sum(stats.rays) == 0
    assert paths.ray_count == 0 and paths.slots == 5
    assert paths.position.shape == (0, 5, 3) and paths.opl.shape == (0, 5)
    assert paths.count.shape == (0,) and paths.event_surfaces.shape == (4,)


def test_record_input_checks(reference_dir: Path) -> None:
    cs = singlet(reference_dir)
    rays = rt.trace.RayBatch(5)
    for selection, match in (([1, 3, 1], "twice"), ([0, 5], "not a ray index"),
                             ([-1], "negative"), ([], "empty")):
        with pytest.raises(ValueError, match=match):
            rt.trace.trace(cs, rays, record_path=True, record_rays=np.array(selection, dtype=int))
    with pytest.raises(ValueError, match="empty"):
        rt.trace.trace(cs, rays, record_path=True, record_rays=[])
    with pytest.raises(ValueError, match=r"2\*\*63"):
        rt.trace.trace(cs, rays, record_path=True,
                       record_rays=np.array([2**63], dtype=np.uint64))
    with pytest.raises(ValueError, match="integer"):
        rt.trace.trace(cs, rays, record_path=True, record_rays=np.array([0.5]))
    with pytest.raises(ValueError, match="integer"):
        rt.trace.trace(cs, rays, record_path=True, record_rays=np.zeros((2, 2), dtype=int))
    with pytest.raises(ValueError, match="record_path"):
        rt.trace.trace(cs, rays, record_rays=np.array([0]))
    with pytest.raises(ValueError, match="record_rays"):
        rt.trace.trace(cs, rays, record_path=True, max_recorded_rays=4)
    # Nothing was traced by the rejected calls.
    assert np.all(np.array(rays.last_surface) == rt.trace.NO_SURFACE)
    _, paths = rt.trace.trace(cs, rays, record_path=True, record_rays=[4, 2], max_recorded_rays=4)
    assert paths.ray_indices.tolist() == [4, 2]
    assert rt.trace.DEFAULT_MAX_RECORDED_RAYS == 10000
