"""Uniaxial crystals (#132, ADR 0026) traced from Python are bitwise equal to C++ (#134).

rtt_py_reference (crystal_cases.cpp) traces the same start rays through the calcite plate of
tests/reference/m4/calcite_walkoff.rtt.json and writes the traced columns including wave_x/y/z
and mode_index, the PRT, the trace statistics and the path transmission as named .npy arrays.
Both sides run with 1 and with 4 threads. Besides the paths "o" and "e" to the detector, the
cases use two paths that end right after the entry ("o entry", "e entry"): there the rays are
inside the crystal, so wave differs from dir (e-mode) and mode_index > 0.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille.trace import RayBatch

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

NONE = np.iinfo(np.uint64).max
Arrays = dict[str, npt.NDArray[Any]]
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "wave_x", "wave_y", "wave_z",
           "mode_index", "opl", "weight", "last_surface", "status")


def start_rays(dx: float, dy: float, dz: float) -> RayBatch:
    """As start_rays() in crystal_cases.cpp: a square grid with a pitch of 0.5 mm inside
    r <= 1 mm at z = 0 (the plate starts at z = 10), direction (dx, dy, dz). Only dir is set;
    wave and mode_index keep their defaults (reading rule, ADR 0026, point 3)."""
    points = [(0.5 * i, 0.5 * j) for i in range(-2, 3) for j in range(-2, 3)
              if (0.5 * i) ** 2 + (0.5 * j) ** 2 <= 1.0]
    rays = RayBatch(len(points))
    for k, (x, y) in enumerate(points):
        rays.pos_x[k] = x
        rays.pos_y[k] = y
        rays.dir_x[k] = dx
        rays.dir_y[k] = dy
        rays.dir_z[k] = dz
    return rays


def normal_rays() -> RayBatch:
    return start_rays(0.0, 0.0, 1.0)


def oblique_rays() -> RayBatch:
    # (1, 1, 5) / sqrt(27) without trigonometry, as in C++; sqrt is correctly rounded.
    norm = math.sqrt(27.0)
    return start_rays(1.0 / norm, 1.0 / norm, 5.0 / norm)


def plate() -> rt.CompiledSystem:
    """As with_entry_paths() in crystal_cases.cpp: the calcite plate with the paths "o entry"
    and "e entry", which end at P.S1 inside the crystal."""
    s = rt.load(REFERENCE_DIR / "m4" / "calcite_walkoff.rtt.json")
    s = rt.apply_patch(s, [
        {"op": "add", "path": "/paths/-", "value": {"name": name, "events": [
            {"surface": "P.S1", "kind": kind}]}}
        for name, kind in (("o entry", "ordinary"), ("e entry", "extraordinary"))])
    return rt.compile(s, rt.MaterialLibrary())


def flatten(system: rt.CompiledSystem, path: str, start: Callable[[], RayBatch],
            threads: int) -> Arrays:
    """Names and order as write_rays() and write_transmission() in crystal_cases.cpp."""
    rays = start()
    stats = rt.trace.trace(system, rays, path=path, threads=threads)
    results: Arrays = {name: np.array(getattr(rays, name)) for name in COLUMNS}
    for row in range(3):
        for col in range(3):
            results[f"prt{row}{col}"] = np.array(rays.prt(row, col))
    results["stats"] = np.array(stats.rays, dtype=np.uint64)
    t = an.path_transmission(system, path, start=start(), threads=threads)
    losses = t.losses
    worst = NONE if losses.worst_surface is None else losses.worst_surface
    results["t_weight"] = np.array(t.rays.weight, dtype=np.float64)
    results["t_status"] = np.array(t.rays.status, dtype=np.uint64)
    results["t_scalars"] = np.array([t.mean, t.min, t.max], dtype=np.float64)
    results["t_losses"] = np.array(
        [losses.launched, *losses.by_status, worst, losses.worst_surface_count], dtype=np.uint64)
    results["t_ints"] = np.array([t.rays_launched, t.rays_arrived, len(t.warnings)],
                                 dtype=np.uint64)
    return results


#: Case name -> (path, start rays), as `cases` in crystal_cases.cpp.
CASES: dict[str, tuple[str, Callable[[], RayBatch]]] = {
    "calcite_o_normal": ("o", normal_rays),
    "calcite_e_normal": ("e", normal_rays),
    "calcite_o_oblique": ("o", oblique_rays),
    "calcite_e_oblique": ("e", oblique_rays),
    "calcite_o_entry_oblique": ("o entry", oblique_rays),
    "calcite_e_entry_normal": ("e entry", normal_rays),
    "calcite_e_entry_oblique": ("e entry", oblique_rays),
}


@pytest.fixture(scope="module")
def compiled_plate() -> rt.CompiledSystem:
    return plate()


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_crystal_equals_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
                                    compiled_plate: rt.CompiledSystem) -> None:
    path, start = CASES[name]
    results = flatten(compiled_plate, path, start, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_cases_cover_what_they_claim(compiled_plate: rt.CompiledSystem) -> None:
    # Sanity of the cases themselves (a bitwise comparison alone stays green if both sides
    # compute the same wrong thing, e.g. rays that miss the plate): every ray arrives Alive; at
    # the detector the rays are back in vacuum (mode_index 0, wave = dir); after the entry
    # they are in the crystal (mode_index > 0), and the e-mode has wave != dir (walk-off).
    alive = int(rt.trace.RayStatus.ALIVE)
    for name, (path, start) in CASES.items():
        r = flatten(compiled_plate, path, start, 1)
        assert set(r["status"].tolist()) == {alive}, name
        wave = np.stack([r["wave_x"], r["wave_y"], r["wave_z"]])
        dir_ = np.stack([r["dir_x"], r["dir_y"], r["dir_z"]])
        if "entry" in path:
            assert (r["mode_index"] > 1.4).all(), name
            walk_off = np.linalg.norm(wave - dir_, axis=0)
            assert ((walk_off > 0.01) if path == "e entry" else (walk_off < 1e-12)).all(), name
        else:
            assert not r["mode_index"].any(), name
            assert (wave == dir_).all(), name
