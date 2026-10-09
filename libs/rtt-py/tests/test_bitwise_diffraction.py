"""Diffraction orders (#127, ADR 0025) traced from Python are bitwise equal to C++ (#134).

rtt_py_reference (diffraction_cases.cpp) traces the same start rays through the grating bench
tests/reference/m4/grating_transmission.rtt.json and two variants of it, and writes the traced
columns, the PRT, the trace statistics and the path transmission as named .npy arrays. Both
sides run with 1 and with 4 threads. The order +6 is evanescent at normal incidence (status
Evanescent at the grating) and propagates off the screen at oblique incidence (Vignetted); the
efficiency variant gives order -1 the efficiency 0 (weight 0, status Alive).
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
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "weight",
           "last_surface", "status")


def start_rays(dx: float, dz: float) -> RayBatch:
    """As start_rays() in diffraction_cases.cpp: a square grid with a pitch of 0.5 mm inside
    r <= 2 mm at z = -5 mm (in front of the stop), direction (dx, 0, dz)."""
    points = [(0.5 * i, 0.5 * j) for i in range(-4, 5) for j in range(-4, 5)
              if (0.5 * i) ** 2 + (0.5 * j) ** 2 <= 4.0]
    rays = RayBatch(len(points))
    for k, (x, y) in enumerate(points):
        rays.pos_x[k] = x
        rays.pos_y[k] = y
        rays.pos_z[k] = -5.0
        rays.dir_x[k] = dx
        rays.dir_z[k] = dz
    return rays


def normal_rays() -> RayBatch:
    return start_rays(0.0, 1.0)


def oblique_rays() -> RayBatch:
    # (-1, 0, 5) / sqrt(26) without trigonometry, as in C++; sqrt is correctly rounded.
    norm = math.sqrt(26.0)
    return start_rays(-1.0 / norm, 5.0 / norm)


def bench() -> rt.System:
    return rt.load(REFERENCE_DIR / "m4" / "grating_transmission.rtt.json")


def grating_pointer(system: rt.System) -> str:
    """JSON pointer of the grating surface G in the edit form."""
    pointer = system.locate_surface("G")
    assert pointer is not None
    return pointer


def reflection_grating() -> rt.System:
    """As reflection_grating() in diffraction_cases.cpp: the grating as an ideal mirror and the
    path "reflect +1" back through the stop."""
    s = bench()
    grating = grating_pointer(s)
    return rt.apply_patch(s, [
        {"op": "replace", "path": grating + "/interaction", "value": {"type": "ideal_mirror"}},
        {"op": "add", "path": "/paths/-", "value": {"name": "reflect +1", "events": [
            {"surface": "STO", "kind": "transmit"},
            {"surface": "G", "kind": "reflect", "order": 1},
            {"surface": "STO", "kind": "transmit"}]}},
    ])


def with_efficiencies() -> rt.System:
    """As with_efficiencies() in diffraction_cases.cpp: order 0 0.3, order +1 0.4, all others
    0."""
    s = bench()
    return rt.apply_patch(s, [{
        "op": "add", "path": grating_pointer(s) + "/diffraction_efficiency",
        "value": [{"order": 0, "efficiency": 0.3}, {"order": 1, "efficiency": 0.4}]}])


Systems = dict[str, rt.CompiledSystem]


def systems() -> Systems:
    lib = rt.MaterialLibrary()
    return {
        "transmission": rt.compile(bench(), lib),
        "reflection": rt.compile(reflection_grating(), lib),
        "efficiency": rt.compile(with_efficiencies(), lib),
    }


def flatten(system: rt.CompiledSystem, path: str, start: Callable[[], RayBatch],
            threads: int) -> Arrays:
    """Names and order as write_rays() and write_transmission() in diffraction_cases.cpp."""
    path_id = system.find_path(path)
    assert path_id is not None
    rays = start()
    stats = rt.trace.trace(system, rays, path=path_id, threads=threads)
    results: Arrays = {name: np.array(getattr(rays, name)) for name in COLUMNS}
    for row in range(3):
        for col in range(3):
            results[f"prt{row}{col}"] = np.array(rays.prt(row, col))
    results["stats"] = np.array(stats.rays, dtype=np.uint64)
    t = an.path_transmission(system, path_id, start=start(), threads=threads)
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


Case = tuple[str, str, Callable[[], RayBatch]]

#: Case name -> (system, path, start rays), as `cases` in diffraction_cases.cpp.
CASES: dict[str, Case] = {
    "grating_m1_normal": ("transmission", "order -1", normal_rays),
    "grating_0_normal": ("transmission", "order 0", normal_rays),
    "grating_p1_normal": ("transmission", "order +1", normal_rays),
    "grating_p6_normal": ("transmission", "order +6", normal_rays),
    "grating_p1_oblique": ("transmission", "order +1", oblique_rays),
    "grating_p6_oblique": ("transmission", "order +6", oblique_rays),
    "reflection_p1": ("reflection", "reflect +1", normal_rays),
    "efficiency_m1": ("efficiency", "order -1", normal_rays),
    "efficiency_0": ("efficiency", "order 0", normal_rays),
    "efficiency_p1": ("efficiency", "order +1", normal_rays),
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_diffraction_equals_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
                                        compiled_systems: Systems) -> None:
    system, path, start = CASES[name]
    results = flatten(compiled_systems[system], path, start, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
def test_cases_cover_the_statuses(compiled_systems: Systems) -> None:
    # Sanity of the cases themselves: each one tests what its name says.
    status = {name: set(flatten(compiled_systems[s], p, start, 1)["status"].tolist())
              for name, (s, p, start) in CASES.items()}
    alive, vignetted, evanescent = (int(rt.trace.RayStatus.ALIVE),
                                    int(rt.trace.RayStatus.VIGNETTED),
                                    int(rt.trace.RayStatus.EVANESCENT))
    expected = {name: {alive} for name in CASES}
    expected["grating_p6_normal"] = {evanescent}
    expected["grating_p6_oblique"] = {vignetted}
    assert status == expected
    weight = flatten(compiled_systems["efficiency"], "order -1", normal_rays, 1)["weight"]
    assert not weight.any()
