"""Configurations and node frames from Python are bitwise equal to C++ (#169).

rtt_py_reference (configuration_cases.cpp) compiles tests/reference/m5/zoom.rtt.json for both
configurations (by index and by name), traces the same start rays and writes the columns and
the frames of every node as named .npy arrays. Both sides run with 1 and with 4 threads.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import layout
from raytatouille.trace import RayBatch

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_y", "dir_z", "opl", "status")


def start_rays() -> RayBatch:
    """As start_rays() in configuration_cases.cpp."""
    rays = RayBatch(4)
    for k in range(3):
        rays.pos_y[k] = k - 1.0
        rays.pos_z[k] = -5.0
    rays.pos_z[3] = -5.0
    rays.dir_y[3] = 0.0625
    rays.dir_z[3] = math.sqrt(1.0 - 0.0625 * 0.0625)
    return rays


def flatten(cs: rt.CompiledSystem, threads: int) -> Arrays:
    rays = start_rays()
    stats = rt.trace.trace(cs, rays, path=0, threads=threads)
    out: Arrays = {name: np.array(getattr(rays, name)) for name in COLUMNS}
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    out["configuration"] = np.array([cs.configuration], dtype=np.uint64)
    # 3 x 4 rows of each 4 x 4 matrix, reference then to_global, node by node.
    out["frames"] = np.array([m[:3, :].ravel() for f in layout.node_frames(cs)
                              for m in (f.reference, f.to_global)], dtype=np.float64).ravel()
    return out


CASES = {"zoom_wide": 0, "zoom_tele": "tele"}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_configurations_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    system = rt.load(REFERENCE_DIR / "m5" / "zoom.rtt.json")
    results = flatten(rt.compile(system, configuration=CASES[name]), threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_differ() -> None:
    # Sanity: the two configurations give different geometry, and the rays arrive.
    system = rt.load(REFERENCE_DIR / "m5" / "zoom.rtt.json")
    wide = flatten(rt.compile(system), 1)
    tele = flatten(rt.compile(system, configuration=1), 1)
    assert not np.array_equal(wide["frames"], tele["frames"])
    alive = int(rt.trace.RayStatus.ALIVE)
    assert set(wide["status"].tolist()) == {alive} and set(tele["status"].tolist()) == {alive}
