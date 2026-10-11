"""Q and P per slot traced from Python are bitwise equal to C++ (#207, ADR 0034).

rtt_py_reference (q_cases.cpp) traces tests/reference/m3/metal_periscope.rtt.json along both
paths with the same three start rays, Q requested and P/Q recorded per slot, and writes the
columns as named .npy arrays (per-slot arrays flattened in C order). Both sides run with 1 and
with 4 threads.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille.trace import RayBatch

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]


def start_rays() -> RayBatch:
    """As start_rays() in q_cases.cpp."""
    rays = RayBatch(3)
    rays.pos_x[1] = 0.5
    rays.dir_x[2], rays.dir_z[2] = 1.0, 0.0
    rays.enable_q()
    return rays


def flatten(path: str, threads: int) -> Arrays:
    cs = rt.compile(rt.load(REFERENCE_DIR / "m3" / "metal_periscope.rtt.json"))
    rays = start_rays()
    stats, paths = rt.trace.trace(cs, rays, path=path, threads=threads, record_path=True,
                                  polarization=True)
    out: Arrays = {name: np.array(getattr(rays, name))
                   for name in ("dir_x", "dir_y", "dir_z", "weight", "status")}
    prt = rays.prt_matrices()
    q = rays.q_matrices()
    assert q is not None
    for row in range(3):
        for col in range(3):
            out[f"prt{row}{col}"] = np.ascontiguousarray(prt[:, row, col])
            out[f"q{row}{col}"] = np.ascontiguousarray(q[:, row, col])
    assert paths.prt is not None and paths.q is not None
    out["paths_prt"] = np.ascontiguousarray(paths.prt).ravel()
    out["paths_q"] = np.ascontiguousarray(paths.q).ravel()
    out["paths_count"] = np.array(paths.count)
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    return out


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("path", ["single", "periscope"])
def test_q_equal_cpp_bitwise(path: str, threads: int, cpp_dir: Path) -> None:
    name = f"q_{path}"
    results = flatten(path, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_have_content() -> None:
    """Rays 0 and 1 arrive, ray 2 is Missed at the first mirror; Q differs from I after the
    single reflection, the per-slot Q of slot 1 too, and the unused slots of ray 2 are NaN."""
    alive = int(rt.trace.RayStatus.ALIVE)
    missed = int(rt.trace.RayStatus.MISSED)
    single = flatten("single", 1)
    assert single["status"].tolist() == [alive, alive, missed]
    assert abs(single["q11"][0] - 1.0) > 0.5  # Householder: q11 = 1 - 2 sin^2 45 = 0
    periscope = flatten("periscope", 1)
    slots = periscope["paths_q"].size // (3 * 9)
    q = periscope["paths_q"].reshape(3, slots, 3, 3)
    assert np.max(np.abs(q[0, 1] - np.eye(3))) > 0.5
    assert np.all(np.isnan(q[2, 2:]))
