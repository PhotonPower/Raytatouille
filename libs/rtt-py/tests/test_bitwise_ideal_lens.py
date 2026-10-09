"""The ideal lens and the ideal cylinder lens traced from Python are bitwise equal to C++ (#178).

rtt_py_reference (ideal_lens_cases.cpp) traces tests/reference/r2/ideal_lens.rtt.json along
both paths with the same start rays and writes the columns as named .npy arrays. Both sides run
with 1 and with 4 threads.
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
from raytatouille.trace import RayBatch

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "weight", "status")


def start_rays(x0: float, y0: float) -> RayBatch:
    """As start_rays() in ideal_lens_cases.cpp: from (x0, y0, -75) towards (x, y, -1) of the
    stop, x, y in {-8, -4, 0, 4, 8}, normalised with sqrt(dx^2 + dy^2 + dz^2) in this order."""
    rays = RayBatch(25)
    i = 0
    for j in range(-2, 3):
        for k in range(-2, 3):
            dx = 4.0 * k - x0
            dy = 4.0 * j - y0
            dz = 74.0
            n = math.sqrt(dx * dx + dy * dy + dz * dz)
            rays.pos_x[i], rays.pos_y[i], rays.pos_z[i] = x0, y0, -75.0
            rays.dir_x[i], rays.dir_y[i], rays.dir_z[i] = dx / n, dy / n, dz / n
            i += 1
    return rays


def flatten(path: str, x0: float, y0: float, threads: int) -> Arrays:
    cs = rt.compile(rt.load(REFERENCE_DIR / "r2" / "ideal_lens.rtt.json"))
    rays = start_rays(x0, y0)
    stats = rt.trace.trace(cs, rays, path=path, threads=threads)
    out: Arrays = {name: np.array(getattr(rays, name)) for name in COLUMNS}
    prt = rays.prt_matrices()
    for row in range(3):
        for col in range(3):
            out[f"prt{row}{col}"] = np.ascontiguousarray(prt[:, row, col])
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    return out


CASES: dict[str, tuple[str, float, float]] = {
    "ideal_lens_image": ("ideal lens", 3.0, -4.0),
    "cylinder_lens_image": ("cylinder lens", 3.0, 2.0),
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_ideal_lens_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = flatten(*CASES[name], threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_have_content() -> None:
    """Every ray arrives and the lens bends it: the image of (3, -4) is (-6, 8) at z = 150, the
    line of the cylinder lens y = -4 there (m = -2)."""
    alive = int(rt.trace.RayStatus.ALIVE)
    ideal = flatten("ideal lens", 3.0, -4.0, 1)
    assert set(ideal["status"].tolist()) == {alive}
    assert np.allclose(ideal["pos_x"], -6.0, atol=1e-9) and np.allclose(ideal["pos_y"], 8.0,
                                                                        atol=1e-9)
    cyl = flatten("cylinder lens", 3.0, 2.0, 1)
    assert set(cyl["status"].tolist()) == {alive}
    assert np.allclose(cyl["pos_y"], -4.0, atol=1e-9)
    assert np.ptp(cyl["pos_x"]) > 1.0  # x is not focused
