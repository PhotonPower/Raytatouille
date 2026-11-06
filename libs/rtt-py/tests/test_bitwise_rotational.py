"""Odd asphere and axicon traced from Python are bitwise equal to C++ (#314, ADR 0043,
acceptance 12).

rtt_py_reference (rotational_cases.cpp) loads the examples under examples/rotational. The odd
asphere singlet is aimed with a hexapolar pupil at FIELDS, real and paraxial, and traced; the
axicon has no paraxial model at its apex (ADR 0043, F6), so the rays of axicon_rays() are traced.
The start rays and the traced columns are written as named .npy arrays. Both sides run with 1
and with 4 threads.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_EXE, REPO_ROOT

import raytatouille as rt
from raytatouille.trace import Aiming

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]
EXAMPLE_DIR = REPO_ROOT / "examples" / "rotational"
#: Field angles in degree, as in rotational_cases.cpp.
FIELDS = ((0.0, 0.0), (0.0, 2.0), (1.5, -1.0))
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "field",
           "last_surface", "status")
START = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z")


def axicon_rays() -> rt.trace.RayBatch:
    """As axicon_rays() in rotational_cases.cpp: 41 rays from z = -1 at heights y = 0.1 j mm,
    x = 0 for j <= 20 and x = 0.03 (j - 20) mm above; every fourth ray (j % 4 == 3) tilted along
    (0.01, -0.02, sqrt(1 - 0.0005)). Ray 0 runs along the axis through the apex."""
    rays = rt.trace.RayBatch(41)
    for j in range(41):
        rays.pos_x[j] = 0.03 * (j - 20.0) if j > 20 else 0.0
        rays.pos_y[j] = 0.1 * j
        rays.pos_z[j] = -1.0
        if j % 4 == 3:
            rays.dir_x[j] = 0.01
            rays.dir_y[j] = -0.02
            rays.dir_z[j] = math.sqrt(1.0 - 0.0005)
        else:
            rays.dir_x[j] = 0.0
            rays.dir_y[j] = 0.0
            rays.dir_z[j] = 1.0
    return rays


def odd_asphere() -> rt.CompiledSystem:
    data = json.loads((EXAMPLE_DIR / "odd_asphere_singlet.rtt.json").read_text(encoding="utf-8"))
    data["fields"]["points"] = [{"x": x, "y": y} for x, y in FIELDS]
    return rt.compile(rt.System.from_json(json.dumps(data)))


def flatten(name: str, threads: int) -> Arrays:
    if name == "rotational_axicon":
        cs = rt.compile(rt.load(EXAMPLE_DIR / "axicon_planoconvex.rtt.json"))
        rays = axicon_rays()
    else:
        aiming = Aiming.REAL if name.endswith("_real") else Aiming.PARAXIAL
        cs = odd_asphere()
        rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(3), aiming=aiming,
                                  threads=threads)
    out: Arrays = {f"start_{column}": np.array(getattr(rays, column)) for column in START}
    stats = rt.trace.trace(cs, rays, threads=threads)
    out.update({column: np.array(getattr(rays, column)) for column in COLUMNS})
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    return out


CASES = ("rotational_odd_asphere_real", "rotational_odd_asphere_paraxial", "rotational_axicon")


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", CASES)
def test_rotational_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = flatten(name, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_odd_asphere_cases_have_content() -> None:
    """Every ray is aimed and arrives, the fields are tilted, real aiming moves the start rays
    away from the paraxial ones, and the bundle is not trivial."""
    alive = int(rt.trace.RayStatus.ALIVE)
    real = flatten("rotational_odd_asphere_real", 1)
    paraxial = flatten("rotational_odd_asphere_paraxial", 1)
    for r in (real, paraxial):
        assert r["status"].size == 3 * 37
        assert set(r["status"].tolist()) == {alive}
        assert np.ptp(r["pos_x"]) > 0.1 and np.ptp(r["pos_y"]) > 0.1
    assert np.all(np.abs(real["start_dir_y"][real["field"] == 1]) > 0.02)
    assert np.max(np.abs(real["start_pos_y"] - paraxial["start_pos_y"])) > 1e-9


def test_the_axicon_case_has_content() -> None:
    """Every ray arrives; the axial ray stays on the axis, the others are bent by the cone
    towards the axis (ADR 0043, acceptance 8), so the case is not a plate."""
    alive = int(rt.trace.RayStatus.ALIVE)
    r = flatten("rotational_axicon", 1)
    assert set(r["status"].tolist()) == {alive}
    assert (r["dir_x"][0], r["dir_y"][0], r["dir_z"][0]) == (0.0, 0.0, 1.0)
    parallel = np.array([j % 4 != 3 and 0 < j <= 20 for j in range(41)])
    assert np.all(r["dir_y"][parallel] < -0.04)
    assert np.all(r["dir_x"][parallel] == 0.0)
