"""Q-con and Q-bfs traced from Python are bitwise equal to C++ (#314, ADR 0043, acceptance 12).

rtt_py_reference (forbes_cases.cpp) loads qcon_singlet and qbfs_singlet under examples/forbes,
aims a hexapolar pupil at FIELDS, real and paraxial, and traces it. The start rays and the traced
columns are written as named .npy arrays. Both sides run with 1 and with 4 threads.
"""

from __future__ import annotations

import json
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
EXAMPLE_DIR = REPO_ROOT / "examples" / "forbes"
#: Field angles in degree, as in forbes_cases.cpp.
FIELDS = ((0.0, 0.0), (0.0, 2.0), (1.5, -1.0))
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "field",
           "last_surface", "status")
START = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z")


def compiled(stem: str) -> rt.CompiledSystem:
    data = json.loads((EXAMPLE_DIR / f"{stem}.rtt.json").read_text(encoding="utf-8"))
    data["fields"]["points"] = [{"x": x, "y": y} for x, y in FIELDS]
    return rt.compile(rt.System.from_json(json.dumps(data)))


def flatten(name: str, threads: int) -> Arrays:
    stem, label = CASES[name]
    aiming = Aiming.REAL if label == "real" else Aiming.PARAXIAL
    cs = compiled(stem)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(3), aiming=aiming, threads=threads)
    out: Arrays = {f"start_{column}": np.array(getattr(rays, column)) for column in START}
    stats = rt.trace.trace(cs, rays, threads=threads)
    out.update({column: np.array(getattr(rays, column)) for column in COLUMNS})
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    return out


CASES: dict[str, tuple[str, str]] = {
    f"forbes_{stem}_{label}": (stem, label)
    for stem in ("qcon_singlet", "qbfs_singlet") for label in ("real", "paraxial")
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_forbes_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = flatten(name, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


@pytest.mark.parametrize("stem", ["qcon_singlet", "qbfs_singlet"])
def test_the_forbes_cases_have_content(stem: str) -> None:
    """Every ray is aimed and arrives, the fields are tilted in y and in x, and the bundle is
    not trivial; the polynomial part matters (the two systems end in different places)."""
    alive = int(rt.trace.RayStatus.ALIVE)
    for label in ("real", "paraxial"):
        r = flatten(f"forbes_{stem}_{label}", 1)
        assert r["status"].size == 3 * 37
        assert set(r["status"].tolist()) == {alive}
        assert np.ptp(r["pos_x"]) > 0.1 and np.ptp(r["pos_y"]) > 0.1
        assert np.all(np.abs(r["start_dir_y"][r["field"] == 1]) > 0.02)
        assert np.all(np.abs(r["start_dir_x"][r["field"] == 2]) > 0.02)
    a = flatten("forbes_qcon_singlet_real", 1)
    b = flatten("forbes_qbfs_singlet_real", 1)
    assert np.max(np.abs(a["pos_y"] - b["pos_y"])) > 1e-3
