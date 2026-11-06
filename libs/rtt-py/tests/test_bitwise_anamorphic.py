"""Anamorphic systems aimed and traced from Python are bitwise equal to C++ (#315).

rtt_py_reference (anamorphic_cases.cpp) loads the examples under examples/anamorphic (biconic,
cylinder, crossed cylinders), replaces the fields by FIELDS, aims a hexapolar pupil per section
(ADR 0042, point 3), real and paraxial, traces the rays and writes the start rays and the traced
columns as named .npy arrays. Both sides run with 1 and with 4 threads.
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
EXAMPLE_DIR = REPO_ROOT / "examples" / "anamorphic"
STEMS = ("biconic_singlet", "crossed_cylinders", "cylinder_planoconvex")
#: Built in code from the cylinder example (toroid_before_stop()), as in anamorphic_cases.cpp.
TOROID = "toroid_before_stop"
#: Field angles in degree, as in anamorphic_cases.cpp.
FIELDS = ((0.0, 0.0), (1.5, 0.0), (0.0, 1.5), (-1.0, 0.5))
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "field",
           "last_surface", "status")
START = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z")


def toroid_before_stop(data: dict[str, Any]) -> None:
    """As make_toroid_before_stop() in anamorphic_cases.cpp: the cylinder example with a toroid
    (profile R = 51.68, radius of rotation 77.52) as L1.S1 and the stop (r = 4) moved behind the
    lens to z = 20, system aperture stop_size; the entrance pupil differs per section."""
    children = data["root"]["children"]
    children[0], children[1] = children[1], children[0]  # lens, stop, detector
    children[0]["surfaces"][0]["shape"]["base"] = {
        "type": "toroid", "radius": 51.68, "radius_of_rotation": 77.52}
    children[1]["pose"] = {"position": [0.0, 0.0, 20.0]}
    children[1]["surfaces"][0]["aperture"]["radius"] = 4.0
    data["aperture"] = {"type": "stop_size"}


def compiled(stem: str) -> rt.CompiledSystem:
    file = "cylinder_planoconvex" if stem == TOROID else stem
    data = json.loads((EXAMPLE_DIR / f"{file}.rtt.json").read_text(encoding="utf-8"))
    data["fields"]["points"] = [{"x": x, "y": y} for x, y in FIELDS]
    if stem == TOROID:
        toroid_before_stop(data)
    return rt.compile(rt.System.from_json(json.dumps(data)))


def flatten(stem: str, aiming: Aiming, threads: int) -> Arrays:
    cs = compiled(stem)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(3), aiming=aiming, threads=threads)
    out: Arrays = {f"start_{name}": np.array(getattr(rays, name)) for name in START}
    stats = rt.trace.trace(cs, rays, threads=threads)
    out.update({name: np.array(getattr(rays, name)) for name in COLUMNS})
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    return out


CASES: dict[str, tuple[str, Aiming]] = {
    f"anamorphic_{stem}_{label}": (stem, aiming)
    for stem in (*STEMS, TOROID)
    for label, aiming in (("real", Aiming.REAL), ("paraxial", Aiming.PARAXIAL))
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_anamorphic_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = flatten(*CASES[name], threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


@pytest.mark.parametrize("stem", STEMS)
def test_the_cases_have_content(stem: str) -> None:
    """Every ray is aimed and arrives, the fields are tilted in x and in y, and the bundles are
    not trivial (else both sides could agree on something wrong)."""
    alive = int(rt.trace.RayStatus.ALIVE)
    for aiming in (Aiming.REAL, Aiming.PARAXIAL):
        r = flatten(stem, aiming, 1)
        assert r["status"].size == 4 * 37
        assert set(r["status"].tolist()) == {alive}
        # Field 1 is tilted in x only, field 2 in y only.
        f1 = r["field"] == 1
        f2 = r["field"] == 2
        assert np.all(np.abs(r["start_dir_x"][f1]) > 0.02)
        assert np.all(r["start_dir_y"][f1] == 0.0)
        assert np.all(np.abs(r["start_dir_y"][f2]) > 0.02)
        assert np.all(r["start_dir_x"][f2] == 0.0)
        assert np.ptp(r["pos_x"]) > 0.1 and np.ptp(r["pos_y"]) > 0.1


def test_the_toroid_case_aims_per_section() -> None:
    """Toroid in front of the stop: real aiming reaches every target (ALIVE) and moves the rays
    away from the paraxial start rays (separate pupil planes, elliptic image of the round EP in
    the stop, Newton at work), so the bitwise case covers more than the paraxial line."""
    alive = int(rt.trace.RayStatus.ALIVE)
    real = flatten(TOROID, Aiming.REAL, 1)
    paraxial = flatten(TOROID, Aiming.PARAXIAL, 1)
    assert real["status"].size == 4 * 37
    assert set(real["status"].tolist()) == {alive}
    moved = (np.abs(real["start_pos_x"] - paraxial["start_pos_x"])
             + np.abs(real["start_pos_y"] - paraxial["start_pos_y"])
             + np.abs(real["start_dir_x"] - paraxial["start_dir_x"])
             + np.abs(real["start_dir_y"] - paraxial["start_dir_y"]))
    assert np.max(moved) > 1e-6
