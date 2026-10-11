"""The Jones pupil from Python is bitwise equal to C++ (#208 B2; ADR 0034, point 8).

rtt_py_reference (jones_cases.cpp) runs the same calls and writes every result as named .npy
arrays; flatten_jones() below builds the same arrays from the Python results. Both sides run
with 1 and with 4 threads.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import analysis as an

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None, reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)"
)

NONE = np.iinfo(np.uint64).max
Arrays = dict[str, npt.NDArray[Any]]
Systems = dict[str, rt.CompiledSystem]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def losses(lo: an.RayLosses, warnings: int) -> npt.NDArray[np.uint64]:
    worst = NONE if lo.worst_surface is None else lo.worst_surface
    return u8([lo.launched, *list(lo.by_status), worst, lo.worst_surface_count, warnings])


# Flattening; names and order as write_jones() in jones_cases.cpp.

def flatten_jones(p: an.JonesPupil) -> Arrays:
    frames = [x for v in (p.x_in, p.y_in, p.k_in, p.x_out, p.y_out, p.k_out) for x in v.tolist()]
    return {
        "j": np.array(p.j).ravel(),
        "px": np.array(p.px), "py": np.array(p.py), "w": np.array(p.w),
        "transmission": np.array(p.transmission),
        "diattenuation": np.array(p.diattenuation),
        "diattenuation_orientation": np.array(p.diattenuation_orientation),
        "diattenuation_ellipticity": np.array(p.diattenuation_ellipticity),
        "retardance": np.array(p.retardance),
        "retardance_orientation": np.array(p.retardance_orientation),
        "retardance_ellipticity": np.array(p.retardance_ellipticity),
        "status": np.array(p.status),
        "frames": f8(frames),
        "sphere": f8([*p.sphere.centre.tolist(), p.sphere.radius]),
        "scalars": f8([p.wavelength_um, p.reference_wavelength_um, p.global_phase]),
        "ints": u8([p.path, p.field, p.wavelength, int(p.basis), int(p.remove_global_phase),
                    int(p.field_amplitude)]),
        "losses": losses(p.losses, len(p.warnings)),
    }


def parallel_rays() -> rt.trace.RayBatch:
    """As parallel_rays() in jones_cases.cpp: 5 x 5 rays along +z from z = 0, chief index 12."""
    rays = rt.trace.RayBatch(25)
    for i in range(25):
        x = -3.0 + 1.5 * (i % 5)
        y = -3.0 + 1.5 * (i // 5)
        rays.pos_x[i] = x
        rays.pos_y[i] = y
        rays.pupil_x[i] = x / 3.0
        rays.pupil_y[i] = y / 3.0
    return rays


def systems() -> Systems:
    return {
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
        "periscope": rt.compile(rt.load(REFERENCE_DIR / "m3" / "metal_periscope.rtt.json")),
    }


CASES: dict[str, Callable[[Systems, int], Arrays]] = {
    "singlet_jones": lambda s, t: flatten_jones(an.jones_pupil(
        s["singlet"], 0, field=1, rays="grid:9", threads=t)),
    "singlet_jones_dipole": lambda s, t: flatten_jones(an.jones_pupil(
        s["singlet"], 0, field=1, wavelength=0, rays="grid:9", basis=an.JonesBasis.DIPOLE,
        remove_global_phase=True, threads=t)),
    "singlet_jones_field": lambda s, t: flatten_jones(an.jones_pupil(
        s["singlet"], 0, field=1, rays="grid:9", field_amplitude=True, threads=t)),
    "periscope_jones": lambda s, t: flatten_jones(an.jones_pupil(
        s["periscope"], "periscope", start=parallel_rays(), chief=12, remove_global_phase=True,
        threads=t)),
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_jones_results_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
                                         compiled_systems: Systems) -> None:
    results = CASES[name](compiled_systems, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_have_content(compiled_systems: Systems) -> None:
    """The compared arrays carry content (AGENTS.md): arrived points, J that is not trivial and
    differs between the cases, diattenuation from the refractions and the metal mirrors."""
    alive = int(rt.trace.RayStatus.ALIVE)
    cases = {name: CASES[name](compiled_systems, 1) for name in CASES}
    for name, a in cases.items():
        assert np.count_nonzero(a["status"] == alive) > 20, name
        assert np.any(np.isfinite(a["j"])) and np.any(a["j"] != 0.0), name
    assert np.nanmax(cases["singlet_jones"]["diattenuation"]) > 1e-6
    assert np.min(cases["periscope_jones"]["diattenuation"]) > 1e-2
    assert np.isnan(cases["periscope_jones"]["sphere"][3])  # the form with start rays
    assert cases["periscope_jones"]["scalars"][2] != 0.0  # a global phase was removed
    for other in ("singlet_jones_dipole", "singlet_jones_field"):
        assert cases[other]["j"].tobytes() != cases["singlet_jones"]["j"].tobytes(), other
    assert cases["singlet_jones_dipole"]["scalars"][0] != cases["singlet_jones"]["scalars"][0]
