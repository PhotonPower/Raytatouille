"""M6 acceptance (#212): PSF, MTF, encircled energy, through-focus Strehl, Zernike fit and Jones
pupil of the acceptance systems from Python are bitwise equal to C++.

rtt_py_reference (m6_cases.cpp) runs the same calls on tests/reference/m6/ and r2/ideal_lens and
writes every result as named .npy arrays; the flatten_* functions below build the same arrays
from the Python results. Both sides run with 1 and with 4 threads.
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

pytestmark = [
    pytest.mark.skipif(REFERENCE_EXE is None,
                       reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)"),
    pytest.mark.filterwarnings(
        r"ignore:warning \[pupil\.mapping_nonlinear\]:raytatouille.errors.RaytatouilleWarning"),
]

Arrays = dict[str, npt.NDArray[Any]]
Systems = dict[str, rt.CompiledSystem]
RADII = [0.002, 0.004, 0.008]
SHIFTS = [0.0, 0.5, 1.0]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64).ravel()


# Flattening; names and order as the write_* functions in m6_cases.cpp.

def flatten_psf(p: an.Psf) -> Arrays:
    b = p.lattice
    return {"intensity": f8(p.intensity),
            "scalars": f8([p.strehl, p.edge_energy_deficit, b[0, 0], b[0, 1], b[1, 0], b[1, 1]])}


def flatten_mtf(m: an.Mtf) -> Arrays:
    return {"t_frequency": f8(m.tangential.frequency),
            "t_modulation": f8(m.tangential.modulation),
            "s_modulation": f8(m.sagittal.modulation)}


def flatten_jones(j: an.JonesPupil) -> Arrays:
    return {"j": np.array(j.j, dtype=np.complex128).ravel(), "w": f8(j.w),
            "diattenuation": f8(j.diattenuation), "retardance": f8(j.retardance)}


def parallel_rays() -> rt.trace.RayBatch:
    """As parallel_rays() in m6_cases.cpp: 5 x 5 rays along +z from z = 0, chief index 12."""
    rays = rt.trace.RayBatch(25)
    for i in range(25):
        x = -3.0 + 1.5 * (i % 5)
        y = -3.0 + 1.5 * (i // 5)
        rays.pos_x[i] = x
        rays.pos_y[i] = y
        rays.pupil_x[i] = x / 3.0
        rays.pupil_y[i] = y / 3.0
    return rays


def converging_rays() -> tuple[rt.trace.RayBatch, int]:
    """As converging_rays() in m6_cases.cpp: pitch 1.25 mm inside r = 5 from z = -1."""
    points = [(1.25 * ix, 1.25 * iy) for iy in range(-4, 5) for ix in range(-4, 5)
              if (1.25 * ix) ** 2 + (1.25 * iy) ** 2 <= 25.0]
    rays = rt.trace.RayBatch(len(points))
    for i, (x, y) in enumerate(points):
        rays.pos_x[i] = x
        rays.pos_y[i] = y
        rays.pos_z[i] = -1.0
        rays.pupil_x[i] = x / 5.0
        rays.pupil_y[i] = y / 5.0
    return rays, points.index((0.0, 0.0))


def systems() -> Systems:
    m6 = REFERENCE_DIR / "m6"
    return {
        "f10": rt.compile(rt.load(m6 / "paraboloid_f10.rtt.json")),
        "f50": rt.compile(rt.load(m6 / "paraboloid_f50.rtt.json")),
        "plano": rt.compile(rt.load(m6 / "planoconvex_f10.rtt.json")),
        "metal": rt.compile(rt.load(m6 / "fold45_metal.rtt.json")),
        "focus": rt.compile(rt.load(m6 / "fold45_focus.rtt.json")),
        "ideal": rt.compile(rt.load(REFERENCE_DIR / "r2" / "ideal_lens.rtt.json")),
    }


def focus_jones(s: Systems, t: int) -> an.JonesPupil:
    start, chief = converging_rays()
    return an.jones_pupil(s["focus"], start=start, chief=chief, threads=t)


# Same calls as run_m6_cases() in m6_cases.cpp.
CASES: dict[str, Callable[[Systems, int], Arrays]] = {
    "m6_f10_psf": lambda s, t: flatten_psf(an.psf(s["f10"], grid=32, threads=t)),
    "m6_f10_mtf": lambda s, t: flatten_mtf(an.mtf(s["f10"], grid=32, threads=t)),
    "m6_f10_ee": lambda s, t: {"energy": f8(an.encircled_energy(
        s["f10"], RADII, grid=32, threads=t).energy)},
    "m6_f50_through_focus": lambda s, t: {"strehl": f8(an.through_focus(
        s["f50"], shifts=SHIFTS, grid=32, threads=t).strehl)},
    "m6_planoconvex_zernike": lambda s, t: {"coefficients": f8(an.zernike_fit(
        an.pupil_function(s["plano"], grid=32, threads=t)).coefficients)},
    "m6_fold_metal_jones": lambda s, t: flatten_jones(an.jones_pupil(
        s["metal"], start=parallel_rays(), chief=12, threads=t)),
    "m6_fold_focus_jones": lambda s, t: flatten_jones(focus_jones(s, t)),
    "m6_ideal_lens_psf": lambda s, t: flatten_psf(an.psf(s["ideal"], grid=32, threads=t)),
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_m6_results_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
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
    """The compared arrays carry the acceptance content (AGENTS.md): the aberration-free PSFs
    peak with Strehl 1, the MTF falls, the EE grows, defocus lowers the Strehl ratio, the
    plano-convex lens has spherical aberration, the metal fold polarizes and the converging
    fold has an OPD and a varying diattenuation."""
    s = compiled_systems
    for name in ("m6_f10_psf", "m6_ideal_lens_psf"):
        psf = CASES[name](s, 1)
        assert psf["scalars"][0] == pytest.approx(1.0, abs=1e-9), name
        assert np.count_nonzero(psf["intensity"] > 1e-6) > 10, name
    mtf = CASES["m6_f10_mtf"](s, 1)["t_modulation"]
    assert mtf[0] == 1.0 and mtf[-1] == 0.0 and np.all(np.diff(mtf) <= 0.0)
    ee = CASES["m6_f10_ee"](s, 1)["energy"]
    assert np.all(np.diff(ee) > 0.0) and ee[-1] < 1.0
    tf = CASES["m6_f50_through_focus"](s, 1)["strehl"]
    assert tf[0] > tf[1] > tf[2]
    c = CASES["m6_planoconvex_zernike"](s, 1)["coefficients"]
    assert abs(c[10]) > 0.01  # Z11, spherical aberration
    metal = CASES["m6_fold_metal_jones"](s, 1)
    assert np.min(metal["diattenuation"]) > 1e-2 and np.min(metal["retardance"]) > 1e-2
    focus = CASES["m6_fold_focus_jones"](s, 1)
    assert np.max(np.abs(focus["w"])) > 1.0
    assert np.ptp(focus["diattenuation"]) > 1e-4
