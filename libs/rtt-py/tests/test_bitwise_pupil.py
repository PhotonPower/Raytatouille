"""Pupil function, refocus, Zernike fit and OPD at pupil points from Python are bitwise equal to
C++ (#209, P1a).

rtt_py_reference (pupil_cases.cpp) runs the same calls and writes every result as named .npy
arrays; the flatten_* functions below build the same arrays from the Python results. Both sides
run with 1 and with 4 threads.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import analysis as an

pytestmark = [
    pytest.mark.skipif(REFERENCE_EXE is None,
                       reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)"),
    pytest.mark.filterwarnings(
        r"ignore:warning \[pupil\.mapping_nonlinear\]:raytatouille.errors.RaytatouilleWarning"),
]

NONE = np.iinfo(np.uint64).max
Arrays = dict[str, npt.NDArray[Any]]
Systems = dict[str, rt.CompiledSystem]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def sphere(s: an.ReferenceSphere) -> npt.NDArray[np.float64]:
    return f8([*s.centre.tolist(), s.radius])


def losses(lo: an.RayLosses, warnings: int) -> npt.NDArray[np.uint64]:
    worst = NONE if lo.worst_surface is None else lo.worst_surface
    return u8([lo.launched, *list(lo.by_status), worst, lo.worst_surface_count, warnings])


# Flattening; names and order as the write_* functions in pupil_cases.cpp.

def flatten_pupil(p: an.PupilFunction) -> Arrays:
    m = p.mapping
    mapping = (f8([np.nan] * 18) if m is None else
               f8([*m.m.ravel().tolist(), *m.s0.tolist(), *m.centre.tolist(), *m.e_x.tolist(),
                   *m.e_y.tolist(), *m.d_c.tolist(), m.delta_max, m.delta_rms]))
    return {
        "fraction": np.array(p.fraction), "px": np.array(p.px), "py": np.array(p.py),
        "amplitude": np.array(p.amplitude), "w": np.array(p.w), "status": np.array(p.status),
        "prt": np.array(p.prt).ravel(), "sphere": sphere(p.sphere), "mapping": mapping,
        "scalars": f8([p.wavelength_um, p.reference_wavelength_um, p.focus_shift,
                       p.apodization, p.arrived_fraction]),
        "losses": losses(p.losses, len(p.warnings)),
    }


def flatten_zernike(f: an.ZernikeFit) -> Arrays:
    return {
        "coefficients": np.array(f.coefficients),
        "scalars": f8([f.residual_rms, f.sigma, f.sigma_bf, f.strehl, f.arrived_fraction]),
        "ints": u8([f.cells, len(f.warnings)]),
    }


def flatten_points(p: an.OpdPupilPoints) -> Arrays:
    return {
        "px": np.array(p.points.px), "py": np.array(p.points.py), "w": np.array(p.points.w),
        "status": np.array(p.points.status), "sphere": sphere(p.sphere),
        "losses": losses(p.losses, len(p.warnings)),
    }


def systems() -> Systems:
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG_DIR / "m2" / "schott.agf")
    return {
        "paraboloid": rt.compile(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json")),
        "cooke": rt.compile(rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json"), schott),
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
    }


def cooke_pupil(s: Systems, t: int) -> an.PupilFunction:
    return an.pupil_function(s["cooke"], 0, 2, 0, grid=24, with_prt=True, focus_shift=0.3,
                             threads=t)


CASES: dict[str, Callable[[Systems, int], Arrays]] = {
    "paraboloid_pupil": lambda s, t: flatten_pupil(an.pupil_function(
        s["paraboloid"], grid=24, with_prt=True, apodization=0.5, threads=t)),
    "cooke_pupil": lambda s, t: flatten_pupil(cooke_pupil(s, t)),
    "cooke_refocus": lambda s, t: flatten_pupil(an.refocus(s["cooke"], cooke_pupil(s, t), -0.2)),
    "cooke_zernike": lambda s, t: flatten_zernike(an.zernike_fit(cooke_pupil(s, t))),
    "singlet_points": lambda s, t: flatten_points(an.opd_points(
        s["singlet"], 0, 1, 0, rays=rt.trace.HexapolarPupil(4), focus_shift=0.2, threads=t)),
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_pupil_results_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
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
    """The compared arrays carry content: mixed statuses (cells outside the circle), W and PRT
    that are not trivially zero or identity, coefficients that are not all zero."""
    alive = int(rt.trace.RayStatus.ALIVE)
    for name in ("paraboloid_pupil", "cooke_pupil", "cooke_refocus"):
        a = CASES[name](compiled_systems, 1)
        assert np.any(a["status"] == alive) and np.any(a["status"] != alive), name
        assert np.any(a["prt"] != 0.0), name
        assert not np.isnan(a["mapping"][0]), name
    cooke = CASES["cooke_pupil"](compiled_systems, 1)
    assert np.max(np.abs(cooke["w"])) > 0.1  # field 20 deg, defocused: many waves
    assert cooke["scalars"][0] != cooke["scalars"][1]  # lambda != lambda_ref
    para = CASES["paraboloid_pupil"](compiled_systems, 1)
    arrived = para["status"] == alive
    assert np.min(para["amplitude"][arrived]) < 0.8  # the apodisation
    assert np.any(CASES["cooke_refocus"](compiled_systems, 1)["w"] != cooke["w"])
    z = CASES["cooke_zernike"](compiled_systems, 1)
    assert np.count_nonzero(np.abs(z["coefficients"]) > 1e-3) > 5
    points = CASES["singlet_points"](compiled_systems, 1)
    assert np.any(points["w"] != 0.0) and len(points["w"]) == 1 + 6 * (1 + 2 + 3 + 4)  # 4 rings
