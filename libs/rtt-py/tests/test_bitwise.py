"""Reference case of #32: results from Python are bitwise equal to C++.

The program rtt_py_reference (tests/rtt_py_reference.cpp) runs the same cases in C++ and
writes every RayBatch column after the trace, the trace statistics and the first-order values
as .npy files. Program and extension module come from the same CMake tree with the same flags,
so equality is exact (including NaN bit patterns and signed zeros). Both sides run with 1 and
with 4 threads; the result must not depend on it (ADR 0004, rule 7).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_EXE

import raytatouille as rt
from raytatouille.trace import Aiming, PupilSampling, RayStatus

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

COLUMNS = [
    "pos_x",
    "pos_y",
    "pos_z",
    "dir_x",
    "dir_y",
    "dir_z",
    "wl",
    "opl",
    "weight",
    "field",
    "pupil_x",
    "pupil_y",
    "last_surface",
    "status",
]

# Same order as first_order_values() in rtt_py_reference.cpp.
FIRST_ORDER_FIELDS = [
    "object_index",
    "image_index",
    "image_direction",
    "power",
    "efl",
    "front_focal_length",
    "rear_focal_length",
    "ffl",
    "bfl",
    "front_focal_z",
    "rear_focal_z",
    "front_principal_z",
    "rear_principal_z",
    "image_z",
    "lateral_magnification",
    "angular_magnification",
]


@dataclass(frozen=True)
class Case:
    name: str
    file: str
    catalog: bool
    sampling: PupilSampling | None  # None: hand_filled_rays()
    wavelength: int | None
    aiming: Aiming


# Same cases as cases() in rtt_py_reference.cpp.
CASES = [
    Case("singlet_hex", "m1/singlet_const.rtt.json", False, rt.trace.HexapolarPupil(6), None,
         Aiming.REAL),
    Case("paraboloid_random", "m2/paraboloid_stop.rtt.json", False,
         rt.trace.RandomPupil(500, 42), None, Aiming.PARAXIAL),
    Case("achromat_grid", "m2/achromat.rtt.json", True, rt.trace.GridPupil(15), 0, Aiming.REAL),
    Case("singlet_hand_filled", "m1/singlet_const.rtt.json", False, None, None, Aiming.REAL),
]


def hand_filled_rays() -> rt.trace.RayBatch:
    """Rays set through the NumPy views, as hand_filled_rays() in rtt_py_reference.cpp, so that
    status and last_surface vary: ALIVE at IMG, VIGNETTED at the stop for heights above its
    radius of 10 mm, MISSED without any surface for rays travelling towards -z."""
    rays = rt.trace.RayBatch(48)
    i = np.arange(48, dtype=np.float64)
    rays.pos_z[:] = -10.0
    # Parallel to the axis at y = 0 ... 19.5 mm, all three wavelengths.
    rays.pos_y[:40] = 0.5 * i[:40]
    rays.wl[:40] = np.arange(40) % 3
    # Oblique in the x-z plane.
    rays.pos_y[40:44] = 2.0 * (i[40:44] - 40.0)
    rays.dir_x[40:44] = 0.6
    rays.dir_z[40:44] = 0.8
    rays.wl[40:44] = 1
    # Away from the system.
    rays.pos_y[44:] = i[44:] - 44.0
    rays.dir_z[44:] = -1.0
    rays.wl[44:] = 1
    return rays


def first_order_values(fo: rt.paraxial.FirstOrder) -> npt.NDArray[np.float64]:
    def value(v: float | None) -> float:
        return math.nan if v is None else float(v)

    values = [value(getattr(fo, name)) for name in FIRST_ORDER_FIELDS]
    for pupil in (fo.entrance_pupil, fo.exit_pupil):
        values += [math.nan, math.nan] if pupil is None else [value(pupil.z),
                                                              value(pupil.diameter)]
    return np.array(values, dtype=np.float64)


def python_results(case: Case, reference_dir: Path, catalog_dir: Path,
                   threads: int) -> dict[str, npt.NDArray[np.generic]]:
    lib = rt.MaterialLibrary()
    if case.catalog:
        lib.add_catalog(catalog_dir / "schott.agf")
    cs = rt.compile(rt.load(reference_dir / case.file), lib)
    if case.sampling is None:
        rays = hand_filled_rays()
    else:
        rays = rt.trace.make_rays(cs, case.sampling, path=0, wavelength=case.wavelength,
                                  aiming=case.aiming)
    stats = rt.trace.trace(cs, rays, path=0, threads=threads)
    results: dict[str, npt.NDArray[np.generic]] = {
        name: np.array(getattr(rays, name)) for name in COLUMNS
    }
    for row in range(3):
        for col in range(3):
            results[f"prt{row}{col}"] = np.array(rays.prt(row, col))
    results["stats"] = np.array(stats.rays, dtype=np.uint64)
    wl = cs.reference_wavelength if case.wavelength is None else case.wavelength
    results["first_order"] = first_order_values(rt.paraxial.first_order(cs, 0, wl))
    return results


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_python_equals_cpp_bitwise(case: Case, threads: int, cpp_dir: Path,
                                   reference_dir: Path, catalog_dir: Path) -> None:
    results = python_results(case, reference_dir, catalog_dir, threads)
    assert len(results["pos_x"]) > 0
    if case.sampling is None:
        # The case must separate status and last_surface (review of #32); weight is constant
        # before M3.
        assert np.unique(results["status"]).size > 1
        assert np.unique(results["last_surface"]).size > 1
        assert set(np.unique(results["status"])) >= {
            int(RayStatus.ALIVE), int(RayStatus.VIGNETTED), int(RayStatus.MISSED)}
    for name, actual in results.items():
        expected = np.load(cpp_dir / f"{case.name}.{name}.npy")
        assert actual.dtype == expected.dtype, name
        assert actual.shape == expected.shape, name
        assert actual.tobytes() == expected.tobytes(), name
