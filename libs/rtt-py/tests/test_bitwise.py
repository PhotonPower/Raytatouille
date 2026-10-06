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
    sampling: PupilSampling | None  # None: rays set by hand, see `hand`
    wavelength: int | None
    aiming: Aiming
    path: int = 0
    coatings: bool = False  # load <catalog dir>/coatings (ADR 0019)
    hand: str = "singlet"  # hand-filled rays when sampling is None: singlet, plate, michelson
    paraxial: bool = True  # compare first-order values (not for tilted or stop-less systems)


# Same cases as cases() in rtt_py_reference.cpp. Since #61 every case has a non-trivial P and
# weight (Fresnel, coatings, ideal elements, absorption); the last five are the cases of #62.
CASES = [
    Case("singlet_hex", "m1/singlet_const.rtt.json", False, rt.trace.HexapolarPupil(6), None,
         Aiming.REAL),
    Case("paraboloid_random", "m2/paraboloid_stop.rtt.json", False,
         rt.trace.RandomPupil(500, 42), None, Aiming.PARAXIAL),
    Case("achromat_grid", "m2/achromat.rtt.json", True, rt.trace.GridPupil(15), 0, Aiming.REAL),
    Case("singlet_hand_filled", "m1/singlet_const.rtt.json", False, None, None, Aiming.REAL),
    Case("ar_singlet_hex", "m3/ar_singlet.rtt.json", False, rt.trace.HexapolarPupil(6), 0,
         Aiming.REAL, coatings=True),
    Case("absorbing_ar_plate_angles", "m3/absorbing_ar_plate.rtt.json", False, None, None,
         Aiming.REAL, coatings=True, hand="plate", paraxial=False),
    Case("michelson_reference_arm", "m0/michelson.rtt.json", False, None, None, Aiming.REAL,
         path=0, hand="michelson", paraxial=False),
    Case("michelson_test_arm", "m0/michelson.rtt.json", False, None, None, Aiming.REAL,
         path=1, hand="michelson", paraxial=False),
    Case("polarizer_qwp_hex", "m3/polarizer_qwp.rtt.json", False, rt.trace.HexapolarPupil(4),
         None, Aiming.REAL),
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


def plate_rays() -> rt.trace.RayBatch:
    """As plate_rays() in rtt_py_reference.cpp: 9 rays from the origin plane in the y-z plane
    at direction sines -0.4 ... 0.4, so that path length and absorption in the plate vary.
    sqrt is correctly rounded in IEEE 754 on both sides."""
    rays = rt.trace.RayBatch(9)
    for j in range(9):
        s = (j - 4) / 10.0
        rays.pos_y[j] = 0.25 * (j - 4)
        rays.dir_y[j] = s
        rays.dir_z[j] = math.sqrt(1.0 - s * s)
    return rays


def michelson_rays() -> rt.trace.RayBatch:
    """As michelson_rays() in rtt_py_reference.cpp: 8 rays at x = -1.75 ... 1.75 mm, tilted in
    the x-z plane by direction sines -0.02 ... 0.015."""
    rays = rt.trace.RayBatch(8)
    for j in range(8):
        s = (j - 4) / 200.0
        rays.pos_x[j] = 0.5 * j - 1.75
        rays.pos_y[j] = 0.125 * j
        rays.dir_x[j] = s
        rays.dir_z[j] = math.sqrt(1.0 - s * s)
    return rays


def path_results(paths: rt.trace.RayPaths) -> dict[str, npt.NDArray[np.generic]]:
    """RayPaths of the recorded trace, as write_paths() in rtt_py_reference.cpp."""
    return {
        "path_ray_indices": np.array(paths.ray_indices),
        "path_event_surfaces": np.array(paths.event_surfaces),
        "path_position": np.array(paths.position),
        "path_direction": np.array(paths.direction),
        "path_opl": np.array(paths.opl),
        "path_weight": np.array(paths.weight),
        "path_status": np.array(paths.status),
        "path_count": np.array(paths.count),
        "path_lost_at": np.array(paths.lost_at),
    }


HAND_RAYS = {"singlet": hand_filled_rays, "plate": plate_rays, "michelson": michelson_rays}


def polar_results(rays: rt.trace.RayBatch) -> dict[str, npt.NDArray[np.generic]]:
    """raytatouille.polar on the traced rays, as polar_results() in rtt_py_reference.cpp."""
    x = rt.polar.transverse_polarization(rays, np.array([1.0, 0.0, 0.0]))
    d = rt.polar.diattenuation(rays)
    r = rt.polar.retardance(rays)
    return {
        "k0": rt.polar.initial_directions(rays),
        "e_x": x,
        "transmission": rt.polar.transmission(rays),
        "transmission_x": rt.polar.transmission(rays, x),
        "d_value": d.value,
        "d_maximum": d.maximum,
        "d_minimum": d.minimum,
        "d_axis": d.axis,
        "r_value": r.value,
        "r_fast_axis": r.fast_axis,
        "stokes_x": rt.polar.stokes(rays, x, np.array([1.0, 0.0, 0.0])),
    }


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
    coatings = None
    if case.coatings:
        coatings = rt.CoatingLibrary()
        coatings.add_catalog(catalog_dir / "coatings")
    cs = rt.compile(rt.load(reference_dir / case.file), lib, coatings)
    if case.sampling is None:
        rays = HAND_RAYS[case.hand]()
    else:
        rays = rt.trace.make_rays(cs, case.sampling, path=case.path, wavelength=case.wavelength,
                                  aiming=case.aiming)
    # Recorded trace (#80): the columns must still equal the plain C++ trace bitwise.
    stats, paths = rt.trace.trace(cs, rays, path=case.path, threads=threads, record_path=True)
    results: dict[str, npt.NDArray[np.generic]] = {
        name: np.array(getattr(rays, name)) for name in COLUMNS
    }
    for row in range(3):
        for col in range(3):
            results[f"prt{row}{col}"] = np.array(rays.prt(row, col))
    results["stats"] = np.array(stats.rays, dtype=np.uint64)
    if case.paraxial:
        wl = cs.reference_wavelength if case.wavelength is None else case.wavelength
        results["first_order"] = first_order_values(rt.paraxial.first_order(cs, case.path, wl))
    results.update(polar_results(rays))
    results.update(path_results(paths))
    return results


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_python_equals_cpp_bitwise(case: Case, threads: int, cpp_dir: Path,
                                   reference_dir: Path, catalog_dir: Path) -> None:
    results = python_results(case, reference_dir, catalog_dir, threads)
    assert len(results["pos_x"]) > 0
    if case.sampling is None and case.hand == "singlet":
        # The case must separate status and last_surface (review of #32); weight carries the
        # interaction losses since #61 and is compared bit for bit like every other column.
        assert np.unique(results["status"]).size > 1
        assert np.unique(results["last_surface"]).size > 1
        assert set(np.unique(results["status"])) >= {
            int(RayStatus.ALIVE), int(RayStatus.VIGNETTED), int(RayStatus.MISSED)}
    for name, actual in results.items():
        expected = np.load(cpp_dir / f"{case.name}.{name}.npy")
        assert actual.dtype == expected.dtype, name
        assert actual.shape == expected.shape, name
        assert actual.tobytes() == expected.tobytes(), name
