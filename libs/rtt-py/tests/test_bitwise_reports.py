"""Reports from Python are bitwise equal to C++ (#177).

rtt_py_reference (reports_cases.cpp) runs the same calls and writes every report as named .npy
arrays; the flatten_* functions below build the same arrays from the Python reports. Both sides
run with 1 and with 4 threads (the raytrace report traces; the other reports do not).
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

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

NONE = np.iinfo(np.uint64).max
Arrays = dict[str, npt.NDArray[Any]]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def start_rays() -> rt.trace.RayBatch:
    """As start_rays() in reports_cases.cpp: an axial ray, a ray at y = 5 mm and a steep ray from
    y = 15 mm (lost: VIGNETTED at the singlet's stop, MISSED in the Cooke triplet), all from
    z = -1 (in front of the first surface at z = 0)."""
    rays = rt.trace.RayBatch(3)
    for k in range(3):
        rays.pos_z[k] = -1.0
    rays.pos_y[1] = 5.0
    rays.pos_y[2] = 15.0
    rays.dir_y[2] = 0.8
    rays.dir_z[2] = 0.6
    return rays


# Flattening; names and order as the write_* functions in reports_cases.cpp.

RAYTRACE_ARRAYS = ("ray", "slot", "surface", "x", "y", "z", "dx", "dy", "dz", "local_x",
                   "local_y", "local_z", "local_dx", "local_dy", "local_dz", "opl", "weight",
                   "status")
DIMENSION_ARRAYS = ("element", "first_surface", "coaxial", "centre_thickness",
                    "semi_diameter_first", "semi_diameter_second", "aperture_first",
                    "aperture_second", "edge_thickness", "diameter")


def flatten_raytrace(r: an.RaytraceReport) -> Arrays:
    out: Arrays = {name: np.array(getattr(r.rows, name)) for name in RAYTRACE_ARRAYS}
    out["ints"] = u8([r.path, r.rays, r.slots])
    return out


def _value(v: float | None) -> float:
    return np.nan if v is None else v


def flatten_system(r: an.SystemReport) -> Arrays:
    p = r.prescription
    out: Arrays = {
        "ints": u8([r.path, r.wavelength, r.reference_wavelength, r.field_count, r.surface_count,
                    r.event_count, NONE if r.stop is None else r.stop, int(p is not None),
                    len(r.warnings)]),
        "wavelengths_um": f8(r.wavelengths_um),
    }
    if p is None:
        return out
    s = p.surfaces
    out.update({
        "p_surface": u8(s.surface), "p_z": f8(s.z), "p_n": f8(s.n), "p_y": f8(s.y),
        "p_u": f8(s.u), "p_i": f8(s.i), "p_y_bar": f8(s.y_bar), "p_u_bar": f8(s.u_bar),
        "p_i_bar": f8(s.i_bar), "p_lagrange": f8(s.lagrange),
        "p_scalars": f8([p.total_track, _value(p.object_distance),
                         _value(p.paraxial_working_f_number), _value(p.paraxial_image_na),
                         _value(p.lagrange_invariant), _value(p.first_order.efl),
                         _value(p.first_order.bfl)]),
    })
    return out


def flatten_dimensions(d: an.DimensionReport) -> Arrays:
    return {name: np.array(getattr(d.segments, name)) for name in DIMENSION_ARRAYS}


Systems = dict[str, rt.CompiledSystem]


def systems() -> Systems:
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG_DIR / "m2" / "schott.agf")  # with N-LAK9 and N-SF5
    return {
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
        "cooke": rt.compile(rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json"),
                            materials=schott),
        "grating": rt.compile(rt.load(REFERENCE_DIR / "m4" / "grating_transmission.rtt.json")),
        # relative placement and configurations (#169 A): the gap of the configuration "tele"
        "zoom_tele": rt.compile(rt.load(REFERENCE_DIR / "m5" / "zoom.rtt.json"),
                                configuration="tele"),
    }


Case = Callable[[Systems, int], Arrays]

# Same calls as run_reports_cases() in reports_cases.cpp.
CASES: dict[str, Case] = {
    "singlet_raytrace": lambda s, t: flatten_raytrace(an.raytrace_report(
        s["singlet"], 0, start=start_rays(), threads=t)),
    "cooke_raytrace": lambda s, t: flatten_raytrace(an.raytrace_report(
        s["cooke"], 0, start=start_rays(), threads=t)),
    "singlet_system": lambda s, t: flatten_system(an.system_report(s["singlet"], 0, 1)),
    "cooke_system": lambda s, t: flatten_system(an.system_report(s["cooke"], 0)),
    "grating_system": lambda s, t: flatten_system(an.system_report(s["grating"], "order +1", 0)),
    "singlet_dimensions": lambda s, t: flatten_dimensions(an.dimension_report(s["singlet"])),
    "cooke_dimensions": lambda s, t: flatten_dimensions(an.dimension_report(s["cooke"])),
    "zoom_tele_dimensions": lambda s, t: flatten_dimensions(an.dimension_report(s["zoom_tele"])),
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


# The grating order has no paraxial data on purpose: only that case filters the warning.
NO_PARAXIAL = pytest.mark.filterwarnings(
    r"ignore:warning \[report\.paraxial_unavailable\]:raytatouille.errors.RaytatouilleWarning")
CASE_PARAMS = [pytest.param(name, marks=NO_PARAXIAL) if name == "grating_system" else name
               for name in CASES]


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", CASE_PARAMS)
def test_reports_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
                                   compiled_systems: Systems) -> None:
    results = CASES[name](compiled_systems, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


@NO_PARAXIAL
def test_the_cases_have_content(compiled_systems: Systems) -> None:
    """The compared arrays carry content (review of #192): a comparison of rays that all end at
    the first surface, or of NaN dimensions, would also be bitwise equal."""
    status = rt.trace.RayStatus
    lost = {"singlet_raytrace": status.VIGNETTED, "cooke_raytrace": status.MISSED}
    for name, end_of_ray_2 in lost.items():
        a = CASES[name](compiled_systems, 1)
        slots = int(a["ints"][2])
        # The last row of each ray (rows are ray-major): rays 0 and 1 arrive ALIVE in the last
        # slot, ray 2 ends lost as start_rays() says.
        last = {int(r): i for i, r in enumerate(a["ray"])}
        assert sorted(last) == [0, 1, 2], name
        for ray in (0, 1):
            assert (int(a["slot"][last[ray]]), int(a["status"][last[ray]])) == (
                slots - 1, int(status.ALIVE)), (name, ray)
        assert int(a["status"][last[2]]) == int(end_of_ray_2), name
        assert int(a["slot"][last[2]]) < slots - 1, name
    # The prescriptions exist except for the grating order.
    assert CASES["singlet_system"](compiled_systems, 1)["ints"][7] == 1
    assert CASES["cooke_system"](compiled_systems, 1)["ints"][7] == 1
    assert CASES["grating_system"](compiled_systems, 1)["ints"][7] == 0
    # Every segment coaxial with a finite positive centre thickness. The singlet and the zoom
    # have apertures, so finite positive edge thicknesses and diameters; the Cooke triplet has
    # none, so its edge thicknesses and diameters are NaN by definition.
    for name, count in (("singlet_dimensions", 1), ("cooke_dimensions", 3),
                        ("zoom_tele_dimensions", 2)):
        d = CASES[name](compiled_systems, 1)
        assert len(d["element"]) == count, name
        assert bool(np.all(d["coaxial"])), name
        positive = ["centre_thickness"]
        if name != "cooke_dimensions":
            positive += ["edge_thickness", "diameter"]
        for key in positive:
            assert bool(np.all(np.isfinite(d[key]) & (d[key] > 0.0))), (name, key)
    assert bool(np.all(np.isnan(CASES["cooke_dimensions"](compiled_systems, 1)["diameter"])))
