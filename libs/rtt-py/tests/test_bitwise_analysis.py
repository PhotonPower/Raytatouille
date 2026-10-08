"""Acceptance criterion of #33: every analysis from Python is bitwise equal to C++.

rtt_py_reference (analysis_cases.cpp) runs the same calls and writes every result as named
.npy arrays; flatten() below builds the same arrays from the Python results. Both sides run
with 1 and with 4 threads. Three cases use the paraboloid with its mirror aperture reduced to
21 mm (beam radius 30 mm), so vignetted rays and statuses other than ALIVE are compared too.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille.trace import Aiming, RayStatus

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


# Flattening; names and order as the write_* functions in analysis_cases.cpp.


def flatten_spot(s: an.SpotDiagram) -> Arrays:
    st = s.stats
    return {
        "x": f8(s.x),
        "y": f8(s.y),
        "wavelengths": u8(s.wavelengths),
        "weight": f8(s.weight),
        "scalars": f8([s.chief.x, s.chief.y, st.centroid.x, st.centroid.y, st.rms_centroid,
                       st.rms_chief, st.geo_centroid, st.geo_chief, s.vignetted_fraction]),
        "ints": u8([s.field, NONE if s.wavelength is None else s.wavelength, s.image_surface,
                    s.rays_launched, s.rays_arrived]),
    }


def flatten_fan_points(prefix: str, p: an.FanPoints) -> Arrays:
    return {prefix + "p": f8(p.p), prefix + "ex": f8(p.ex), prefix + "ey": f8(p.ey),
            prefix + "status": u8(p.status)}


def flatten_fan(f: an.RayFan) -> Arrays:
    return {**flatten_fan_points("t_", f.tangential), **flatten_fan_points("s_", f.sagittal),
            "scalars": f8([f.chief.x, f.chief.y]),
            "ints": u8([f.field, f.wavelength, f.image_surface])}


def flatten_opd_points(prefix: str, p: an.OpdPoints) -> Arrays:
    return {prefix + "px": f8(p.px), prefix + "py": f8(p.py), prefix + "w": f8(p.w),
            prefix + "status": u8(p.status)}


def sphere_values(s: an.ReferenceSphere) -> list[float]:
    return [*(float(v) for v in s.centre), s.radius]


def flatten_opd_map(m: an.OpdMap) -> Arrays:
    return {**flatten_opd_points("", m.points),
            "scalars": f8([*sphere_values(m.sphere), m.rms, m.pv]),
            "ints": u8([m.field, m.wavelength, m.arrived, m.vignetted])}


def flatten_opd_fan(f: an.OpdFan) -> Arrays:
    return {**flatten_opd_points("t_", f.tangential), **flatten_opd_points("s_", f.sagittal),
            "scalars": f8(sphere_values(f.sphere)), "ints": u8([f.field, f.wavelength])}


def flatten_longitudinal(lc: an.LongitudinalColour) -> Arrays:
    return {"f_wavelength": u8(lc.foci.wavelength), "f_paraxial_z": f8(lc.foci.paraxial_z),
            "f_real_z": f8(lc.foci.real_z), "scalars": f8([lc.paraxial, lc.real]),
            "ints": u8([lc.pair.first, lc.pair.second])}


def flatten_lateral(lc: an.LateralColour) -> Arrays:
    return {"chief_x": f8(lc.chief.x), "chief_y": f8(lc.chief.y),
            "offset_x": f8(lc.offset.x), "offset_y": f8(lc.offset.y), "ints": u8([lc.field])}


def distortion_values(p: an.DistortionPoint) -> list[float]:
    return [p.fraction, p.field.x, p.field.y, p.real_height, p.paraxial_height, p.percent]


def field_curvature_values(p: an.FieldCurvaturePoint) -> list[float]:
    return [p.fraction, p.field.x, p.field.y, p.tangential, p.sagittal, p.astigmatism]


def flatten_distortion(d: an.DistortionSweep) -> Arrays:
    columns = [d.fraction, d.field_x, d.field_y, d.real_height, d.paraxial_height, d.percent]
    return {f"q{i}": f8(c) for i, c in enumerate(columns)}


def flatten_field_curvature(d: an.FieldCurvatureSweep) -> Arrays:
    columns = [d.fraction, d.field_x, d.field_y, d.tangential, d.sagittal, d.astigmatism]
    return {f"q{i}": f8(c) for i, c in enumerate(columns)}


def flatten_seidel(s: rt.paraxial.Seidel) -> Arrays:
    q = s.surfaces
    terms = s.sum
    pair = s.chromatic
    return {
        "surface": u8(q.surface), "y": f8(q.y), "y_bar": f8(q.y_bar), "a": f8(q.a),
        "a_bar": f8(q.a_bar), "lagrange": f8(q.lagrange), "s1": f8(q.s1), "s2": f8(q.s2),
        "s3": f8(q.s3), "s4": f8(q.s4), "s5": f8(q.s5), "c_l": f8(q.c_l), "c_t": f8(q.c_t),
        "scalars": f8([terms.s1, terms.s2, terms.s3, terms.s4, terms.s5, terms.c_l, terms.c_t,
                       s.lagrange, s.marginal.z, s.marginal.y, s.marginal.u, s.chief.z,
                       s.chief.y, s.chief.u]),
        "ints": u8([NONE, NONE] if pair is None else [pair.first, pair.second]),
    }


def flatten_prescription(p: rt.paraxial.Prescription) -> Arrays:
    q = p.surfaces

    def value(v: float | None) -> float:
        return float("nan") if v is None else v

    return {
        "surface": u8(q.surface), "z": f8(q.z), "n": f8(q.n), "y": f8(q.y), "u": f8(q.u),
        "i": f8(q.i), "y_bar": f8(q.y_bar), "u_bar": f8(q.u_bar), "i_bar": f8(q.i_bar),
        "lagrange": f8(q.lagrange),
        "scalars": f8([p.total_track, value(p.object_distance),
                       value(p.paraxial_working_f_number), value(p.paraxial_image_na),
                       value(p.lagrange_invariant)]),
    }


def systems() -> dict[str, rt.CompiledSystem]:
    """The systems of load_systems() in analysis_cases.cpp."""
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG_DIR / "schott.agf")
    data = json.loads(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json").to_json())
    data["root"]["children"][1]["surfaces"][0]["aperture"] = {"type": "circular", "radius": 21.0}
    return {
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
        "achromat": rt.compile(rt.load(REFERENCE_DIR / "m2" / "achromat.rtt.json"), schott),
        "paraboloid21": rt.compile(rt.System.from_json(json.dumps(data))),
    }


Case = Callable[[dict[str, rt.CompiledSystem], int], Arrays]

# Same calls as run_analysis_cases() in analysis_cases.cpp (wavelength None = reference,
# except for spot).
CASES: dict[str, Case] = {
    "spot_mono": lambda s, t: flatten_spot(
        an.spot(s["singlet"], 0, 2, 0, rays="hexapolar:6", threads=t)),
    "spot_poly": lambda s, t: flatten_spot(
        an.spot(s["achromat"], 0, 2, None, rays="random:300:7", threads=t)),
    "spot_vignetted": lambda s, t: flatten_spot(
        an.spot(s["paraboloid21"], 0, 0, None, rays="grid:21", threads=t)),
    "fan": lambda s, t: flatten_fan(an.ray_fan(s["singlet"], 0, 1, 0, points=21, threads=t)),
    "fan_vignetted": lambda s, t: flatten_fan(
        an.ray_fan(s["paraboloid21"], 0, 0, points=41, aiming=Aiming.PARAXIAL, threads=t)),
    "opd_map": lambda s, t: flatten_opd_map(an.opd_map(s["singlet"], 0, 2, grid=17, threads=t)),
    "opd_map_vignetted": lambda s, t: flatten_opd_map(
        an.opd_map(s["paraboloid21"], 0, 0, grid=21, threads=t)),
    "opd_fan": lambda s, t: flatten_opd_fan(an.opd_fan(s["achromat"], 0, 1, 2, threads=t)),
    "longitudinal": lambda s, t: flatten_longitudinal(
        an.longitudinal_colour(s["achromat"], zone=0.7, threads=t)),
    "lateral": lambda s, t: flatten_lateral(an.lateral_colour(s["achromat"], 0, 2, threads=t)),
    "distortion": lambda s, t: flatten_distortion(
        an.distortion(s["singlet"], samples=11, threads=t)),
    "distortion_at": lambda s, t: {"scalars": f8(distortion_values(
        an.distortion_at(s["singlet"], (0.0, 4.0), wavelength=0, threads=t)))},
    "field_curvature": lambda s, t: flatten_field_curvature(
        an.field_curvature(s["achromat"], samples=7, delta=1e-3, threads=t)),
    "field_curvature_at": lambda s, t: {"scalars": f8(field_curvature_values(
        an.field_curvature_at(s["achromat"], rt.Field(1.0, 3.0), delta=2e-3, threads=t)))},
    "seidel": lambda s, t: flatten_seidel(an.seidel(s["achromat"], pair=(0, 2))),
    "prescription": lambda s, t: flatten_prescription(rt.paraxial.prescription(s["achromat"])),
}


@pytest.fixture(scope="module")
def compiled_systems() -> dict[str, rt.CompiledSystem]:
    return systems()


# vignetted on purpose: the rays.lost warning is expected
@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_analysis_equals_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
                                     compiled_systems: dict[str, rt.CompiledSystem]) -> None:
    results = CASES[name](compiled_systems, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


# vignetted on purpose: the rays.lost warning is expected
@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
def test_vignetting_cases_separate_status(compiled_systems: dict[str, rt.CompiledSystem]) -> None:
    # The comparison must see rays that do not arrive (lesson of #51).
    spot = CASES["spot_vignetted"](compiled_systems, 1)
    assert 0.0 < spot["scalars"][-1] < 1.0  # vignetted fraction
    fan = CASES["fan_vignetted"](compiled_systems, 1)
    statuses = set(np.unique(np.concatenate([fan["t_status"], fan["s_status"]])).tolist())
    assert {int(RayStatus.ALIVE), int(RayStatus.VIGNETTED)} <= statuses
    opd = CASES["opd_map_vignetted"](compiled_systems, 1)
    assert len(np.unique(opd["status"])) > 1 and opd["ints"][3] > 0
