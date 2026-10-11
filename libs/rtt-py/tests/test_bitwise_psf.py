"""PSF, through-focus Strehl, MTF and encircled energy from Python are bitwise equal to C++
(#209, P1b).

rtt_py_reference (psf_cases.cpp) runs the same calls and writes every result as named .npy
arrays; the flatten_* functions below build the same arrays from the Python results. Both sides
run with 1 and with 4 threads. The systems: m2/paraboloid_stop at D = 10 mm (F/10) and 2 mm
(F/50) with the fields 0 and +2 deg (coma), and F/10 with the wavelengths 0.5876 and 0.7345 um.
"""

from __future__ import annotations

import json
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
    pytest.mark.filterwarnings(
        r"ignore:warning \[psf\.pupil_undersampled\]:raytatouille.errors.RaytatouilleWarning"),
]

NONE = np.iinfo(np.uint64).max
Arrays = dict[str, npt.NDArray[Any]]
Systems = dict[str, rt.CompiledSystem]
RADII = [0.002, 0.004, 0.008]
POINTS = [[0.0, 0.0], [0.001, 0.0], [0.0, 0.002], [-0.003, 0.001]]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64).ravel()


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64).ravel()


def c16(values: Any) -> npt.NDArray[np.complex128]:
    return np.array(values, dtype=np.complex128).ravel()


def losses(lo: an.RayLosses, warnings: int) -> npt.NDArray[np.uint64]:
    worst = NONE if lo.worst_surface is None else lo.worst_surface
    return u8([lo.launched, *list(lo.by_status), worst, lo.worst_surface_count, warnings])


def frame(r: Any) -> npt.NDArray[np.float64]:
    return f8([*r.centre.tolist(), *r.e_x.tolist(), *r.e_y.tolist()])


def wavelengths(w: an.PsfWavelengths) -> Arrays:
    return {
        "wavelengths": f8(np.column_stack([w.weight, w.energy, w.outside_field])),
        "wavelength_ints": u8(np.column_stack([w.wavelength, w.length])),
    }


# Flattening; names and order as the write_* functions in psf_cases.cpp.

def flatten_psf(p: an.Psf) -> Arrays:
    return {
        "intensity": f8(p.intensity), "frame": frame(p), "lattice": f8(p.lattice),
        "scalars": f8([p.strehl, p.edge_energy_deficit, p.polychromatic_scale_error,
                       p.peak.value, p.peak.k1, p.peak.k2, p.peak.x, p.peak.y]),
        "ints": u8([p.field, p.grid, p.length]), **wavelengths(p.wavelengths),
        "losses": losses(p.losses, len(p.warnings)),
    }


def flatten_through_focus(t: an.ThroughFocus) -> Arrays:
    return {
        "shifts": f8(t.shifts), "strehl": f8(t.strehl), "arrived_fraction": f8(t.arrived_fraction),
        "ints": u8([t.field, t.wavelength]), "losses": losses(t.losses, len(t.warnings)),
    }


def flatten_curve(prefix: str, c: an.MtfCurve) -> Arrays:
    return {f"{prefix}_direction": f8(c.direction), f"{prefix}_frequency": f8(c.frequency),
            f"{prefix}_modulation": f8(c.modulation), f"{prefix}_otf": c16(c.otf)}


def flatten_mtf(m: an.Mtf) -> Arrays:
    return {
        **flatten_curve("t", m.tangential), **flatten_curve("s", m.sagittal),
        "frequency_matrix": f8(m.frequency_matrix), "map": f8(m.map), "otf_map": c16(m.otf_map),
        "frame": frame(m), "scalars": f8([m.edge_energy_deficit]), "ints": u8([m.field, m.grid]),
        **wavelengths(m.wavelengths), "losses": losses(m.losses, len(m.warnings)),
    }


def flatten_ee(e: an.EncircledEnergy) -> Arrays:
    return {
        "radii": f8(e.radii), "energy": f8(e.energy), "frame": frame(e),
        "scalars": f8([e.offset.x, e.offset.y, e.radius_max, e.edge_energy_deficit]),
        "ints": u8([e.field, e.grid]), **wavelengths(e.wavelengths),
        "losses": losses(e.losses, len(e.warnings)),
    }


def paraboloid(epd: float, two_wavelengths: bool = False) -> rt.CompiledSystem:
    """As paraboloid() in psf_cases.cpp."""
    data = json.loads(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json").to_json())
    data["aperture"]["value"] = epd
    data["fields"]["points"].append({"y": 2.0})
    if two_wavelengths:
        data["wavelengths"] = [{"um": 0.5876, "reference": True}, {"um": 0.7345, "weight": 3.0}]
    return rt.compile(rt.System.from_json(json.dumps(data)))


def systems() -> Systems:
    return {"f10": paraboloid(10.0), "f50": paraboloid(2.0), "poly": paraboloid(10.0, True)}


def pupil(s: Systems, t: int) -> an.PupilFunction:
    return an.pupil_function(s["f10"], 0, 1, 0, grid=32, threads=t)


CASES: dict[str, Callable[[Systems, int], Arrays]] = {
    "f10_psf": lambda s, t: flatten_psf(an.psf(s["f10"], 0, 1, grid=32, threads=t)),
    "poly_psf": lambda s, t: flatten_psf(an.psf(s["poly"], grid=16, threads=t)),
    "f50_through_focus": lambda s, t: flatten_through_focus(an.through_focus(
        s["f50"], 0, 0, 0, shifts=[0.0, 0.5, 1.0], grid=32, threads=t)),
    "f10_mtf": lambda s, t: flatten_mtf(an.mtf(
        s["f10"], 0, 1, grid=16, map=True, complex_otf=True, threads=t)),
    "poly_mtf": lambda s, t: flatten_mtf(an.mtf(s["poly"], grid=16, threads=t)),
    "f10_ee_reference": lambda s, t: flatten_ee(an.encircled_energy(
        s["f10"], RADII, 0, 1, grid=32, threads=t)),
    "f10_ee_linear_phase": lambda s, t: flatten_ee(an.encircled_energy(
        s["f10"], RADII, 0, 1, grid=32, centre=an.EeCentre.LINEAR_PHASE, threads=t)),
    "refocus_mtf": lambda s, t: flatten_mtf(an.mtf(
        an.refocus(s["f10"], pupil(s, t), 0.05), complex_otf=True)),
    "refocus_ee": lambda s, t: flatten_ee(an.encircled_energy(
        an.refocus(s["f10"], pupil(s, t), 0.05), RADII, centre=an.EeCentre.LINEAR_PHASE)),
    "psf_at": lambda s, t: {"values": f8(an.psf_at(pupil(s, t), POINTS))},
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_psf_results_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
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
    """The compared arrays carry content (AGENTS.md): an aberrated PSF off axis, a longer
    wavelength outside the reference field, defocus lowering the Strehl ratio, an OTF phase,
    different tangential and sagittal curves, a moved LINEAR_PHASE centre and increasing EE."""
    psf = CASES["f10_psf"](compiled_systems, 1)
    assert 0.0 < psf["scalars"][0] < 0.999  # Strehl at +2 deg below 1: coma
    poly = CASES["poly_psf"](compiled_systems, 1)
    assert poly["wavelengths"][5] > 0.0  # outside_field of 0.7345 um
    assert poly["wavelength_ints"][1] != poly["wavelength_ints"][3]  # L_i differ
    tf = CASES["f50_through_focus"](compiled_systems, 1)
    assert tf["strehl"][0] > tf["strehl"][1] > tf["strehl"][2]
    mtf = CASES["f10_mtf"](compiled_systems, 1)
    assert np.max(np.abs(np.angle(mtf["t_otf"]))) > 1e-3
    assert np.max(np.abs(mtf["t_modulation"] - mtf["s_modulation"])) > 1e-4
    assert mtf["map"].size == 31 * 31
    ref = CASES["f10_ee_reference"](compiled_systems, 1)
    lin = CASES["f10_ee_linear_phase"](compiled_systems, 1)
    assert np.all(np.diff(ref["energy"]) > 0.0) and ref["energy"][-1] < 1.0
    assert lin["scalars"][1] != 0.0 and np.any(lin["energy"] != ref["energy"])
    moved = CASES["refocus_mtf"](compiled_systems, 1)
    unmoved = an.mtf(pupil(compiled_systems, 1)).tangential.modulation  # same grid 32
    assert np.max(np.abs(moved["t_modulation"] - unmoved)) > 1e-4
    values = CASES["psf_at"](compiled_systems, 1)["values"]
    assert 0.0 < values[1] < values[0] < 1.0
