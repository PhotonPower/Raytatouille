"""Path evaluation (#122), ghost generator (#123) and ghost ranking (#124) from Python are
bitwise equal to C++ (#133).

rtt_py_reference (paths_cases.cpp) runs the same calls and writes every result as named .npy
arrays; the flatten_* functions below build the same arrays from the Python results. Both sides
run with 1 and with 4 threads. The Michelson cases use the start rays of the main form (a folded
path, which make_rays cannot aim); "michelson_small" has a test mirror of radius 1 mm, so rays
end vignetted and the statuses and losses are compared too.
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


# Flattening; names and order as the write_* functions in paths_cases.cpp.


def losses_values(losses: rt.analysis.RayLosses) -> list[int]:
    worst = NONE if losses.worst_surface is None else losses.worst_surface
    return [losses.launched, *losses.by_status, worst, losses.worst_surface_count]


def flatten_transmission(t: an.PathTransmission) -> Arrays:
    return {
        "px": f8(t.rays.px),
        "py": f8(t.rays.py),
        "weight": f8(t.rays.weight),
        "status": u8(t.rays.status),
        "scalars": f8([t.mean, t.min, t.max]),
        "ints": u8([t.path, t.image_surface, t.rays_launched, t.rays_arrived, len(t.warnings)]),
        "losses": u8(losses_values(t.losses)),
    }


def flatten_difference(d: an.PathOplDifference) -> Arrays:
    return {
        "px": f8(d.points.px),
        "py": f8(d.points.py),
        "delta": f8(d.points.delta),
        "status": u8(d.points.status),
        "scalars": f8([np.nan if d.chief is None else d.chief]),
        "ints": u8([d.path_a, d.path_b, d.image_surface, len(d.warnings)]),
        "losses_a": u8(losses_values(d.losses_a)),
        "losses_b": u8(losses_values(d.losses_b)),
    }


def flatten_ghosts(g: rt.GhostSystem) -> Arrays:
    return {
        "path": u8(g.ghosts.path),
        "base": u8(g.ghosts.base),
        "surface_j": u8(g.ghosts.surface_j),
        "surface_i": u8(g.ghosts.surface_i),
        "event_j": u8(g.ghosts.event_j),
        "event_i": u8(g.ghosts.event_i),
        "ints": u8([len(g.system.path_names)]),
    }


def flatten_ranking(r: an.GhostRanking) -> Arrays:
    e = r.entries
    return {
        "path": u8(e.path),
        "surface_j": u8(e.surface_j),
        "surface_i": u8(e.surface_i),
        "power": f8(e.power),
        "relative_power": f8(e.relative_power),
        "rms_radius": f8(e.rms_radius),
        "rays_arrived": u8(e.rays_arrived),
        "relative_irradiance": f8(e.relative_irradiance),
        "focus_offset": f8(e.focus_offset),
        "paraxial_blur_radius": f8(e.paraxial_blur_radius),
        "losses": u8([n for losses in e.losses for n in losses_values(losses)]),
        "scalars": f8([r.base_power, r.base_rms_radius, r.resolution_radius]),
        "ints": u8([r.base, r.field, r.wavelength, len(r.warnings)]),
    }


def collimated_bundle() -> rt.trace.RayBatch:
    """As collimated_bundle() in paths_cases.cpp: the axis ray first, then a square grid with a
    pitch of 0.5 mm inside r <= 2.25 mm, along +z at z = 0; pupil coordinates r / 2.5."""
    points = [(0.0, 0.0)]
    for i in range(-4, 5):
        for j in range(-4, 5):
            x, y = 0.5 * i, 0.5 * j
            if (i, j) != (0, 0) and x * x + y * y <= 2.25 * 2.25:
                points.append((x, y))
    rays = rt.trace.RayBatch(len(points))
    for k, (x, y) in enumerate(points):
        rays.pos_x[k] = x
        rays.pos_y[k] = y
        rays.pupil_x[k] = x / 2.5
        rays.pupil_y[k] = y / 2.5
    return rays


Systems = dict[str, Any]


def systems() -> Systems:
    """The systems of load_systems() in paths_cases.cpp."""
    data = json.loads(rt.load(REFERENCE_DIR / "m4" / "michelson_offset.rtt.json").to_json())
    data["root"]["children"][3]["surfaces"][0]["aperture"] = {"type": "circular", "radius": 1.0}
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG_DIR / "m2" / "schott.agf")
    schott_bk7 = rt.MaterialLibrary()  # N-BK7 of the feature tour
    schott_bk7.add_catalog(CATALOG_DIR / "schott.agf")
    # The feature tour shows every format element and does not compile as it is: without the
    # Zernike term of A.S1 (M8) and with a bare Fresnel A.S1 instead of the coating "AR_VIS"
    # (no catalog), as in paths_cases.cpp.
    tour = json.loads(rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json").to_json())
    a_s1 = tour["root"]["children"][2]["surfaces"][0]
    del a_s1["shape"]["terms"]
    a_s1["interaction"] = {"type": "fresnel"}
    return {
        "michelson": rt.compile(rt.load(REFERENCE_DIR / "m4" / "michelson_offset.rtt.json")),
        "michelson_small": rt.compile(rt.System.from_json(json.dumps(data))),
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
        "plate": rt.compile_with_ghosts(rt.load(REFERENCE_DIR / "m3" / "fresnel_bk7.rtt.json"),
                                        "main"),
        "cooke": rt.compile_with_ghosts(rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json"),
                                        "main", materials=schott),
        # M4 acceptance (#135)
        "ghost_plates": rt.compile_with_ghosts(
            rt.load(REFERENCE_DIR / "m4" / "ghost_plates.rtt.json"), "main"),
        "ghost_singlet": rt.compile_with_ghosts(
            rt.load(REFERENCE_DIR / "m4" / "ghost_singlet.rtt.json"), "main"),
        "tour": rt.compile(rt.System.from_json(json.dumps(tour)), materials=schott_bk7),
    }


ARMS = ("reference arm", "test arm", "reference arm, return", "test arm, return")

Case = Callable[[Systems, int], Arrays]


def _arm(k: int) -> Case:
    return lambda s, t: flatten_transmission(
        an.path_transmission(s["michelson"], ARMS[k], start=collimated_bundle(), threads=t))


# Same calls as run_paths_cases() in paths_cases.cpp.
CASES: dict[str, Case] = {
    **{f"michelson_t{k}": _arm(k) for k in range(4)},
    "michelson_opl": lambda s, t: flatten_difference(an.opl_difference(
        s["michelson"], "reference arm", "test arm", start=collimated_bundle(), threads=t)),
    "michelson_small_t": lambda s, t: flatten_transmission(an.path_transmission(
        s["michelson_small"], "test arm", start=collimated_bundle(), threads=t)),
    "michelson_small_opl": lambda s, t: flatten_difference(an.opl_difference(
        s["michelson_small"], "test arm", "reference arm", start=collimated_bundle(),
        threads=t)),
    "singlet_t": lambda s, t: flatten_transmission(an.path_transmission(
        s["singlet"], 0, field=1, wavelength=0, rays="hexapolar:6", threads=t)),
    "singlet_opl": lambda s, t: flatten_difference(an.opl_difference(
        s["singlet"], 0, 0, field=2, wavelength=0, rays="grid:9", aiming=Aiming.PARAXIAL,
        threads=t)),
    "plate_ghosts": lambda s, t: flatten_ghosts(s["plate"]),
    "plate_ranking": lambda s, t: flatten_ranking(an.ghost_ranking(
        s["plate"], 0, rays="hexapolar:4", threads=t)),
    "cooke_ghosts": lambda s, t: flatten_ghosts(s["cooke"]),
    "cooke_ranking": lambda s, t: flatten_ranking(an.ghost_ranking(
        s["cooke"], 1, resolution_radius=0.01, threads=t)),
    # M4 acceptance (#135)
    "ghost_plates_ghosts": lambda s, t: flatten_ghosts(s["ghost_plates"]),
    "ghost_plates_ranking": lambda s, t: flatten_ranking(an.ghost_ranking(
        s["ghost_plates"], 0, threads=t)),
    "ghost_singlet_ghosts": lambda s, t: flatten_ghosts(s["ghost_singlet"]),
    "ghost_singlet_ranking": lambda s, t: flatten_ranking(an.ghost_ranking(
        s["ghost_singlet"], 0, threads=t)),
    "tour_first_order": lambda s, t: flatten_transmission(an.path_transmission(
        s["tour"], "first order", start=collimated_bundle(), threads=t)),
}


@pytest.fixture(scope="module")
def compiled_systems() -> Systems:
    return systems()


# vignetted on purpose: the rays.lost warning is expected; the annular stop of the feature tour
# clips the collimated bundle (stop.clips_beam, #135)
@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
@pytest.mark.filterwarnings(
    r"ignore:warning \[stop\.clips_beam\]:raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_paths_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path,
                                 compiled_systems: Systems) -> None:
    results = CASES[name](compiled_systems, threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


# vignetted on purpose: the rays.lost warning is expected; the annular stop of the feature tour
# clips the collimated bundle (stop.clips_beam, #135)
@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
@pytest.mark.filterwarnings(
    r"ignore:warning \[stop\.clips_beam\]:raytatouille.errors.RaytatouilleWarning")
def test_the_cases_separate_status_and_values(compiled_systems: Systems) -> None:
    # The comparison must see lost rays, non-trivial deltas and several ghosts (lesson of #51).
    small = CASES["michelson_small_t"](compiled_systems, 1)
    statuses = set(np.unique(small["status"]).tolist())
    assert {int(RayStatus.ALIVE), int(RayStatus.VIGNETTED)} <= statuses
    opl = CASES["michelson_small_opl"](compiled_systems, 1)
    assert set(np.unique(opl["status"]).tolist()) == statuses
    # -15 mm (to rounding) where the ray arrived on both paths, 0 where it was lost.
    assert np.unique(np.round(opl["delta"], 9)).tolist() == [-15.0, 0.0]
    cooke = CASES["cooke_ranking"](compiled_systems, 1)
    assert cooke["path"].size == 15
    assert np.unique(cooke["relative_irradiance"]).size > 1
    # M4 (#135): six plate ghosts with different rho, a singlet ghost with a paraxial focus, and
    # the feature tour with lost rays (annular stop) next to arrived ones.
    plates = CASES["ghost_plates_ranking"](compiled_systems, 1)
    assert plates["path"].size == 6 and np.unique(plates["relative_irradiance"]).size == 6
    singlet = CASES["ghost_singlet_ranking"](compiled_systems, 1)
    assert np.isfinite(singlet["focus_offset"]).all()
    tour = CASES["tour_first_order"](compiled_systems, 1)
    assert {int(RayStatus.ALIVE), int(RayStatus.VIGNETTED)} <= set(np.unique(tour["status"]).tolist())
    assert np.unique(tour["weight"][tour["status"] == int(RayStatus.ALIVE)]).size > 1

