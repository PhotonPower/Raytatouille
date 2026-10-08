"""Path evaluation (#122), ghost generator (#123) and ghost ranking (#124) from Python (#133):
the two forms of the path functions, ghost names and events, warnings, errors and run control.
Bitwise equality with C++ is in test_bitwise_paths.py; the reference values here are analytic.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from conftest import REFERENCE_DIR
from test_bitwise_paths import collimated_bundle

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille.model import EventKind
from raytatouille.trace import Aiming, RayStatus

#: Fresnel at normal incidence, n = 1.5168 in VACUUM (Byrnes, arXiv:1603.02720v5, Eq. (6);
#: T = 1 - R for a lossless interface, Eqs. (21)-(23); docs/quellen.md).
N_PLATE = 1.5168
R_PLATE = ((N_PLATE - 1.0) / (N_PLATE + 1.0)) ** 2
T_PLATE = 1.0 - R_PLATE


@pytest.fixture(scope="module")
def michelson() -> rt.CompiledSystem:
    return rt.compile(rt.load(REFERENCE_DIR / "m4" / "michelson_offset.rtt.json"))


@pytest.fixture(scope="module")
def plate() -> rt.System:
    return rt.load(REFERENCE_DIR / "m3" / "fresnel_bk7.rtt.json")


def test_michelson_main_form(michelson: rt.CompiledSystem) -> None:
    # #122: the four paths carry R T, T R, T T, R R = 0.25 each (ideal splitter, R = 0.5), so
    # every start ray adds up to (R + T)^2 = 1; the arms differ by 2 Delta = 15 mm (VACUUM).
    start = collimated_bundle()
    columns = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "weight", "status",
               "last_surface")
    before = {name: np.array(getattr(start, name)) for name in columns}
    total = np.zeros(len(start))
    for name in ("reference arm", "test arm", "reference arm, return", "test arm, return"):
        t = an.path_transmission(michelson, name, start=start)
        assert (t.rays_launched, t.rays_arrived, len(t.rays)) == (len(start),) * 3
        assert t.path == michelson.find_path(name)
        assert math.isclose(t.mean, 0.25, abs_tol=1e-12)
        assert np.all(t.rays.status == int(RayStatus.ALIVE))
        np.testing.assert_array_equal(t.rays.px, start.pupil_x)
        total += t.rays.weight
    np.testing.assert_allclose(total, 1.0, rtol=0.0, atol=1e-12)
    d = an.opl_difference(michelson, "reference arm", "test arm", start=start)
    assert d.chief is not None and math.isclose(d.chief, 15.0, abs_tol=1e-10)
    np.testing.assert_allclose(d.points.delta, 15.0, rtol=0.0, atol=1e-10)
    for name, values in before.items():  # the start rays are copied, not traced
        np.testing.assert_array_equal(getattr(start, name), values)


def test_the_forms_do_not_mix(michelson: rt.CompiledSystem) -> None:
    start = collimated_bundle()
    for extra in ({"field": 0}, {"wavelength": 0}, {"rays": "hexapolar:2"},
                  {"aiming": Aiming.REAL}):
        with pytest.raises(ValueError, match="only apply without start rays"):
            an.path_transmission(michelson, 0, start=start, **extra)
        with pytest.raises(ValueError, match="only apply without start rays"):
            an.opl_difference(michelson, 0, 1, start=start, **extra)


def test_convenience_form_and_its_limits(reference_dir: Path,
                                         michelson: rt.CompiledSystem) -> None:
    # make_rays on a symmetric path: the same rays as the main form with make_rays (#122).
    singlet = rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))
    convenient = an.path_transmission(singlet, "main", field=1, rays="hexapolar:3")
    start = rt.trace.make_rays(singlet, rt.trace.HexapolarPupil(rings=3), fields=[1])
    main = an.path_transmission(singlet, "main", start=start)
    assert convenient.rays.weight.tobytes() == main.rays.weight.tobytes()
    assert convenient.mean == main.mean
    # The folded Michelson cannot be aimed paraxially: use start rays there.
    with pytest.raises(rt.ParaxialError):
        an.path_transmission(michelson, "test arm")
    # Different image surfaces: no OPL difference.
    with pytest.raises(ValueError):
        an.opl_difference(michelson, "reference arm", "test arm, return",
                          start=collimated_bundle())


def test_lost_rays_warn_and_count(reference_dir: Path) -> None:
    data = json.loads(rt.load(reference_dir / "m4" / "michelson_offset.rtt.json").to_json())
    data["root"]["children"][3]["surfaces"][0]["aperture"] = {"type": "circular", "radius": 1.0}
    small = rt.compile(rt.System.from_json(json.dumps(data)))
    with pytest.warns(rt.RaytatouilleWarning, match=r"rays\.lost"):
        t = an.path_transmission(small, "test arm", start=collimated_bundle())
    assert [w.code for w in t.warnings] == ["rays.lost"]
    assert 0 < t.rays_arrived < t.rays_launched
    assert t.losses.count(RayStatus.VIGNETTED) == t.rays_launched - t.rays_arrived
    lost = t.rays.status != int(RayStatus.ALIVE)
    assert np.all(t.rays.weight[lost] == 0.0)
    assert t.min == pytest.approx(0.25, abs=1e-12) and t.max == pytest.approx(0.25, abs=1e-12)


def test_ghost_paths_of_the_plate(plate: rt.System) -> None:
    # ADR 0027: one ghost, base[0..j-1], Reflect at P.S2, Reflect at P.S1, base[i+1..end].
    compiled = rt.compile(plate)
    [ghost] = rt.ghost_paths(compiled, "main")
    assert ghost.name == "main ghost P.S2/P.S1"
    assert not ghost.automatic
    assert [(e.surface, e.kind) for e in ghost.events] == [
        ("STO", EventKind.TRANSMIT), ("P.S1", EventKind.REFRACT), ("P.S2", EventKind.REFLECT),
        ("P.S1", EventKind.REFLECT), ("P.S2", EventKind.REFRACT), ("IMG", EventKind.TRANSMIT)]
    assert rt.ghost_paths(compiled, 0) == [ghost]
    with pytest.raises(ValueError, match="max_paths"):
        rt.ghost_paths(compiled, "main", max_paths=0)


def test_compile_with_ghosts(plate: rt.System) -> None:
    before = plate.to_json()
    g = rt.compile_with_ghosts(plate, "main")
    assert plate.to_json() == before  # the System itself is not changed
    names = g.system.path_names
    model_paths = [p.name for p in plate.paths]
    assert names == [*model_paths, "main ghost P.S2/P.S1"]  # the model's paths, then the ghost
    assert len(g.ghosts) == 1
    assert g.ghosts.path.tolist() == [len(names) - 1]
    assert g.ghosts.base.tolist() == [0]
    ids = g.system.surface_ids
    assert [ids[g.ghosts.surface_j[0]], ids[g.ghosts.surface_i[0]]] == ["P.S2", "P.S1"]
    assert (g.ghosts.event_j.tolist(), g.ghosts.event_i.tolist()) == ([2], [1])
    by_index = rt.compile_with_ghosts(plate, 0)
    assert by_index.system.path_names == names
    with pytest.raises(ValueError, match="does not exist"):
        rt.compile_with_ghosts(plate, 99)
    with pytest.raises(ValueError):
        rt.compile_with_ghosts(plate, "no such path")


def test_ghost_ranking_of_the_plate(plate: rt.System) -> None:
    # Axis ray only: P_b = T^2, P_g = T R R T, both spots of radius 0, so
    # rho = (P_g / P_b) (0 + r0^2) / (0 + r0^2) = R^2 (ADR 0027, addendum #124).
    g = rt.compile_with_ghosts(plate, "main")
    r = an.ghost_ranking(g, 0, rays="single:0,0")
    assert len(r.entries) == 1
    assert math.isclose(r.base_power, T_PLATE**2, rel_tol=1e-12)
    assert math.isclose(r.entries.power[0], T_PLATE * R_PLATE * R_PLATE * T_PLATE,
                        rel_tol=1e-12)
    assert math.isclose(r.entries.relative_irradiance[0], R_PLATE**2, rel_tol=1e-12)
    assert r.resolution_radius == 0.005
    # A plane plate in a collimated beam: the ghost leaves collimated, no paraxial focus.
    assert math.isnan(r.entries.focus_offset[0])
    assert len(r.entries.losses) == 1 and r.entries.losses[0].launched == 1
    with pytest.raises(ValueError):
        an.ghost_ranking(g, 0, resolution_radius=0.0)


def test_cooke_ranking_is_sorted(reference_dir: Path, catalog_dir: Path) -> None:
    schott = rt.MaterialLibrary()
    schott.add_catalog(catalog_dir / "m2" / "schott.agf")
    g = rt.compile_with_ghosts(rt.load(reference_dir / "m2" / "cooke_triplet.rtt.json"), "main",
                               materials=schott)
    r = an.ghost_ranking(g, 0, rays="hexapolar:3")
    rho = r.entries.relative_irradiance
    assert len(rho) == 15
    assert np.all(np.diff(rho) <= 0.0)  # descending
    assert sorted(r.entries.path.tolist()) == g.ghosts.path.tolist()


def test_run_control(reference_dir: Path, michelson: rt.CompiledSystem) -> None:
    singlet = rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))
    stages: list[str] = []
    an.path_transmission(singlet, 0, field=1, rays="hexapolar:2",
                         progress=lambda done, total, stage: stages.append(stage))
    assert set(stages) == {"aim", "trace"}
    token = rt.CancelToken()
    token.cancel()
    with pytest.raises(rt.Cancelled):
        an.opl_difference(michelson, "reference arm", "test arm", start=collimated_bundle(),
                          cancel=token)


def test_ghost_ranking_options_and_run_control() -> None:
    # Every keyword reaches C++ at its place: resolution_radius shows in the result and changes
    # rho, lost_warning_fraction is checked, cancel and progress act (stages aim and trace). The
    # singlet's ghost is defocused (r_g != r_b), so rho depends on r0; for the plate it would not.
    g = rt.compile_with_ghosts(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json"), "main")
    default = an.ghost_ranking(g, 0, rays="hexapolar:2", threads=1)
    wide = an.ghost_ranking(g, 0, rays="hexapolar:2", resolution_radius=0.02, threads=1)
    assert (default.resolution_radius, wide.resolution_radius) == (0.005, 0.02)
    assert wide.entries.relative_irradiance[0] != default.entries.relative_irradiance[0]
    assert wide.entries.relative_power[0] == default.entries.relative_power[0]
    with pytest.raises(ValueError):
        an.ghost_ranking(g, 0, lost_warning_fraction=1.5)
    stages: list[str] = []
    an.ghost_ranking(g, 0, rays="hexapolar:2",
                     progress=lambda done, total, stage: stages.append(stage))
    assert set(stages) == {"aim", "trace"}
    token = rt.CancelToken()
    token.cancel()
    with pytest.raises(rt.Cancelled):
        an.ghost_ranking(g, 0, cancel=token)
