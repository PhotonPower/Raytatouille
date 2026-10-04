"""rt.analysis beyond the bitwise comparison: arguments, sampling shorthand, results."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille.trace import (
    FanXPupil,
    FanYPupil,
    GridPupil,
    HexapolarPupil,
    RandomPupil,
    RayStatus,
    SinglePupilPoint,
)


def singlet(reference_dir: Path) -> rt.System:
    return rt.load(reference_dir / "m1" / "singlet_const.rtt.json")


def test_sampling_shorthand() -> None:
    hexapolar = an.sampling("hexapolar:12")
    assert isinstance(hexapolar, HexapolarPupil) and hexapolar.rings == 12
    grid = an.sampling("grid:33")
    assert isinstance(grid, GridPupil) and grid.n == 33
    assert isinstance(an.sampling("fan_x:5"), FanXPupil)
    assert isinstance(an.sampling("fan_y:7"), FanYPupil)
    random = an.sampling("random:500")
    assert isinstance(random, RandomPupil) and (random.count, random.seed) == (500, 0)
    random = an.sampling("random:500:42")
    assert isinstance(random, RandomPupil) and (random.count, random.seed) == (500, 42)
    single = an.sampling("single:0.5,-0.25")
    assert isinstance(single, SinglePupilPoint) and (single.px, single.py) == (0.5, -0.25)
    single = an.sampling("single:-1e-1,.5")
    assert isinstance(single, SinglePupilPoint) and (single.px, single.py) == (-0.1, 0.5)
    hexapolar = an.sampling("hexapolar:0")  # centre ray only
    assert isinstance(hexapolar, HexapolarPupil) and hexapolar.rings == 0
    random = an.sampling("random:0:18446744073709551615")
    assert isinstance(random, RandomPupil) and random.seed == 2**64 - 1
    objects = GridPupil(n=3)
    assert an.sampling(objects) is objects


@pytest.mark.parametrize(
    "bad",
    ["hexapolar", "hexapolar:", "hexapolar:x", "hexapolar:1:2", "grid:1.5", "spiral:3",
     "random:-1", "random:10:-2", "random:1:2:3", "single:1", "single:a,b", "",
     # Ranges and overflow (review of #33), lenient int()/float() forms, non-finite values.
     "grid:0", "fan_x:0", "fan_y:0", "hexapolar:-1", "hexapolar:99999999999",
     "random:99999999999999999999999", "random:1:18446744073709551616", "hexapolar:1_0",
     "hexapolar: 3", "grid:+3", "single:nan,0", "single:0,inf", "single:1e999,0",
     "single:1,2,3", "single:0x1,0"],
)
def test_malformed_shorthand_names_the_string(bad: str) -> None:
    with pytest.raises(ValueError) as info:
        an.sampling(bad)
    assert repr(bad) in str(info.value)


def test_system_and_compiled_system_give_the_same_result(reference_dir: Path) -> None:
    system = singlet(reference_dir)
    compiled = rt.compile(system)
    from_system = an.spot(system, "main", 1, rays="hexapolar:4")
    from_compiled = an.spot(compiled, "main", 1, rays="hexapolar:4")
    assert from_system.x.tobytes() == from_compiled.x.tobytes()
    assert from_system.y.tobytes() == from_compiled.y.tobytes()
    assert from_system.stats.rms_centroid == from_compiled.stats.rms_centroid
    with pytest.raises(ValueError, match="materials only applies to a System"):
        an.spot(compiled, materials=rt.MaterialLibrary())


def test_system_with_materials(reference_dir: Path, catalog_dir: Path) -> None:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    achromat = rt.load(reference_dir / "m2" / "achromat.rtt.json")
    with pytest.raises(rt.CompileError):
        an.lateral_colour(achromat)  # catalogue glass without the catalogue
    lateral = an.lateral_colour(achromat, field=2, materials=lib)
    assert len(lateral.chief) == 3 and len(lateral.offset) == 3


def test_seidel_alias_is_the_same_object(reference_dir: Path) -> None:
    assert an.seidel is rt.paraxial.seidel
    s = rt.paraxial.seidel(singlet(reference_dir))
    assert s.chromatic is None and s.sum.c_l == 0.0
    compiled = rt.compile(singlet(reference_dir))
    with_pair = rt.paraxial.seidel(compiled, pair=rt.paraxial.ChromaticPair(0, 2))
    pair = with_pair.chromatic
    assert pair is not None and (pair.first, pair.second) == (0, 2)
    assert len(s.surfaces) == len(s.surfaces.s1)


def test_result_arrays_are_read_only_copies(reference_dir: Path) -> None:
    spot = an.spot(singlet(reference_dir), rays="hexapolar:3")
    x = spot.x
    assert not x.flags.writeable
    with pytest.raises(ValueError):
        x[0] = 1.0
    assert not np.shares_memory(x, spot.x)  # every access is a new copy
    fan = an.ray_fan(singlet(reference_dir))
    assert not fan.tangential.ey.flags.writeable
    assert fan.tangential.status.dtype == np.uint8


def test_wavelength_none_is_polychromatic_only_for_spot(reference_dir: Path) -> None:
    compiled = rt.compile(singlet(reference_dir))
    poly = an.spot(compiled, field=1)
    assert poly.wavelength is None
    assert set(np.unique(poly.wavelengths).tolist()) == {0, 1, 2}
    mono = an.spot(compiled, field=1, wavelength=compiled.reference_wavelength)
    assert mono.wavelength == compiled.reference_wavelength
    fan = an.ray_fan(compiled, field=1)
    assert fan.wavelength == compiled.reference_wavelength


def test_physics_sanity(reference_dir: Path) -> None:
    # Paraboloid on axis: perfect focus (architecture, validation; C++ test_spot.cpp).
    paraboloid = rt.load(reference_dir / "m2" / "paraboloid_stop.rtt.json")
    assert an.spot(paraboloid, rays="hexapolar:6").stats.rms_centroid < 1e-9
    assert an.opd_map(paraboloid).rms < 1e-6
    # Distortion is 0 on axis by definition (field.hpp).
    d = an.distortion(singlet(reference_dir), samples=3)
    assert len(d) == 3 and d.fraction.tolist() == [0.0, 0.5, 1.0] and d.percent[0] == 0.0
    assert d[2].percent == d.percent[2] and d[-1].percent == d.percent[-1]
    with pytest.raises(IndexError):
        d[3]
    with pytest.raises(IndexError):
        d[-4]
    point = an.distortion_at(singlet(reference_dir), (0.0, 4.0))
    point.field.x = 99.0  # nested structs are copies: the result stays unchanged
    assert point.field.x == 0.0
    curvature = an.field_curvature_at(singlet(reference_dir), (0.0, 0.0))
    assert curvature.astigmatism == curvature.tangential - curvature.sagittal


def test_analysis_error(reference_dir: Path) -> None:
    # An annular mirror blocks the pupil centre: the chief ray never reaches the image.
    data = json.loads(rt.load(reference_dir / "m2" / "paraboloid_stop.rtt.json").to_json())
    data["root"]["children"][1]["surfaces"][0]["aperture"] = {
        "type": "circular", "radius": 31.0, "inner_radius": 5.0}
    with pytest.raises(rt.AnalysisError) as info:
        an.spot(rt.System.from_json(json.dumps(data)))
    assert isinstance(info.value, rt.RaytatouilleError)
    assert "chief ray" in str(info.value)


def test_vignetted_rays_keep_their_status(reference_dir: Path) -> None:
    data = json.loads(rt.load(reference_dir / "m2" / "paraboloid_stop.rtt.json").to_json())
    data["root"]["children"][1]["surfaces"][0]["aperture"] = {"type": "circular", "radius": 21.0}
    opd = an.opd_map(rt.System.from_json(json.dumps(data)), grid=21)
    status = opd.points.status
    assert opd.vignetted == np.count_nonzero(status != int(RayStatus.ALIVE)) > 0
    assert np.all(opd.points.w[status != int(RayStatus.ALIVE)] == 0.0)


def test_threads_argument(reference_dir: Path) -> None:
    compiled = rt.compile(singlet(reference_dir))
    one = an.spot(compiled, field=2, rays="random:400:1", threads=1)
    four = an.spot(compiled, field=2, rays="random:400:1", threads=4)
    assert one.x.tobytes() == four.x.tobytes()
    with pytest.raises(ValueError):
        an.spot(compiled, threads=0)
