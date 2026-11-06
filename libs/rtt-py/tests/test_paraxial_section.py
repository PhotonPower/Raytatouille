"""First-order data per section x-z and y-z in Python (ADR 0042, point 2; #315):
rt.paraxial.Section, first_order_section, trace_ray_section and is_anamorphic.

Oracle: the examples under examples/anamorphic (VACUUM, n = 1.5168, t = 4 mm, stop at z = 0,
lens front at z = 5). A curved side of radius R in front of a plane rear side has
efl = R / (n - 1) and bfl = efl - t / n (Greivenkamp, OPTI-201/202, Sec. 7, p. 7-4 with
P2 = 0, and Sec. 9, p. 9-12; docs/quellen.md): R = 51.68 gives 100 mm, R = 77.52 gives 150 mm.
The crossed cylinders focus both sections at z = 106.36286919831225 (example file). A section
without power is afocal: efl None. The C++ tests (test_paraxial_section.cpp) check the values
against an independent matrix computation; here the bindings are checked.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from conftest import REPO_ROOT

import raytatouille as rt
from raytatouille import paraxial

EXAMPLE_DIR = REPO_ROOT / "examples" / "anamorphic"
N = 1.5168


def example(name: str) -> rt.CompiledSystem:
    return rt.compile(rt.load(EXAMPLE_DIR / f"{name}.rtt.json"))


def test_section_enum() -> None:
    assert {s.name for s in paraxial.Section} == {"X", "Y"}
    assert paraxial.Section.X != paraxial.Section.Y


def test_cylinder_lens_sections() -> None:
    cs = example("cylinder_planoconvex")
    y = paraxial.first_order_section(cs, paraxial.Section.Y)
    assert y.efl is not None and y.bfl is not None
    assert y.efl == pytest.approx(100.0, rel=1e-12)
    assert y.bfl == pytest.approx(100.0 - 4.0 / N, rel=1e-12)
    x = paraxial.first_order_section(cs, paraxial.Section.X, path="main", wavelength=None)
    assert x.efl is None and x.bfl is None
    assert x.entrance_pupil is not None and x.entrance_pupil.z == 0.0
    assert x.entrance_pupil.diameter == 10.0
    assert paraxial.is_anamorphic(cs)
    assert paraxial.is_anamorphic(cs, "main")


def test_biconic_and_crossed_cylinders() -> None:
    cs = example("biconic_singlet")
    x = paraxial.first_order_section(cs, paraxial.Section.X, 0, 0)
    y = paraxial.first_order_section(cs, paraxial.Section.Y, 0, 0)
    assert x.efl == pytest.approx(150.0, rel=1e-12)
    assert y.efl == pytest.approx(100.0, rel=1e-12)
    crossed = example("crossed_cylinders")
    for s in (paraxial.Section.X, paraxial.Section.Y):
        fo = paraxial.first_order_section(crossed, s)
        assert fo.rear_focal_z == pytest.approx(106.36286919831225, rel=1e-12)
    fx = paraxial.first_order_section(crossed, paraxial.Section.X)
    fy = paraxial.first_order_section(crossed, paraxial.Section.Y)
    # Content guard: the sections have different focal lengths (75 mm lens at z = 30 in x,
    # 100 mm lens at z = 5 in y).
    assert fx.efl is not None and fy.efl is not None and abs(fx.efl - fy.efl) > 1.0


def test_first_order_still_rejects_anamorphic_paths() -> None:
    cs = example("biconic_singlet")
    with pytest.raises(rt.ParaxialError) as error:
        paraxial.first_order(cs)
    assert "not rotationally symmetric in second order" in str(error.value)
    assert error.value.surface == "L1.S1"
    assert error.value.location == "/root/children/1/surfaces/0"


def test_rotationally_symmetric_path_bitwise(reference_dir: Path) -> None:
    cs = rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))
    assert not paraxial.is_anamorphic(cs)
    for w in range(len(cs.wavelengths_um)):
        fo = paraxial.first_order(cs, 0, w)
        assert fo.efl is not None  # content guard
        for s in (paraxial.Section.X, paraxial.Section.Y):
            fs = paraxial.first_order_section(cs, s, 0, w)
            assert fs.efl == fo.efl
            assert fs.bfl == fo.bfl
            assert fs.rear_principal_z == fo.rear_principal_z
            assert fs.entrance_pupil is not None and fo.entrance_pupil is not None
            assert fs.entrance_pupil.z == fo.entrance_pupil.z
            assert fs.exit_pupil is not None and fo.exit_pupil is not None
            assert fs.exit_pupil.z == fo.exit_pupil.z


def test_trace_ray_section() -> None:
    # Ray parallel to the axis at height 2 from z = -1: the x section of the cylinder lens has
    # no power (the ray stays parallel, u = 0 after every event); the y section bends it
    # towards the axis and meets it at the rear focus.
    cs = example("cylinder_planoconvex")
    x = paraxial.trace_ray_section(cs, paraxial.Section.X, -1.0, 2.0, 0.0)
    assert [r.z for r in x] == pytest.approx([0.0, 5.0, 9.0, 106.36286919831225], rel=1e-12)
    assert all(r.y == 2.0 and r.u == 0.0 for r in x)
    assert [r.n for r in x] == [1.0, N, 1.0, 1.0]
    y = paraxial.trace_ray_section(cs, paraxial.Section.Y, -1.0, 2.0, 0.0, path="main")
    assert len(y) == 4
    assert y[-1].u < 0.0
    # The image plane of the example is the paraxial y focus: the height there is 0.
    assert abs(y[-1].y) < 1e-12
    assert isinstance(y[0], paraxial.RayAtEvent)


def test_principal_axes_not_aligned() -> None:
    system = rt.load(EXAMPLE_DIR / "biconic_singlet.rtt.json")
    rotated = rt.apply_patch(system, [{
        "op": "add", "path": "/root/children/1/surfaces/0/pose",
        "value": {"position": [0.0, 0.0, 0.0], "rotation_deg": [0.0, 0.0, 45.0]}}])
    cs = rt.compile(rotated)
    for s in (paraxial.Section.X, paraxial.Section.Y):
        with pytest.raises(rt.ParaxialError) as error:
            paraxial.first_order_section(cs, s)
        assert "principal axes not aligned with x and y" in str(error.value)
        assert error.value.surface == "L1.S1"
        assert error.value.location == "/root/children/1/surfaces/0"
