"""Configurations of the parameter table and node frames in Python (#169; ADR 0028, 0029).

The zoom of tests/reference/m5/zoom.rtt.json: L2 is placed relative to L1.S2 at the air gap
G (wide 20 mm, tele 5 mm), the image relative to L2.S2 at B = TOTAL - G with TOTAL = 40 mm;
L1 at z = 10 with L1.S2 3 mm behind, L2.S2 2 mm behind L2.S1. So L2.S1 lies at z = 33 (wide) or
18 (tele) and the image at z = 55 in both.
"""

from __future__ import annotations

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import layout


def zoom() -> rt.System:
    return rt.load(REFERENCE_DIR / "m5" / "zoom.rtt.json")


def z_of(cs: rt.CompiledSystem, surface: str) -> float:
    return float(layout.surfaces(cs)[cs.surface_ids.index(surface)].translation[2])


def test_compile_by_index_and_by_name() -> None:
    s = zoom()
    for configuration, name, l2 in ((0, "wide", 33.0), (1, "tele", 18.0), ("tele", "tele", 18.0)):
        cs = rt.compile(s, configuration=configuration)
        assert cs.configuration == (0 if name == "wide" else 1)
        assert cs.configuration_name == name
        assert z_of(cs, "L2.S1") == l2
        assert z_of(cs, "IMG") == 55.0
    nominal = rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json"))
    assert (nominal.configuration, nominal.configuration_name) == (0, "")


@pytest.mark.parametrize("configuration", [2, "Tele", "nominal"])
def test_unknown_configuration_is_a_compile_error(configuration: int | str) -> None:
    with pytest.raises(rt.CompileError) as info:
        rt.compile(zoom(), configuration=configuration)
    codes = [d.code for d in info.value.diagnostics]
    assert codes == ["config.unknown"]
    assert info.value.diagnostics[0].location == "/configurations"


def test_compile_with_ghosts_of_a_configuration() -> None:
    g = rt.compile_with_ghosts(zoom(), "main", configuration="tele")
    tele = rt.compile(zoom(), configuration="tele")
    assert g.system.configuration == 1
    assert z_of(g.system, "L2.S1") == z_of(tele, "L2.S1") == 18.0
    assert len(g.ghosts) > 0


def test_resolved_replaces_bound_params_by_values() -> None:
    s = zoom()
    l2 = s.root.children[2]
    assert (l2.pose.position[2].param, l2.pose.position[2].value) == ("G", None)
    for configuration, g, b in ((0, 20.0, 20.0), ("tele", 5.0, 35.0)):
        r = s.resolved(configuration)
        z_l2 = r.root.children[2].pose.position[2]
        z_image = r.root.children[3].pose.position[2]
        assert (z_l2.param, z_l2.value, z_image.value) == (None, g, b)
        assert r.parameters == s.parameters and r.configurations == s.configurations
        # A resolved system compiles to the same geometry as the column of the original.
        assert z_of(rt.compile(r), "L2.S1") == z_of(rt.compile(s, configuration=configuration),
                                                     "L2.S1")
    # The default is configuration 0; the original is unchanged.
    assert s.resolved().root.children[2].pose.position[2].value == 20.0
    assert s.root.children[2].pose.position[2].value is None


def test_resolved_errors() -> None:
    s = zoom()
    with pytest.raises(ValueError, match="no configuration 'x'"):
        s.resolved("x")
    with pytest.raises(ValueError, match="configuration 2 does not exist"):
        s.resolved(2)
    broken = rt.apply_patch(s, [{"op": "replace", "path": "/parameters/2/expression",
                                 "value": "TOTAL / (G - 5)"}], check="structure_only")
    with pytest.raises(ValueError, match="parameters.not_finite"):
        broken.resolved(0)


def test_locate_parameter() -> None:
    s = zoom()
    assert s.locate_parameter("G") == "/parameters/1"
    assert s.json_at("/parameters/1/name") == "G"
    assert s.locate_parameter("nope") is None


def test_node_frames_and_reference_frame() -> None:
    s = zoom()
    cs = rt.compile(s, configuration="tele")
    frames = layout.node_frames(cs)
    assert frames[0].location == "/root"
    np.testing.assert_array_equal(frames[0].reference, np.eye(4))
    by_location = {f.location: f for f in frames}
    assert len(by_location) == len(frames)
    # Surfaces: the same global transform as the layout, bit for bit.
    for surface in layout.surfaces(cs):
        frame = by_location[s.locate_surface(surface.id) or ""]
        assert frame.to_global.tobytes() == surface.to_global.tobytes()
    # L2 is relative to L1.S2 (z = 13), the image to L2.S2 (z = 20 in tele).
    l2 = s.locate_node("L2")
    image = s.locate_node("image")
    assert l2 is not None and image is not None
    assert layout.reference_frame(cs, l2)[2, 3] == 13.0
    assert layout.reference_frame(cs, image)[2, 3] == 20.0
    # An absolute node: its parent (the root assembly at the origin).
    stop = s.locate_node("stop")
    assert stop is not None
    np.testing.assert_array_equal(layout.reference_frame(cs, stop), np.eye(4))
    with pytest.raises(KeyError):
        layout.reference_frame(cs, "/root/children/9")
