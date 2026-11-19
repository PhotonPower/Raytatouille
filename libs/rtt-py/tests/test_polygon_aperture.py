"""Regular polygon and polygon apertures in Python (#320, ADR 0050, point 2): the model classes,
the layout fields aperture and aperture_vertices, the dimension report and the result format
0.1.21."""

from __future__ import annotations

import json
import math

import numpy as np
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import analysis, layout, model, results

MIRROR = REFERENCE_DIR / "m8" / "hexagon_mirror.rtt.json"


def test_model_classes_of_the_reference_file() -> None:
    s = rt.load(MIRROR)
    mirror = s.root.children[1]
    detector = s.root.children[2]
    assert isinstance(mirror, model.Element) and isinstance(detector, model.Element)
    hexagon = mirror.surfaces[0].aperture
    assert isinstance(hexagon, model.RegularPolygonAperture)
    assert (hexagon.sides, hexagon.circumradius.value, hexagon.rotation_deg) == (6, 34.0, 30.0)
    assert hexagon.circumradius.variable
    square = detector.surfaces[0].aperture
    assert isinstance(square, model.PolygonAperture)
    assert square.vertices == [[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0], [-2.0, 2.0]]


def test_layout_vertices_and_parameters() -> None:
    surfaces = layout.surfaces(rt.compile(rt.load(MIRROR)))
    stop, mirror, detector = surfaces
    assert stop.aperture_vertices is None
    assert mirror.aperture == ("regular_polygon",
                               {"sides": 6.0, "circumradius": 34.0, "rotation_deg": 30.0})
    v = mirror.aperture_vertices
    assert v is not None and v.shape == (6, 2)
    # Vertex k at 30 + 60 k degree, counterclockwise from +x (ADR 0050, point 2).
    for k in range(6):
        phi = math.radians(30.0 + 60.0 * k)
        assert np.allclose(v[k], [34.0 * math.cos(phi), 34.0 * math.sin(phi)], atol=1e-13)
    assert detector.aperture == ("polygon", {})
    assert detector.aperture_vertices is not None
    assert detector.aperture_vertices.tolist() == [[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0],
                                                   [-2.0, 2.0]]


def test_result_format_carries_the_vertices() -> None:
    surfaces = layout.surfaces(rt.compile(rt.load(MIRROR)))
    text = surfaces[1].to_json()
    doc = json.loads(text)
    assert doc["schema_version"] == results.SCHEMA_VERSION == "0.1.23"
    loaded = results.load_json(text).data
    vertices = surfaces[1].aperture_vertices
    assert vertices is not None
    assert np.array_equal(loaded["aperture_vertices"], vertices)
    assert json.loads(surfaces[0].to_json())["data"]["aperture_vertices"] is None
    # An older file has no key: no vertices.
    old = dict(doc)
    old["schema_version"] = "0.1.20"
    old["data"] = {k: v for k, v in doc["data"].items() if k != "aperture_vertices"}
    assert results.load_json(json.dumps(old)).data["aperture_vertices"] is None


def test_dimension_report_kinds() -> None:
    assert analysis.ApertureKind.REGULAR_POLYGON.value == 4
    assert analysis.ApertureKind.POLYGON.value == 5
    data = json.loads(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json").to_json())
    lens = data["root"]["children"][1]["surfaces"]
    lens[0]["aperture"] = {"type": "regular_polygon", "sides": 6, "circumradius": 12.7}
    lens[1]["aperture"] = {"type": "polygon",
                           "vertices": [[-9.0, -9.0], [9.0, -9.0], [9.0, 9.0], [-9.0, 9.0]]}
    segments = analysis.dimension_report(rt.compile(rt.System.from_json(json.dumps(data))))
    s = segments.segments
    assert s.aperture_first[0] == analysis.ApertureKind.REGULAR_POLYGON
    assert s.aperture_second[0] == analysis.ApertureKind.POLYGON
    assert s.semi_diameter_first[0] == 12.7
    assert math.isclose(s.semi_diameter_second[0], 9.0 * math.sqrt(2.0), rel_tol=1e-15)


def test_the_hexagon_vignettes_across_its_flats_not_towards_its_vertex() -> None:
    """End to end through m8/hexagon_mirror (ADR 0050, point 2).

    EPD 60 mm, object at infinity: pupil coordinate p maps to the height 30 p mm, and the
    collimated rays meet the paraboloid at the same (x, y). The hexagon has R = 34 mm and
    rotation 30 deg: vertices at 30, 90, ..., 330 deg, edge normals at psi = 60, 120, 180, ...,
    360 deg, inradius r = R cos 30 deg = 29.445 mm.
    - (0.995, 0) hits x = 29.85 mm; along psi = 0 that is beyond r: VIGNETTED at the mirror.
    - (0, 0.995) hits y = 29.85 mm towards the vertex at 90 deg; its largest projection is
      29.85 sin 60 deg = 25.85 mm < r: it passes, reaches the focus inside the +-2 mm square of
      the detector and stays ALIVE.
    The stop (radius 31 mm) passes both."""
    cs = rt.compile(rt.load(MIRROR))
    rays = rt.trace.make_rays(cs, rt.trace.PointListPupil([0.995, 0.0], [0.0, 0.995]))
    rt.trace.trace(cs, rays)
    assert rays.status[0] == rt.trace.RayStatus.VIGNETTED
    assert rays.status[1] == rt.trace.RayStatus.ALIVE
    ids = cs.surface_ids
    assert ids[int(rays.last_surface[0])] == "M"
    assert ids[int(rays.last_surface[1])] == "F"
    # The hand numbers.
    assert 0.995 * 30.0 > 34.0 * math.cos(math.pi / 6.0) > 0.995 * 30.0 * math.sin(math.pi / 3.0)
