"""The grid sag term in Python (#316, ADR 0051, point 5, PR 4): the model class with its lists,
the layout field terms, the sag at the grid points and between them, the domain and the missing
paraxial model."""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from conftest import REFERENCE_DIR

import raytatouille as rt
from raytatouille import layout, model

MIRROR = REFERENCE_DIR / "m8" / "grid_sag_mirror.rtt.json"


def _astig(x: float, y: float) -> float:
    return 2e-5 * (x * x - y * y) / 484.0


def _plane_with_grid() -> rt.CompiledSystem:
    data = json.loads(MIRROR.read_text(encoding="utf-8"))
    shape = data["root"]["children"][1]["surfaces"][0]["shape"]
    shape["base"] = {"type": "plane"}
    return rt.compile(rt.System.from_json(json.dumps(data)))


def test_model_class_of_the_reference_file() -> None:
    mirror = rt.load(MIRROR).root.children[1]
    assert isinstance(mirror, model.Element)
    (grid,) = mirror.surfaces[0].shape.terms
    assert isinstance(grid, model.GridSag)
    assert (grid.nx, grid.ny, grid.spacing_x, grid.spacing_y) == (9, 9, 6.0, 6.0)
    assert (grid.decenter_x, grid.decenter_y) == (0.0, 0.0)
    assert len(grid.sag) == 81 and grid.dz_dx == [] and grid.d2z_dxdy == []
    assert math.isclose(grid.sag[9 * 4 + 8], _astig(24.0, 0.0), rel_tol=1e-12)  # i = 8, j = 4


def test_layout_lists_the_grid() -> None:
    stop, mirror, detector = layout.surfaces(rt.compile(rt.load(MIRROR)))
    assert mirror.shape == "conic"
    assert mirror.terms == [("grid_sag", {"nx": 9.0, "ny": 9.0, "spacing_x": 6.0,
                                          "spacing_y": 6.0, "decenter_x": 0.0,
                                          "decenter_y": 0.0})]
    assert stop.terms == [] and detector.terms == []
    assert json.loads(mirror.to_json())["data"]["terms"][0][0] == "grid_sag"


def test_sag_at_and_between_the_grid_points() -> None:
    """At a grid point the value of the file; between them the interpolant of the quadratic.
    Finite differences of second order are exact for a quadratic, so the bicubic interpolant is
    the quadratic itself up to rounding."""
    cs = _plane_with_grid()
    xs = np.array([0.0, 6.0, -18.0, 3.7, -11.2, 23.9])
    ys = np.array([0.0, -12.0, 24.0, 1.3, 7.75, -5.5])
    z = rt.layout.sag(cs, 1, xs, ys)
    for k in range(len(xs)):
        assert math.isclose(float(z[k]), _astig(float(xs[k]), float(ys[k])), rel_tol=1e-12,
                            abs_tol=1e-18)
    assert float(z[2]) != 0.0  # content guard


def test_outside_the_grid_is_nan_and_no_paraxial_model() -> None:
    cs = _plane_with_grid()
    assert math.isnan(float(rt.layout.sag(cs, 1, np.array([24.5]), np.array([0.0]))[0]))
    with pytest.raises(rt.errors.ParaxialError):
        rt.paraxial.first_order(cs, path=0)
