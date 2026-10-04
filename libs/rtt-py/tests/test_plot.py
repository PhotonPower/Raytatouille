"""raytatouille.plot: every plot draws the result data unchanged (matplotlib, backend Agg)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille import plot

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


@pytest.fixture
def compiled(reference_dir: Path, catalog_dir: Path) -> rt.CompiledSystem:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    return rt.compile(rt.load(reference_dir / "m2" / "achromat.rtt.json"), lib)


@pytest.fixture(autouse=True)
def close_figures() -> object:
    yield None
    plt.close("all")


def test_spot(compiled: rt.CompiledSystem) -> None:
    spot = an.spot(compiled, field=2, rays="hexapolar:3")
    ax = plot.spot(spot)
    drawn = np.concatenate([c.get_offsets() for c in ax.collections])
    expected = np.concatenate(
        [np.column_stack([spot.x, spot.y])[spot.wavelengths == wl]
         for wl in np.unique(spot.wavelengths)])
    assert np.array_equal(drawn, expected)
    _, given = plt.subplots()
    assert plot.spot(spot, ax=given) is given


def test_ray_fan_and_opd(compiled: rt.CompiledSystem) -> None:
    fan = an.ray_fan(compiled, field=1)
    lines = plot.ray_fan(fan).get_lines()
    assert np.array_equal(lines[0].get_ydata(), fan.tangential.ey)
    assert np.array_equal(lines[1].get_ydata(), fan.sagittal.ex)
    opd = an.opd_map(compiled, field=1, grid=9)
    points = plot.opd_map(opd).collections[0]
    assert np.array_equal(np.asarray(points.get_array()), opd.points.w)
    opd_fan = an.opd_fan(compiled, field=1)
    lines = plot.opd_fan(opd_fan).get_lines()
    assert np.array_equal(lines[0].get_ydata(), opd_fan.tangential.w)


def test_field_and_colour(compiled: rt.CompiledSystem) -> None:
    d = an.distortion(compiled, samples=5)
    line = plot.distortion(d).get_lines()[0]
    assert np.array_equal(line.get_xdata(), d.percent)
    fc = an.field_curvature(compiled, samples=5)
    lines = plot.field_curvature(fc).get_lines()
    assert np.array_equal(lines[0].get_xdata(), fc.tangential)
    assert np.array_equal(lines[1].get_xdata(), fc.sagittal)
    lc = an.longitudinal_colour(compiled)
    lines = plot.longitudinal_colour(lc).get_lines()
    assert np.array_equal(lines[0].get_xdata(), lc.foci.paraxial_z)


def test_missing_matplotlib_gives_a_hint(monkeypatch: pytest.MonkeyPatch,
                                         compiled: rt.CompiledSystem) -> None:
    spot = an.spot(compiled, rays="hexapolar:2")
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    with pytest.raises(ImportError, match=r"raytatouille\[plot\]"):
        plot.spot(spot)
