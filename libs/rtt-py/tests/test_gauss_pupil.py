"""The Gaussian pupil sampling from Python (#168): rt.trace.GaussPupil, its quadrature weights
and the analyses that reject it (coordinator's decision (a))."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt


def compiled_singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def test_defaults_and_weights() -> None:
    g = rt.trace.GaussPupil()
    assert (g.rings, g.arms) == (3, 6)
    q = rt.trace.gauss_pupil_weights(rt.trace.GaussPupil(rings=2, arms=4))
    assert len(q) == 8
    assert math.isclose(sum(q), 1.0, abs_tol=1e-14)


def test_make_rays_with_the_gauss_sampling(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    g = rt.trace.GaussPupil(rings=3, arms=6)
    rays = rt.trace.make_rays(cs, g, fields=[1])
    assert len(rays) == 18
    q = np.asarray(rt.trace.gauss_pupil_weights(g))
    # Mean of u = px^2 + py^2 over the pupil with the quadrature weights: 1/2 (exact).
    u = rays.pupil_x**2 + rays.pupil_y**2
    assert abs(float(np.sum(q * u)) - 0.5) <= 1e-14
    # The ray weight stays the power (ADR 0021).
    assert np.all(rays.weight == 1.0)


def test_the_averaging_analyses_reject_it(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    with pytest.raises(ValueError, match="gauss"):
        rt.analysis.spot(cs, rays=rt.trace.GaussPupil())
    with pytest.raises(ValueError, match="gauss"):
        rt.analysis.path_transmission(cs, rays=rt.trace.GaussPupil())


def test_input_errors(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    with pytest.raises(ValueError):
        rt.trace.make_rays(cs, rt.trace.GaussPupil(rings=0))
    with pytest.raises(ValueError):
        rt.trace.gauss_pupil_weights(rt.trace.GaussPupil(arms=0))
