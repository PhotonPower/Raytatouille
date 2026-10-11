"""Explicit pupil point lists from Python (ADR 0032, point 2; Refs #202): rt.trace.PointListPupil
aims every point exactly like rt.trace.SinglePupilPoint."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt


def compiled_singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def test_columns_are_read_only_copies() -> None:
    p = rt.trace.PointListPupil(px=[0.25, 1.5, -0.0], py=np.array([-0.5, 0.0, 0.75]))
    assert len(p) == 3
    assert p.px.tolist() == [0.25, 1.5, -0.0]
    assert math.copysign(1.0, float(p.px[2])) == -1.0  # the sign bit of -0.0 is kept
    assert p.py.tolist() == [-0.5, 0.0, 0.75]
    assert not p.px.flags.writeable and not p.py.flags.writeable


def test_make_rays_equals_single_points(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    px = [0.0, 0.3, -0.6, 0.0]
    py = [0.0, -0.7, 0.2, 1.0]
    rays = rt.trace.make_rays(cs, rt.trace.PointListPupil(px, py), fields=[0, 1])
    assert len(rays) == 2 * len(px)
    for f in (0, 1):
        for k in range(len(px)):
            i = f * len(px) + k
            single = rt.trace.make_rays(cs, rt.trace.SinglePupilPoint(px[k], py[k]), fields=[f])
            for column in ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z"):
                assert getattr(rays, column)[i] == getattr(single, column)[0], (f, k, column)
            assert rays.pupil_x[i] == px[k] and rays.pupil_y[i] == py[k]
            assert rays.field[i] == f
    # Content guard: aimed rays starting at different points (object at infinity: one direction).
    assert np.all(rays.status == int(rt.trace.RayStatus.ALIVE))
    assert len(set(np.round(rays.pos_y[: len(px)], 12).tolist())) > 1


def test_input_errors(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    with pytest.raises(ValueError, match="same length"):
        rt.trace.PointListPupil([0.0, 0.1], [0.0])
    with pytest.raises(ValueError, match="empty"):
        rt.trace.make_rays(cs, rt.trace.PointListPupil([], []))
    with pytest.raises(ValueError, match="point 1"):
        rt.trace.make_rays(cs, rt.trace.PointListPupil([0.0, math.nan], [0.0, 0.0]))
