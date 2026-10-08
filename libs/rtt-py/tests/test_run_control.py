"""Cancellation and progress (#83): CancelToken, cancel= and progress= on trace, make_rays and
the bundle and sweep analyses, rt.errors.Cancelled."""

from __future__ import annotations

import enum
import struct
import threading
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import raytatouille as rt
from raytatouille import analysis as an

COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "wl", "opl", "weight", "field",
           "pupil_x", "pupil_y", "last_surface", "status")


def compiled_singlet(reference_dir: Path) -> rt.CompiledSystem:
    # Three fields (0, 3.5 and 5 deg) and three wavelengths.
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


def same(a: object, b: object, depth: int = 0) -> bool:
    """Bitwise equality of two analysis results: arrays and floats by their bytes, enums and
    integers by value, other objects attribute by attribute."""
    assert depth < 6, "result nested deeper than expected"
    if isinstance(a, np.ndarray):
        return (isinstance(b, np.ndarray) and a.dtype == b.dtype and a.shape == b.shape
                and a.tobytes() == b.tobytes())
    if isinstance(a, float):
        return isinstance(b, float) and struct.pack("<d", a) == struct.pack("<d", b)
    if a is None or isinstance(a, (bool, int, str, enum.Enum)):
        return bool(a == b)
    if isinstance(a, (list, tuple)):
        return (isinstance(b, (list, tuple)) and len(a) == len(b)
                and all(same(x, y, depth + 1) for x, y in zip(a, b)))
    names = [n for n in dir(a) if not n.startswith("_")]
    for name in names:
        value = getattr(a, name)
        if callable(value):
            continue
        if not same(value, getattr(b, name), depth + 1):
            return False
    return True


class Recorder:
    """Progress callback that stores every report."""

    def __init__(self) -> None:
        self.reports: list[tuple[int, int, str]] = []

    def __call__(self, done: int, total: int, stage: str) -> None:
        self.reports.append((done, total, stage))

    def stages(self) -> set[str]:
        return {stage for _, _, stage in self.reports}


def analyses(cs: rt.CompiledSystem) -> dict[str, Callable[..., Any]]:
    """Every analysis with cancel= and progress=, called on the outer field of the singlet."""
    return {
        "spot": lambda **kw: an.spot(cs, field=2, rays="hexapolar:6", **kw),
        "ray_fan": lambda **kw: an.ray_fan(cs, field=2, points=11, **kw),
        "opd_map": lambda **kw: an.opd_map(cs, field=2, grid=9, **kw),
        "opd_fan": lambda **kw: an.opd_fan(cs, field=2, points=11, **kw),
        "distortion": lambda **kw: an.distortion(cs, samples=5, **kw),
        "field_curvature": lambda **kw: an.field_curvature(cs, samples=5, **kw),
        "longitudinal_colour": lambda **kw: an.longitudinal_colour(cs, **kw),
        "lateral_colour": lambda **kw: an.lateral_colour(cs, field=2, **kw),
    }


def test_cancel_token() -> None:
    token = rt.CancelToken()
    assert rt.CancelToken is rt.trace.CancelToken and not token.cancelled
    token.cancel()
    assert token.cancelled
    token.cancel()  # idempotent
    assert token.cancelled


def test_cancelled_is_a_raytatouille_error() -> None:
    assert issubclass(rt.errors.Cancelled, rt.errors.RaytatouilleError)
    assert rt.Cancelled is rt.errors.Cancelled


def test_trace_with_control_is_bitwise_the_same(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    sampling = rt.trace.RandomPupil(count=3000, seed=5)
    plain = rt.trace.make_rays(cs, sampling)
    controlled = rt.trace.make_rays(cs, sampling, cancel=rt.CancelToken(), progress=Recorder())
    for name in COLUMNS:
        assert getattr(plain, name).tobytes() == getattr(controlled, name).tobytes(), name
    stats = rt.trace.trace(cs, plain)
    for threads in (1, 4):
        rays = rt.trace.make_rays(cs, sampling)
        controlled_stats = rt.trace.trace(cs, rays, threads=threads, cancel=rt.CancelToken(),
                                          progress=Recorder())
        assert list(controlled_stats.rays) == list(stats.rays)
        for name in COLUMNS:
            assert getattr(plain, name).tobytes() == getattr(rays, name).tobytes(), name
        assert plain.prt_matrices().tobytes() == rays.prt_matrices().tobytes()


def test_analyses_with_control_are_bitwise_the_same(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    for name, run in analyses(cs).items():
        assert same(run(), run(cancel=rt.CancelToken(), progress=Recorder())), name


def test_progress_reports_of_trace_and_make_rays(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    aim = Recorder()
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=20), progress=aim)
    n = len(rays)
    assert aim.stages() == {"aim"}
    assert aim.reports[-1] == (n, n, "aim")
    traced = Recorder()
    rt.trace.trace(cs, rays, progress=traced)
    assert traced.stages() == {"trace"}
    assert traced.reports[-1] == (n, n, "trace")
    done = [d for d, _, _ in traced.reports]
    assert done == sorted(done) and all(total == n for _, total, _ in traced.reports)


def test_progress_stages_of_the_analyses(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    expected = {
        "spot": {"aim", "trace"},
        "ray_fan": {"aim", "trace"},
        "opd_map": {"aim", "trace"},
        "opd_fan": {"aim", "trace"},
        "distortion": {"field"},
        "field_curvature": {"field"},
        "longitudinal_colour": {"wavelength"},
        "lateral_colour": {"wavelength"},
    }
    for name, run in analyses(cs).items():
        progress = Recorder()
        run(progress=progress)
        assert progress.stages() == expected[name], name
        done, total, _ = progress.reports[-1]
        assert done == total, name


def test_cancel_before_start(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    token = rt.CancelToken()
    token.cancel()
    with pytest.raises(rt.errors.Cancelled):
        rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=6), cancel=token)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=6))
    with pytest.raises(rt.errors.Cancelled):
        rt.trace.trace(cs, rays, cancel=token)
    for name, run in analyses(cs).items():
        with pytest.raises(rt.errors.Cancelled):
            run(cancel=token)
            pytest.fail(name)


def test_cancel_from_another_python_thread(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=100),
                              aiming=rt.trace.Aiming.PARAXIAL)
    token = rt.CancelToken()
    started = threading.Event()
    requested = threading.Event()

    def cancel_when_started() -> None:
        if started.wait(timeout=30.0):
            token.cancel()
            requested.set()

    def progress(done: int, total: int, stage: str) -> None:
        # The first report blocks until the other thread has cancelled: the request certainly
        # falls into the running trace. Waiting releases the GIL for the other thread.
        if not started.is_set():
            started.set()
            assert requested.wait(timeout=30.0)

    canceller = threading.Thread(target=cancel_when_started)
    canceller.start()
    try:
        with pytest.raises(rt.errors.Cancelled):
            rt.trace.trace(cs, rays, cancel=token, progress=progress)
    finally:
        started.set()  # never leave the other thread waiting
        canceller.join()
    assert token.cancelled


def test_exception_of_the_callback_is_raised(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=20))
    calls: list[int] = []

    class Stop(Exception):
        pass

    def progress(done: int, total: int, stage: str) -> None:
        calls.append(done)
        raise Stop("from the callback")

    with pytest.raises(Stop, match="from the callback"):
        rt.trace.trace(cs, rays, progress=progress)
    assert len(calls) == 1
    with pytest.raises(Stop):
        an.distortion(cs, samples=5, progress=progress)
    with pytest.raises(Stop):
        rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=6), progress=progress)


def test_invalid_control_arguments(reference_dir: Path) -> None:
    cs = compiled_singlet(reference_dir)
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=2))
    with pytest.raises(TypeError):
        rt.trace.trace(cs, rays, progress=42)  # type: ignore[call-overload]
    with pytest.raises(TypeError):
        rt.trace.trace(cs, rays, cancel=True)  # type: ignore[call-overload]
    # Recording a path has no cancellation (it is meant for a few rays of a drawing).
    with pytest.raises(ValueError, match="record_path"):
        rt.trace.trace(cs, rays, record_path=True, cancel=rt.CancelToken())
    with pytest.raises(ValueError, match="record_path"):
        rt.trace.trace(cs, rays, record_path=True, progress=Recorder())
