"""The Python examples under examples/python run without errors."""

from __future__ import annotations

import runpy
from pathlib import Path

import pytest
from conftest import REPO_ROOT


def test_singlet_example(capsys: pytest.CaptureFixture[str]) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "singlet.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "EFL" in out and "field 2" in out


def test_polarization_example(capsys: pytest.CaptureFixture[str]) -> None:
    # main() returns 1 if any value deviates from its analytic value by more than 1e-12.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "polarization.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "Malus's law" in out and "Stokes after the quarter-wave plate" in out


def test_materials_example(capsys: pytest.CaptureFixture[str]) -> None:
    # main() returns 1 if the vectorized dispersion curve differs from the scalar calls.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "materials.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "Catalogues: SCHOTT, SCHOTT_M2" in out and "SCHOTT_M2:N-LAK9" in out


def test_ray_paths_example(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    # main() returns 1 if a last recorded slot differs from the final ray state.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "ray_paths.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "recorded with 5 slots" in out and "STO, L1.S1, L1.S2, IMG" in out
    pytest.importorskip("matplotlib")
    png = tmp_path / "ray_paths.png"
    assert namespace["main"](str(png)) == 0
    assert png.stat().st_size > 0


def test_layout_example(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    # main() returns 1 if an outline is not closed.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "layout.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "3 lens outlines in the y-z plane" in out and "L2.S1" in out
    pytest.importorskip("matplotlib")
    png = tmp_path / "layout.png"
    assert namespace["main"](str(png)) == 0
    assert png.stat().st_size > 0


def test_analysis_example(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "analysis.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "Spot field 2" in out and "Longitudinal colour" in out
    pytest.importorskip("matplotlib")
    png = tmp_path / "analysis.png"
    assert namespace["main"](str(png)) == 0
    assert png.stat().st_size > 0
