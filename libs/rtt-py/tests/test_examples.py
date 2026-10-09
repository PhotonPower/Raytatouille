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
    assert "paraxial working F/#" in out and "L1.S1" in out


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


def test_results_example(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    # main() returns 1 if the JSON round trip changes the data.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "results.py"))
    assert namespace["main"](str(tmp_path)) == 0
    out = capsys.readouterr().out
    assert "SpotDiagram written to" in out and "EFL " in out and "arrived" in out
    assert (tmp_path / "spot.result.json").is_file()


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
    assert "2 lens outlines in the y-z plane" in out and "L1.S2" in out
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


def test_model_tree_example(capsys: pytest.CaptureFixture[str]) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "model_tree.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "lens 'L2' at z 10.691 mm, SCHOTT:N-SF5" in out
    assert "L1.S1: R 23.713 mm" in out and "path 'main'" in out


def test_cancel_progress_example(capsys: pytest.CaptureFixture[str]) -> None:
    # main() returns 1 unless every stage of the complete run ended with done == total and the
    # second run was cancelled.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "cancel_progress.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "spot complete" in out and "second run: cancelled" in out


def test_edit_undo_example(capsys: pytest.CaptureFixture[str]) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "edit_undo.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "rejected: [surface.id_duplicate]" in out
    assert "after two commands: radius 50.0 mm" in out and "same system: True" in out


def test_ghosts_example(capsys: pytest.CaptureFixture[str]) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "ghosts.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "15 ghosts of 'main'" in out and "main ghost " in out
    assert "reference arm: transmission 0.250 (3 of 3 rays)" in out
    assert "OPL(test arm) - OPL(reference arm) = 15.000000 mm" in out


def test_grating_orders_example(capsys: pytest.CaptureFixture[str]) -> None:
    # main() returns 1 if a direction or a transmission deviates from its closed form by more
    # than 1e-12.
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "grating_orders.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "lambda0 G = 0.17628" in out and "order +6  EVANESCENT" in out
    assert "order +1  transmission 0.400" in out
