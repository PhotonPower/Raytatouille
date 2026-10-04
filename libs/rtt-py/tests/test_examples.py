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


def test_analysis_example(capsys: pytest.CaptureFixture[str], tmp_path: Path) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "analysis.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "Spot field 2" in out and "Longitudinal colour" in out
    pytest.importorskip("matplotlib")
    png = tmp_path / "analysis.png"
    assert namespace["main"](str(png)) == 0
    assert png.stat().st_size > 0
