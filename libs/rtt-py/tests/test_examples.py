"""The Python examples under examples/python run without errors."""

from __future__ import annotations

import runpy

import pytest
from conftest import REPO_ROOT


def test_singlet_example(capsys: pytest.CaptureFixture[str]) -> None:
    namespace = runpy.run_path(str(REPO_ROOT / "examples" / "python" / "singlet.py"))
    assert namespace["main"]() == 0
    out = capsys.readouterr().out
    assert "EFL" in out and "field 2" in out
