"""Shared fixtures of the rtt-py tests.

ctest sets RTT_REFERENCE_DIR, RTT_CATALOG_DIR and RTT_PY_REFERENCE_EXE; without them (tests
against an installed package from a source checkout) the directories of the repository are
used and the bitwise comparison with C++ is skipped.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]


def pytest_configure(config: pytest.Config) -> None:
    # A RaytatouilleWarning that a test does not expect (pytest.warns) or filter on purpose
    # (filterwarnings mark) fails the test, so that a new, unexpected warning shows (#86).
    config.addinivalue_line("filterwarnings", "error::raytatouille.errors.RaytatouilleWarning")


def _directory(variable: str, default: Path) -> Path:
    value = os.environ.get(variable)
    return Path(value) if value else default


REFERENCE_DIR = _directory("RTT_REFERENCE_DIR", REPO_ROOT / "tests" / "reference")
CATALOG_DIR = _directory("RTT_CATALOG_DIR", REPO_ROOT / "tests" / "catalogs")
#: The C++ program rtt_py_reference (ctest only); None skips the bitwise tests.
REFERENCE_EXE = os.environ.get("RTT_PY_REFERENCE_EXE") or None


def reference_files() -> list[Path]:
    """All reference systems (*.rtt.json), sorted like the C++ round-trip test."""
    return sorted(REFERENCE_DIR.rglob("*.rtt.json"))


@pytest.fixture
def reference_dir() -> Path:
    return REFERENCE_DIR


@pytest.fixture
def catalog_dir() -> Path:
    return CATALOG_DIR


@pytest.fixture(scope="session", params=[1, 4], ids=lambda t: f"cpp{t}threads")
def cpp_dir(request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Output directory of rtt_py_reference, run once per number of threads (1 and 4)."""
    assert REFERENCE_EXE is not None
    out = tmp_path_factory.mktemp(f"cpp{request.param}")
    subprocess.run(
        [REFERENCE_EXE, str(REFERENCE_DIR), str(CATALOG_DIR), str(out), str(request.param)],
        check=True,
    )
    return out
