"""Shared fixtures of the rtt-py tests.

ctest sets RTT_REFERENCE_DIR, RTT_CATALOG_DIR and RTT_PY_REFERENCE_EXE; without them (tests
against an installed package from a source checkout) the directories of the repository are
used and the bitwise comparison with C++ is skipped.

Every test runs under a watchdog (#35): a test that takes longer than RTT_PY_TEST_TIMEOUT
seconds (default 300) prints the tracebacks of all threads and ends the process with status 1.
"""

from __future__ import annotations

import faulthandler
import os
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]

#: Watchdog per test in seconds (#35): generous, a hang of the whole run is what it catches.
TEST_TIMEOUT = float(os.environ.get("RTT_PY_TEST_TIMEOUT") or 300.0)
#: Copy of the stderr file descriptor from before the output capture (as pytest's faulthandler
#: plugin keeps one): the watchdog writes there, and the capture is lost when it ends the process.
_WATCHDOG_FD = pytest.StashKey[int]()


def pytest_configure(config: pytest.Config) -> None:
    # A RaytatouilleWarning that a test does not expect (pytest.warns) or filter on purpose
    # (filterwarnings mark) fails the test, so that a new, unexpected warning shows (#86).
    config.addinivalue_line("filterwarnings", "error::raytatouille.errors.RaytatouilleWarning")
    config.stash[_WATCHDOG_FD] = os.dup(sys.stderr.fileno())


def pytest_unconfigure(config: pytest.Config) -> None:
    if _WATCHDOG_FD in config.stash:
        os.close(config.stash[_WATCHDOG_FD])
        del config.stash[_WATCHDOG_FD]


@pytest.fixture(autouse=True)
def watchdog(request: pytest.FixtureRequest) -> Iterator[None]:
    """Ends a hanging test run with the tracebacks of all threads and exit status 1 (#35).

    faulthandler's watchdog is a C thread: it fires also while another thread holds the GIL in
    C code, e.g. a binding that waits for a worker which waits for the GIL (#83). pytest-timeout
    would not: its timer thread and its SIGALRM handler both run Python code, which needs the
    GIL.
    """
    faulthandler.dump_traceback_later(
        TEST_TIMEOUT, exit=True, file=request.config.stash[_WATCHDOG_FD]
    )
    yield
    faulthandler.cancel_dump_traceback_later()


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
