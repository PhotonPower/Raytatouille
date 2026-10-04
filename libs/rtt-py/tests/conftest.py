"""Shared fixtures of the rtt-py tests.

ctest sets RTT_REFERENCE_DIR, RTT_CATALOG_DIR and RTT_PY_REFERENCE_EXE; without them (tests
against an installed package from a source checkout) the directories of the repository are
used and the bitwise comparison with C++ is skipped.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]


def _directory(variable: str, default: Path) -> Path:
    value = os.environ.get(variable)
    return Path(value) if value else default


REFERENCE_DIR = _directory("RTT_REFERENCE_DIR", REPO_ROOT / "tests" / "reference")
CATALOG_DIR = _directory("RTT_CATALOG_DIR", REPO_ROOT / "tests" / "catalogs")


def reference_files() -> list[Path]:
    """All reference systems (*.rtt.json), sorted like the C++ round-trip test."""
    return sorted(REFERENCE_DIR.rglob("*.rtt.json"))


@pytest.fixture
def reference_dir() -> Path:
    return REFERENCE_DIR


@pytest.fixture
def catalog_dir() -> Path:
    return CATALOG_DIR
