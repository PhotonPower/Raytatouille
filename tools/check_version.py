"""Checks that every place that states the version of Raytatouille states the same one (#336).

The version lives in four files without a common source: the CMake project, the Python package
metadata (pyproject.toml), the vcpkg manifest and raytatouille.__version__. This check runs in
the CI (job "JSON schema") and from the repository root:

    python tools/check_version.py

Exit status 0 if all agree, 1 with a list of the places otherwise.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _match(path: Path, pattern: str) -> str | None:
    found = re.search(pattern, path.read_text(encoding="utf-8"), re.MULTILINE)
    return found.group(1) if found else None


def versions() -> dict[str, str | None]:
    """The version stated by each place, None where it is not found."""
    return {
        "CMakeLists.txt (project VERSION)": _match(
            ROOT / "CMakeLists.txt", r"^project\(Raytatouille\s+VERSION\s+([0-9][0-9.]*)"
        ),
        "pyproject.toml ([project] version)": _match(
            ROOT / "pyproject.toml", r'^version\s*=\s*"([^"]+)"'
        ),
        "vcpkg.json (version-semver)": json.loads(
            (ROOT / "vcpkg.json").read_text(encoding="utf-8")
        ).get("version-semver"),
        "raytatouille/__init__.py (__version__)": _match(
            ROOT / "libs" / "rtt-py" / "python" / "raytatouille" / "__init__.py",
            r'^__version__\s*=\s*"([^"]+)"',
        ),
    }


def main() -> int:
    found = versions()
    if None not in found.values() and len(set(found.values())) == 1:
        print(f"version {next(iter(found.values()))} in all {len(found)} places")
        return 0
    print("the version differs between its places:", file=sys.stderr)
    for place, version in found.items():
        print(f"  {place}: {version}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
