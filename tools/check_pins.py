"""Pins of dependencies outside vcpkg must equal the vcpkg port (ADR 0033, point 1; #204).

cmake/RttPocketfft.cmake downloads pocketfft with FetchContent when find_path does not find the
header (presets with system packages). It must load exactly the version of the vcpkg port in
the builtin-baseline. This script compares

- RTT_POCKETFFT_REF and RTT_POCKETFFT_SHA512 in cmake/RttPocketfft.cmake with REF and SHA512 in
  <vcpkg root>/ports/pocketfft/portfile.cmake (CI checks vcpkg out at the builtin-baseline), and
- VCPKG_COMMIT in .github/workflows/ci.yml with "builtin-baseline" in vcpkg.json,

and exits with 1 on any difference, a missing file or a value it cannot read; it is never
silently green.

    python tools/check_pins.py <vcpkg root>

Only the standard library is used.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CMAKE_PIN = ROOT / "cmake" / "RttPocketfft.cmake"
VCPKG_JSON = ROOT / "vcpkg.json"
CI_YML = ROOT / ".github" / "workflows" / "ci.yml"


class PinError(Exception):
    """A file or value is missing, or the pins differ."""


def read(path: Path) -> str:
    if not path.is_file():
        raise PinError(f"{path}: file not found")
    return path.read_text(encoding="utf-8")


def one(pattern: str, text: str, what: str) -> str:
    """The single capture group of `pattern` in `text`; exactly one match, else PinError."""
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if len(matches) != 1:
        raise PinError(f"{what}: expected exactly one match of {pattern!r}, found {len(matches)}")
    return matches[0]


def check(vcpkg_root: Path) -> list[str]:
    """Differences as messages; raises PinError if something cannot be read."""
    cmake = read(CMAKE_PIN)
    ours_ref = one(r'^set\(RTT_POCKETFFT_REF "([0-9a-f]{40})"\)', cmake, f"{CMAKE_PIN} REF")
    ours_sha = one(r'^set\(RTT_POCKETFFT_SHA512 "([0-9a-f]{128})"\)', cmake,
                   f"{CMAKE_PIN} SHA512")

    portfile = vcpkg_root / "ports" / "pocketfft" / "portfile.cmake"
    port = read(portfile)
    port_ref = one(r"^\s*REF\s+([0-9a-f]{40})\s*$", port, f"{portfile} REF")
    port_sha = one(r"^\s*SHA512\s+([0-9a-f]{128})\s*$", port, f"{portfile} SHA512")

    try:
        baseline = json.loads(read(VCPKG_JSON))["builtin-baseline"]
    except (json.JSONDecodeError, KeyError) as e:
        raise PinError(f"{VCPKG_JSON}: no builtin-baseline ({e})") from e
    ci_commit = one(r"^\s*VCPKG_COMMIT:\s*([0-9a-f]{40})\s*$", read(CI_YML), f"{CI_YML} VCPKG_COMMIT")

    diffs = []
    if ours_ref != port_ref:
        diffs.append(f"pocketfft REF: {CMAKE_PIN.name} {ours_ref} != port {port_ref}")
    if ours_sha != port_sha:
        diffs.append(f"pocketfft SHA512: {CMAKE_PIN.name} {ours_sha} != port {port_sha}")
    if baseline != ci_commit:
        diffs.append(f"vcpkg: builtin-baseline {baseline} != VCPKG_COMMIT {ci_commit} in ci.yml")
    return diffs


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    try:
        diffs = check(Path(argv[1]))
    except PinError as e:
        print(f"check_pins: error: {e}", file=sys.stderr)
        return 1
    for d in diffs:
        print(f"check_pins: {d}", file=sys.stderr)
    if diffs:
        return 1
    print("check_pins: pocketfft pin equals the vcpkg port; vcpkg baseline equals VCPKG_COMMIT")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
