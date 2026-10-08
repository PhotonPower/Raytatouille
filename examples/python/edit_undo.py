"""Change a system with JSON Patch, undo and redo, and replay the history (ADR 0024).

Run from the repository root after building the package (pip install .):

    python examples/python/edit_undo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]


def main(file: Path = REPO_ROOT / "tests" / "reference" / "m0" / "singlet.rtt.json") -> int:
    system = rt.load(file)
    ed = rt.Editor(system)

    # Stable anchors: the pointer of a surface by its id, not by an index.
    s1 = ed.system.locate_surface("L1.S1")
    assert s1 is not None
    radius = s1 + "/shape/base/radius"
    print(f"{s1}: radius {ed.system.json_at(radius)}")

    ed.set(radius, 50.0)  # a Param: only its value changes, "variable" stays
    ed.insert("/wavelengths/-", {"um": 0.55})
    print(f"after two commands: radius {ed.system.json_at(radius)['value']} mm, "
          f"{len(ed.system.wavelengths)} wavelengths")

    # A command that would add an error is rejected and changes nothing.
    try:
        ed.set(s1 + "/id", "L1.S2")
    except rt.EditError as e:
        print(f"rejected: [{e.code}] {e.location}")

    ed.undo()
    print(f"undo: {len(ed.system.wavelengths)} wavelengths, redo possible: {ed.can_redo}")

    # The history is plain JSON: the base file and the patches (script commands).
    replay = rt.Editor.from_history(ed.export_history())
    same = replay.system.to_json() == ed.system.to_json()
    print(f"replayed {len(replay.history)} patch(es), same system: {same}")
    for patch in ed.history:
        print("  ", patch)
    return 0 if same else 1


if __name__ == "__main__":
    sys.exit(main())
