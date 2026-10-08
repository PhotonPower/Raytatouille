"""Read the model tree of a system: assemblies, elements, surfaces and paths (ADR 0024).

The objects are immutable copies (raytatouille.model); a Param holds the value and whether the
optimizer may vary it. Run from the repository root after building the package (pip install .):

    python examples/python/model_tree.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt
from raytatouille import model

REPO_ROOT = Path(__file__).resolve().parents[2]


def param(p: model.Param) -> str:
    """A Param as text; a trailing V marks a variable."""
    return f"{p.value:g}{' V' if p.variable else ''}"


def shape(s: model.ShapeStack) -> str:
    base = s.base
    if isinstance(base, model.Plane):
        text = "plane"
    elif isinstance(base, model.Conic):
        text = f"R {param(base.radius)} mm, k {param(base.conic)}"
    else:
        text = f"asphere R {param(base.radius)} mm, {len(base.coefficients)} coefficients"
    return text + (f" + {len(s.terms)} Zernike term(s)" if s.terms else "")


def show(node: model.Assembly | model.Element, depth: int) -> None:
    indent = "  " * depth
    z = param(node.pose.position[2])
    if isinstance(node, model.Assembly):
        print(f"{indent}assembly {node.name!r} at z {z} mm")
        for child in node.children:
            show(child, depth + 1)
        return
    material = node.material or ", ".join(node.segment_materials) or "no material"
    print(f"{indent}{node.kind.name.lower()} {node.name!r} at z {z} mm, {material}")
    for s in node.surfaces:
        aperture = s.aperture
        clear = f", r {aperture.radius:g} mm" if isinstance(aperture, model.CircularAperture) else ""
        print(f"{indent}  {s.id}: {shape(s.shape)}{clear}, {type(s.interaction).__name__}")


def main(file: Path = REPO_ROOT / "tests" / "reference" / "m2" / "cooke_triplet.rtt.json") -> int:
    system = rt.load(file)
    print(system.name)
    a = system.aperture
    print(f"  aperture {a.type.name.lower()} {param(a.value)}, "
          f"{len(system.fields.points)} field point(s) ({system.fields.type.name.lower()})")
    show(system.root, 1)
    for p in system.paths:
        events = "auto" if p.automatic else " -> ".join(e.surface for e in p.events)
        print(f"  path {p.name!r}: {events}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
