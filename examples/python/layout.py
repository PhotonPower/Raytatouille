"""Draw a lens section from the geometry export (raytatouille.layout, #81).

Prints the surfaces of the Cooke triplet with their global vertex, shape and media, and, with
an output file name, draws the y-z section: glass outlines of the lenses and the profiles of
the other surfaces (the stop; the image plane has no aperture and is skipped).

Run from the repository root after building the package (pip install .):

    python examples/python/layout.py [out.png]
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
FILE = REPO_ROOT / "tests" / "reference" / "m2" / "cooke_triplet.rtt.json"
CATALOG = REPO_ROOT / "tests" / "catalogs" / "m2" / "schott.agf"


def main(png: str | None = None, file: Path = FILE, catalog: Path = CATALOG) -> int:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog)
    compiled = rt.compile(rt.load(file), lib)
    media = [m.reference for m in compiled.media]
    for s in rt.layout.surfaces(compiled):
        radius = f"R = {1.0 / s.curvature:9.3f}" if s.curvature != 0.0 else "plane      "
        print(f"{s.id:6s} {s.kind:9s} z = {s.translation[2]:8.3f} mm  {radius}  "
              f"{media[s.medium_front]} | {media[s.medium_back]}")

    elements = rt.layout.elements(compiled)
    outlines = [o for i in range(len(elements)) for o in rt.layout.outlines(compiled, i)]
    print(f"{len(outlines)} lens outlines in the y-z plane")
    if not all((o[0] == o[-1]).all() for o in outlines):
        return 1

    if png is not None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 3))
        for o in outlines:
            ax.fill(o[:, 2], o[:, 1], color="lightsteelblue", edgecolor="k", linewidth=0.8)
        for i, e in enumerate(elements):
            if not e.segmented:
                for s in range(e.first_surface, e.first_surface + e.surface_count):
                    try:
                        pieces = rt.layout.profile(compiled, s)
                    except ValueError:
                        continue  # unbounded surface (no aperture)
                    for p in pieces:
                        ax.plot(p[:, 2], p[:, 1], color="k", linewidth=1.2)
        ax.set_aspect("equal")
        ax.set_xlabel("z in mm")
        ax.set_ylabel("y in mm")
        fig.savefig(png, dpi=120, bbox_inches="tight")
        plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
