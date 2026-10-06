"""Record the paths of rays for drawing them in a layout (trace(..., record_path=True), #80).

Traces a hexapolar bundle per field through the singlet, records every 4th ray (an explicit
selection, not the first rays of the bundle, which would be the inner rings only) and prints
the y-z points of each recorded ray, with the surface where it was lost, if any. Optionally
draws the y-z section with matplotlib.

Run from the repository root after building the package (pip install .):

    python examples/python/ray_paths.py [out.png]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
FILE = REPO_ROOT / "tests" / "reference" / "m1" / "singlet_const.rtt.json"


def main(png: str | None = None, file: Path = FILE) -> int:
    compiled = rt.compile(rt.load(file))
    rays = rt.trace.make_rays(compiled, rt.trace.HexapolarPupil(rings=3))
    selection = np.arange(0, len(rays), 4)
    stats, paths = rt.trace.trace(compiled, rays, record_path=True, record_rays=selection)
    names = [compiled.surface_ids[int(s)] for s in paths.event_surfaces]
    print(f"{len(rays)} rays traced, {paths.ray_count} recorded with {paths.slots} slots "
          f"(start, then {', '.join(names)})")

    ok = True
    for r in range(paths.ray_count):
        i = int(paths.ray_indices[r])
        last = int(paths.count[r]) - 1
        points = " -> ".join(f"({y:.3f}, {z:.3f})" for _, y, z in paths.position[r, : last + 1])
        lost = int(paths.lost_at[r])
        where = f", lost at {names[lost]}" if lost >= 0 else ""
        print(f"  ray {i} (field {rays.field[i]}): {points}{where}")
        # The last recorded slot is the final state of the ray in the batch.
        ok &= bool(paths.position[r, last, 1] == rays.pos_y[i])

    if png is not None:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 3))
        for r in range(paths.ray_count):
            last = int(paths.count[r])
            ax.plot(paths.position[r, :last, 2], paths.position[r, :last, 1],
                    color=f"C{rays.field[int(paths.ray_indices[r])]}", linewidth=0.8)
        ax.set_xlabel("z in mm")
        ax.set_ylabel("y in mm")
        fig.savefig(png, dpi=120, bbox_inches="tight")
        plt.close(fig)
    print(f"alive {stats.count(rt.trace.RayStatus.ALIVE)}, "
          f"vignetted {stats.count(rt.trace.RayStatus.VIGNETTED)}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
