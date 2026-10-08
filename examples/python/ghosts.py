"""Ghosts of the Cooke triplet and the arms of a Michelson interferometer (#122 to #124).

First derives the two-reflection ghosts of the triplet's path "main" (raytatouille
.compile_with_ghosts, ADR 0027) and ranks them by their irradiance at the image relative to the
useful image (raytatouille.analysis.ghost_ranking). The rank value depends on the resolution
radius r0 of the detector, a model choice; the table shows the three strongest ghosts.

Then evaluates the Michelson interferometer of tests/reference/m4 with start rays of its own
(a folded path, which make_rays cannot aim): the transmission of both arms to the camera and
their optical path difference, twice the 7.5 mm by which the test arm is longer.

Run from the repository root after building the package (pip install .):

    python examples/python/ghosts.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE = REPO_ROOT / "tests" / "reference"
CATALOG = REPO_ROOT / "tests" / "catalogs" / "m2" / "schott.agf"


def ghosts() -> None:
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG)
    g = rt.compile_with_ghosts(rt.load(REFERENCE / "m2" / "cooke_triplet.rtt.json"), "main",
                               materials=schott)
    ranking = rt.analysis.ghost_ranking(g, field=0, resolution_radius=0.01)
    names = g.system.path_names
    e = ranking.entries
    print(f"{len(e)} ghosts of 'main', r0 = {ranking.resolution_radius * 1e3:.0f} um, "
          f"useful image: P = {ranking.base_power:.4f}, "
          f"r = {ranking.base_rms_radius * 1e3:.1f} um")
    print("  ghost                         rel. irradiance  rel. power  RMS radius")
    for k in range(3):
        print(f"  {names[e.path[k]]:<30}{e.relative_irradiance[k]:>12.3e}"
              f"{e.relative_power[k]:>12.3e}{e.rms_radius[k]:>10.3f} mm")


def michelson() -> None:
    compiled = rt.compile(rt.load(REFERENCE / "m4" / "michelson_offset.rtt.json"))
    start = rt.trace.RayBatch(3)  # three rays along +z: on the axis and 1 mm off it in x, y
    start.pos_x[1] = 1.0
    start.pos_y[2] = 1.0
    for arm in ("reference arm", "test arm"):
        t = rt.analysis.path_transmission(compiled, arm, start=start)
        print(f"{arm}: transmission {t.mean:.3f} ({t.rays_arrived} of {t.rays_launched} rays)")
    d = rt.analysis.opl_difference(compiled, "reference arm", "test arm", start=start)
    print(f"OPL(test arm) - OPL(reference arm) = {d.chief:.6f} mm")


def main() -> int:
    ghosts()
    michelson()
    return 0


if __name__ == "__main__":
    sys.exit(main())
