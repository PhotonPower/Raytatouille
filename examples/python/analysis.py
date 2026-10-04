"""Analyses of the cemented achromat N-BK7/F2 with the SCHOTT test catalogue, optionally plotted.

Run from the repository root after building the package (pip install .[plot]):

    python examples/python/analysis.py            # prints the results
    python examples/python/analysis.py out.png    # also saves the plots (needs matplotlib)
"""

from __future__ import annotations

import sys
from pathlib import Path

import raytatouille as rt
from raytatouille import analysis as an

REPO_ROOT = Path(__file__).resolve().parents[2]


def main(plot_file: str | None = None) -> int:
    lib = rt.MaterialLibrary()
    lib.add_catalog(REPO_ROOT / "tests" / "catalogs" / "schott.agf")
    # Compile once: every analysis below works on the same CompiledSystem.
    achromat = rt.load(REPO_ROOT / "tests" / "reference" / "m2" / "achromat.rtt.json")
    system = rt.compile(achromat, lib)

    fo = rt.paraxial.first_order(system)
    print(f"EFL {fo.efl:.4f} mm (absolute indices, AIR = Ciddor air)")
    seidel = an.seidel(system, pair=(0, 2))
    print(f"Seidel S_I {seidel.sum.s1:.3e} mm, C_L {seidel.sum.c_l:.3e} mm (F - C)")

    spot = an.spot(system, path="main", field=2, rays="hexapolar:12")  # polychromatic
    print(f"Spot field 2: RMS {1000 * spot.stats.rms_centroid:.2f} um, "
          f"{spot.rays_arrived} of {spot.rays_launched} rays")
    colour = an.longitudinal_colour(system)
    print(f"Longitudinal colour F - C: paraxial {colour.paraxial * 1000:.2f} um, "
          f"real {colour.real * 1000:.2f} um")
    distortion = an.distortion(system, samples=5)
    print(f"Distortion at full field {distortion.percent[-1]:.4f} %")
    curvature = an.field_curvature(system, samples=5)
    print(f"Astigmatism at full field {curvature.astigmatism[-1] * 1000:.1f} um")
    opd = an.opd_map(system, field=0, grid=21)
    print(f"OPD on axis: RMS {opd.rms:.4f} waves, PV {opd.pv:.4f} waves")

    if plot_file is not None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from raytatouille import plot

        fig, axes = plt.subplots(2, 3, figsize=(15, 9))
        plot.spot(spot, axes[0, 0])
        plot.ray_fan(an.ray_fan(system, field=2), axes[0, 1])
        plot.opd_map(opd, axes[0, 2])
        plot.distortion(distortion, axes[1, 0])
        plot.field_curvature(curvature, axes[1, 1])
        plot.longitudinal_colour(colour, axes[1, 2])
        fig.tight_layout()
        fig.savefig(plot_file)
        plt.close(fig)
        print(f"Plots in {plot_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
