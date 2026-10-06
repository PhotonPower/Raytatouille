"""Material library from Python: catalogue alias, glass listing, glass map, dispersion curve.

Loads the two test catalogues that both are called SCHOTT (tests/catalogs/schott.agf and
tests/catalogs/m2/schott.agf) side by side with an alias, lists their glasses as written in the
AGF files, collects the n_d/v_d values for a glass map and evaluates a dispersion curve over a
NumPy array of wavelengths (#85). Manufacturer catalogues load the same way.

Run from the repository root after building the package (pip install .):

    python examples/python/materials.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOGS = REPO_ROOT / "tests" / "catalogs"


def main(catalogs: Path = CATALOGS) -> int:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalogs / "schott.agf")  # catalogue SCHOTT (file name)
    lib.add_catalog(catalogs / "m2" / "schott.agf", name="SCHOTT_M2")  # same stem: alias
    print("Catalogues:", ", ".join(lib.catalogs()))

    for catalog in lib.catalogs():
        for g in lib.glasses(catalog):
            low, high = g.wavelength_range_um or (float("nan"), float("nan"))
            print(f"  {g.reference:18s} formula {g.formula}  n_d {g.nd:.5f}  v_d {g.vd:6.2f}"
                  f"  {low:.2f}-{high:.2f} um  density {g.density_g_per_cm3} g/cm^3")

    chart = lib.glass_map()
    print("Glass map:", list(zip(chart.reference, chart.vd.tolist(), chart.nd.tolist())))

    # Dispersion curve of N-BK7 at 20 degC and 1 atm: one call for the whole array, element by
    # element identical with the scalar call.
    wavelengths = np.linspace(0.4, 0.7, 7)
    curve = lib.index("SCHOTT:N-BK7", wavelengths)
    scalar = np.array([lib.index("SCHOTT:N-BK7", float(wl)) for wl in wavelengths])
    print("N-BK7 n(lambda):", " ".join(f"{n.real:.5f}" for n in curve))
    return 0 if curve.tobytes() == scalar.tobytes() else 1


if __name__ == "__main__":
    sys.exit(main())
