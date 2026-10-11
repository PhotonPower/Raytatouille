"""Wavefront of the Cooke triplet: pupil function, Zernike fit and Strehl approximation per
field (ADR 0032; #202, #203, Python #209).

For each field of tests/reference/m2/cooke_triplet.rtt.json at the reference wavelength the
pupil function on 48 x 48 cells is fitted with 37 Noll terms. The table shows defocus (Z4),
astigmatism (Z5, Z6), coma (Z7, Z8) and spherical aberration (Z11) in waves, the RMS about the
diffraction focus and the Strehl approximation after Marechal, exp(-(2 pi sigma_bf)^2) (good
down to about 0.1). A through-focus line for the axis follows: refocus() moves the reference
point along the chief ray without tracing again.

main() returns 1 if a value is not finite or a Strehl value lies outside (0, 1].

Run from the repository root after building the package (pip install .):

    python examples/python/wavefront.py
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import raytatouille as rt
from raytatouille import analysis as an

REPO_ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    library = rt.MaterialLibrary()
    library.add_catalog(REPO_ROOT / "tests" / "catalogs" / "m2" / "schott.agf")
    cs = rt.compile(rt.load(REPO_ROOT / "tests" / "reference" / "m2" / "cooke_triplet.rtt.json"),
                    library)
    ok = True
    print("field     Z4      Z5      Z6      Z7      Z8     Z11   sigma_bf  Strehl")
    for field in range(cs.field_count):
        with warnings.catch_warnings():
            # The mapping of the wide field is nonlinear (pupil.mapping_nonlinear); it does not
            # matter for the Zernike fit.
            warnings.simplefilter("ignore", rt.errors.RaytatouilleWarning)
            pupil = an.pupil_function(cs, field=field, grid=48)
        fit = an.zernike_fit(pupil)
        c = fit.coefficients
        print(f"{field:5d} " + " ".join(f"{c[j - 1]:7.3f}" for j in (4, 5, 6, 7, 8, 11))
              + f"  {fit.sigma_bf:8.4f}  {fit.strehl:6.3f}")
        ok = ok and all(math.isfinite(float(x)) for x in c) and 0.0 < fit.strehl <= 1.0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", rt.errors.RaytatouilleWarning)
        axis = an.pupil_function(cs, field=0, grid=48)
    print("focus shift (mm)   Z4 (waves)   Strehl")
    for shift in (-0.05, 0.0, 0.05):
        fit = an.zernike_fit(an.refocus(cs, axis, shift))
        print(f"{shift:16.3f}   {fit.coefficients[3]:10.4f}   {fit.strehl:6.3f}")
        ok = ok and 0.0 < fit.strehl <= 1.0
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
