"""Diffraction of an aberration-free paraboloid: PSF, MTF, encircled energy and the Strehl ratio
through focus (ADR 0033; #205, #206, Python #209).

tests/reference/m2/paraboloid_stop.rtt.json (f = 100 mm) with the entrance pupil diameter set
to 10 mm (F/10) images the axial point without aberration, so the results are those of the
circular pupil: the PSF is the Airy pattern with its first zero at 1.2197 lambda F#_w, the MTF
is (2 / pi) [acos x - x sqrt(1 - x^2)] for x = nu / nu_c, and the encircled energy inside the
first zero is 83.8 %. F#_w = 1 / (2 |M_11|) comes from the pupil mapping. The through-focus
line uses D = 2 mm (F/50), where the Strehl ratio of a defocus of w waves is sinc^2(w).

main() returns 1 if a value lies farther than 1e-3 (relative) from the circular pupil or the
Strehl ratio at the focus is not 1. With a directory argument it also saves the plots (needs
matplotlib).

Run from the repository root after building the package (pip install .):

    python examples/python/diffraction.py [plot_dir]
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

import raytatouille as rt
from raytatouille import analysis as an

REPO_ROOT = Path(__file__).resolve().parents[2]


def paraboloid(epd: float) -> rt.CompiledSystem:
    data = json.loads(rt.load(REPO_ROOT / "tests" / "reference" / "m2" /
                              "paraboloid_stop.rtt.json").to_json())
    data["aperture"]["value"] = epd
    return rt.compile(rt.System.from_json(json.dumps(data)))


def mtf_circle(x: float) -> float:
    return 2.0 / math.pi * (math.acos(x) - x * math.sqrt(1.0 - x * x))


def main(plot_dir: str | None = None) -> int:
    cs = paraboloid(10.0)
    pupil = an.pupil_function(cs, grid=128)
    assert pupil.mapping is not None
    lambda_mm = pupil.wavelength_um * 1e-3
    unit = lambda_mm / (2.0 * abs(pupil.mapping.m[0, 0]))  # lambda F#_w, mm
    ok = True

    psf = an.psf(cs, grid=64)
    r = np.linspace(1.0, 1.45, 451) * unit
    profile = an.psf_at(pupil, np.column_stack([r, np.zeros_like(r)]))
    zero = float(r[int(np.argmin(profile))]) / unit
    print(f"PSF: L = {psf.length}, Strehl {psf.strehl:.12f}, first zero {zero:.4f} lambda F#_w")
    ok = ok and abs(psf.strehl - 1.0) <= 1e-9 and abs(zero - 1.2197) <= 1e-3

    mtf = an.mtf(cs, grid=128)
    nu_c = 1.0 / unit
    print("x = nu / nu_c   MTF tangential   circular pupil")
    for x in (0.25, 0.5, 0.75):
        value = float(np.interp(x * nu_c, mtf.tangential.frequency, mtf.tangential.modulation))
        print(f"{x:12.2f}   {value:16.6f}   {mtf_circle(x):14.6f}")
        ok = ok and abs(value / mtf_circle(x) - 1.0) <= 1e-3

    r0 = 1.2197 * unit
    ee = an.encircled_energy(cs, [0.5 * unit, r0, 2.0 * unit], grid=128)
    print("encircled energy at 0.5, 1.22 and 2 lambda F#_w: "
          + ", ".join(f"{e:.4f}" for e in ee.energy))
    ok = ok and 0.83 < float(ee.energy[1]) < 0.84

    slow = paraboloid(2.0)
    tf = an.through_focus(slow, shifts=[0.0, 0.5, 1.0, 2.0], grid=64)
    print("focus shift (mm)   Strehl")
    for shift, strehl in zip(tf.shifts, tf.strehl):
        print(f"{shift:16.2f}   {strehl:6.4f}")
    ok = ok and abs(float(tf.strehl[0]) - 1.0) <= 1e-9 and float(tf.strehl[-1]) < 1.0

    if plot_dir is not None:
        from raytatouille import plot

        out = Path(plot_dir)
        for name, draw in (("psf", lambda: plot.psf(psf)), ("mtf", lambda: plot.mtf(mtf)),
                           ("encircled_energy", lambda: plot.encircled_energy(ee)),
                           ("through_focus", lambda: plot.through_focus(tf))):
            ax = draw()
            ax.figure.savefig(out / f"paraboloid_{name}.png")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
