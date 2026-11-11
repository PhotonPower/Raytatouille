#!/usr/bin/env python3
"""Writes the ZMX test cases of the import (ADR 0044) and their expected results.

Each case is described once below, as an optical prescription. From it the script writes
- <case>.zmx: the file as OpticStudio would print it (lengths with 13 significant digits,
  curvatures with 19, in the unit of the case);
- <case>.rtt.json: the system the import must produce, built here from the prescription by the
  rules of ADR 0044 and NOT with the C++ code (acceptance 1);
- first_order.json: EFL, BFL and the rear focal length from an own paraxial y-nu trace
  (Greivenkamp, OPTI-201/202, p. 9-2) with own indices: Sellmeier (AGF formula 2, OpticStudio
  User Guide) from the catalogues in tests/catalogs and air after Ciddor (NIST, Eqs. (A21) to
  (A41), dry air, 450 umol/mol CO2 as the AIR material of rtt-material), absolute and relative
  to air (acceptance 3).

Only the Python standard library. Run from anywhere; it rewrites the files next to it:
    python libs/rtt-io/tests/zmx/make_zmx_cases.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CATALOGS = HERE.parents[3] / "tests" / "catalogs"

# ----------------------------------------------------------------------------- air and glass


def ciddor_air_index(wavelength_um: float, temperature_c: float, pressure_atm: float) -> float:
    """Phase index of dry air after Ciddor (NIST "Engineering Metrology Toolbox", Appendix
    A-III, Eqs. (A21)-(A41) with x_v = 0) at a vacuum wavelength in um."""
    p = pressure_atm * 101325.0  # (A27) 1 atm in Pa
    k0, k1, k2, k3 = 238.0185, 5792105.0, 57.362, 167917.0  # (A22)
    a0, a1, a2 = 1.58123e-6, -2.9331e-8, 1.1043e-10  # (A23)
    d = 1.83e-11  # (A26)
    p_r1, t_r1, z_a, r_gas = 101325.0, 288.15, 0.9995922115, 8.314472  # (A27), (A28), (A30)
    s = 1.0 / (wavelength_um * wavelength_um)  # (A31)
    r_as = 1e-8 * (k1 / (k0 - s) + k3 / (k2 - s))  # (A32)
    x_co2 = 450.0
    m_a = 0.0289635 + 1.2011e-8 * (x_co2 - 400.0)  # (A34)
    r_axs = r_as * (1.0 + 5.34e-7 * (x_co2 - 450.0))  # (A35)
    t = temperature_c
    t_k = t + 273.15  # (A36)
    p_over_t = p / t_k
    z_m = 1.0 - p_over_t * (a0 + a1 * t + a2 * t * t) + p_over_t * p_over_t * d  # (A37)
    rho_axs = p_r1 * m_a / (z_a * r_gas * t_r1)  # (A38)
    rho_a = p * m_a / (z_m * r_gas * t_k)  # (A40)
    return 1.0 + (rho_a / rho_axs) * r_axs  # (A41)


def n_air(vacuum_um: float) -> float:
    return ciddor_air_index(vacuum_um, 20.0, 1.0)


def vacuum_wavelength(air_um: float) -> float:
    """lambda_vac = lambda_air n_air(lambda_vac) at 20 degC and 1 atm, fixed point from
    lambda_air (ADR 0044, point 4)."""
    v = air_um
    for _ in range(10):
        nxt = air_um * n_air(v)
        if nxt == v:
            break
        v = nxt
    return v


def read_sellmeier(catalog: Path) -> dict[str, tuple[list[float], float]]:
    """Glasses of an AGF file with formula 2: name -> (K1 L1 K2 L2 K3 L3, T_ref in degC)."""
    out: dict[str, tuple[list[float], float]] = {}
    name = ""
    formula = 0
    for line in catalog.read_text(encoding="latin-1").splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "NM":
            name, formula = f[1], int(float(f[2]))
        elif f[0] == "CD" and formula == 2:
            out[name] = ([float(x) for x in f[1:7]], 20.0)
        elif f[0] == "TD" and name in out:
            out[name] = (out[name][0], float(f[7]))
    return out


def n_rel(coeffs: list[float], lam_um: float) -> float:
    """Sellmeier 1 (AGF formula 2): n^2 - 1 = sum K_i lam^2 / (lam^2 - L_i), relative to air."""
    l2 = lam_um * lam_um
    s = 0.0
    for k, lc in zip(coeffs[0::2], coeffs[1::2]):
        s += k * l2 / (l2 - lc)
    return math.sqrt(1.0 + s)


# ---------------------------------------------------------------------------- the cases

UNITS = {"MM": 1.0, "CM": 10.0, "IN": 25.4, "METER": 1000.0}


def surf(**kw):
    base = dict(type="STANDARD", curv=0.0, coni=0.0, parm=[], xdat=[], disz=0.0, glas="",
                diam=(0.0, 0, 0), clap=None, flap=None, sqap=None, elap=None, stop=False)
    base.update(kw)
    return base


def singlet(unit: str) -> dict:
    return dict(
        name="Singlet", unit=unit, catalog="schott.agf", aperture=("ENPD", 10.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 5.0, 1.0), (0.0, 7.0, 0.5)],
        waves=[(0.4861327, 1.0), (0.5875618, 1.0), (0.6562725, 1.0)], pwav=2, object=None,
        surfaces=[
            surf(disz=5.0, diam=(5.0, 0, 0), stop=True),
            surf(curv=1.0 / 50.0, disz=6.0, glas="N-BK7", diam=(8.0, 1, 0), flap=8.0),
            surf(curv=-1.0 / 50.0, disz=45.0, diam=(8.0, 2, 2), flap=8.0),
            surf(),
        ],
        merit=True, thic=(3, 45.0))


def asphere(unit: str) -> dict:
    # Even asphere with c != 0 (even_asphere) and a plane asphere (odd_asphere, c = 0); the stop
    # is the front of the plate with a fixed semi-diameter and a FLAP; aperture float by stop.
    return dict(
        name="Aspheres", unit=unit, catalog="schott.agf", aperture=("FLOA",), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 3.0, 1.0), (0.0, 5.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(type="EVENASPH", parm=[0.0, 2.0e-6, -1.0e-9], disz=2.0, glas="N-BK7",
                 diam=(6.0, 1, 0), flap=6.0, stop=True),
            surf(disz=10.0, diam=(6.0, 0, 0)),
            surf(type="EVENASPH", curv=1.0 / 30.0, coni=-0.8, parm=[0.0, -1.5e-6, 3.0e-10, -2e-13],
                 disz=4.0, glas="N-BK7", diam=(8.0, 0, 0)),
            surf(curv=-1.0 / 120.0, disz=40.0, diam=(8.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None)


CASES = {
    "singlet": singlet("MM"),
    "singlet_in": singlet("IN"),
    "singlet_m": singlet("METER"),
    # Cemented achromat N-BK7 / F2 behind a dropped dummy, a CLAP dummy (thin element) and a stop
    # with a fixed semi-diameter; F/# aperture; the first lens surface has a FLAP next to an
    # automatic DIAM (dropped with a warning).
    "doublet": dict(
        name="Achromat", unit="MM", catalog="schott.agf", aperture=("FNUM", 5.0), ftyp=2,
        fields=[(0.0, 0.0, 1.0), (0.0, 3.0, 1.0), (1.5, 4.0, 1.0)],
        waves=[(0.4861327, 1.0), (0.5875618, 2.0), (0.6562725, 1.0)], pwav=2, object=None,
        surfaces=[
            surf(disz=2.0),
            surf(disz=3.0, clap=(0.0, 12.0)),
            surf(disz=4.0, diam=(10.5, 1, 0), stop=True),
            surf(curv=1.0 / 61.47, disz=6.0, glas="N-BK7", diam=(12.7, 0, 0), flap=12.7),
            surf(curv=-1.0 / 44.64, disz=2.5, glas="F2", diam=(12.7, 0, 0)),
            surf(curv=-1.0 / 129.94, disz=95.0, diam=(12.7, 0, 0)),
            surf(),
        ],
        merit=False, thic=None),
    # Parabolic mirror with the stop on the mirror (aperture.stop, ADR 0039); thickness -100
    # after the mirror (as tests/reference/m1/paraboloid_mirror.rtt.json).
    "mirror": dict(
        name="Paraboloid", unit="MM", catalog="schott.agf", aperture=("ENPD", 20.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 0.5, 1.0)], waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(curv=-0.005, coni=-1.0, disz=-100.0, glas="MIRROR", clap=(0.0, 10.5), stop=True),
            surf(),
        ],
        merit=False, thic=None),
    # Finite object 100 mm in front of surface 1, object-space NA, object heights.
    "finite": dict(
        name="Finite NA", unit="MM", catalog="schott.agf", aperture=("OBNA", 0.05), ftyp=1,
        fields=[(0.0, 0.0, 1.0), (0.0, 2.0, 1.0), (0.0, 3.0, 1.0)],
        waves=[(0.5875618, 1.0), (0.6562725, 1.0)], pwav=1, object=100.0,
        surfaces=[
            surf(disz=2.0, diam=(5.5, 0, 0), stop=True),
            surf(curv=1.0 / 40.0, disz=5.0, glas="N-BK7", diam=(9.0, 0, 0)),
            surf(curv=-1.0 / 40.0, disz=60.0, diam=(9.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None),
    "asphere": asphere("MM"),
    "asphere_in": asphere("IN"),
    "asphere_cm": asphere("CM"),
    # Object 10 mm inside N-BK7 (object medium, ADR 0037), surface 1 starts the lens; the stop is
    # the back surface with a semi-diameter picked up from surface 1; object heights in cm.
    "immersion": dict(
        name="Immersion", unit="CM", catalog="schott.agf", aperture=("OBNA", 0.1), ftyp=1,
        fields=[(0.0, 0.0, 1.0), (0.0, 0.5, 1.0), (0.0, 1.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=10.0, object_glass="N-BK7",
        surfaces=[
            surf(disz=4.0, glas="F2", diam=(2.0, 1, 0)),
            surf(curv=-1.0 / 4.0, disz=22.0, diam=(2.0, 2, 1), stop=True),
            surf(),
        ],
        merit=False, thic=None),
    # Cooke triplet after Geiser (Sasian, OPTI 517, L20; docs/quellen.md), stop on L2.S1 with
    # radius 10 as tests/reference/m2/cooke_triplet_surface_stop.rtt.json (acceptance 4).
    "sasian": dict(
        name="Sasian triplet", unit="MM", catalog="m2/schott.agf", aperture=("ENPD", 12.5), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 14.0, 1.0), (0.0, 20.0, 1.0)],
        waves=[(0.4861327, 1.0), (0.5875618, 1.0), (0.6562725, 1.0)], pwav=2, object=None,
        surfaces=[
            surf(curv=1.0 / 23.713, disz=4.831, glas="N-LAK9", diam=(10.0, 0, 0)),
            surf(curv=1.0 / 7331.288, disz=5.86, diam=(10.0, 0, 0)),
            surf(curv=-1.0 / 24.456, disz=0.975, glas="N-SF5", diam=(10.0, 1, 0),
                 clap=(0.0, 10.0), stop=True),
            surf(curv=1.0 / 21.896, disz=4.822, diam=(10.0, 0, 0)),
            surf(curv=1.0 / 86.759, disz=3.127, glas="N-LAK9", diam=(10.0, 0, 0)),
            surf(curv=-1.0 / 20.4942, disz=41.10346, diam=(10.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None),
    # PR C: ideal lens f = 100 as the stop (PARAXIAL, OPD mode 2), image in its focal plane.
    "paraxial": dict(
        name="Paraxial lens", unit="MM", catalog="schott.agf", aperture=("ENPD", 10.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 2.0, 1.0), (0.0, 4.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(type="PARAXIAL", parm=[100.0, 2.0], disz=100.0, diam=(5.0, 1, 0), stop=True),
            surf(),
        ],
        merit=False, thic=None),
    # PR C: rectangular and elliptical apertures on dummy surfaces in front of a singlet, in
    # inches (the half widths scale with the unit).
    "apertures": dict(
        name="Apertures", unit="IN", catalog="schott.agf", aperture=("ENPD", 10.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 3.0, 1.0), (0.0, 5.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=2.0, sqap=(7.0, 6.0)),
            surf(disz=3.0, elap=(7.0, 6.0)),
            surf(disz=5.0, diam=(8.0, 1, 0), stop=True),
            surf(curv=1.0 / 50.0, disz=6.0, glas="N-BK7", diam=(9.0, 1, 0), flap=9.0),
            surf(curv=-1.0 / 50.0, disz=45.0, diam=(9.0, 2, 4), flap=9.0),
            surf(),
        ],
        merit=False, thic=None),
}


def odd_asphere(unit: str) -> dict:
    # #360: an odd asphere with c != 0 and the terms beta2 ... beta8 (UG "Odd Asphere"); beta2
    # adds 2 beta2 to the vertex curvature, beta1 = 0 keeps a paraxial model.
    return dict(
        name="Odd asphere", unit=unit, catalog="schott.agf", aperture=("ENPD", 8.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 2.0, 1.0), (0.0, 4.0, 1.0)],
        waves=[(0.4861327, 1.0), (0.5875618, 1.0), (0.6562725, 1.0)], pwav=2, object=None,
        surfaces=[
            surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
            surf(type="ODDASPHE", curv=1.0 / 45.0, coni=-0.6,
                 parm=[0.0, 2e-4, 1e-6, -3e-7, 2e-9, 1e-10, -5e-12, 1e-13], disz=6.0,
                 glas="N-BK7", diam=(8.0, 1, 0), flap=8.0),
            surf(curv=-1.0 / 80.0, disz=50.0, diam=(8.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None)


def toroid(unit: str) -> dict:
    # #360: a toroid (UG "Toroidal": radius of rotation in parameter 1, alpha1 ... alpha7 in
    # 2 ... 8) and a cylinder (radius of rotation 0), with fields in x and y.
    return dict(
        name="Toroid", unit=unit, catalog="schott.agf", aperture=("ENPD", 8.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 3.0, 1.0), (2.0, 0.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
            surf(type="TOROIDAL", curv=1.0 / 50.0, parm=[80.0, 0.0, -1e-6, 2e-9], disz=5.0,
                 glas="N-BK7", diam=(7.0, 1, 0), flap=7.0),
            surf(disz=4.0, diam=(7.0, 0, 0)),
            surf(type="TOROIDAL", curv=1.0 / 60.0, coni=-0.3, parm=[0.0, 0.0, 5e-7], disz=4.0,
                 glas="N-BK7", diam=(7.0, 0, 0)),
            surf(disz=60.0, diam=(7.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None)


CASES["oddasph"] = odd_asphere("MM")
CASES["oddasph_in"] = odd_asphere("IN")
# A cone at the vertex (beta1 != 0) on a plane base with a conic, which has no effect at c = 0
# and is dropped: imported, without a paraxial model.
CASES["oddasph_cone"] = dict(
    odd_asphere("MM"), name="Odd asphere cone",
    surfaces=[
        surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
        surf(type="ODDASPHE", coni=-2.0, parm=[0.02, 0.0, 1e-5], disz=6.0, glas="N-BK7",
             diam=(8.0, 1, 0), flap=8.0),
        surf(curv=-1.0 / 30.0, disz=40.0, diam=(8.0, 0, 0)),
        surf(),
    ])
CASES["toroid"] = toroid("MM")
CASES["toroid_in"] = toroid("IN")


def biconic(unit: str) -> dict:
    # #360, PR B: a biconic (UG "Biconic": parameter 1 = Rx, 0 = infinite; 2 = kx; Ry and ky in
    # the radius and conic columns, decision H3) and one with Rx = 0, a cylinder in y.
    return dict(
        name="Biconic", unit=unit, catalog="schott.agf", aperture=("ENPD", 8.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 3.0, 1.0), (2.0, 0.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
            surf(type="BICONICX", curv=1.0 / 50.0, coni=-0.5, parm=[70.0, 0.3], disz=5.0,
                 glas="N-BK7", diam=(7.0, 1, 0), flap=7.0),
            surf(disz=4.0, diam=(7.0, 0, 0)),
            surf(type="BICONICX", curv=1.0 / 60.0, coni=-0.2, parm=[0.0, 0.0], disz=4.0,
                 glas="N-BK7", diam=(7.0, 0, 0)),
            surf(disz=60.0, diam=(7.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None)


def qtype(unit: str, kind: int, curv: float, coni: float, coeffs: list[float]) -> dict:
    # #360, PR B: a Q-type asphere (UG "Q-Type Asphere (sequential)": extra data 13 = type,
    # 14 = maximum term number, 15 = normalization radius, 16 ff. = A0, A1, ...; written as
    # XDAT 1, 2, ..., decision H4) on the front of a singlet.
    return dict(
        name="Q type " + ("con" if kind else "bfs"), unit=unit, catalog="schott.agf",
        aperture=("ENPD", 8.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 2.0, 1.0), (0.0, 4.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
            surf(type="QED_TYPE", curv=curv, coni=coni,
                 xdat=[float(kind), float(len(coeffs)), 8.0] + coeffs, disz=6.0, glas="N-BK7",
                 diam=(8.0, 1, 0), flap=8.0),
            surf(curv=-1.0 / 80.0, disz=50.0, diam=(8.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None)


CASES["biconic"] = biconic("MM")
CASES["biconic_in"] = biconic("IN")
CASES["qcon"] = qtype("MM", 1, 1.0 / 40.0, -0.2, [1e-3, -2e-4, 5e-5])
CASES["qcon_in"] = qtype("IN", 1, 1.0 / 40.0, -0.2, [1e-3, -2e-4, 5e-5])
CASES["qbfs"] = qtype("MM", 0, 1.0 / 45.0, 0.0, [2e-3, -1e-4])

# #360, decision H2: a toroid with an infinite y-z radius has power in x only; it is imported as
# a cylinder turned by 90 deg about z, with a turned frame after it. Its expected data are the
# first order per section, not an expected file (the pose of the untilt rounds).
ROTATED = {
    # a biconic with an infinite Ry: power in x only, the cylinder with Rx and kx turned
    "xbiconic": dict(
        name="Biconic in x", unit="MM", catalog="schott.agf", aperture=("ENPD", 8.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (2.0, 0.0, 1.0), (0.0, 2.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
            surf(type="BICONICX", parm=[60.0, -0.4], disz=5.0, glas="N-BK7", diam=(7.0, 1, 0),
                 sqap=(7.0, 5.0)),
            surf(curv=-1.0 / 200.0, disz=90.0, diam=(7.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None),
    "xcylinder": dict(
        name="Cylinder in x", unit="MM", catalog="schott.agf", aperture=("ENPD", 8.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (2.0, 0.0, 1.0), (0.0, 2.0, 1.0)],
        waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=3.0, diam=(4.0, 1, 0), stop=True),
            surf(type="TOROIDAL", parm=[70.0], disz=5.0, glas="N-BK7", diam=(7.0, 1, 0),
                 elap=(7.0, 5.0)),
            surf(curv=-1.0 / 200.0, disz=90.0, diam=(7.0, 0, 0)),
            surf(),
        ],
        merit=False, thic=None),
}


def cb(dx=0.0, dy=0.0, a=0.0, b=0.0, c=0.0, order=0, disz=0.0, **kw):
    """A coordinate break (UG "Coordinate Break"): decenters in mm, tilts in degrees."""
    return surf(type="COORDBRK", parm=[dx, dy, a, b, c, float(order)], disz=disz, **kw)


# PR C, folded and tilted systems: the expected data are the frames of OpticStudio (below), not
# an expected file, because the poses of a decomposed transform round.
FOLDED = {
    # tests/reference/m6/fold45_focus.rtt.json as ZMX: ideal lens f = 50 (the stop), fold +45
    # deg about x at z = 20, detector 32 behind the fold (2 mm behind the focus) at (0, 32, 20).
    "fold": dict(
        name="Fold 45", unit="MM", catalog="schott.agf", aperture=("ENPD", 10.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 1.0, 1.0)], waves=[(0.55, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(type="PARAXIAL", parm=[50.0, 2.0], disz=20.0, diam=(5.0, 1, 0), stop=True),
            cb(a=45.0),
            surf(glas="MIRROR", disz=0.0),
            cb(a=45.0, disz=-32.0),
            surf(),
        ],
        merit=False, thic=None),
    # Off-axis parabola: decenter -10 in y (flag 0), the parabola R = -200, then three breaks in
    # a row (a compound tilt, its reversal with flag 1, the decenter undone with flag 1) and the
    # image in the plane of the focus at z = 0.
    "oap": dict(
        name="Off-axis parabola", unit="MM", catalog="schott.agf", aperture=("ENPD", 6.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 0.2, 1.0)], waves=[(0.55, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=100.0, diam=(5.0, 1, 0), stop=True),
            cb(dy=-10.0),
            surf(curv=-0.005, coni=-1.0, glas="MIRROR", disz=0.0),
            cb(a=10.0, b=5.0),
            cb(a=-10.0, b=-5.0, order=1),
            cb(dy=10.0, order=1, disz=-100.0),
            surf(),
        ],
        merit=False, thic=None),
    # Pivoted lens after the OISC file of the same name: a break of 5 deg at surface 1, a
    # cemented lens, the way back to the pivot with a negative thickness and the break undone;
    # finite object, object heights.
    "pivot": dict(
        name="Pivoted lens", unit="MM", catalog="schott.agf", aperture=("ENPD", 8.0), ftyp=1,
        fields=[(0.0, 0.0, 1.0), (0.0, 5.0, 1.0)], waves=[(0.5875618, 1.0)], pwav=1,
        object=200.0,
        surfaces=[
            cb(a=5.0),
            surf(curv=1.0 / 40.0, disz=6.0, glas="N-BK7", diam=(8.0, 1, 0), flap=8.0, stop=True),
            surf(curv=-1.0 / 30.0, disz=3.0, glas="F2", diam=(8.0, 1, 0), flap=8.0),
            surf(curv=-1.0 / 100.0, disz=-9.0, diam=(8.0, 1, 0), flap=8.0),
            cb(a=-5.0, disz=150.0),
            surf(),
        ],
        merit=False, thic=None),
    # A singlet with a tilted image plane (Scheimpflug, TILTSURF at the image).
    "tilt_image": dict(
        name="Tilted image", unit="MM", catalog="schott.agf", aperture=("ENPD", 10.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 3.0, 1.0)], waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=5.0, diam=(5.0, 1, 0), stop=True),
            surf(curv=1.0 / 50.0, disz=6.0, glas="N-BK7", diam=(8.0, 1, 0), flap=8.0),
            surf(curv=-1.0 / 50.0, disz=45.0, diam=(8.0, 1, 0), flap=8.0),
            surf(type="TILTSURF", parm=[0.1, -0.3]),
        ],
        merit=False, thic=None),
    # A plate with a tilted entrance face (TILTSURF, the first surface of the lens) and a plane
    # exit face; the exit face takes the tilt out exactly (rotate_first, PR C point 4).
    "tilt_plate": dict(
        name="Tilted plate", unit="MM", catalog="schott.agf", aperture=("ENPD", 6.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 1.0, 1.0)], waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=5.0, diam=(5.0, 1, 0), stop=True),
            surf(type="TILTSURF", parm=[0.05, 0.2], disz=5.0, glas="N-BK7"),
            surf(disz=20.0),
            surf(),
        ],
        merit=False, thic=None),
    # A prism: tilted entrance face, then a break of -3 deg inside the glass (it takes the medium
    # of the surface before it), the exit face in the new frame.
    "tilt_prism": dict(
        name="Tilted prism", unit="MM", catalog="schott.agf", aperture=("ENPD", 6.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 1.0, 1.0)], waves=[(0.5875618, 1.0)], pwav=1, object=None,
        surfaces=[
            surf(disz=5.0, diam=(5.0, 1, 0), stop=True),
            surf(type="TILTSURF", parm=[0.05, 0.2], disz=5.0, glas="N-BK7"),
            cb(a=-3.0),
            surf(disz=20.0),
            surf(),
        ],
        merit=False, thic=None),
}


def chained_breaks(name: str, breaks: list) -> dict:
    """Breaks in a row (multiplied out and decomposed, PR C point 3, rule 3) in front of a thin
    element with a CLAP, whose frame then checks the decomposition, and the image behind it."""
    return dict(
        name=name, unit="MM", catalog="schott.agf", aperture=("ENPD", 4.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0)], waves=[(0.55, 1.0)], pwav=1, object=None,
        surfaces=[surf(disz=5.0, diam=(2.0, 1, 0), stop=True)] + breaks + [
            surf(disz=10.0, clap=(0.0, 20.0)),
            surf(),
        ],
        merit=False, thic=None)


def random_breaks() -> list:
    """Three breaks with angles from a fixed seed (Python's Mersenne Twister, seed 317), flag 0
    and flag 1 alternating, printed with 13 digits like any other value."""
    import random
    rng = random.Random(317)
    out = []
    for i in range(3):
        a, b, c = (round(rng.uniform(-80.0, 80.0), 6) for _ in range(3))
        dx, dy = (round(rng.uniform(-3.0, 3.0), 6) for _ in range(2))
        out.append(cb(dx=dx, dy=dy, a=a, b=b, c=c, order=i % 2, disz=1.0))
    return out


# Ry(45) Ry(45) = Ry(90) and Rx(20) Ry(45) Ry(45): b = 90 deg, where only a + c is defined, built
# by composition as R02 = 1 +- 1 ulp (the case of the review of the proposal).
FOLDED["gimbal"] = chained_breaks("Gimbal", [cb(b=45.0), cb(b=45.0, disz=3.0)])
FOLDED["gimbal_x"] = chained_breaks("Gimbal x", [cb(a=20.0), cb(b=45.0), cb(b=45.0, disz=3.0)])
FOLDED["random_breaks"] = chained_breaks("Random breaks", random_breaks())

# --------------------------------------------------------------------- frames of OpticStudio


def mat_mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mat_vec(a, v):
    return [sum(a[i][k] * v[k] for k in range(3)) for i in range(3)]


def rot(axis: int, deg: float):
    """Right-handed rotation about the x, y or z axis (UG "Coordinate Break": "right-handed with
    respect to the positive axes")."""
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    i, j = [(1, 2), (2, 0), (0, 1)][axis]
    m[i][i], m[i][j], m[j][i], m[j][j] = c, -s, s, c
    return m


class Frame:
    """Origin and axes (columns of r) of the current coordinate system of OpticStudio, global."""

    def __init__(self):
        self.r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self.o = [0.0, 0.0, 0.0]

    def shift(self, v):  # translate along the current axes
        d = mat_vec(self.r, v)
        self.o = [self.o[i] + d[i] for i in range(3)]

    def turn(self, m):  # rotate the current axes (intrinsic)
        self.r = mat_mul(self.r, m)


def frames(c: dict) -> dict:
    """The frame of every surface after UG "Coordinate Break" (flag 0: decenter, tilt about x, the
    new y, the new z; flag != 0: tilt about z, the new y, the new x, then decenter; the thickness
    last, along the new z), "Tilted" (a plane z = x tan_x + y tan_y that keeps the frame: only its
    normal is given) and "Thickness" (a mirror does not turn the frame). Surface 1 is the origin."""
    s = UNITS[c["unit"]]
    f = Frame()
    out = {}
    for k, sf in enumerate(c["surfaces"], start=1):
        p = [float(x) for x in sf["parm_s"]]
        if sf["type"] == "COORDBRK":
            dx, dy, a, b, cz, order = p[0] * s, p[1] * s, p[2], p[3], p[4], p[5]
            if order == 0.0:
                f.shift([dx, dy, 0.0])
                f.turn(rot(0, a))
                f.turn(rot(1, b))
                f.turn(rot(2, cz))
            else:
                f.turn(rot(2, cz))
                f.turn(rot(1, b))
                f.turn(rot(0, a))
                f.shift([dx, dy, 0.0])
        else:
            entry = {"origin": list(f.o)}
            if sf["type"] == "TILTSURF":
                n = [-p[0], -p[1], 1.0]
                norm = math.sqrt(sum(x * x for x in n))
                entry["normal"] = mat_vec(f.r, [x / norm for x in n])
            else:
                entry["axes"] = [list(row) for row in f.r]  # rows of the matrix, axes in columns
            out[f"S{k}"] = entry
        f.shift([0.0, 0.0, float(sf["disz_s"]) * s])
    return out


def refract(d, n, n1, n2):
    """Vector Snell's law for a unit direction d at a unit normal n of either orientation:
    t = mu d + (cos_t - mu cos_i) m, m the normal on the side d points to, mu = n1 / n2."""
    cos_i = sum(d[i] * n[i] for i in range(3))
    m = n if cos_i > 0.0 else [-x for x in n]
    cos_i = abs(cos_i)
    mu = n1 / n2
    cos_t = math.sqrt(1.0 - mu * mu * (1.0 - cos_i * cos_i))
    return [mu * d[i] + (cos_t - mu * cos_i) * m[i] for i in range(3)]


def to_plane(p, d, origin, normal):
    t = (sum((origin[i] - p[i]) * normal[i] for i in range(3)) /
         sum(d[i] * normal[i] for i in range(3)))
    return [p[i] + t * d[i] for i in range(3)]


def plate_ray(c: dict, fr: dict, glasses: dict) -> dict:
    """A ray from (0, 1, -1) along +z through the tilted entrance face, the glass and the exit
    face to the image plane, with the own absolute indices above (Sellmeier, Ciddor)."""
    lam = vacuum_wavelength(float(e(c["waves"][0][0])))
    na = n_air(lam)
    ng = n_rel(glasses["N-BK7"][0], lam / na) * na
    entry, exit_, img = (fr[k] for k in sorted(fr, key=lambda k: int(k[1:]))[1:])
    p, d = [0.0, 1.0, -1.0], [0.0, 0.0, 1.0]
    p = to_plane(p, d, entry["origin"], entry["normal"])
    d = refract(d, entry["normal"], na, ng)
    z_exit = [row[2] for row in exit_["axes"]]
    p = to_plane(p, d, exit_["origin"], z_exit)
    d = refract(d, z_exit, ng, na)
    p = to_plane(p, d, img["origin"], [row[2] for row in img["axes"]])
    return {"start": [0.0, 1.0, -1.0], "hit": p}


# ------------------------------------------------------------------------------ printing


def e(x: float, digits: int = 12) -> str:
    """OpticStudio style: mantissa with `digits` decimals, exponent with three digits."""
    m, ex = f"{x:.{digits}E}".split("E")
    sign = "-" if ex[0] == "-" else "+"
    return f"{m}E{sign}{abs(int(ex)):03d}"


def printed(case: dict) -> dict:
    """The case in its unit, every number as printed in the file and read back as a double."""
    s = UNITS[case["unit"]]
    out = dict(case)
    out["surfaces"] = []
    for sf in case["surfaces"]:
        p = dict(sf)
        p["curv_s"] = e(sf["curv"] * s, 18)
        p["disz_s"] = e(sf["disz"] / s)
        p["diam_s"] = e(sf["diam"][0] / s)
        p["clap_s"] = None if sf["clap"] is None else (e(sf["clap"][0] / s), e(sf["clap"][1] / s))
        p["flap_s"] = None if sf["flap"] is None else e(sf["flap"] / s)
        if sf["type"] == "COORDBRK":
            # decenters in lens units, tilts in degrees, the order flag (UG "Coordinate Break")
            p["parm_s"] = [e(v / s) if i < 2 else e(v) for i, v in enumerate(sf["parm"])]
        elif sf["type"] == "PARAXIAL":
            # focal length in lens units, OPD mode (UG "Paraxial")
            p["parm_s"] = [e(sf["parm"][0] / s), e(sf["parm"][1])]
        elif sf["type"] == "TILTSURF":
            p["parm_s"] = [e(v) for v in sf["parm"]]  # tangents (UG "Tilted")
        elif sf["type"] == "ODDASPHE":
            # beta_n in L^(1 - n): the value in the unit of the file is beta_mm * s^(n - 1)
            p["parm_s"] = [e(b * s ** i) for i, b in enumerate(sf["parm"])]
        elif sf["type"] == "BICONICX":
            # Rx in lens units, kx (UG "Biconic")
            p["parm_s"] = [e(sf["parm"][0] / s), e(sf["parm"][1])]
        elif sf["type"] == "TOROIDAL":
            # radius of rotation in lens units, then alpha_i in L^(1 - 2i) as EVENASPH
            p["parm_s"] = [e(sf["parm"][0] / s)] + [
                e(a * s ** (2 * i - 1)) for i, a in enumerate(sf["parm"][1:], start=1)]
        else:
            # alpha_i in L^(1 - 2i): the value in the unit of the file is alpha_mm * s^(2i - 1)
            p["parm_s"] = [e(a * s ** (2 * i + 1)) for i, a in enumerate(sf["parm"])]
        for k in ("sqap", "elap"):
            p[k + "_s"] = None if sf[k] is None else (e(sf[k][0] / s), e(sf[k][1] / s))
        # QED_TYPE extra data: type, number of terms, normalization radius (lens units), then
        # the coefficients A0, A1, ... (lens units)
        x = sf["xdat"]
        p["xdat_s"] = [] if not x else [e(x[0]), e(x[1]), e(x[2] / s)] + [e(a / s) for a in x[3:]]
        out["surfaces"].append(p)
    ap = case["aperture"]
    out["aperture_s"] = (ap[0],) if ap[0] == "FLOA" else (
        ap[0], e(ap[1] / s) if ap[0] == "ENPD" else e(ap[1]))
    out["fields_s"] = [(e(x / s) if case["ftyp"] in (1, 2) else e(x),
                        e(y / s) if case["ftyp"] in (1, 2) else e(y), e(w)) for x, y, w in case["fields"]]
    out["object_s"] = None if case["object"] is None else e(case["object"] / s)
    return out


def write_zmx(name: str, c: dict) -> None:
    n_f, n_w = len(c["fields"]), len(c["waves"])
    pad = lambda vals, fill: vals + [fill] * (12 - len(vals))  # noqa: E731
    L = ["VERS 60801 67 21410", "MODE SEQ", f"NAME {c['name']}", "NOTE 0 written by make_zmx_cases.py",
         f"UNIT {c['unit']} X W X CM MR"]
    ap = c["aperture_s"]
    L.append({"ENPD": lambda: f"ENPD {ap[1]}", "FNUM": lambda: f"FNUM {ap[1]} 0",
              "OBNA": lambda: f"OBNA {ap[1]} 0", "FLOA": lambda: "FLOA"}[ap[0]]())
    L += ["PUPD 0.000000000000E+000 0.000000000000E+000",
          f"ENVD {e(20.0)} {e(1.0)} 0", f"GFAC {e(0.0)} 0", "GCAT SCHOTT",
          "RAIM 1.0E-8 1 1 1 0 0 0 0", f"FTYP {c['ftyp']} 0 {n_f} {n_w} 0 0 0", "ROPD 2", "PICB 1"]
    zero = e(0.0)
    L.append("XFLN " + " ".join(pad([f[0] for f in c["fields_s"]], zero)))
    L.append("YFLN " + " ".join(pad([f[1] for f in c["fields_s"]], zero)))
    L.append("FWGN " + " ".join(pad([f[2] for f in c["fields_s"]], e(1.0))))
    for k in ("VDXN", "VDYN", "VCXN", "VCYN", "VANN"):
        L.append(f"{k} " + " ".join([zero] * 12))
    for i in range(24):
        wl, w = c["waves"][i] if i < n_w else (0.55, 1.0)
        L.append(f"WAVM {i + 1} {e(wl)} {e(w)}")
    L += [f"PWAV {c['pwav']}", "POLS 1 0 1 0 0 1 0",
          "COFN COATING.DAT SCATTER_PROFILE.DAT ABG_DATA.DAT PROFILE.GRD"]
    surfaces = [surf()] + c["surfaces"]
    obj = dict(surfaces[0], curv_s=e(0.0, 18), diam_s=zero, clap_s=None, flap_s=None, parm_s=[],
               xdat_s=[],
               sqap_s=None, elap_s=None, glas=c.get("object_glass", ""))
    obj["disz_s"] = "INFINITY" if c["object_s"] is None else c["object_s"]
    surfaces[0] = obj
    for k, sf in enumerate(surfaces):
        L.append(f"SURF {k}")
        if sf["stop"]:
            L.append("  STOP")
        L.append(f"  TYPE {sf['type']}")
        L.append(f"  CURV {sf['curv_s']} 0 {zero} {zero} 0")
        for i, a in enumerate(sf["parm_s"]):
            L.append(f"  PARM {i + 1} {a}")
        for i, a in enumerate(sf.get("xdat_s", [])):
            L.append(f'  XDAT {i + 1} {a} 0 0 1.000000000000E+00 0.000000000000E+00 0 ""')
        L += ["  HIDE 0 0 0 0 0 0 0 0", "  MIRR 2 0.000000000E+000", f"  SLAB {k + 1}"]
        if sf["coni"] != 0.0:
            L.append(f"  CONI {e(sf['coni'])}")
        L.append(f"  DISZ {sf['disz_s']}")
        if sf["glas"] == "MIRROR":
            L.append("  GLAS MIRROR 0 0 1.5 0 0 0 0 0 0 0")
        elif sf["glas"]:
            L.append(f"  GLAS {sf['glas']} 0 0 1.50000000 40.00000000 0 0 0 0 0 0")
        d = sf["diam"]
        L.append(f"  DIAM {sf['diam_s']} {d[1]} {d[2]} 0 1.000000000000E+000")
        L.append("  POPS 0 0 0 0 0 0 0 0 1 1 1.000000000000E+000 1.000000000000E+000 0 0 0")
        if sf["clap_s"] is not None:
            L.append(f"  CLAP {sf['clap_s'][0]} {sf['clap_s'][1]} 0")
        if sf["flap_s"] is not None:
            L.append(f"  FLAP {zero} {sf['flap_s']} 0")
        for k in ("sqap", "elap"):
            if sf.get(k + "_s") is not None:
                L.append(f"  {k.upper()} {sf[k + '_s'][0]} {sf[k + '_s'][1]} 0")
    if c["merit"]:
        L += ["EFFL 0 2 0 0 0 0 0 0 1.0E+002", "BLNK Default merit function",
              "DMFS 0 0 0 0 0 0 0 0", "TRAC 0 1 0 0 0.3 0.4 0 0.1", "OPDX 0 1 0 0 0.3 0.4 0 0.1"]
    L += ["TOL TOFF 0 0 0 0 0 0 0", "TOL TWAV 0 0 0.55 0 0 0 0", "MNUM 1",
          'MOFF 0 1 "" 0 0 0 1 1 0 0.0 ']
    if c["thic"] is not None:
        surf_no, value = c["thic"]
        L.append(f"THIC {surf_no} 1 {e(value / UNITS[c['unit']])} 0 0 0 1 1 0")
    (HERE / f"{name}.zmx").write_bytes(("\r\n".join(L) + "\r\n").encode("ascii"))


# --------------------------------------------------------------------- expected system


def build_expected(c: dict, glasses: dict) -> tuple[dict, dict]:
    """The system of ADR 0044 for the case, from the printed numbers, and the data of the
    y-nu trace (indices per space)."""
    s = UNITS[c["unit"]]
    inv_s = 1.0 / s
    waves = []
    for i, (wl_s, w) in enumerate(c["waves"]):
        wl = float(e(wl_s))
        waves.append({"um": vacuum_wavelength(wl), "weight": float(e(w))} |
                     ({"reference": True} if i + 1 == c["pwav"] else {}))
    lam_ref = waves[c["pwav"] - 1]["um"]
    nair = n_air(lam_ref)

    sys: dict = {"schema_version": "0.10.0", "name": c["name"],
                 "units": {"length": "mm", "wavelength": "um"},
                 "environment": {"temperature_c": 20.0, "pressure_atm": 1.0, "medium": "AIR"}}
    if c["object_s"] is not None:
        sys["object"] = {"at_infinity": False, "distance": float(c["object_s"]) * s}
        if c.get("object_glass"):
            sys["object"]["medium"] = f"SCHOTT:{c['object_glass']}"
    sys["wavelengths"] = waves
    ap = c["aperture_s"]
    if ap[0] == "ENPD":
        sys["aperture"] = {"type": "epd", "value": float(ap[1]) * s}
    elif ap[0] == "FNUM":
        sys["aperture"] = {"type": "image_fnumber", "value": float(ap[1]) / nair}
    elif ap[0] == "OBNA":
        sys["aperture"] = {"type": "object_na", "value": nair * float(ap[1]), "reading": "real"}
    else:
        sys["aperture"] = {"type": "stop_size"}
    ftype = {0: "angle_deg", 1: "object_height", 2: "paraxial_image_height"}[c["ftyp"]]
    pts = []
    for x, y, w in c["fields_s"]:
        k = s if c["ftyp"] in (1, 2) else 1.0
        pts.append({"x": float(x) * k, "y": float(y) * k, "weight": float(w)})
    sys["fields"] = {"type": ftype, "points": pts}

    def scaled(v: float, n: int) -> float:
        """v times 1/s, n times (as the import does it, no pow)."""
        for _ in range(n):
            v *= inv_s
        return v

    def shape(sf) -> dict:
        cv = float(sf["curv_s"]) / s
        if sf["type"] == "ODDASPHE":
            # beta_n in L^(1 - n): (n - 1) multiplications with 1/s
            betas = [scaled(float(b), i) for i, b in enumerate(sf["parm_s"])]
            while betas and betas[-1] == 0.0:
                betas.pop()
            if betas:
                # the conic has no effect at c = 0 and is dropped there
                return {"type": "odd_asphere", "curvature": cv,
                        "conic": float(e(sf["coni"])) if cv != 0.0 else 0.0,
                        "coefficients": betas}
        if sf["type"] == "BICONICX":
            rx = float(sf["parm_s"][0]) * s
            kx = float(sf["parm_s"][1])
            assert cv != 0.0  # power in x only: the ROTATED cases
            if rx == 0.0:
                return {"type": "cylinder", "radius": 1.0 / cv, "conic": float(e(sf["coni"])),
                        "coefficients": []}
            return {"type": "biconic", "radius_x": rx, "radius_y": 1.0 / cv, "conic_x": kx,
                    "conic_y": float(e(sf["coni"]))}
        if sf["type"] == "QED_TYPE":
            x = [float(v) for v in sf["xdat_s"]]
            n = int(x[1])
            coeffs = [v * s for v in x[3:3 + n]]
            while coeffs and coeffs[-1] == 0.0:
                coeffs.pop()
            return {"type": "qcon" if x[0] == 1.0 else "qbfs", "curvature": cv,
                    "conic": float(e(sf["coni"])), "normalization_radius": x[2] * s,
                    "coefficients": coeffs}
        if sf["type"] == "TOROIDAL":
            rot = float(sf["parm_s"][0]) * s
            alphas = [scaled(float(a), 2 * i - 1) for i, a in enumerate(sf["parm_s"][1:], start=1)]
            assert not alphas or alphas[0] == 0.0  # alpha1 (y^2) is an error, not a case here
            coeffs = alphas[1:]
            while coeffs and coeffs[-1] == 0.0:
                coeffs.pop()
            assert cv != 0.0  # power in x only: the ROTATED cases
            profile = {"radius": 1.0 / cv, "conic": float(e(sf["coni"])), "coefficients": coeffs}
            if rot == 0.0:
                return {"type": "cylinder"} | profile
            return {"type": "toroid"} | profile | {"radius_of_rotation": rot}
        if sf["type"] == "EVENASPH" and any(float(a) != 0.0 for a in sf["parm_s"]):
            alphas = []
            for i, a in enumerate(sf["parm_s"]):
                x = float(a)
                for _ in range(2 * (i + 1) - 1):
                    x *= inv_s
                alphas.append(x)
            while alphas and alphas[-1] == 0.0:
                alphas.pop()
            assert alphas[0] == 0.0  # alpha1 (r^2) is an error, not a case here
            if cv != 0.0:
                return {"type": "even_asphere", "radius": 1.0 / cv, "conic": float(e(sf["coni"])),
                        "coefficients": alphas[1:]}
            coeffs = []
            for a in alphas[1:]:
                coeffs += [0.0, a]
            return {"type": "odd_asphere", "curvature": 0.0, "conic": 0.0,
                    "coefficients": [0.0, 0.0, 0.0] + coeffs[1:]}
        if cv == 0.0:
            return {"type": "plane"}
        return {"type": "conic", "radius": 1.0 / cv, "conic": float(e(sf["coni"]))}

    def aperture(sf) -> dict | None:
        if sf["sqap_s"] is not None:
            return {"type": "rectangular", "half_width_x": float(sf["sqap_s"][0]) * s,
                    "half_width_y": float(sf["sqap_s"][1]) * s}
        if sf["elap_s"] is not None:
            return {"type": "elliptical", "semi_axis_x": float(sf["elap_s"][0]) * s,
                    "semi_axis_y": float(sf["elap_s"][1]) * s}
        if sf["clap_s"] is not None:
            a = {"type": "circular", "radius": float(sf["clap_s"][1]) * s}
            if float(sf["clap_s"][0]) != 0.0:
                a["inner_radius"] = float(sf["clap_s"][0]) * s
            return a
        if sf["flap_s"] is not None and sf["diam"][1] in (1, 2):
            return {"type": "circular", "radius": float(sf["flap_s"]) * s}
        return None

    def stop_radius(sf) -> float:
        a = aperture(sf)
        if a is not None:
            return a["radius"]
        v = float(sf["diam_s"]) * s
        if c["aperture"][0] == "FLOA" or sf["diam"][1] in (1, 2):
            return v
        return 1.1 * v

    children = []
    counters = {"L": 0, "M": 0, "T": 0}
    surfaces = c["surfaces"]
    pending = 0.0
    first = True
    spaces = []  # (thickness mm, medium) for the y-nu trace, plus surface data
    lens = None

    def place(node: dict, z: float) -> None:
        nonlocal first
        if first:
            node["pose"] = {"position": [0.0, 0.0, z]}
            first = False
        else:
            node["pose"] = {"reference": "relative_to_preceding", "position": [0.0, 0.0, z]}

    def surface_node(k: int, sf) -> dict:
        node = {"id": f"S{k}", "shape": {"base": shape(sf)}}
        a = aperture(sf)
        if sf["stop"] and a is None:
            a = {"type": "circular", "radius": stop_radius(sf)}
        if a is not None:
            node["aperture"] = a
        return node

    before = "AIR"
    for k, sf in enumerate(surfaces, start=1):
        is_image = k == len(surfaces)
        after = "AIR" if (sf["glas"] in ("", "MIRROR")) else sf["glas"]
        disz = float(sf["disz_s"]) * s
        if is_image:
            node = {"type": "detector", "name": "image", "surfaces": [surface_node(k, sf)]}
            place(node, pending)
            children.append(node)
            break
        if lens is not None:
            sn = surface_node(k, sf)
            sn["pose"] = {"reference": "relative_to_preceding", "position": [0.0, 0.0, pending]}
            lens["surfaces"].append(sn)
            if sf["stop"]:
                sys["aperture"]["stop"] = {"surface": f"S{k}"}
            if after == "AIR":
                mats = lens.pop("_mats")
                if all(m == mats[0] for m in mats):
                    lens["material"] = mats[0]
                else:
                    lens["material"] = mats
                lens = None
            else:
                lens["_mats"].append(f"SCHOTT:{after}")
        elif sf["type"] == "PARAXIAL" and float(sf["parm_s"][0]) != 0.0:
            # ideal lens on a thin element (ADR 0031), f in mm; the stop sits on its surface
            counters["T"] += 1
            sn = surface_node(k, sf)
            sn["interaction"] = {"type": "ideal_lens", "focal_length": float(sf["parm_s"][0]) * s}
            node = {"type": "thin_element", "name": f"T{counters['T']}", "surfaces": [sn]}
            place(node, pending)
            children.append(node)
            if sf["stop"]:
                sys["aperture"]["stop"] = {"surface": f"S{k}"}
        elif sf["glas"] == "MIRROR":
            counters["M"] += 1
            sn = surface_node(k, sf)
            sn["interaction"] = {"type": "ideal_mirror"}
            node = {"type": "mirror", "name": f"M{counters['M']}", "surfaces": [sn]}
            place(node, pending)
            children.append(node)
            if sf["stop"]:
                sys["aperture"]["stop"] = {"surface": f"S{k}"}
        elif after != "AIR":
            counters["L"] += 1
            lens = {"type": "lens", "name": f"L{counters['L']}", "surfaces": [surface_node(k, sf)],
                    "_mats": [f"SCHOTT:{after}"]}
            place(lens, pending)
            children.append(lens)
            if sf["stop"]:
                sys["aperture"]["stop"] = {"surface": f"S{k}"}
        elif sf["stop"]:
            node = {"type": "stop", "name": "stop", "surfaces": [surface_node(k, sf)]}
            place(node, pending)
            children.append(node)
        elif aperture(sf) is not None:
            counters["T"] += 1
            node = {"type": "thin_element", "name": f"T{counters['T']}",
                    "surfaces": [surface_node(k, sf)]}
            place(node, pending)
            children.append(node)
        else:
            pending += disz
            continue
        pending = disz
    sys["root"] = {"type": "assembly", "name": "system", "children": children}
    sys["paths"] = [{"name": "main", "events": "auto"}]
    return sys, {"lam_ref": lam_ref, "nair": nair}


# ------------------------------------------------------------------------------ y-nu trace


def vertex_curvature(sf: dict, s: float, section: str | None) -> float:
    """Paraxial curvature of a surface at its vertex, 1/mm, in the y-z or x-z section.
    ODDASPHE: c + 2 beta2 (beta1 = 0). TOROIDAL: c in y-z; in x-z the circle of radius R about
    the axis at z = R, so 1/R, and 0 for R = 0 (cylinder)."""
    cv = float(sf["curv_s"]) / s
    if sf["type"] == "ODDASPHE":
        betas = [float(b) * s ** (-i) for i, b in enumerate(sf["parm_s"])]
        assert betas[0] == 0.0
        return cv + 2.0 * (betas[1] if len(betas) > 1 else 0.0)
    if sf["type"] == "BICONICX" and section == "x":
        rx = float(sf["parm_s"][0]) * s
        return 0.0 if rx == 0.0 else 1.0 / rx
    if sf["type"] == "QED_TYPE":
        # Q-con: the polynomial term starts at u^4, the vertex curvature is c (ADR 0043)
        assert float(sf["xdat_s"][0]) == 1.0
        return cv
    if sf["type"] == "TOROIDAL" and section == "x":
        rot = float(sf["parm_s"][0]) * s
        return 0.0 if rot == 0.0 else 1.0 / rot
    return cv


def first_order(c: dict, glasses: dict, lam_ref: float, nair: float,
                section: str | None = None) -> dict:
    """EFL, BFL and rear focal length from a paraxial y-nu trace (Greivenkamp p. 9-2):
    n'u' = nu - y c (n' - n), y' = y + t u'; with absolute indices (as rtt) and relative to air
    (air = 1, glass n_rel at the air wavelength). Signed indices after a mirror (n' = -n)."""
    s = UNITS[c["unit"]]
    lam_air = lam_ref / nair

    def trace(absolute: bool):
        def index(medium: str) -> float:
            if medium == "AIR":
                return nair if absolute else 1.0
            coeffs, t_ref = glasses[medium]
            assert t_ref == 20.0
            if absolute:
                return n_rel(coeffs, lam_ref / n_air(lam_ref)) * n_air(lam_ref)
            return n_rel(coeffs, lam_air)

        y, nu = 1.0, 0.0
        n = index(c.get("object_glass") or "AIR")
        sum_abs = 0.0
        surfaces = c["surfaces"][:-1]
        z = 0.0
        z_last = 0.0
        y1 = None
        for k, sf in enumerate(surfaces):
            cv = vertex_curvature(sf, s, section)
            if sf["type"] == "PARAXIAL":
                # ideal thin lens n'u' = nu - y phi with phi = |n| / f (UG "Paraxial" in air,
                # ADR 0031 with the absolute index); the medium does not change
                phi = abs(n) / (float(sf["parm_s"][0]) * s)
                if y1 is None:
                    y1 = y
                z_last = z
                sum_abs += abs(y * phi)
                nu = nu - y * phi
                t = float(sf["disz_s"]) * s
                y = y + t * nu / n
                z += t
                continue
            if sf["glas"] == "MIRROR":
                n_next = -n
            else:
                n_next = (1.0 if n > 0 else -1.0) * index(sf["glas"] or "AIR")
            phi = cv * (n_next - n)
            if phi != 0.0 or n_next != n:
                if y1 is None:
                    y1 = y
                z_last = z
            sum_abs += abs(y * phi)
            nu = nu - y * phi
            n = n_next
            t = float(sf["disz_s"]) * s
            y = y + t * nu / n
            z += t
        # y now at the image plane; back to the last vertex
        t_img = float(surfaces[-1]["disz_s"]) * s
        u = nu / n
        y_last = y - t_img * u
        power = -nu / y1
        z_focus = (z - t_img) - y_last / u
        direction = 1.0 if n > 0 else -1.0
        return {"efl": 1.0 / power, "rear_focal_length": abs(n) / power,
                "bfl": (z_focus - z_last) * direction, "condition": sum_abs / abs(y1 * power)}

    a = trace(True)
    r = trace(False)
    out = {"n_air": nair, "efl": a["efl"], "bfl": a["bfl"],
           "rear_focal_length": a["rear_focal_length"], "efl_relative": r["efl"],
           "bfl_relative": r["bfl"], "condition": a["condition"]}
    ap = c["aperture_s"]
    if ap[0] == "FNUM":
        # EPD = |EFL| / F# in OpticStudio, with its EFL relative to air (acceptance 3)
        out["epd_relative"] = abs(r["efl"]) / float(ap[1])
    if ap[0] == "OBNA":
        # paraxial marginal slope u = NA / n with the relative index of object space
        glass = c.get("object_glass")
        n_obj = n_rel(glasses[glass][0], lam_air) if glass else 1.0
        out["na_over_n_relative"] = float(ap[1]) / n_obj
    return out


def main() -> None:
    results = {}
    for name, case in CASES.items():
        glasses = read_sellmeier(CATALOGS / case["catalog"])
        c = printed(case)
        write_zmx(name, c)
        system, data = build_expected(c, glasses)
        (HERE / f"{name}.rtt.json").write_text(json.dumps(system, indent=2) + "\n",
                                               encoding="utf-8", newline="\n")
        if name in ("oddasph_cone", "qbfs"):
            continue  # no first order from here: a cone, a Q-bfs vertex curvature c + 2 S(0)/rho^2
        if name.startswith(("toroid", "biconic")):
            results[name] = {sec: first_order(c, glasses, data["lam_ref"], data["nair"], sec)
                             for sec in ("x", "y")}
        else:
            results[name] = first_order(c, glasses, data["lam_ref"], data["nair"])
    for name, case in ROTATED.items():
        glasses = read_sellmeier(CATALOGS / case["catalog"])
        c = printed(case)
        write_zmx(name, c)
        lam_ref = vacuum_wavelength(float(e(case["waves"][case["pwav"] - 1][0])))
        results[name] = {sec: first_order(c, glasses, lam_ref, n_air(lam_ref), sec)
                         for sec in ("x", "y")}
    (HERE / "first_order.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8",
                                           newline="\n")
    folded = {}
    for name, case in FOLDED.items():
        glasses = read_sellmeier(CATALOGS / case["catalog"])
        c = printed(case)
        write_zmx(name, c)
        fr = frames(c)
        folded[name] = {"frames": fr}
        if name in ("tilt_plate", "tilt_prism"):
            folded[name]["ray"] = plate_ray(c, fr, glasses)
    (HERE / "frames.json").write_text(json.dumps(folded, indent=2) + "\n", encoding="utf-8",
                                      newline="\n")


if __name__ == "__main__":
    main()
