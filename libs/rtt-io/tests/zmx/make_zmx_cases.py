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
    base = dict(type="STANDARD", curv=0.0, coni=0.0, parm=[], disz=0.0, glas="", diam=(0.0, 0, 0),
                clap=None, flap=None, stop=False)
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
        name="Achromat", unit="MM", catalog="schott.agf", aperture=("FNUM", 5.0), ftyp=0,
        fields=[(0.0, 0.0, 1.0), (0.0, 2.0, 1.0), (1.0, 3.0, 1.0)],
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
}

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
        # alpha_i in L^(1 - 2i): the value in the unit of the file is alpha_mm * s^(2i - 1)
        p["parm_s"] = [e(a * s ** (2 * i + 1)) for i, a in enumerate(sf["parm"])]
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
    obj = dict(surfaces[0], curv_s=e(0.0, 18), diam_s=zero, clap_s=None, flap_s=None, parm_s=[])
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

    sys: dict = {"schema_version": "0.9.0", "name": c["name"],
                 "units": {"length": "mm", "wavelength": "um"},
                 "environment": {"temperature_c": 20.0, "pressure_atm": 1.0, "medium": "AIR"}}
    if c["object_s"] is not None:
        sys["object"] = {"at_infinity": False, "distance": float(c["object_s"]) * s}
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

    def shape(sf) -> dict:
        cv = float(sf["curv_s"]) / s
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


def first_order(c: dict, glasses: dict, lam_ref: float, nair: float) -> dict:
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
        n = index("AIR")
        sum_abs = 0.0
        surfaces = c["surfaces"][:-1]
        z = 0.0
        z_last = 0.0
        y1 = None
        for k, sf in enumerate(surfaces):
            cv = float(sf["curv_s"]) / s
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
    return {"n_air": nair, "efl": a["efl"], "bfl": a["bfl"],
            "rear_focal_length": a["rear_focal_length"], "efl_relative": r["efl"],
            "bfl_relative": r["bfl"], "condition": a["condition"]}


def main() -> None:
    results = {}
    for name, case in CASES.items():
        glasses = read_sellmeier(CATALOGS / case["catalog"])
        c = printed(case)
        write_zmx(name, c)
        system, data = build_expected(c, glasses)
        (HERE / f"{name}.rtt.json").write_text(json.dumps(system, indent=2) + "\n",
                                               encoding="utf-8", newline="\n")
        results[name] = first_order(c, glasses, data["lam_ref"], data["nair"])
    (HERE / "first_order.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8",
                                           newline="\n")


if __name__ == "__main__":
    main()
