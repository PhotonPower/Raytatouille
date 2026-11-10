#!/usr/bin/env python3
"""Writes the cases of the ZMX export and the surface table (ADR 0045, acceptance 4; addendum
#366, acceptance 4) and the expected files.

The case is described once below, as an optical prescription. From it the script writes
- export_singlet.rtt.json: the system that the test exports;
- export_singlet.zmx: the file the export must write, built here from the prescription by the
  rules of ADR 0045, point 4, and NOT with the C++ code. Air wavelengths and the F/# relative to
  air come from the own Ciddor of make_zmx_cases.py (NIST, Eqs. (A21) to (A41));
- export_singlet.csv: the surface table of the same system by the rules of point 7 (CRLF, a 0 in
  conic, a4 ... a16 and the inner radius left empty);
- export_fold.rtt.json and export_fold.zmx: the singlet folded by a plane mirror at 45 deg about x
  behind the lens (ADR 0045, addendum #366). The expected file follows the pattern of
  "Add Fold Mirror" in the Ansys Knowledge Base, KA-01573: a coordinate break with the tilt
  before the mirror, MIRROR, a break with the same tilt after it, and the thicknesses after the
  mirror with the opposite sign (KA-01588). The fold turns +z into +y (ADR 0028, point 4).

The test compares the two files word for word, numbers as doubles: bitwise, except the air
wavelengths and FNUM (two Ciddor implementations, <= 1e-15 relative). Python prints the shortest
representation (repr) like std::to_chars, only in another style ("50.0" for "50").

Only the Python standard library. Run from anywhere; it rewrites the files next to it:
    python libs/rtt-io/tests/zmx/make_export_cases.py
"""

from __future__ import annotations

import json
from pathlib import Path

from make_zmx_cases import n_air

HERE = Path(__file__).resolve().parent

# ----------------------------------------------------------------------------- prescription

# Vacuum wavelengths (um), weights; the second one is the reference.
WAVELENGTHS = [(0.4862, 1.0), (0.5877, 1.0), (0.6565, 0.5)]
REFERENCE = 1
FNUMBER = 4.0  # image_fnumber, absolute (model)
FIELDS = [(0.0, 0.0, 1.0), (0.0, 5.0, 1.0), (0.0, 7.0, 0.5)]  # angle_deg x, y, weight

# Surfaces in path order. "t" is the thickness to the next surface (mm).
STOP = {"radius": 5.0, "inner": 1.0, "t": 5.0}
FRONT = {"radius": 50.0, "conic": -1.0, "a": [1e-6, -2e-9], "aperture": 8.0, "t": 6.0}
BACK = {"radius": -50.0, "conic": 0.5, "aperture": 8.0, "t": 45.0}
IMAGE = {"half_x": 10.0, "half_y": 8.0}
GLASS = ("SCHOTT", "N-BK7")
# The fold: a plane mirror 20 mm behind L1.S2, turned by 45 deg about x, the image 25 mm along the
# folded axis (+y).
FOLD = {"to_mirror": 20.0, "angle": 45.0, "aperture": 10.0, "to_image": 25.0}


def system() -> dict:
    def pose(t: float, first: bool = False) -> dict:
        p = {"position": [0.0, 0.0, t]}
        if not first:
            p = {"reference": "relative_to_preceding", **p}
        return p

    stop = {
        "type": "stop",
        "name": "stop",
        "surfaces": [{
            "id": "STO",
            "shape": {"base": {"type": "plane"}},
            "aperture": {"type": "circular", "radius": STOP["radius"],
                         "inner_radius": STOP["inner"]},
        }],
        "pose": pose(0.0, first=True),
    }
    lens = {
        "type": "lens",
        "name": "L1",
        "surfaces": [
            {
                "id": "L1.S1",
                "shape": {"base": {"type": "even_asphere", "radius": FRONT["radius"],
                                   "conic": FRONT["conic"], "coefficients": FRONT["a"]}},
                "aperture": {"type": "circular", "radius": FRONT["aperture"]},
            },
            {
                "id": "L1.S2",
                "shape": {"base": {"type": "conic", "radius": BACK["radius"],
                                   "conic": BACK["conic"]}},
                "aperture": {"type": "circular", "radius": BACK["aperture"]},
                "pose": pose(FRONT["t"]),
            },
        ],
        "pose": pose(STOP["t"]),
        "material": f"{GLASS[0]}:{GLASS[1]}",
    }
    image = {
        "type": "detector",
        "name": "image",
        "surfaces": [{
            "id": "IMG",
            "shape": {"base": {"type": "plane"}},
            "aperture": {"type": "rectangular", "half_width_x": IMAGE["half_x"],
                         "half_width_y": IMAGE["half_y"]},
        }],
        "pose": pose(BACK["t"]),
    }
    return {
        "schema_version": "0.9.0",
        "name": "Export singlet",
        "units": {"length": "mm", "wavelength": "um"},
        "environment": {"temperature_c": 20.0, "pressure_atm": 1.0, "medium": "AIR"},
        "wavelengths": [
            {"um": um, "weight": w, **({"reference": True} if i == REFERENCE else {})}
            for i, (um, w) in enumerate(WAVELENGTHS)
        ],
        "aperture": {"type": "image_fnumber", "value": FNUMBER},
        "fields": {"type": "angle_deg",
                   "points": [{"x": x, "y": y, "weight": w} for x, y, w in FIELDS]},
        "root": {"type": "assembly", "name": "system", "children": [stop, lens, image]},
        "paths": [{"name": "main", "events": "auto"}],
    }


def expected_zmx(fold: bool = False) -> str:
    """ADR 0045, point 4: header, then one SURF per event, indented, CRLF."""
    ref_um = WAVELENGTHS[REFERENCE][0]
    lines = [
        "VERS 60801 67 21410",
        "MODE SEQ",
        "NAME Export fold" if fold else "NAME Export singlet",
        "UNIT MM",
        f"FNUM {FNUMBER * n_air(ref_um)!r} 0",  # F/# relative to air
        "ENVD 20 1 0",
        f"GCAT {GLASS[0]}",
        f"FTYP 0 0 {len(FIELDS)} {len(WAVELENGTHS)} 0 0 0",
        "XFLN " + " ".join(repr(f[0]) for f in FIELDS),
        "YFLN " + " ".join(repr(f[1]) for f in FIELDS),
        "FWGN " + " ".join(repr(f[2]) for f in FIELDS),
    ]
    for i, (um, w) in enumerate(WAVELENGTHS):
        lines.append(f"WAVM {i + 1} {um / n_air(um)!r} {w!r}")  # lambda_air = lambda_vac / n_air
    lines.append(f"PWAV {REFERENCE + 1}")

    def surf(number: int, body: list[str]) -> None:
        lines.append(f"SURF {number}")
        lines.extend("  " + b for b in body)

    surf(0, ["TYPE STANDARD", "CURV 0", "DISZ INFINITY", "DIAM 0 0"])
    surf(1, ["STOP", "TYPE STANDARD", "CURV 0", f"DISZ {STOP['t']!r}",
             f"DIAM {STOP['radius']!r} 1", f"CLAP {STOP['inner']!r} {STOP['radius']!r}"])
    a = FRONT["a"] + [0.0] * (7 - len(FRONT["a"]))
    surf(2, ["TYPE EVENASPH", f"CURV {1.0 / FRONT['radius']!r}", f"CONI {FRONT['conic']!r}",
             "PARM 1 0"] + [f"PARM {i + 2} {v!r}" for i, v in enumerate(a)] +
         [f"DISZ {FRONT['t']!r}", f"GLAS {GLASS[1]} 0", "DIAM 0 0",
          f"CLAP 0 {FRONT['aperture']!r}"])
    back_t = FOLD["to_mirror"] if fold else BACK["t"]
    surf(3, ["TYPE STANDARD", f"CURV {1.0 / BACK['radius']!r}", f"CONI {BACK['conic']!r}",
             f"DISZ {back_t!r}", "DIAM 0 0", f"CLAP 0 {BACK['aperture']!r}"])
    image = 4
    if fold:
        def coordinate_break(number: int, tilt_x: float, thickness: float) -> None:
            # UG "Coordinate Break": PARM 1, 2 decenter, 3 to 5 tilts about x, y, z, 6 order flag
            surf(number, ["TYPE COORDBRK", "PARM 1 0", "PARM 2 0", f"PARM 3 {tilt_x!r}",
                          "PARM 4 0", "PARM 5 0", "PARM 6 0", f"DISZ {thickness!r}"])

        half = FOLD["angle"]  # the mirror at 45 deg folds by 90 deg: twice the tilt (KA-01573)
        coordinate_break(4, half, 0.0)
        surf(5, ["TYPE STANDARD", "CURV 0", "DISZ 0", "GLAS MIRROR 0", "DIAM 0 0",
                 f"CLAP 0 {FOLD['aperture']!r}"])
        coordinate_break(6, half, -FOLD["to_image"])  # one mirror: negative (KA-01588)
        image = 7
    surf(image, ["TYPE STANDARD", "CURV 0", "DISZ 0", "DIAM 0 0",
                 f"SQAP {IMAGE['half_x']!r} {IMAGE['half_y']!r}"])
    lines.append("MNUM 1")
    return "".join(line + "\r\n" for line in lines)


COLUMNS = ["row", "surface", "element", "kind", "shape", "radius", "conic", "a4", "a6", "a8",
           "a10", "a12", "a14", "a16", "thickness", "material", "aperture", "aperture_1",
           "aperture_2", "stop", "focal_length", "object_distance"]


def expected_csv() -> str:
    """ADR 0045, point 7: the header, the object row, one row per surface in path order."""
    def row(**kw) -> list[str]:
        return [kw.get(c, "") for c in COLUMNS]

    def opt(v: float) -> str:  # 0 is written as an empty field
        return "" if v == 0.0 else repr(v)

    a = FRONT["a"] + [0.0] * (7 - len(FRONT["a"]))
    asph = {f"a{4 + 2 * i}": opt(v) for i, v in enumerate(a)}
    rows = [
        COLUMNS,
        row(row="object", thickness="Infinity"),
        row(row="surface", surface="STO", element="stop", kind="stop", shape="plane",
            thickness=repr(STOP["t"]), aperture="circular", aperture_1=repr(STOP["radius"]),
            aperture_2=opt(STOP["inner"]), stop="true"),
        row(row="surface", surface="L1.S1", element="L1", kind="lens", shape="even_asphere",
            radius=repr(FRONT["radius"]), conic=opt(FRONT["conic"]), thickness=repr(FRONT["t"]),
            material=f"{GLASS[0]}:{GLASS[1]}", aperture="circular",
            aperture_1=repr(FRONT["aperture"]), **asph),
        row(row="surface", surface="L1.S2", element="L1", kind="lens", shape="conic",
            radius=repr(BACK["radius"]), conic=opt(BACK["conic"]), thickness=repr(BACK["t"]),
            aperture="circular", aperture_1=repr(BACK["aperture"])),
        row(row="surface", surface="IMG", element="image", kind="detector", shape="plane",
            aperture="rectangular", aperture_1=repr(IMAGE["half_x"]),
            aperture_2=repr(IMAGE["half_y"])),
    ]
    return "".join(",".join(r) + "\r\n" for r in rows)


def fold_system() -> dict:
    """The singlet with the fold mirror: the mirror relative to L1.S2, turned in its own frame;
    the image with rotate_first and (-135, 0, 0), so that its z lies along the folded axis
    (ADR 0028, point 4)."""
    s = system()
    s["name"] = "Export fold"
    stop, lens, image = s["root"]["children"]
    mirror = {
        "type": "mirror",
        "name": "M",
        "surfaces": [{
            "id": "M.S",
            "shape": {"base": {"type": "plane"}},
            "aperture": {"type": "circular", "radius": FOLD["aperture"]},
            "interaction": {"type": "ideal_mirror"},
        }],
        "pose": {"reference": "relative_to_preceding", "position": [0.0, 0.0, FOLD["to_mirror"]],
                 "rotation_deg": [FOLD["angle"], 0.0, 0.0]},
    }
    image["pose"] = {"reference": "relative_to_preceding", "order": "rotate_first",
                     "position": [0.0, 0.0, FOLD["to_image"]],
                     "rotation_deg": [FOLD["angle"] - 180.0, 0.0, 0.0]}
    s["root"]["children"] = [stop, lens, mirror, image]
    return s


def main() -> None:
    (HERE / "export_singlet.rtt.json").write_text(json.dumps(system(), indent=2) + "\n",
                                                  encoding="utf-8", newline="\n")
    (HERE / "export_singlet.zmx").write_bytes(expected_zmx().encode("utf-8"))
    (HERE / "export_singlet.csv").write_bytes(expected_csv().encode("utf-8"))
    (HERE / "export_fold.rtt.json").write_text(json.dumps(fold_system(), indent=2) + "\n",
                                               encoding="utf-8", newline="\n")
    (HERE / "export_fold.zmx").write_bytes(expected_zmx(fold=True).encode("utf-8"))


if __name__ == "__main__":
    main()
