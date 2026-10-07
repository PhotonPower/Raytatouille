"""Results as data: save an analysis as JSON and read it in another process (#86, ADR 0023).

Computes a spot diagram and the first-order data of the reference singlet, writes both in the
format raytatouille-result, reads them back as plain data (dicts and NumPy arrays, no C++
objects) and checks that the round trip is bit for bit. The spot also says where its rays
went (losses by status, worst loss surface).

Run from the repository root after building the package (pip install .):

    python examples/python/results.py [output directory]
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
FILE = REPO_ROOT / "tests" / "reference" / "m1" / "singlet_const.rtt.json"


def main(out_dir: str | None = None, file: Path = FILE) -> int:
    if out_dir is None:
        with tempfile.TemporaryDirectory() as tmp:
            return run(Path(tmp), file)
    return run(Path(out_dir), file)


def run(directory: Path, file: Path) -> int:
    system = rt.compile(rt.load(file))
    spot = rt.analysis.spot(system, field=1, rays="hexapolar:4")
    first = rt.paraxial.first_order(system)

    for name, result in (("spot", spot), ("first_order", first)):
        path = directory / f"{name}.result.json"
        path.write_text(result.to_json(indent=1), encoding="utf-8")
        print(f"{type(result).__name__} written to {path} ({path.stat().st_size} bytes)")

    # Another process (a GUI) reads the data back without the C++ objects.
    data = rt.results.load_json((directory / "spot.result.json").read_text(encoding="utf-8")).data
    losses = data["losses"]
    print(f"Spot of field {data['field']}: {losses['launched']} rays, "
          f"{losses['by_status'][rt.trace.RayStatus.ALIVE.value]} arrived, "
          f"RMS radius {data['stats']['rms_centroid'] * 1e3:.3f} um")
    efl = rt.results.load_json(
        (directory / "first_order.result.json").read_text(encoding="utf-8")).data["efl"]
    print(f"EFL {efl:.6f} mm")

    if not (np.array_equal(data["x"], spot.x) and data["x"].tobytes() == spot.x.tobytes()
            and efl == first.efl):
        print("round trip changed the data", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else None))
