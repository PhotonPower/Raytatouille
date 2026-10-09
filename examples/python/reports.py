"""Reports of the Cooke triplet as data and as CSV (#177).

Three reports from raytatouille.analysis:

- dimension_report: centre thickness, edge thickness and diameter of every lens;
- system_report: system settings and the paraxial prescription (EFL, F/#, total track);
- raytrace_report: chief and marginal rays of the outer field, surface by surface.

Each report has to_dict()/to_json() (result format, ADR 0023) and to_csv(). The CSV writes
floats in the shortest form that reads back to the same double; main() checks that for the
centre thicknesses and returns 1 if a value does not come back bit for bit.

Run from the repository root after building the package (pip install .):

    python examples/python/reports.py [output directory for the CSV files]
"""

from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE = REPO_ROOT / "tests" / "reference"
CATALOG = REPO_ROOT / "tests" / "catalogs" / "m2" / "schott.agf"


def main(out_dir: str | None = None) -> int:
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG)
    cs = rt.compile(rt.load(REFERENCE / "m2" / "cooke_triplet.rtt.json"), materials=schott)
    names = rt.layout.elements(cs)

    dims = rt.analysis.dimension_report(cs)
    s = dims.segments
    print("Lens  centre thickness  edge thickness  diameter   (mm)")
    for k in range(len(s)):
        print(f"{names[s.element[k]].name:<5} {s.centre_thickness[k]:16.4f} "
              f"{s.edge_thickness[k]:15.4f} {s.diameter[k]:9.3f}")

    system = rt.analysis.system_report(cs, "main")
    p = system.prescription
    assert p is not None and p.first_order.efl is not None
    print(f"\nSystem: {len(system.wavelengths_um)} wavelengths, {system.field_count} fields, "
          f"{system.event_count} events, stop {cs.surface_ids[system.stop or 0]}")
    print(f"EFL {p.first_order.efl:.4f} mm, paraxial working F/# "
          f"{p.paraxial_working_f_number:.3f}, total track {p.total_track:.3f} mm")

    # Lower marginal, chief and upper marginal ray of the outer field.
    start = rt.trace.make_rays(cs, rt.trace.FanYPupil(3), path="main", fields=[2])
    trace = rt.analysis.raytrace_report(cs, "main", start=start)
    rows = trace.rows
    print(f"\nRaytrace report: {trace.rays} rays x {trace.slots} slots, {len(rows)} rows")
    print("ray  surface   y (mm)      local y (mm)")
    for i in range(len(rows)):
        if rows.slot[i] == 0:
            continue
        surface = cs.surface_ids[rows.surface[i]]
        print(f"{rows.ray[i]:3}  {surface:<8} {rows.y[i]:10.5f}  {rows.local_y[i]:10.5f}")

    reports = {"dimensions": dims, "system": system, "raytrace": trace}
    if out_dir is not None:
        for name, report in reports.items():
            path = Path(out_dir) / f"cooke_{name}.csv"
            path.write_text(report.to_csv(), encoding="utf-8", newline="")
            print(f"{name} written to {path}")

    # The CSV keeps every bit of the doubles.
    lines = list(csv.DictReader(io.StringIO(dims.to_csv())))
    return 0 if all(float(line["centre_thickness"]) == s.centre_thickness[k]
                    for k, line in enumerate(lines)) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else None))
