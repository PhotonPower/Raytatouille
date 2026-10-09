"""Reports as data from Python (#177): raytrace, system data and dimension reports, CSV."""

from __future__ import annotations

import csv
import io
import math
import struct
from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille.analysis import ApertureKind
from raytatouille.trace import RayStatus


@pytest.fixture
def singlet(reference_dir: Path) -> rt.CompiledSystem:
    """m1/singlet_const: stop at z = 0, plano-convex lens CONST:1.5168 (L1.S1 R = 51.68 mm,
    L1.S2 plane 4 mm behind), both apertures 12.7 mm."""
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


@pytest.fixture
def grating(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m4" / "grating_transmission.rtt.json"))


def bits(x: float) -> bytes:
    return struct.pack("<d", x)


def rows_of(text: str) -> list[list[str]]:
    return list(csv.reader(io.StringIO(text)))


def start_with_lost_ray(cs: rt.CompiledSystem) -> rt.trace.RayBatch:
    """As the C++ test: hexapolar rings 3 over all fields, ray 1 turned steeply off axis."""
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=3))
    rays.dir_x[1], rays.dir_y[1], rays.dir_z[1] = 0.0, 0.8, 0.6
    return rays


# ------------------------------------------------------------------------- raytrace report ---


def test_raytrace_report_is_bitwise_ray_paths(singlet: rt.CompiledSystem) -> None:
    start = start_with_lost_ray(singlet)
    before = {k: np.array(v) for k, v in start.to_dict().items()}
    report = rt.analysis.raytrace_report(singlet, start=start)
    _, paths = rt.trace.trace(singlet, start_with_lost_ray(singlet), record_path=True)
    rows = report.rows
    assert (report.path, report.rays, report.slots) == (0, len(start), paths.slots)
    assert len(rows) == int(np.sum(paths.count))
    assert paths.count[1] < paths.slots  # the off-axis ray stopped
    i = 0
    for r in range(paths.ray_count):
        for s in range(int(paths.count[r])):
            assert (rows.ray[i], rows.slot[i]) == (paths.ray_indices[r], s)
            for k, axis in enumerate("xyz"):
                assert bits(getattr(rows, axis)[i]) == bits(paths.position[r, s, k])
                assert bits(getattr(rows, "d" + axis)[i]) == bits(paths.direction[r, s, k])
            assert bits(rows.opl[i]) == bits(paths.opl[r, s])
            assert bits(rows.weight[i]) == bits(paths.weight[r, s])
            assert rows.status[i] == paths.status[r, s]
            if s == 0:
                assert rows.surface[i] == rt.trace.NO_SURFACE
                assert math.isnan(rows.local_x[i]) and math.isnan(rows.local_dz[i])
            else:
                assert rows.surface[i] == paths.event_surfaces[s - 1]
            i += 1
    after = start.to_dict()  # the start rays are not changed
    assert all(np.array(after[k]).tobytes() == v.tobytes() for k, v in before.items())


def test_raytrace_report_axial_ray_meets_the_vertices(singlet: rt.CompiledSystem) -> None:
    start = rt.trace.RayBatch(1)  # on the axis along +z, in front of the stop at z = 0
    start.pos_z[0] = -1.0
    report = rt.analysis.raytrace_report(singlet, 0, start=start)
    rows = report.rows
    assert len(rows) == report.slots  # the ray arrives
    layout = rt.layout.surfaces(singlet)
    for s in range(1, report.slots):
        vertex = layout[int(rows.surface[s])].translation
        assert abs(rows.x[s] - vertex[0]) <= 1e-12
        assert abs(rows.y[s] - vertex[1]) <= 1e-12
        assert abs(rows.z[s] - vertex[2]) <= 1e-12
        assert abs(rows.local_x[s]) <= 1e-12 and abs(rows.local_z[s]) <= 1e-12
        assert abs(rows.local_dz[s] - 1.0) <= 1e-12
        assert rows.status[s] == RayStatus.ALIVE


def test_raytrace_report_columns_are_read_only_copies(singlet: rt.CompiledSystem) -> None:
    report = rt.analysis.raytrace_report(singlet, start=rt.trace.RayBatch(2))
    dtypes = {"ray": np.uint64, "slot": np.uint64, "surface": np.uint32, "x": np.float64,
              "local_dz": np.float64, "opl": np.float64, "status": np.uint8}
    for name, dtype in dtypes.items():
        column = getattr(report.rows, name)
        assert column.dtype == dtype, name
        assert not column.flags.writeable, name


def test_raytrace_report_errors(singlet: rt.CompiledSystem) -> None:
    with pytest.raises(ValueError, match="no start rays"):
        rt.analysis.raytrace_report(singlet, start=rt.trace.RayBatch(0))
    with pytest.raises(ValueError):
        rt.analysis.raytrace_report(singlet, 7, start=rt.trace.RayBatch(1))
    with pytest.raises(ValueError):
        rt.analysis.raytrace_report(singlet, "no such path", start=rt.trace.RayBatch(1))


# --------------------------------------------------------------------- system data report ---


def test_system_report_settings_and_prescription(singlet: rt.CompiledSystem) -> None:
    r = rt.analysis.system_report(singlet, "main", 1)
    assert (r.path, r.wavelength) == (0, 1)
    assert list(r.wavelengths_um) == list(singlet.wavelengths_um)
    assert r.reference_wavelength == singlet.reference_wavelength
    assert (r.field_count, r.surface_count) == (3, len(singlet.surface_ids))
    assert r.stop is not None and singlet.surface_ids[r.stop] == "STO"
    assert r.warnings == []
    assert r.prescription is not None
    # Bitwise the prescription of raytatouille.paraxial (the JSON text keeps every bit).
    assert r.prescription.to_json() == rt.paraxial.prescription(singlet, 0, 1).to_json()
    # None means the reference wavelength.
    assert rt.analysis.system_report(singlet).wavelength == singlet.reference_wavelength


def test_system_report_without_paraxial_data_warns(grating: rt.CompiledSystem) -> None:
    with pytest.warns(rt.RaytatouilleWarning) as caught:
        r = rt.analysis.system_report(grating, "order +1", 0)
    [w] = caught
    assert w.filename == __file__
    assert isinstance(w.message, rt.RaytatouilleWarning)
    assert w.message.code == "report.paraxial_unavailable"
    assert r.prescription is None
    [d] = r.warnings
    assert (d.code, d.severity) == ("report.paraxial_unavailable", rt.Severity.WARNING)
    assert r.event_count == 3


def test_system_report_errors(singlet: rt.CompiledSystem) -> None:
    with pytest.raises(ValueError):
        rt.analysis.system_report(singlet, 0, 9)
    with pytest.raises(ValueError):
        rt.analysis.system_report(singlet, 3)


# ----------------------------------------------------------------------- dimension report ---


def test_dimension_report_plano_convex_singlet(singlet: rt.CompiledSystem) -> None:
    d = rt.analysis.dimension_report(singlet)
    s = d.segments
    assert len(s) == 1  # only the lens has two surfaces
    assert singlet.surface_ids[s.first_surface[0]] == "L1.S1"
    assert rt.layout.elements(singlet)[s.element[0]].name == "L1"
    assert s.coaxial.dtype == np.bool_ and bool(s.coaxial[0])
    assert abs(s.centre_thickness[0] - 4.0) <= 1e-12
    assert s.semi_diameter_first[0] == s.semi_diameter_second[0] == 12.7
    assert s.aperture_first[0] == s.aperture_second[0] == ApertureKind.CIRCULAR
    assert s.diameter[0] == 25.4
    # Edge thickness at h = 12.7: t - (R - sqrt(R^2 - h^2)) for the sphere R = 51.68 and the
    # plane (k = 0).
    r, h = 51.68, 12.7
    assert abs(s.edge_thickness[0] - (4.0 - (r - math.sqrt(r * r - h * h)))) <= 1e-12


# -------------------------------------------------------------------------------- CSV ---


def test_raytrace_csv(singlet: rt.CompiledSystem) -> None:
    report = rt.analysis.raytrace_report(singlet, start=start_with_lost_ray(singlet))
    text = report.to_csv()
    assert text == rt.reports.to_csv(report)
    assert text.endswith("\n") and "\r" not in text
    lines = rows_of(text)
    assert tuple(lines[0]) == rt.reports.RAYTRACE_COLUMNS
    assert len(lines) == len(report.rows) + 1
    rows = report.rows
    for i, line in enumerate(lines[1:]):
        field = dict(zip(lines[0], line))
        assert int(field["ray"]) == rows.ray[i] and int(field["slot"]) == rows.slot[i]
        assert field["status"] == RayStatus(int(rows.status[i])).name
        if rows.slot[i] == 0:
            assert field["surface"] == "" and field["local_x"] == ""  # NO_SURFACE, NaN
        else:
            assert int(field["surface"]) == rows.surface[i]
        for name in ("x", "dy", "local_z", "opl", "weight"):
            value = getattr(rows, name)[i]
            if math.isnan(value):
                assert field[name] == "", name
            else:
                assert bits(float(field[name])) == bits(value), name  # shortest round trip
    statuses = {line[-1] for line in lines[1:]}
    assert "ALIVE" in statuses and len(statuses) > 1  # the off-axis ray ends lost


def test_dimension_csv(singlet: rt.CompiledSystem) -> None:
    report = rt.analysis.dimension_report(singlet)
    lines = rows_of(report.to_csv())
    assert tuple(lines[0]) == rt.reports.DIMENSION_COLUMNS
    [line] = lines[1:]
    field = dict(zip(lines[0], line))
    assert field["coaxial"] == "true"
    assert field["aperture_first"] == field["aperture_second"] == "CIRCULAR"
    assert field["semi_diameter_first"] == "12.7" and field["diameter"] == "25.4"
    seg = report.segments
    assert bits(float(field["edge_thickness"])) == bits(seg.edge_thickness[0])
    assert bits(float(field["centre_thickness"])) == bits(seg.centre_thickness[0])


def test_system_csv_keys(singlet: rt.CompiledSystem) -> None:
    report = rt.analysis.system_report(singlet)
    lines = rows_of(report.to_csv())
    assert lines[0] == ["key", "value"]
    values = dict((k, v) for k, v in lines[1:])
    assert len(values) == len(lines) - 1  # keys are unique
    p = report.prescription
    assert p is not None and p.first_order.efl is not None
    assert values["path"] == "0"
    assert bits(float(values["wavelengths_um[0]"])) == bits(report.wavelengths_um[0])
    assert values["stop"] == str(report.stop)
    assert bits(float(values["prescription.first_order.efl"])) == bits(p.first_order.efl)
    assert bits(float(values["prescription.surfaces.z[3]"])) == bits(p.surfaces.z[3])
    assert "prescription.first_order.entrance_pupil.diameter" in values
    assert not any(k.startswith("warnings") for k in values)  # an empty list has no line
    # One line per scalar of to_dict, in its order.
    assert [k for k, _ in lines[1:4]] == ["path", "wavelength", "wavelengths_um[0]"]


def test_system_csv_without_prescription(grating: rt.CompiledSystem) -> None:
    with pytest.warns(rt.RaytatouilleWarning):
        report = rt.analysis.system_report(grating, "order +1", 0)
    values = dict((k, v) for k, v in rows_of(report.to_csv())[1:])
    assert values["prescription"] == ""
    assert values["warnings[0].code"] == "report.paraxial_unavailable"
    assert values["warnings[0].severity"] == "WARNING"
    assert values["warnings[0].message"] == report.warnings[0].message


def test_csv_special_values() -> None:
    field = rt.reports._field
    assert [field(math.nan), field(None), field(math.inf), field(-math.inf)] == [
        "", "", "Infinity", "-Infinity"]
    assert [field(True), field(np.bool_(False)), field(np.uint8(3)), field(-0.0)] == [
        "true", "false", "3", "-0.0"]
    assert field(0.1 + 0.2) == "0.30000000000000004"


def test_csv_rejects_other_objects(singlet: rt.CompiledSystem) -> None:
    with pytest.raises(TypeError, match="not a report"):
        rt.reports.to_csv(rt.paraxial.prescription(singlet))  # type: ignore[arg-type]
