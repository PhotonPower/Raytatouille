"""Result format raytatouille-result (ADR 0023, #86): to_dict, to_json, load_json."""

from __future__ import annotations

import json
import struct
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from conftest import REFERENCE_DIR
from optim_systems import gap_merit, singlet_merit

import raytatouille as rt


def same(a: Any, b: Any, where: str = "data") -> None:
    """Asserts a == b bit for bit: arrays with dtype and shape, floats including NaN and -0.0."""
    assert type(a) is type(b), f"{where}: {type(a).__name__} != {type(b).__name__}"
    if isinstance(a, np.ndarray):
        assert a.dtype == b.dtype and a.shape == b.shape, where
        assert np.ascontiguousarray(a).tobytes() == np.ascontiguousarray(b).tobytes(), where
    elif isinstance(a, float):
        assert struct.pack("<d", a) == struct.pack("<d", b), f"{where}: {a!r} != {b!r}"
    elif isinstance(a, complex):
        same(a.real, b.real, where + ".re")
        same(a.imag, b.imag, where + ".im")
    elif isinstance(a, dict):
        assert list(a) == list(b), where
        for key in a:
            same(a[key], b[key], f"{where}.{key}")
    elif isinstance(a, list):
        assert len(a) == len(b), where
        for i, (x, y) in enumerate(zip(a, b)):
            same(x, y, f"{where}[{i}]")
    else:
        assert a == b, f"{where}: {a!r} != {b!r}"


def no_constants(name: str) -> None:
    raise AssertionError(f"non-standard JSON constant {name}")


@pytest.fixture
def singlet(reference_dir: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(reference_dir / "m1" / "singlet_const.rtt.json"))


@pytest.fixture
def library(catalog_dir: Path) -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    return lib


def traced(cs: rt.CompiledSystem) -> tuple[rt.trace.RayBatch, rt.trace.TraceStats,
                                           rt.trace.RayPaths]:
    rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=2), fields=[2])
    # Ray 0 turned steeply off axis: it misses a surface, so its later slots are NaN.
    rays.dir_x[0], rays.dir_y[0], rays.dir_z[0] = 0.0, 0.8, 0.6
    stats, paths = rt.trace.trace(cs, rays, record_path=True)
    return rays, stats, paths


def preamble_warning() -> rt.LoadWarning:
    """The warning of a catalogue excerpt with the RadiantZemax header line (#71)."""
    lib = rt.MaterialLibrary()
    text = ("Reproduced here by permission of RadiantZemax (www.radiantzemax.com).\n"
            "CC excerpt\nNM A 2 517642 1.5168 64.17\nCD 1.04 0.006 0.23 0.02 1.01 103.6\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lib.add_catalog_text(text, "Z", "z.agf")
    [w] = lib.load_warnings
    return w


def singlet_ghosts() -> rt.GhostSystem:
    """The singlet with the ghost of its path "main" (#123, #133)."""
    return rt.compile_with_ghosts(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json"), "main")


def merit_at_start() -> rt.optim.MeritEvaluation:
    """The merit function of the M5 singlet at its start values (#169)."""
    merit = rt.optim.MeritFunction(singlet_merit())
    return merit.evaluate(merit.start())


PRODUCERS: dict[str, Callable[[rt.CompiledSystem, rt.MaterialLibrary], Any]] = {
    "SpotDiagram": lambda cs, lib: rt.analysis.spot(cs, field=1, rays="hexapolar:3"),
    "RayFan": lambda cs, lib: rt.analysis.ray_fan(cs, field=2, points=7),
    "OpdMap": lambda cs, lib: rt.analysis.opd_map(cs, field=1, grid=5),
    "OpdFan": lambda cs, lib: rt.analysis.opd_fan(cs, field=1, points=5),
    "LongitudinalColour": lambda cs, lib: rt.analysis.longitudinal_colour(cs, pair=(0, 2)),
    "LateralColour": lambda cs, lib: rt.analysis.lateral_colour(cs, field=2),
    "DistortionSweep": lambda cs, lib: rt.analysis.distortion(cs, samples=3),
    "DistortionPoint": lambda cs, lib: rt.analysis.distortion_at(cs, (0.0, 4.0)),
    "FieldCurvatureSweep": lambda cs, lib: rt.analysis.field_curvature(cs, samples=3),
    "FieldCurvaturePoint": lambda cs, lib: rt.analysis.field_curvature_at(cs, (0.0, 4.0)),
    "PathTransmission": lambda cs, lib: rt.analysis.path_transmission(cs, field=1,
                                                                      rays="hexapolar:2"),
    "PathOplDifference": lambda cs, lib: rt.analysis.opl_difference(cs, 0, 0, field=1,
                                                                    rays="fan_y:3"),
    "GhostRanking": lambda cs, lib: rt.analysis.ghost_ranking(singlet_ghosts(), rays="hexapolar:2"),
    # reports (#177): a y fan of three rays at field 2
    "RaytraceReport": lambda cs, lib: rt.analysis.raytrace_report(
        cs, start=rt.trace.make_rays(cs, rt.trace.FanYPupil(3), fields=[2])),
    "SystemReport": lambda cs, lib: rt.analysis.system_report(cs),
    "DimensionReport": lambda cs, lib: rt.analysis.dimension_report(cs),
    # optimization (#169): the M5 gap system with EFL and focus, the singlet's merit at start
    "OptimResult": lambda cs, lib: rt.optim.optimize(gap_merit()),
    "MeritEvaluation": lambda cs, lib: merit_at_start(),
    "FirstOrder": lambda cs, lib: rt.paraxial.first_order(cs),
    "Seidel": lambda cs, lib: rt.paraxial.seidel(cs, pair=(0, 2)),
    "Prescription": lambda cs, lib: rt.paraxial.prescription(cs),
    "TraceStats": lambda cs, lib: traced(cs)[1],
    "RayBatch": lambda cs, lib: traced(cs)[0],
    "RayPaths": lambda cs, lib: traced(cs)[2],
    "Diattenuation": lambda cs, lib: rt.polar.diattenuation(
        np.diag([1.0, 0.5, 1.0]).astype(complex), [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
    "Diattenuations": lambda cs, lib: rt.polar.diattenuation(traced(cs)[0]),
    "Retardance": lambda cs, lib: rt.polar.retardance(
        np.diag([1.0, np.exp(0.3j), 1.0]), [0.0, 0.0, 1.0]),
    "Retardances": lambda cs, lib: rt.polar.retardance(traced(cs)[0]),
    "GlassInfo": lambda cs, lib: lib.glass("SCHOTT:N-BK7"),
    "GlassMap": lambda cs, lib: lib.glass_map(["SCHOTT"]),
    "SurfaceLayout": lambda cs, lib: rt.layout.surfaces(cs)[1],
    "CompiledElement": lambda cs, lib: rt.layout.elements(cs)[1],
    "LoadWarning": lambda cs, lib: preamble_warning(),
    "Diagnostic": lambda cs, lib: rt.compile(
        rt.load(Path(__file__).parent / "data" / "bare_asphere.rtt.json")).diagnostics[0],
}


def test_every_result_type_has_a_case() -> None:
    assert set(rt.results.TYPES) == set(PRODUCERS)


@pytest.mark.filterwarnings(r"ignore:warning \[shape\.asphere_without_coefficients\]"
                            r":raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("name", sorted(PRODUCERS))
def test_round_trip_is_bit_identical(name: str, singlet: rt.CompiledSystem,
                                     library: rt.MaterialLibrary) -> None:
    obj = PRODUCERS[name](singlet, library)
    data = obj.to_dict()
    same(rt.results.to_dict(obj), data)
    text = obj.to_json()
    assert text == rt.results.to_json(obj)
    envelope = json.loads(text, parse_constant=no_constants)  # standard JSON only
    assert list(envelope) == ["format", "schema_version", "type", "data"]
    assert envelope["format"] == "raytatouille-result"
    assert envelope["schema_version"] == rt.results.SCHEMA_VERSION == "0.1.6"
    assert envelope["type"] == name
    loaded = rt.results.load_json(text)
    assert (loaded.type, loaded.schema_version) == (name, "0.1.6")
    same(loaded.data, data)


def test_encoding_of_arrays_floats_enums_and_options(singlet: rt.CompiledSystem) -> None:
    spot = rt.analysis.spot(singlet, field=1, rays="hexapolar:3")
    data = spot.to_dict()
    same(data["x"], spot.x)
    assert data["wavelength"] is None  # polychromatic
    assert set(data["chief"]) == {"x", "y"}
    raw = json.loads(spot.to_json())["data"]
    assert raw["x"] == {"dtype": "float64", "shape": [len(spot.x)],
                        "values": [float(v) for v in spot.x]}
    fan = json.loads(rt.analysis.ray_fan(singlet, field=2, points=3).to_json())["data"]
    assert fan["tangential"]["status"]["dtype"] == "uint8"
    first = rt.paraxial.first_order(singlet).to_dict()
    assert first["entrance_pupil"] is not None and set(first["entrance_pupil"]) == {"z",
                                                                                    "diameter"}


def small_lens(reference_dir: Path, tmp_path: Path,
               radius: float = 6.0) -> rt.CompiledSystem:
    """m1/singlet_const with the aperture of L1.S1 shrunk from 12.7 mm to `radius`."""
    text = (reference_dir / "m1" / "singlet_const.rtt.json").read_text(encoding="utf-8")
    old = '"base": {"type": "conic", "radius": 51.68}'
    i = text.index(old)
    j = text.index('"radius": 12.7', i)
    path = tmp_path / "small_lens.rtt.json"
    path.write_text(text[:j] + f'"radius": {radius}' + text[j + len('"radius": 12.7'):],
                    encoding="utf-8")
    return rt.compile(rt.load(path))


# vignetted on purpose: the rays.lost warning is expected
@pytest.mark.filterwarnings(
    r"ignore:warning \[rays\.lost\]:raytatouille.errors.RaytatouilleWarning")
def test_ray_losses_are_data(reference_dir: Path, tmp_path: Path) -> None:
    cs = small_lens(reference_dir, tmp_path)
    spot = rt.analysis.spot(cs, field=0)
    losses = spot.losses
    assert losses.launched == spot.rays_launched == sum(losses.by_status)
    assert losses.count(rt.trace.RayStatus.ALIVE) == spot.rays_arrived
    lost = spot.rays_launched - spot.rays_arrived
    assert lost > 0 and losses.count(rt.trace.RayStatus.VIGNETTED) == lost
    assert losses.worst_surface is not None
    assert cs.surface_ids[losses.worst_surface] == "L1.S1"
    assert losses.worst_surface_count == lost
    assert spot.to_dict()["losses"] == {"launched": losses.launched,
                                        "by_status": list(losses.by_status),
                                        "worst_surface": 1, "worst_surface_count": lost}


def test_analyses_warn_about_lost_rays(reference_dir: Path, tmp_path: Path) -> None:
    # Lens aperture 3 mm, EPD 20 mm on axis: most rays end Vignetted at L1.S1.
    cs = small_lens(reference_dir, tmp_path, 3.0)
    lens = "/root/children/1/surfaces/0"
    calls: dict[str, Callable[[float], Any]] = {
        "spot": lambda f: rt.analysis.spot(cs, field=0, lost_warning_fraction=f),
        "ray_fan": lambda f: rt.analysis.ray_fan(cs, field=0, lost_warning_fraction=f),
        "opd_map": lambda f: rt.analysis.opd_map(cs, field=0, lost_warning_fraction=f),
        "opd_fan": lambda f: rt.analysis.opd_fan(cs, field=0, lost_warning_fraction=f)}
    for name, call in calls.items():
        with pytest.warns(rt.errors.RaytatouilleWarning) as caught:
            result = call(0.5)
        [w] = caught
        assert w.filename == __file__, name  # points at the line that called the analysis
        assert isinstance(w.message, rt.errors.RaytatouilleWarning)
        assert (w.message.code, w.message.location) == ("rays.lost", lens), name
        [d] = result.warnings
        assert (d.code, d.severity, d.location) == ("rays.lost", rt.Severity.WARNING, lens)
        assert result.to_dict()["warnings"] == [d.to_dict()]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert call(0.99).warnings == []  # above the threshold: silent
    with pytest.raises(ValueError, match="lost_warning_fraction"):
        rt.analysis.spot(cs, field=0, lost_warning_fraction=1.5)


def test_nan_and_infinity_are_strings(singlet: rt.CompiledSystem) -> None:
    _, _, paths = traced(singlet)
    text = paths.to_json()
    assert '"NaN"' in text and "NaN," not in text.replace('"NaN",', "")
    loaded = rt.results.load_json(text).data["position"]
    assert np.isnan(loaded).any()
    same(loaded, paths.position)


def test_scalar_special_floats_are_tagged() -> None:
    data = {"rays": [], "a": float("nan"), "b": float("inf"), "c": -float("inf"), "d": "NaN",
            "e": -0.0}
    text = rt.results.encode_json("TraceStats", data)
    raw = json.loads(text, parse_constant=no_constants)["data"]
    assert raw == {"rays": [], "a": {"float": "NaN"}, "b": {"float": "Infinity"},
                   "c": {"float": "-Infinity"}, "d": "NaN", "e": -0.0}
    same(rt.results.load_json(text).data, data)


def test_prescription_without_stop(reference_dir: Path) -> None:
    """Without a stop the pupil values are None and the chief-ray arrays NaN (#84); both come
    back bit for bit."""
    cs = rt.compile(rt.load(reference_dir / "m1" / "paraboloid_mirror.rtt.json"))
    p = rt.paraxial.prescription(cs)
    data = p.to_dict()
    assert data["chief_start"] is None
    assert data["lagrange_invariant"] is None
    assert np.isnan(data["surfaces"]["y_bar"]).all()
    text = p.to_json()
    assert '"NaN"' in text
    same(rt.results.load_json(text).data, data)


def test_complex_values(singlet: rt.CompiledSystem) -> None:
    d = rt.polar.diattenuation(traced(singlet)[0])
    raw = json.loads(d.to_json())["data"]["axis"]
    assert raw["dtype"] == "complex128" and raw["shape"] == list(d.axis.shape)
    assert all(isinstance(v, list) and len(v) == 2 for v in raw["values"])
    same(rt.results.load_json(d.to_json()).data["axis"], d.axis)


@pytest.mark.filterwarnings(r"ignore:warning \[shape\.asphere_without_coefficients\]"
                            r":raytatouille.errors.RaytatouilleWarning")
def test_enum_names(reference_dir: Path) -> None:
    d = rt.compile(rt.load(Path(__file__).parent / "data" / "bare_asphere.rtt.json")).diagnostics[0]
    data = d.to_dict()
    assert data == {"severity": "WARNING", "code": "shape.asphere_without_coefficients",
                    "location": "/root/children/1/surfaces/0/shape/base/coefficients",
                    "message": d.message}


def envelope(**changes: Any) -> str:
    e: dict[str, Any] = {"format": "raytatouille-result", "schema_version": "0.1.0",
                         "type": "TraceStats", "data": {"rays": {"dtype": "uint64", "shape": [7],
                                                                 "values": [1, 0, 0, 0, 0, 0, 0]}}}
    e.update(changes)
    return json.dumps(e)


def test_encode_rejects_dtypes_outside_the_format() -> None:
    with pytest.raises(ValueError, match="float16"):
        rt.results.encode_json("TraceStats", {"rays": np.zeros(2, dtype=np.float16)})


def test_load_json_accepts_patch_versions() -> None:
    assert rt.results.load_json(envelope(schema_version="0.1.7")).type == "TraceStats"


@pytest.mark.parametrize("text, match", [
    (envelope(format="raytatouille-system"), "format"),
    (envelope(schema_version="0.2.0"), "schema_version"),
    (envelope(schema_version="1.0.0"), "schema_version"),
    (envelope(type="Spot"), "type"),
    (envelope(data=[]), "data"),
    (envelope(data={"rays": {"dtype": "uint64", "shape": [2], "values": [1]}}), "shape"),
    (envelope(data={"rays": {"dtype": "object", "shape": [1], "values": [1]}}), "dtype"),
    (envelope(data={"rays": {"float": "nan"}}), "float"),
    ('{"format": "raytatouille-result", "schema_version": "0.1.0", "type": "TraceStats", '
     '"data": {"x": NaN}}', "NaN"),
    ("[1, 2]", "object"),
    # second review of #104 (H1): numbers that JSON reads as infinite, too large integers
    ('{"format": "raytatouille-result", "schema_version": "0.1.0", "type": "TraceStats", '
     '"data": {"rays": [], "x": 1e400}}', "finite"),
    ('{"format": "raytatouille-result", "schema_version": "0.1.0", "type": "TraceStats", '
     '"data": {"rays": {"dtype": "float64", "shape": [1], "values": [1e400]}}}', "finite"),
    ('{"format": "raytatouille-result", "schema_version": "0.1.0", "type": "TraceStats", '
     '"data": {"rays": {"dtype": "float64", "shape": [1], "values": [1' + "0" * 400 + ']}}}',
     "finite"),
    # the same mistakes as BROKEN in tests/schema/test_result_schema.py (review of #86 B, P1)
    (envelope(data={"rays": {"dtype": "uint8", "shape": [1], "values": [1.5]}}), "uint8"),
    (envelope(data={"rays": {"dtype": "int32", "shape": [1], "values": ["5"]}}), "int32"),
    (envelope(data={"rays": {"dtype": "int64", "shape": [1], "values": [True]}}), "int64"),
    (envelope(data={"rays": {"dtype": "uint64", "shape": [1], "values": [-1]}}), "uint64"),
    (envelope(data={"rays": {"dtype": "bool", "shape": [1], "values": [2]}}), "bool"),
    (envelope(data={"rays": {"dtype": "float64", "shape": [1], "values": [True]}}), "number"),
    (envelope(data={"rays": {"dtype": "float64", "shape": [True], "values": [1.0]}}), "shape"),
    (envelope(data={"rays": {"dtype": ["float64"], "shape": [1], "values": [1.0]}}), "dtype"),
    (envelope(data={"rays": {"dtype": "float64", "values": [1.0]}}), "reserved"),
    (envelope(data={"rays": {"float": "NaN", "x": 1}}), "reserved"),
    (envelope(data={"rays": {"float": ["NaN"]}}), "float"),
    (envelope(data={}), "misses rays"),
    ('{"format": "raytatouille-result", "format": "raytatouille-result", '
     '"schema_version": "0.1.0", "type": "TraceStats", "data": {"rays": []}}', "duplicate"),
])
def test_load_json_rejects_invalid_input(text: str, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        rt.results.load_json(text)


def test_example_files_load(reference_dir: Path) -> None:
    """tests/reference/results holds one example per type; tests/schema checks them against
    schema/raytatouille-result.schema.json, here they must load."""
    files = sorted((reference_dir / "results").glob("*.result.json"))
    types = {rt.results.load_json(f.read_text(encoding="utf-8")).type for f in files}
    assert types == set(rt.results.TYPES)


@pytest.mark.filterwarnings(r"ignore:warning \[shape\.asphere_without_coefficients\]"
                            r":raytatouille.errors.RaytatouilleWarning")
def test_schema_lists_the_keys_of_every_type(reference_dir: Path, singlet: rt.CompiledSystem,
                                              library: rt.MaterialLibrary) -> None:
    """schema/raytatouille-result.schema.json requires exactly the keys that to_dict writes."""
    schema = json.loads((reference_dir.parent.parent / "schema"
                         / "raytatouille-result.schema.json").read_text(encoding="utf-8"))
    required: dict[str, list[str]] = {}
    for rule in schema["allOf"]:
        required.setdefault(rule["if"]["properties"]["type"]["const"], []).extend(
            rule["then"]["properties"]["data"]["required"])
    assert schema["properties"]["type"]["enum"] == list(rt.results.TYPES)
    for name, make in PRODUCERS.items():
        assert required[name] == list(make(singlet, library).to_dict()), name


def test_trace_stats_carry_the_evanescent_status(singlet: rt.CompiledSystem) -> None:
    """Status Evanescent (#127): the status arrays grow by one entry at the end, a compatible
    addition (ADR 0023), hence patch version 0.1.3; the count survives the round trip."""
    rays = rt.trace.make_rays(singlet, rt.trace.HexapolarPupil(rings=1))
    rays.status[1] = int(rt.trace.RayStatus.EVANESCENT)
    stats = rt.trace.trace(singlet, rays, path="main")
    loaded = rt.results.load_json(stats.to_json())
    assert loaded.schema_version == rt.results.SCHEMA_VERSION == "0.1.6"
    counts = list(loaded.data["rays"])
    assert len(counts) == 8
    assert counts[int(rt.trace.RayStatus.EVANESCENT)] == 1


WAVE_KEYS = ["wave_x", "wave_y", "wave_z", "mode_index"]


def test_older_ray_batch_reads_with_the_reading_rule(singlet: rt.CompiledSystem) -> None:
    """RayBatch gained wave_x/y/z and mode_index in 0.1.5 (#134). A file of an older version
    without them still loads; the data then have wave = dir and mode_index 0 (reading rule,
    ADR 0026, point 3), in the key order of to_dict. A 0.1.5 file must have them."""
    rays = traced(singlet)[0]
    data = rays.to_dict()
    assert list(data)[-4:] == WAVE_KEYS
    envelope = json.loads(rays.to_json())
    for key in WAVE_KEYS:
        del envelope["data"][key]
    envelope["schema_version"] = "0.1.4"
    loaded = rt.results.load_json(json.dumps(envelope))
    assert loaded.schema_version == "0.1.4"
    assert list(loaded.data) == list(data)
    for axis in "xyz":
        same(loaded.data[f"wave_{axis}"], data[f"dir_{axis}"])
    same(loaded.data["mode_index"], np.zeros(len(rays)))
    envelope["schema_version"] = "0.1.5"
    with pytest.raises(ValueError, match="misses wave_x, wave_y, wave_z, mode_index"):
        rt.results.load_json(json.dumps(envelope))


def test_ray_batch_carries_wave_and_mode_index(singlet: rt.CompiledSystem) -> None:
    """The new columns go through the round trip bit for bit, also with values that only a
    crystal gives (wave != dir, mode_index > 0)."""
    rays = traced(singlet)[0]
    rays.wave_x[1], rays.wave_y[1], rays.wave_z[1] = 0.6, 0.0, 0.8
    rays.mode_index[1] = 1.5653568260606650
    loaded = rt.results.load_json(rays.to_json())
    same(loaded.data, rays.to_dict())
    assert loaded.data["mode_index"][1] == 1.5653568260606650
