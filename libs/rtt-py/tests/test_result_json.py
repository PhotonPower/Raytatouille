"""The C++ writer of the result format (ADR 0052, point 1) against to_json in Python.

rtt_py_reference (result_json_cases.cpp) computes the results of CASES and writes each with
rtt::io::to_result_json as json_<case>.result.json; CASES below make the same calls in Python.
The texts must be byte for byte equal (indent=None), and each C++ text must be valid against
schema/raytatouille-result.schema.json. Every case has a content guard: the expected status and
values that are not trivially empty, so that equal wrong results on both sides do not pass. The
cases reach every branch of the writer: Q and polarization, lost rays (status values, NaN slots,
lost_at, a worst surface), with and without a chromatic pair, mono- and polychromatic spots, the
MTF with and without maps, a first order without pupils and an optimization with a variable row
per configuration.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import jsonschema  # type: ignore[import-untyped]  # no py.typed; test tool (ADR 0018)
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille import results

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

SCHEMA = REFERENCE_DIR.parent.parent / "schema" / "raytatouille-result.schema.json"


def _systems() -> dict[str, Any]:
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG_DIR / "schott.agf")
    clipped = json.loads(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json").to_json())
    clipped["root"]["children"][1]["surfaces"][0]["aperture"] = {"type": "circular",
                                                                 "radius": 21.0}
    return {
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
        "achromat": rt.compile(rt.load(REFERENCE_DIR / "m2" / "achromat.rtt.json"), schott),
        "paraboloid": rt.compile(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json")),
        "paraboloid21": rt.compile(rt.System.from_json(json.dumps(clipped))),
        "no_stop": rt.compile(rt.load(REFERENCE_DIR / "m1" / "paraboloid_mirror.rtt.json")),
        "optim_singlet": rt.load(REFERENCE_DIR / "m5" / "singlet_optim.rtt.json"),
        "optim_configurations": rt.load(
            REFERENCE_DIR / "m5" / "singlet_optim_configurations.rtt.json"),
    }


def _traced(s: dict[str, Any], q: bool = False) -> tuple[Any, Any, Any]:
    rays = rt.trace.make_rays(s["singlet"], rt.trace.HexapolarPupil(rings=2))
    if q:
        rays.enable_q()
    stats, paths = rt.trace.trace(s["singlet"], rays, record_path=True, polarization=q)
    return rays, paths, stats


def _lost(s: dict[str, Any]) -> tuple[Any, Any, Any]:
    rays = rt.trace.make_rays(s["paraboloid21"], rt.trace.GridPupil(n=21), fields=[0])
    stats, paths = rt.trace.trace(s["paraboloid21"], rays, record_path=True)
    return rays, paths, stats


Case = Callable[[dict[str, Any]], Any]
Guard = Callable[[dict[str, Any]], bool]

# Same calls as run_result_json_cases() in result_json_cases.cpp, each with its type and a
# content guard on the data of the JSON text.
CASES: dict[str, tuple[str, Case, Guard]] = {
    "rays": ("RayBatch", lambda s: _traced(s)[0],
             lambda d: len(d["pos_x"]["values"]) == 57 and d["q"] is None
             and set(d["status"]["values"]) == {0}),
    "paths": ("RayPaths", lambda s: _traced(s)[1],
              lambda d: d["ray_count"] == 57 and d["prt"] is None and d["slots"] > 2),
    "stats": ("TraceStats", lambda s: _traced(s)[2], lambda d: d["rays"][0] == 57),
    "first_order": ("FirstOrder", lambda s: rt.paraxial.first_order(s["achromat"], path=0),
                    lambda d: d["efl"] > 0 and d["entrance_pupil"] is not None),
    "seidel": ("Seidel", lambda s: rt.paraxial.seidel(s["achromat"], path=0, pair=(0, 2)),
               lambda d: d["chromatic"] == {"first": 0, "second": 2} and d["sum"]["c_l"] != 0),
    "spot": ("SpotDiagram", lambda s: an.spot(s["achromat"], 0, 2, None, rays="random:300:7"),
             lambda d: d["wavelength"] is None and d["rays_arrived"] > 0
             and d["stats"]["rms_centroid"] > 0),
    "opd": ("OpdMap", lambda s: an.opd_map(s["singlet"], 0, 2, grid=17),
            lambda d: d["arrived"] > 0 and d["rms"] > 0),
    "zernike": ("ZernikeFit", lambda s: an.zernike_fit(an.pupil_function(s["singlet"], 0, 2,
                                                                          grid=24)),
                lambda d: any(c != 0.0 for c in d["coefficients"]["values"])),
    "psf": ("Psf", lambda s: an.psf(s["paraboloid"], 0, 0, grid=32),
            lambda d: 0 < d["strehl"] <= 1 and d["warnings"]),
    "mtf": ("Mtf", lambda s: an.mtf(s["paraboloid"], 0, 0, grid=32, map=True, complex_otf=True),
            lambda d: d["map"]["shape"] == [63, 63] and d["otf_map"]["dtype"] == "complex128"
            and d["warnings"]),
    "optimize": ("OptimResult", lambda s: rt.optim.optimize(s["optim_singlet"]),
                 lambda d: d["status"].startswith("CONVERGED_") and d["variables"]),
    # Branches the cases above do not reach.
    "rays_q": ("RayBatch", lambda s: _traced(s, q=True)[0],
               lambda d: d["q"] is not None and d["q"]["shape"] == [57, 3, 3]),
    "paths_q": ("RayPaths", lambda s: _traced(s, q=True)[1],
                lambda d: d["polarization"] and d["prt"] is not None and d["q"] is not None),
    "rays_lost": ("RayBatch", lambda s: _lost(s)[0],
                  lambda d: len(set(d["status"]["values"])) > 1),
    "paths_lost": ("RayPaths", lambda s: _lost(s)[1],
                   lambda d: max(d["lost_at"]["values"]) >= 0
                   and "NaN" in d["position"]["values"]),
    "stats_lost": ("TraceStats", lambda s: _lost(s)[2],
                   lambda d: d["rays"][0] > 0 and sum(d["rays"][1:]) > 0),
    "spot_lost": ("SpotDiagram",
                  lambda s: an.spot(s["paraboloid21"], 0, 0, None, rays="grid:21"),
                  lambda d: d["losses"]["worst_surface"] is not None and d["warnings"]),
    "seidel_plain": ("Seidel", lambda s: rt.paraxial.seidel(s["achromat"], path=0),
                     lambda d: d["chromatic"] is None and d["sum"]["s1"] != 0),
    "spot_mono": ("SpotDiagram", lambda s: an.spot(s["singlet"], 0, 2, 0),
                  lambda d: d["wavelength"] == 0 and d["rays_arrived"] > 0),
    "mtf_plain": ("Mtf", lambda s: an.mtf(s["paraboloid"], 0, 0, grid=32),
                  lambda d: d["map"]["shape"] == [0, 0] and d["tangential"]["otf"]["shape"] == [0]
                  and len(d["tangential"]["modulation"]["values"]) > 1),
    "first_order_no_stop": ("FirstOrder", lambda s: rt.paraxial.first_order(s["no_stop"], path=0),
                            lambda d: d["entrance_pupil"] is None and d["efl"] is not None),
    "optimize_configurations": (
        "OptimResult", lambda s: rt.optim.optimize(s["optim_configurations"]),
        lambda d: {(v["row"], v["configuration"]) for v in d["variables"]} >= {("K", 0),
                                                                               ("K", 1)}),
}


@pytest.fixture(scope="module")
def systems() -> dict[str, Any]:
    return _systems()


@pytest.fixture(scope="module")
def validator() -> jsonschema.protocols.Validator:
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    cls = jsonschema.validators.validator_for(schema)
    cls.check_schema(schema)
    return cls(schema)


# The paraboloid warns (pupil.mapping_nonlinear, psf.aliasing) and the clipped paraboloid loses
# rays (rays.lost): wanted, the warnings fill the field warnings of the results with diagnostics.
@pytest.mark.filterwarnings("ignore::raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("name", list(CASES))
def test_cpp_writer_equals_to_json(name: str, cpp_dir: Path, systems: dict[str, Any],
                                   validator: jsonschema.protocols.Validator) -> None:
    type_name, case, guard = CASES[name]
    result = case(systems)
    assert type(result).__name__ == type_name
    expected = results.to_json(result)
    actual = (cpp_dir / f"json_{name}.result.json").read_bytes().decode("utf-8")
    if actual != expected:  # the first difference, for the report
        i = next((k for k, (a, b) in enumerate(zip(actual, expected)) if a != b),
                 min(len(actual), len(expected)))
        pytest.fail(f"{name}: differs at {i}: C++ {actual[max(i - 60, 0):i + 60]!r} "
                    f"Python {expected[max(i - 60, 0):i + 60]!r}")
    doc = json.loads(actual)
    validator.validate(doc)  # the schema of the result format
    assert doc["type"] == type_name
    assert doc["schema_version"] == results.SCHEMA_VERSION
    assert guard(doc["data"]), f"{name}: content guard"


def test_every_case_is_written(cpp_dir: Path) -> None:
    assert sorted(p.name for p in cpp_dir.glob("json_*.result.json")) == sorted(
        f"json_{name}.result.json" for name in CASES)
    assert {t for t, _, _ in CASES.values()} == {
        "RayBatch", "RayPaths", "TraceStats", "FirstOrder", "Seidel", "SpotDiagram", "OpdMap",
        "ZernikeFit", "Mtf", "Psf", "OptimResult"}
