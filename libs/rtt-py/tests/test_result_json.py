"""The C++ writer of the result format (ADR 0052, point 1) against to_json in Python.

rtt_py_reference (result_json_cases.cpp) computes one result of each of the 11 types of the
command line and writes it with rtt::io::to_result_json as json_<case>.result.json; CASES below
make the same calls in Python. The texts must be byte for byte equal (indent=None). Content
guards: every case has the expected type and values that are not trivially empty.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille import results

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)


def _systems() -> dict[str, Any]:
    schott = rt.MaterialLibrary()
    schott.add_catalog(CATALOG_DIR / "schott.agf")
    return {
        "singlet": rt.compile(rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")),
        "achromat": rt.compile(rt.load(REFERENCE_DIR / "m2" / "achromat.rtt.json"), schott),
        "paraboloid": rt.compile(rt.load(REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json")),
        "optim_singlet": rt.load(REFERENCE_DIR / "m5" / "singlet_optim.rtt.json"),
    }


def _traced(s: dict[str, Any]) -> tuple[Any, Any, Any]:
    rays = rt.trace.make_rays(s["singlet"], rt.trace.HexapolarPupil(rings=2))
    stats, paths = rt.trace.trace(s["singlet"], rays, record_path=True)
    return rays, paths, stats


Case = Callable[[dict[str, Any]], Any]

# Same calls as run_result_json_cases() in result_json_cases.cpp.
CASES: dict[str, tuple[str, Case]] = {
    "rays": ("RayBatch", lambda s: _traced(s)[0]),
    "paths": ("RayPaths", lambda s: _traced(s)[1]),
    "stats": ("TraceStats", lambda s: _traced(s)[2]),
    "first_order": ("FirstOrder", lambda s: rt.paraxial.first_order(s["achromat"], path=0)),
    "seidel": ("Seidel", lambda s: rt.paraxial.seidel(s["achromat"], path=0, pair=(0, 2))),
    "spot": ("SpotDiagram", lambda s: an.spot(s["achromat"], 0, 2, None, rays="random:300:7")),
    "opd": ("OpdMap", lambda s: an.opd_map(s["singlet"], 0, 2, grid=17)),
    "zernike": ("ZernikeFit", lambda s: an.zernike_fit(an.pupil_function(s["singlet"], 0, 2,
                                                                          grid=24))),
    "psf": ("Psf", lambda s: an.psf(s["paraboloid"], 0, 0, grid=32)),
    "mtf": ("Mtf", lambda s: an.mtf(s["paraboloid"], 0, 0, grid=32, map=True, complex_otf=True)),
    "optimize": ("OptimResult", lambda s: rt.optim.optimize(s["optim_singlet"])),
}


@pytest.fixture(scope="module")
def systems() -> dict[str, Any]:
    return _systems()


# The paraboloid and the singlet warn (pupil.mapping_nonlinear): wanted, the warnings fill the
# field warnings of the results with diagnostics.
@pytest.mark.filterwarnings("ignore::raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("name", list(CASES))
def test_cpp_writer_equals_to_json(name: str, cpp_dir: Path, systems: dict[str, Any]) -> None:
    type_name, case = CASES[name]
    result = case(systems)
    assert type(result).__name__ == type_name
    expected = results.to_json(result)
    actual = (cpp_dir / f"json_{name}.result.json").read_bytes().decode("utf-8")
    if actual != expected:  # the first difference, for the report
        i = next((k for k, (a, b) in enumerate(zip(actual, expected)) if a != b),
                 min(len(actual), len(expected)))
        pytest.fail(f"{name}: differs at {i}: C++ {actual[i - 60:i + 60]!r} "
                    f"Python {expected[i - 60:i + 60]!r}")
    # Content guard: the type and data with something in them.
    doc = json.loads(actual)
    assert doc["type"] == type_name
    assert doc["schema_version"] == results.SCHEMA_VERSION
    assert doc["data"] and len(actual) > 100


def test_every_case_is_written(cpp_dir: Path) -> None:
    assert sorted(p.name for p in cpp_dir.glob("json_*.result.json")) == sorted(
        f"json_{name}.result.json" for name in CASES)
    assert {t for t, _ in CASES.values()} == {
        "RayBatch", "RayPaths", "TraceStats", "FirstOrder", "Seidel", "SpotDiagram", "OpdMap",
        "ZernikeFit", "Mtf", "Psf", "OptimResult"}
