"""The command line against the Python API (ADR 0052, point 8): rtt trace and rtt analyze write
byte for byte to_json of the same call with the same parameters, so the options of the command
line keep the defaults of the API. ctest sets RTT_CLI_EXE; without it the tests fail instead of
being skipped (AGENTS.md, bitwise comparisons). Content guards: the expected status and values
that are not trivially empty."""

from __future__ import annotations

import json
import os
import subprocess
from collections.abc import Callable
from typing import Any

import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR

import raytatouille as rt
from raytatouille import analysis as an
from raytatouille import results

CLI = os.environ.get("RTT_CLI_EXE")
SCHOTT = CATALOG_DIR / "schott.agf"
SINGLET = REFERENCE_DIR / "m1" / "singlet_const.rtt.json"
ACHROMAT = REFERENCE_DIR / "m2" / "achromat.rtt.json"
PARABOLOID = REFERENCE_DIR / "m2" / "paraboloid_stop.rtt.json"


def _compiled(path: Any, schott: bool = False) -> rt.CompiledSystem:
    lib = rt.MaterialLibrary()
    if schott:
        lib.add_catalog(SCHOTT)
    return rt.compile(rt.load(path), lib)


def _traced(paths: bool) -> Any:
    cs = _compiled(SINGLET)
    rays = rt.trace.make_rays(cs, rt.trace.GridPupil(n=5))
    if paths:
        return rt.trace.trace(cs, rays, record_path=True)[1]
    rt.trace.trace(cs, rays)
    return rays


Case = tuple[list[str], Callable[[], Any], Callable[[dict[str, Any]], bool]]

# Arguments after "rtt", the Python call with the same parameters (API defaults otherwise),
# a content guard on the data.
CASES: dict[str, Case] = {
    "trace": (["trace", str(SINGLET), "--rays", "grid:5"], lambda: _traced(False),
              lambda d: len(d["pos_x"]["values"]) > 10 and set(d["status"]["values"]) == {0}),
    "trace_paths": (["trace", str(SINGLET), "--rays", "grid:5", "--paths"],
                    lambda: _traced(True),
                    lambda d: d["ray_count"] > 10 and set(d["status"]["values"]) == {0}),
    "first_order": (["analyze", "first-order", str(ACHROMAT), "--catalog", str(SCHOTT)],
                    lambda: rt.paraxial.first_order(_compiled(ACHROMAT, True)),
                    lambda d: d["efl"] > 0),
    "seidel": (["analyze", "seidel", str(ACHROMAT), "--catalog", str(SCHOTT)],
               lambda: rt.paraxial.seidel(_compiled(ACHROMAT, True)),
               lambda d: d["sum"]["s1"] != 0),
    "spot": (["analyze", "spot", str(ACHROMAT), "--catalog", str(SCHOTT), "--field", "2"],
             lambda: an.spot(_compiled(ACHROMAT, True), 0, 2),
             lambda d: d["rays_arrived"] > 0 and d["wavelength"] is None),
    "opd": (["analyze", "opd", str(SINGLET), "--field", "2"],
            lambda: an.opd_map(_compiled(SINGLET), 0, 2), lambda d: d["rms"] > 0),
    "zernike": (["analyze", "zernike", str(SINGLET), "--field", "2", "--grid", "24"],
                lambda: an.zernike_fit(an.pupil_function(_compiled(SINGLET), 0, 2, grid=24)),
                lambda d: any(c != 0.0 for c in d["coefficients"]["values"])),
    "psf": (["analyze", "psf", str(PARABOLOID), "--grid", "32"],
            lambda: an.psf(_compiled(PARABOLOID), 0, 0, grid=32),
            lambda d: 0 < d["strehl"] <= 1),
    "mtf": (["analyze", "mtf", str(PARABOLOID), "--grid", "32"],
            lambda: an.mtf(_compiled(PARABOLOID), 0, 0, grid=32),
            lambda d: len(d["tangential"]["modulation"]["values"]) > 1),
}


def test_the_command_line_is_given() -> None:
    assert CLI, "RTT_CLI_EXE not set: run the Python tests through ctest (ADR 0052, point 8)"


@pytest.mark.skipif(not CLI, reason="RTT_CLI_EXE not set (test_the_command_line_is_given fails)")
@pytest.mark.filterwarnings("ignore::raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("name", list(CASES))
def test_command_line_equals_to_json(name: str) -> None:
    args, call, guard = CASES[name]
    assert CLI is not None
    done = subprocess.run([CLI, *args], capture_output=True, check=False)
    assert done.returncode == 0, done.stderr.decode("utf-8", "replace")
    expected = results.to_json(call())
    assert done.stdout.decode("utf-8") == expected
    assert guard(json.loads(expected)["data"]), f"{name}: content guard"


@pytest.mark.skipif(not CLI, reason="RTT_CLI_EXE not set (test_the_command_line_is_given fails)")
def test_optimize_equals_the_python_api(tmp_path: Any) -> None:
    """rtt optimize: the system on stdout is System.to_json() of rt.optim.optimize, --result is
    to_json of the OptimResult (ADR 0052, point 3)."""
    assert CLI is not None
    file = REFERENCE_DIR / "m5" / "singlet_optim.rtt.json"
    result_file = tmp_path / "optimize.result.json"
    done = subprocess.run([CLI, "optimize", str(file), "--result", str(result_file)],
                          capture_output=True, check=False)
    assert done.returncode == 0, done.stderr.decode("utf-8", "replace")
    r = rt.optim.optimize(rt.load(file))
    assert done.stdout.decode("utf-8") == r.system.to_json()
    assert result_file.read_text(encoding="utf-8") == results.to_json(r)
    assert r.status.name.startswith("CONVERGED_") and r.iterations > 0  # content guard
