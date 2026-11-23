"""The freeform mirrors of S3 and the aperture systems of F1 from Python are bitwise equal to C++
(#324, plan v1, section 3).

rtt_py_reference (m8_cases.cpp) loads zernike_mirror, xy_mirror, hexagon_mirror and
spider_telescope under tests/reference/m8, replaces the fields by FIELDS, aims a 21 x 21 grid
pupil, traces it and writes the start rays, the traced columns, the spot of field 1, first_order
and the Seidel sums as named .npy arrays. The grid-sag mirror has no paraxial model (ADR 0051,
point 5) and the Python API no aiming reference: its start rays come from its nominal design (the
conic without the term) and are traced through the grid mirror. Both sides run with 1 and with
4 threads.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import REFERENCE_DIR, REFERENCE_EXE
from test_bitwise import first_order_values
from test_bitwise_analysis import flatten_spot

import raytatouille as rt
from raytatouille import analysis as an


Arrays = dict[str, npt.NDArray[Any]]
#: Field angles in degree, as in m8_cases.cpp.
FIELDS = ((0.0, 0.0), (0.0, 0.5), (0.3, 0.0))
COLUMNS = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z", "opl", "field",
           "last_surface", "status")
START = ("pos_x", "pos_y", "pos_z", "dir_x", "dir_y", "dir_z")
AIMED = ("zernike_mirror", "xy_mirror", "hexagon_mirror", "spider_telescope")
GRID = "grid_sag_mirror"
SAMPLING = rt.trace.GridPupil(n=21)


def data(stem: str) -> dict[str, Any]:
    d: dict[str, Any] = json.loads(
        (REFERENCE_DIR / "m8" / f"{stem}.rtt.json").read_text(encoding="utf-8"))
    d["fields"]["points"] = [{"x": x, "y": y} for x, y in FIELDS]
    return d


def compiled(d: dict[str, Any]) -> rt.CompiledSystem:
    return rt.compile(rt.System.from_json(json.dumps(d)))


def rays_arrays(start: rt.trace.RayBatch, cs: rt.CompiledSystem, threads: int) -> Arrays:
    out: Arrays = {f"start_{name}": np.array(getattr(start, name)) for name in START}
    stats = rt.trace.trace(cs, start, threads=threads)
    out.update({name: np.array(getattr(start, name)) for name in COLUMNS})
    out["stats"] = np.array(stats.rays, dtype=np.uint64)
    return out


def seidel_sums(s: rt.paraxial.Seidel) -> npt.NDArray[np.float64]:
    """As write_seidel() in m8_cases.cpp: the "scalars" of flatten_seidel()."""
    t = s.sum
    return np.array([t.s1, t.s2, t.s3, t.s4, t.s5, t.c_l, t.c_t, s.lagrange, s.marginal.z,
                     s.marginal.y, s.marginal.u, s.chief.z, s.chief.y, s.chief.u],
                    dtype=np.float64)


def flatten(stem: str, threads: int) -> Arrays:
    if stem == GRID:
        d = data(GRID)
        cs = compiled(d)
        d["root"]["children"][1]["surfaces"][0]["shape"]["terms"] = []
        start = rt.trace.make_rays(compiled(d), SAMPLING, threads=threads)
        return rays_arrays(start, cs, threads)
    cs = compiled(data(stem))
    out = rays_arrays(rt.trace.make_rays(cs, SAMPLING, threads=threads), cs, threads)
    spot = an.spot(cs, 0, 1, None, rays=SAMPLING, threads=threads)
    out.update({f"spot_{k}": v for k, v in flatten_spot(spot).items()})
    out["first_order"] = first_order_values(rt.paraxial.first_order(cs))
    out["seidel"] = seidel_sums(rt.paraxial.seidel(cs))
    return out


CASES = {f"m8_{stem}": stem for stem in (*AIMED, GRID)}


@pytest.mark.skipif(REFERENCE_EXE is None,
                    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)")
@pytest.mark.filterwarnings("ignore::raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_m8_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = flatten(CASES[name], threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


ALIVE = int(rt.trace.RayStatus.ALIVE)
VIGNETTED = int(rt.trace.RayStatus.VIGNETTED)


@pytest.mark.filterwarnings("ignore::raytatouille.errors.RaytatouilleWarning")
@pytest.mark.parametrize("stem", (*AIMED, GRID))
def test_the_cases_have_content(stem: str) -> None:
    """The bundles reach the focus, differ between the fields and are not trivial; the freeform
    mirrors keep every ray, the aperture systems lose some at their apertures (else both sides
    could agree on something wrong)."""
    r = flatten(stem, 1)
    n_grid = 21 * 21
    inside = sum(1 for i in range(21) for j in range(21)
                 if (i - 10) ** 2 + (j - 10) ** 2 <= 100)  # GridPupil: px^2 + py^2 <= 1
    assert r["status"].size == 3 * inside < 3 * n_grid
    status = set(r["status"].tolist())
    if stem in ("hexagon_mirror", "spider_telescope"):
        assert status == {ALIVE, VIGNETTED}
    else:
        assert status == {ALIVE}
    alive = r["status"] == ALIVE
    # Field 1 is tilted in y, field 2 in x: the traced bundles move apart on the detector.
    centre = {f: (np.mean(r["pos_x"][alive & (r["field"] == f)]),
                  np.mean(r["pos_y"][alive & (r["field"] == f)])) for f in (0, 1, 2)}
    assert abs(centre[1][1] - centre[0][1]) > 0.1 and abs(centre[1][0] - centre[0][0]) < 0.1
    assert abs(centre[2][0] - centre[0][0]) > 0.1 and abs(centre[2][1] - centre[0][1]) < 0.1
    assert np.ptp(r["pos_x"][alive]) > 1e-3 and np.ptp(r["pos_y"][alive]) > 1e-3
    assert np.all(np.isfinite(r["opl"][alive]))
    if stem == "spider_telescope":
        # Blocked rays end at the stop (disk, vanes) and at the mirror (decentred ellipse).
        ids = compiled(data(stem)).surface_ids
        lost = set(r["last_surface"][~alive].tolist())
        assert lost == {ids.index("STO"), ids.index("M")}
    if stem != GRID:
        assert np.isfinite(r["first_order"][4])  # efl
        assert r["spot_ints"][4] > 0  # rays arrived
        assert r["spot_scalars"][4] > 0.0  # rms about the centroid


def test_the_grid_moves_the_rays() -> None:
    """The grid sag changes the traced rays against its nominal design, so the bitwise case
    covers the grid term and not only the conic."""
    d = data(GRID)
    grid = flatten(GRID, 1)
    d["root"]["children"][1]["surfaces"][0]["shape"]["terms"] = []
    cs = compiled(d)
    rays = rt.trace.make_rays(cs, SAMPLING, threads=1)
    plain = rays_arrays(rays, cs, 1)
    assert np.max(np.abs(grid["pos_x"] - plain["pos_x"])) > 1e-5
