"""Tolerance classes from Python are bitwise equal to C++ (#261).

rtt_py_reference (tolerance_cases.cpp) expands the tolerance set of the M7 reference achromat
with the SCHOTT catalogue; flatten_items below builds the same arrays from the Python result.
The abbe limits (a share of nu_d of the glass) are only checked exactly here, because nu_d is
not available in Python. expand runs on one thread; the cases are compared for both runs of
the reference program (1 and 4 threads).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REFERENCE_EXE

import raytatouille as rt
from raytatouille import model, tolerance

pytestmark = pytest.mark.skipif(
    REFERENCE_EXE is None,
    reason="RTT_PY_REFERENCE_EXE not set (run through ctest in the build tree)",
)

Arrays = dict[str, npt.NDArray[Any]]


def f8(values: Any) -> npt.NDArray[np.float64]:
    return np.array(values, dtype=np.float64)


def u8(values: Any) -> npt.NDArray[np.uint64]:
    return np.array(values, dtype=np.uint64)


def text(lines: list[str]) -> npt.NDArray[np.uint8]:
    return np.frombuffer("\n".join(lines).encode("utf-8"), dtype=np.uint8)


def optional_index(i: int | None) -> int:
    return 0 if i is None else i + 1


# Flattening; names and order as write_items() in tolerance_cases.cpp.


def flatten_items(r: tolerance.ToleranceItems) -> Arrays:
    return {
        "ids": text([i.id for i in r.items]),
        "limit": f8([i.limit for i in r.items]),
        "kind": u8([i.kind.value for i in r.items]),
        "distribution": u8([i.distribution.value for i in r.items]),
        "origin": u8([n for o in r.origins for n in (o.source.value, o.index)]),
        "class_ids": text([c.id for c in r.classes.items]),
        "class_ints": u8([n for c in r.classes.items for n in (
            c.kind.value, c.level.value, c.tolerance_class.value, c.skip.value,
            optional_index(c.value_override), optional_index(c.replaced_by))]),
        "overrides": u8([n for v in r.classes.value_overrides
                         for n in (v.index, v.effect.value, v.skip.value)]),
    }


def achromat(threads: int) -> Arrays:
    del threads  # expand is serial
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    system = rt.load(REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json")
    return flatten_items(tolerance.expand(system, materials=lib))


Case = Callable[[int], Arrays]

# Same calls as run_tolerance_cases() in tolerance_cases.cpp.
CASES: dict[str, Case] = {
    "achromat_expand": achromat,
}


@pytest.mark.parametrize("threads", [1, 4], ids=lambda t: f"py{t}threads")
@pytest.mark.parametrize("name", list(CASES))
def test_tolerance_equal_cpp_bitwise(name: str, threads: int, cpp_dir: Path) -> None:
    results = CASES[name](threads)
    expected_files = sorted(cpp_dir.glob(f"{name}.*.npy"))
    assert sorted(f"{name}.{k}.npy" for k in results) == [f.name for f in expected_files]
    for array, actual in results.items():
        expected = np.load(cpp_dir / f"{name}.{array}.npy")
        assert actual.dtype == expected.dtype, array
        assert actual.shape == expected.shape, array
        assert actual.tobytes() == expected.tobytes(), array


def test_the_cases_have_content() -> None:
    """26 items (23 generated, 3 explicit); the two abbe limits differ (a crown and a flint);
    at least three different skip reasons; both single values applied."""
    a = CASES["achromat_expand"](1)
    assert len(a["limit"]) == 26
    abbe = a["limit"][a["kind"] == model.ToleranceKind.ABBE.value]
    assert len(abbe) == 2 and abbe[0] != abbe[1] and bool(np.all(abbe > 0.0))
    skips = a["class_ints"].reshape(-1, 6)[:, 3]
    assert len(set(skips.tolist()) - {tolerance.ClassSkip.NONE.value}) >= 3
    applied = a["overrides"].reshape(-1, 3)[:, 1]
    assert applied.tolist() == [tolerance.OverrideEffect.APPLIED.value] * 2
