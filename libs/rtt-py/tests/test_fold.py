"""Paths folded by plane mirrors in Python (#357, ADR 0054, PR B2): FirstOrder.folds and
image_space (UnfoldedFrame), their place in the result format 0.1.24, and the analyses behind a
fold against the unfolded singlet (the C++ tests check the numbers to 1e-12, test_fold_*.cpp)."""
from __future__ import annotations

import json
from pathlib import Path

import jsonschema  # type: ignore[import-untyped]  # no py.typed; test tool (ADR 0018)
import numpy as np
import pytest

from conftest import REFERENCE_DIR, REPO_ROOT

import raytatouille as rt
from raytatouille import results

FOLDED = REFERENCE_DIR / "m6" / "singlet_fold45.rtt.json"
FLAT = REFERENCE_DIR / "m1" / "singlet_const.rtt.json"
RESULT_SCHEMA = json.loads(
    (REPO_ROOT / "schema" / "raytatouille-result.schema.json").read_text(encoding="utf-8"))


def _compiled(path: Path) -> rt.CompiledSystem:
    return rt.compile(rt.load(path), rt.MaterialLibrary())


def test_first_order_has_folds_and_the_image_space_frame() -> None:
    fo = rt.paraxial.first_order(_compiled(FOLDED))
    flat = rt.paraxial.first_order(_compiled(FLAT))
    assert fo.folds == 1
    assert flat.folds == 0
    np.testing.assert_array_equal(flat.image_space.linear, np.eye(3))
    np.testing.assert_array_equal(flat.image_space.offset, np.zeros(3))
    # one fold: a reflection (parity change, determinant -1)
    assert np.linalg.det(fo.image_space.linear) == pytest.approx(-1.0, abs=1e-14)
    assert fo.efl == pytest.approx(flat.efl, rel=1e-12)
    assert fo.image_z == pytest.approx(flat.image_z, rel=1e-12)
    # the unfolded image point lies on the reflected axis: image_z - 49 along +y from the mirror
    assert fo.image_z is not None
    image = fo.image_space.to_global((0.0, 0.0, fo.image_z))
    np.testing.assert_allclose(image, [0.0, fo.image_z - 49.0, 49.0], atol=1e-12)
    np.testing.assert_allclose(fo.image_space.to_unfolded(tuple(image)), [0.0, 0.0, fo.image_z],
                               atol=1e-12)
    assert fo.image_space.linear.flags.writeable is False


def test_the_result_format_carries_the_frame() -> None:
    fo = rt.paraxial.first_order(_compiled(FOLDED))
    text = results.to_json(fo)
    jsonschema.Draft202012Validator(RESULT_SCHEMA).validate(json.loads(text))
    data = results.load_json(text).data
    assert data["folds"] == 1
    np.testing.assert_array_equal(data["image_space"]["linear"], fo.image_space.linear)
    np.testing.assert_array_equal(data["image_space"]["offset"], fo.image_space.offset)
    # a file of 0.1.0 (before 0.1.24) reads with no fold and the global frame
    old = results.load_json((REFERENCE_DIR / "results" / "FirstOrder.result.json").read_text())
    assert old.data["folds"] == 0
    np.testing.assert_array_equal(old.data["image_space"]["linear"], np.eye(3))
    np.testing.assert_array_equal(old.data["image_space"]["offset"], np.zeros(3))


def test_analyses_behind_a_fold_equal_the_unfolded_singlet() -> None:
    folded, flat = _compiled(FOLDED), _compiled(FLAT)
    a = rt.analysis.opd_map(folded, field=2)
    b = rt.analysis.opd_map(flat, field=2)
    assert a.rms == pytest.approx(b.rms, abs=1e-9)
    assert b.rms > 0.01  # content guard: the singlet has aberrations
    pa = rt.paraxial.prescription(folded)
    pb = rt.paraxial.prescription(flat)
    assert len(pa.surfaces) == len(pb.surfaces) + 1  # the fold mirror
    assert pa.total_track == pytest.approx(pb.total_track, rel=1e-12)
