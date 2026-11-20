"""Zernike surface errors as a tolerance kind in Python (#316, ADR 0051, point 7, PR 5): the
model class carries ``term``, the kind ZERNIKE, the sensitivity of a defocus error and the
result format."""

from __future__ import annotations

import json

import raytatouille as rt
from conftest import REFERENCE_DIR
from raytatouille import model, results

MIRROR = REFERENCE_DIR / "m8" / "zernike_mirror.rtt.json"


def _with_zernike_item(term: int) -> rt.System:
    data = json.loads(MIRROR.read_text(encoding="utf-8"))
    data["tolerances"] = {
        "items": [{"id": "z", "surface": "M", "kind": "zernike", "term": term, "limit": 1e-4}],
        "criteria": [{"id": "opd", "operand": {"type": "opd_rms", "path": "main"}, "max": 5.0}]}
    return rt.System.from_json(json.dumps(data))


def test_the_item_has_kind_and_term() -> None:
    (item,) = _with_zernike_item(4).tolerances.items
    assert isinstance(item, model.ToleranceItem)
    assert item.kind == model.ToleranceKind.ZERNIKE and item.term == 4
    assert item.surface == "M" and item.segment is None
    assert not rt.validate(_with_zernike_item(4))
    codes = {(d.code, d.location) for d in rt.validate(_with_zernike_item(1))}
    assert ("tolerances.target_invalid", "/tolerances/items/0/term") in codes


def test_sensitivity_of_a_defocus_error() -> None:
    sens = rt.tolerance.sensitivity(_with_zernike_item(4), confirm=False)
    (criterion,) = sens.criteria
    (contribution,) = criterion.items
    assert contribution.plus is not None and contribution.minus is not None
    assert contribution.plus != contribution.minus  # content guard: the error acts
    assert contribution.plus != criterion.nominal


def test_result_format_carries_the_term() -> None:
    items = rt.tolerance.expand(_with_zernike_item(9))
    doc = json.loads(items.to_json())
    assert doc["schema_version"] == results.SCHEMA_VERSION
    assert doc["data"]["items"][0]["term"] == 9
    assert doc["data"]["items"][0]["kind"] == "ZERNIKE"  # enum names, as all kinds
