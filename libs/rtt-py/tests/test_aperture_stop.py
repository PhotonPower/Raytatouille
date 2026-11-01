"""aperture.stop from Python (#53, ADR 0039, point 10): read, compile, edit."""

from __future__ import annotations

import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR

import raytatouille as rt
from raytatouille import analysis as an

FILE = REFERENCE_DIR / "m2" / "cooke_triplet_surface_stop.rtt.json"


def schott() -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "m2" / "schott.agf")
    return lib


def test_the_stop_surface_is_read() -> None:
    s = rt.load(FILE)
    assert s.aperture.stop is not None and s.aperture.stop.surface == "L2.S1"
    # The Cooke triplet with its stop element has none.
    assert rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json").aperture.stop is None


def test_the_system_report_names_the_stop_surface() -> None:
    cs = rt.compile(rt.load(FILE), schott())
    report = an.system_report(cs, 0)
    assert report.stop is not None
    assert cs.surface_ids[report.stop] == "L2.S1"


def test_the_stop_is_removed_and_added_by_a_patch() -> None:
    s = rt.load(FILE)
    # Without aperture.stop and without a stop element the triplet has no stop: valid, compile
    # warns nothing (stop.not_on_path needs a stop), but the pupils need one.
    none = rt.apply_patch(s, [{"op": "remove", "path": "/aperture/stop"}])
    assert none.aperture.stop is None
    assert rt.compile(none, schott()).diagnostics == []
    back = rt.apply_patch(none, [
        {"op": "add", "path": "/aperture/stop", "value": {"surface": "L2.S1"}}])
    assert back.aperture.stop is not None and back.aperture.stop.surface == "L2.S1"
    assert back == s


@pytest.mark.parametrize(("surface", "code"), [
    ("X", "stop.surface_unknown"),
    ("IMG", "stop.surface_kind"),
])
def test_a_wrong_surface_is_an_error(surface: str, code: str) -> None:
    s = rt.load(FILE)
    bad = rt.apply_patch(s, [
        {"op": "replace", "path": "/aperture/stop", "value": {"surface": surface}}],
        check="structure_only")
    with pytest.raises(rt.CompileError) as error:
        rt.compile(bad, schott())
    assert code in error.value.codes
