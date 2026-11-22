"""Design helpers and lens catalogues (ADR 0053, point 6).

singlet, doublet and triplet reproduce the reference systems byte for byte (to_json); the
catalogue gives the same systems through from_catalog; wrong arguments are ValueError, an
invalid built system is DesignError, a malformed catalogue LensCatalogError.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR

import raytatouille as rt
from raytatouille import design

SINGLET = REFERENCE_DIR / "m1" / "singlet_const.rtt.json"
ACHROMAT = REFERENCE_DIR / "m2" / "achromat.rtt.json"
COOKE = REFERENCE_DIR / "m2" / "cooke_triplet_surface_stop.rtt.json"
DEMO = CATALOG_DIR / "lenses" / "demo.json"

SINGLET_NAME = "Plano-convex singlet f = 100 mm, constant index 1.5168"
ACHROMAT_NAME = "Cemented achromat N-BK7/F2, f about 100 mm"
COOKE_WAVELENGTHS = (0.4862661791313749, 0.587721870012617, 0.6564506940778654)

# The arguments of ADR 0053, points 6.1 to 6.3, without the lens data.
SINGLET_ARGS: dict[str, Any] = {
    "epd": 20.0, "stop_distance": 5.0, "semi_diameter": 12.7, "fields_y": (0.0, 3.5, 5.0),
    "image_distance": 97.363, "name": SINGLET_NAME}
ACHROMAT_ARGS: dict[str, Any] = {
    "epd": 20.0, "stop_distance": 5.0, "semi_diameter": 12.7, "fields_y": (0.0, 3.5, 5.0),
    "image_distance": 96.614, "name": ACHROMAT_NAME}


def _cooke_args() -> dict[str, Any]:
    return {"epd": 12.5, "stop_surface": "L2.S1", "stop_radius": 10.0,
            "fields_y": (0.0, 14.0, 20.0), "wavelengths": COOKE_WAVELENGTHS,
            "image_distance": 41.10346, "name": rt.load(COOKE).name}


def _singlet() -> rt.System:
    return design.singlet((51.68, None), 4.0, "CONST:1.5168", **SINGLET_ARGS)


def _doublet() -> rt.System:
    return design.doublet((61.5, -35.2, -136.9), (6.0, 3.0), ("SCHOTT:N-BK7", "SCHOTT:F2"),
                          **ACHROMAT_ARGS)


def _triplet() -> rt.System:
    return design.triplet((23.713, 7331.288, -24.456, 21.896, 86.759, -20.4942),
                          (4.831, 0.975, 3.127), (5.86, 4.822),
                          ("SCHOTT:N-LAK9", "SCHOTT:N-SF5", "SCHOTT:N-LAK9"), **_cooke_args())


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_singlet_equals_the_reference_byte_for_byte() -> None:
    assert _singlet().to_json() == _text(SINGLET)


def test_doublet_equals_the_achromat_byte_for_byte() -> None:
    # The acceptance of the roadmap (P4, F5).
    assert _doublet().to_json() == _text(ACHROMAT)


def test_triplet_equals_the_cooke_triplet_byte_for_byte() -> None:
    # The detector relative to the last node: 16.488 + (3.127 + 41.10346) == 60.71846, while
    # (16.488 + 3.127) + 41.10346 == 60.71845999999999 (ADR 0053, point 1).
    assert _triplet().to_json() == _text(COOKE)


@pytest.mark.parametrize(
    ("build", "catalog"),
    [(_singlet, None), (_doublet, "schott.agf"), (_triplet, "m2/schott.agf")],
    ids=["singlet", "doublet", "triplet"])
def test_the_built_systems_trace(build: Any, catalog: str | None) -> None:
    # Content guard (ADR 0053, point 6.4): valid, compiles, an on-axis bundle reaches IMG.
    system = build()
    assert rt.validate(system) == []
    lib = rt.MaterialLibrary()
    if catalog is not None:
        lib.add_catalog(CATALOG_DIR / catalog)
    compiled = rt.compile(system, lib)
    rays = rt.trace.make_rays(compiled, rt.trace.HexapolarPupil(rings=3))
    rt.trace.trace(compiled, rays)
    status = np.asarray(rays.status)
    assert status.size > 10
    assert (status == int(rt.trace.RayStatus.ALIVE)).all()


def test_the_default_wavelengths_and_names() -> None:
    s = design.singlet((50.0, -50.0), 4.0, "CONST:1.5", epd=10.0, stop_distance=2.0,
                       image_distance=45.0)
    assert s.name == "Singlet"
    assert [w.um for w in s.wavelengths] == list(design.DEFAULT_WAVELENGTHS)
    assert [w.reference for w in s.wavelengths] == [False, True, False]
    assert s.schema_version == rt.SCHEMA_VERSION
    sto = s.locate_surface("STO")
    assert sto is not None
    assert s.json_at(sto + "/aperture/radius/value") == 5.0  # epd / 2
    d = design.doublet((50.0, -40.0, -100.0), (4.0, 2.0), ("CONST:1.5", "CONST:1.6"), epd=10.0,
                       stop_distance=2.0, image_distance=45.0, reference=0)
    assert d.name == "Doublet"
    assert [w.reference for w in d.wavelengths] == [True, False, False]


def test_the_lens_catalogue_reads_the_demo() -> None:
    catalog = design.load_lens_catalog(DEMO)
    assert catalog.name == "DEMO"
    assert len(catalog) == 3
    assert [e.kind for e in catalog.lenses] == ["singlet", "doublet", "triplet"]
    pcx = catalog["PCX-25-100"]
    assert pcx.radii == [51.68, None]
    assert pcx.diameter == 25.4
    assert catalog["TRP-50"].gaps == [5.86, 4.822]
    assert catalog["TRP-50"].diameter is None
    with pytest.raises(KeyError, match="NONE"):
        catalog["NONE"]
    assert [e.name for e in catalog] == ["PCX-25-100", "ACH-25-100", "TRP-50"]
    assert "TRP-50" in catalog
    assert "NONE" not in catalog


def test_from_catalog_gives_the_reference_systems() -> None:
    catalog = design.load_lens_catalog(DEMO)
    single = {k: v for k, v in SINGLET_ARGS.items() if k != "semi_diameter"}  # from diameter
    assert design.from_catalog(catalog["PCX-25-100"], **single).to_json() == _text(SINGLET)
    assert design.from_catalog(catalog["ACH-25-100"], **ACHROMAT_ARGS).to_json() == _text(ACHROMAT)
    assert design.from_catalog(catalog["TRP-50"], **_cooke_args()).to_json() == _text(COOKE)
    # Without a name: the name of the lens.
    named = design.from_catalog(catalog["PCX-25-100"], epd=20.0, stop_distance=5.0,
                                image_distance=97.363)
    assert named.name == "PCX-25-100"
    # A given semi_diameter wins over diameter / 2 (12.7).
    small = design.from_catalog(catalog["PCX-25-100"], epd=20.0, stop_distance=5.0,
                                image_distance=97.363, semi_diameter=11.0)
    s1 = small.locate_surface("L1.S1")
    assert s1 is not None
    assert small.json_at(s1 + "/aperture/radius/value") == 11.0


def test_a_malformed_catalogue_is_a_lens_catalog_error() -> None:
    text = ('{"format": "raytatouille-lenses", "schema_version": "0.1.0", "catalog": "D", '
            '"lenses": [{"name": "L", "kind": "singlet", "radii": [1.0], "thicknesses": [1.0], '
            '"materials": ["M"]}]}')
    with pytest.raises(rt.LensCatalogError) as info:
        design.parse_lens_catalog(text, "cat.json")
    assert info.value.file == "cat.json"
    assert info.value.pointer == "/lenses/0/radii"
    assert isinstance(info.value, ValueError)
    with pytest.raises(rt.LensCatalogError, match="cannot open"):
        design.load_lens_catalog(CATALOG_DIR / "lenses" / "none.json")


def test_the_lens_schema_is_in_the_package() -> None:
    schema = rt.schema("lenses")
    assert schema["properties"]["format"]["const"] == "raytatouille-lenses"


def test_numpy_and_integer_arguments() -> None:
    # Integers and NumPy scalars are numbers; +-inf (also float32) is plano.
    s = design.singlet((51.68, np.float32("inf")), 4, "CONST:1.5168", epd=20, stop_distance=5,
                       semi_diameter=12.7, fields_y=(0, 3.5, 5), image_distance=97.363,
                       reference=np.int64(1), name=SINGLET_NAME)
    assert s.to_json() == _text(SINGLET)


GOOD: dict[str, Any] = {"epd": 10.0, "stop_distance": 2.0, "image_distance": 45.0}


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: design.singlet((50.0,), 4.0, "M:A", **GOOD), "radii needs 2"),
        (lambda: design.singlet((0.0, None), 4.0, "M:A", **GOOD), "must not be 0"),
        (lambda: design.singlet((50.0, None), 0.0, "M:A", **GOOD), "thickness must be > 0"),
        (lambda: design.singlet((50.0, None), 4.0, "", **GOOD), "material"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", epd=10.0, image_distance=45.0),
         "exactly one"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", stop_surface="L1.S1", **GOOD),
         "exactly one"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", epd=10.0, image_distance=45.0,
                                stop_distance=0.0), "stop_distance must be > 0"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", epd=10.0, image_distance=45.0,
                                stop_surface="L1.S1"), "needs stop_radius"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", epd=10.0, image_distance=45.0,
                                stop_surface="L2.S1", stop_radius=4.0), "none of"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", epd=0.0, stop_distance=2.0,
                                image_distance=45.0), "epd must be > 0"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", reference=3, **GOOD), "reference"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", fields_y=(), **GOOD), "fields_y"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", wavelengths=(), **GOOD),
         "wavelengths"),
        (lambda: design.doublet((50.0, -40.0, -100.0), (4.0,), ("M:A", "M:B"), **GOOD),
         "thicknesses needs 2"),
        (lambda: design.doublet((50.0, -40.0, -100.0), (4.0, 2.0), ("M:A",), **GOOD),
         "materials needs 2"),
        (lambda: design.triplet((1.0, 2.0, 3.0, 4.0, 5.0, 6.0), (1.0, 1.0, 1.0), (1.0, -1.0),
                                ("M:A", "M:B", "M:A"), **GOOD), "gaps"),
        (lambda: design.triplet((1.0, 2.0, 3.0, 4.0, 5.0, 6.0), (1.0, 1.0, 1.0), (1.0,),
                                ("M:A", "M:B", "M:A"), **GOOD), "gaps needs 2"),
        (lambda: design.singlet((50.0, None), True, "M:A", **GOOD), "must be a number"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", reference=True, **GOOD), "reference"),
        (lambda: design.singlet((50.0, None), 4.0, "M:A", name="", **GOOD), "name"),
    ],
)
def test_wrong_arguments_are_value_errors(call: Any, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()


def _diagnostics(change: Any) -> list[rt.Diagnostic]:
    """Real diagnostics of validate for the singlet reference changed by ``change(data)``."""
    data = json.loads(_text(SINGLET))
    change(data)
    found = rt.validate(rt.System.from_json(json.dumps(data)))
    assert found
    return found


def test_an_invalid_built_system_is_a_design_error(monkeypatch: pytest.MonkeyPatch) -> None:
    # The helpers check their arguments first, so validate finds no error for a valid call; the
    # branch is tested with the diagnostics of an invalid system (no field point: fields.empty).
    found = _diagnostics(lambda d: d["fields"].update(points=[]))
    assert all(d.severity == rt.Severity.ERROR for d in found)
    monkeypatch.setattr(design, "validate", lambda system: found)
    with pytest.raises(rt.DesignError) as info:
        _singlet()
    assert info.value.codes == ["fields.empty"]
    assert isinstance(info.value, ValueError)


def test_warnings_of_validate_are_raytatouille_warnings(monkeypatch: pytest.MonkeyPatch) -> None:
    # An object-space NA of 1.2 in air: aperture.na_not_physical, a warning.
    found = _diagnostics(lambda d: d.update(aperture={"type": "object_na", "value": 1.2}))
    assert [d.severity for d in found] == [rt.Severity.WARNING]
    monkeypatch.setattr(design, "validate", lambda system: found)
    with pytest.warns(rt.RaytatouilleWarning) as record:
        system = _singlet()
    assert [w.message.code for w in record if isinstance(w.message, rt.RaytatouilleWarning)] == [
        "aperture.na_not_physical"]
    assert system.to_json() == _text(SINGLET)  # a warning does not stop the helper
    # The warning names the caller, for the helpers and for from_catalog.
    assert all(w.filename == __file__ for w in record)
    catalog = design.load_lens_catalog(DEMO)
    with pytest.warns(rt.RaytatouilleWarning) as record:
        design.from_catalog(catalog["PCX-25-100"], epd=20.0, stop_distance=5.0,
                            image_distance=97.363)
    assert record and all(w.filename == __file__ for w in record)
