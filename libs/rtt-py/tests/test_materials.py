"""Material library for a GUI (#85, G9): listing, alias, loading from memory, vectorized index,
glass map. Expected values are the lines of the test catalogues (AGF records after the Ansys
OpticStudio User Guide 2025 R1; units checked against the SCHOTT N-BK7 data sheet,
docs/quellen.md)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import raytatouille as rt
from raytatouille.materials import ClassRange, GlassInfo, MechanicalData, ThermalData

# Invented catalogue with every listed record and the manufacturers' placeholders.
ALL = (
    "CC listing test\n"
    "NM ALL 2 517642.251 1.5168 64.17 1 3 5\n"
    "GC   free text comment  \n"
    "ED 7.1 8.3 2.51 -0.0009 0\n"
    "CD 1.039612120E+00 6.000698670E-03 2.317923440E-01 2.001791440E-02 1.010469450E+00 "
    "1.035606530E+02 0 0\n"
    "TD 1.86e-6 1.31e-8 -1.37e-11 4.34e-7 6.27e-10 0.17 20\n"
    "MD 82.00 0.21 610 _ 1.11\n"
    "OD 1.0 1-2 0.0 _ 2.3 -1\n"
    "LD 0.3 2.5\n"
    "IT 0.3 0.05 25\n"
    "IT 0.31 0.25 25\n"
    "NM OLD 5 1 1.5 60 0 2 -\n"
    "CD 1 0.01 0.001\n"
    "OD - 3.0 2.0 52.3 4.3 4.3\n"  # as SCHOTT line 6427: "-" relative cost
)


def test_listing_matches_the_agf_lines() -> None:
    lib = rt.MaterialLibrary()
    lib.add_catalog_text(ALL, "ALLCAT", "all.agf")
    assert lib.catalogs() == ["ALLCAT"]
    first, old = lib.glasses("ALLCAT")
    assert isinstance(first, GlassInfo)
    assert first.reference == "ALLCAT:ALL"
    assert first.line == 2
    assert first.formula == 2 and first.supported and first.unsupported_reason is None
    assert first.nd == 1.5168 and first.vd == 64.17
    assert first.exclude_substitution is True
    assert first.status == 3 and first.melt_frequency == 5
    assert first.comment == "free text comment"
    assert first.tce_m30_70 == 7.1 and first.tce_100_300 == 8.3
    assert first.density_g_per_cm3 == 2.51 and first.dpgf == -0.0009
    assert first.ignore_thermal_expansion is False
    assert first.thermal == ThermalData(1.86e-6, 1.31e-8, -1.37e-11, 4.34e-7, 6.27e-10, 0.17, 20.0)
    assert first.mechanical == MechanicalData(82.0, 0.21, 610.0, None, 1.11)
    assert first.relative_cost == 1.0
    assert first.climate_resistance == ClassRange(1.0, 2.0)
    assert first.stain_resistance == ClassRange(0.0, 0.0)
    assert first.acid_resistance is None  # "_"
    assert first.alkali_resistance == ClassRange(2.3, 2.3)
    assert first.phosphate_resistance is None  # -1: not available
    assert first.wavelength_range_um == (0.3, 2.5)
    assert first.transmission.shape == (2, 3)
    assert np.array_equal(first.transmission, [[0.3, 0.05, 25.0], [0.31, 0.25, 25.0]])
    # Conrady (formula 5) is listed, but not supported (#42); "-" melt frequency: not given.
    assert old.formula == 5 and not old.supported
    assert old.unsupported_reason is not None and "#42" in old.unsupported_reason
    assert old.status == 2 and old.melt_frequency is None
    assert old.thermal is None and old.mechanical is None and old.relative_cost is None
    assert old.climate_resistance == ClassRange(3.0, 3.0)
    assert old.transmission.shape == (0, 3)
    with pytest.raises(rt.UnknownMaterial):
        lib.index("ALLCAT:OLD", 0.55)
    assert lib.glass("ALLCAT:ALL").reference == first.reference


def test_listing_of_a_manufacturer_excerpt(catalog_dir: Path) -> None:
    # tests/catalogs/nikon/nikon-hikari.agf (byte copy, #42), NICF-V:
    #   NM NICF-V 6 1 1.433837 95.260792 0 0 0
    #   OD -1.00000 1.00000 0.00000 4.50000 2.30000 1.00000
    #   36 IT lines from 0.157/0.995/10 to 2.5/0.998/10
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "nikon" / "nikon-hikari.agf")
    g = lib.glass("NIKON-HIKARI:NICF-V")
    assert (g.formula, g.nd, g.vd) == (6, 1.433837, 95.260792)
    assert g.status == 0 and g.exclude_substitution is False and g.melt_frequency == 0
    assert g.comment == "TCE value is available for 0 to 25 degrees Celsius."
    assert g.relative_cost is None  # -1
    assert g.climate_resistance == ClassRange(1.0, 1.0)
    assert g.acid_resistance == ClassRange(4.5, 4.5)
    assert g.density_g_per_cm3 == 3.18
    assert g.transmission.shape == (36, 3)
    assert tuple(g.transmission[0]) == (0.157, 0.995, 10.0)
    assert tuple(g.transmission[-1]) == (2.5, 0.998, 10.0)
    assert [x.name for x in lib.glasses("NIKON-HIKARI")] == [
        "NICF-V", "NIFS-V", "E-LAKH1", "E-KZFH1", "J-SFH1", "J-FK5"]
    with pytest.raises(rt.UnknownMaterial):
        lib.glasses("NIKON")
    with pytest.raises(rt.UnknownMaterial):
        lib.glass("NIKON-HIKARI:NOPE")


def test_alias_loads_two_files_of_the_same_name(catalog_dir: Path) -> None:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    with pytest.raises(ValueError):
        lib.add_catalog(catalog_dir / "m2" / "schott.agf")  # unchanged without an alias
    lib.add_catalog(catalog_dir / "m2" / "schott.agf", name="SCHOTT_M2")
    assert lib.catalogs() == ["SCHOTT", "SCHOTT_M2"]
    plain = rt.MaterialLibrary()
    plain.add_catalog(catalog_dir / "m2" / "schott.agf")
    assert lib.index("SCHOTT_M2:N-LAK9", 0.5876) == plain.index("SCHOTT:N-LAK9", 0.5876)
    for bad in ("", "A:B", "A B", "CONST", "SCHOTT"):
        with pytest.raises(ValueError):
            lib.add_catalog(catalog_dir / "m2" / "schott.agf", name=bad)
    with pytest.raises(ValueError):
        lib.add_catalog(catalog_dir / "m2", name="M2")  # a name needs a single file
    assert lib.catalogs() == ["SCHOTT", "SCHOTT_M2"]
    # Case-sensitive like every reference: "schott" is another catalogue.
    lib.add_catalog(catalog_dir / "schott.agf", name="schott")
    assert lib.catalogs() == ["SCHOTT", "SCHOTT_M2", "schott"]


def test_memory_and_file_give_the_same_catalogue(catalog_dir: Path) -> None:
    for file in (catalog_dir / "schott.agf", catalog_dir / "utf16" / "schott.agf"):
        from_file = rt.MaterialLibrary()
        from_file.add_catalog(file, name="CAT")
        from_bytes = rt.MaterialLibrary()
        from_bytes.add_catalog_text(file.read_bytes(), "CAT", file.name)
        a = from_file.glasses("CAT")
        b = from_bytes.glasses("CAT")
        assert [g.name for g in a] == [g.name for g in b]
        assert [g.coefficients for g in a] == [g.coefficients for g in b]
        for g in a:
            assert from_file.index(g.reference, 0.55, 30.0, 0.9) == from_bytes.index(
                g.reference, 0.55, 30.0, 0.9)
    lib = rt.MaterialLibrary()
    with pytest.raises(rt.AgfError, match="pasted:2"):
        lib.add_catalog_text("CC c\nXX 1\n", "MEM", "pasted")
    assert lib.catalogs() == []


def test_vectorized_index_is_bit_identical_with_the_scalar_one(catalog_dir: Path) -> None:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    lib.add_catalog(catalog_dir / "nikon" / "nikon-hikari.agf")
    # dtype given: NumPy < 2.3 types a plain linspace as floating[Any] (mypy on Python 3.10).
    wl = np.linspace(0.4, 0.7, 1000, dtype=np.float64)
    for ref in ("SCHOTT:N-BK7", "SCHOTT:F2", "NIKON-HIKARI:NICF-V", "NIKON-HIKARI:E-LAKH1",
                "NIKON-HIKARI:J-SFH1", "AIR", "CONST:1.5,1e-3"):
        many = lib.index(ref, wl, 23.5, 0.95)
        assert many.dtype == np.complex128 and many.shape == wl.shape
        one = np.array([lib.index(ref, float(x), 23.5, 0.95) for x in wl])
        assert many.tobytes() == one.tobytes(), ref  # bit for bit
    grid = wl.reshape(25, 40)
    assert lib.index("SCHOTT:N-BK7", grid).shape == (25, 40)
    assert np.array_equal(lib.index("SCHOTT:N-BK7", grid).ravel(), lib.index("SCHOTT:N-BK7", wl))
    # 0-d, strided and float32 arrays are converted like any NumPy input; lists are not accepted.
    zero_d = lib.index("SCHOTT:N-BK7", np.array(0.55))
    assert zero_d.shape == () and zero_d == lib.index("SCHOTT:N-BK7", 0.55)
    strided = lib.index("SCHOTT:N-BK7", wl[::200])
    assert strided.tobytes() == lib.index("SCHOTT:N-BK7", np.ascontiguousarray(wl[::200])).tobytes()
    single = lib.index("SCHOTT:N-BK7", np.array([0.55], dtype=np.float32))
    assert single[0] == lib.index("SCHOTT:N-BK7", float(np.float32(0.55)))
    with pytest.raises(TypeError):
        lib.index("SCHOTT:N-BK7", [0.5, 0.6])  # type: ignore[call-overload]


def test_glass_map_uses_the_nm_records(catalog_dir: Path) -> None:
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "schott.agf")
    lib.add_catalog(catalog_dir / "nikon" / "nikon-hikari.agf")
    chart = lib.glass_map()
    infos = lib.glasses("NIKON-HIKARI") + lib.glasses("SCHOTT")  # catalogues in ascending order
    assert chart.reference == [g.reference for g in infos]
    assert np.array_equal(chart.nd, [g.nd for g in infos])
    assert np.array_equal(chart.vd, [g.vd for g in infos])
    assert chart.reference[-2:] == ["SCHOTT:N-BK7", "SCHOTT:F2"]
    assert tuple(chart.nd[-2:]) == (1.5168, 1.62004)
    assert tuple(chart.status[-2:]) == (1, 1)
    only = lib.glass_map(["SCHOTT"])
    assert only.reference == ["SCHOTT:N-BK7", "SCHOTT:F2"]


def test_material_library_still_compiles_systems(catalog_dir: Path, reference_dir: Path) -> None:
    # A system with catalogue glasses (Cooke triplet, #34) through the Python subclass.
    lib = rt.MaterialLibrary()
    lib.add_catalog(catalog_dir / "m2" / "schott.agf")
    assert isinstance(lib, rt._core.MaterialLibrary)
    system = rt.load(reference_dir / "m2" / "cooke_triplet.rtt.json")
    compiled = rt.compile(system, lib)
    n = lib.index("SCHOTT:N-LAK9", 0.5875618, system.environment.temperature_c, 1.0)
    assert n.real > 1.69
    assert compiled is not None


def test_supported_matches_index(catalog_dir: Path) -> None:
    # S1 of the review: supported is False exactly when index() raises, also for a TD record
    # without 7 values (CatalogMaterial rejects it like an unverified formula).
    lib = rt.MaterialLibrary()
    lib.add_catalog_text(ALL + "NM TD6 2 1 1.5 60\nCD 1 0.01 0 0 0 0\nTD 1 2 3 4 5 6\n", "ALLCAT")
    lib.add_catalog(catalog_dir / "nikon" / "nikon-hikari.agf")
    td6 = lib.glass("ALLCAT:TD6")
    assert not td6.supported and td6.unsupported_reason is not None and "TD" in td6.unsupported_reason
    for catalog in lib.catalogs():
        for g in lib.glasses(catalog):
            try:
                lib.index(g.reference, 0.55)
                raised = False
            except rt.UnknownMaterial:
                raised = True
            assert raised == (not g.supported), g.reference


def test_names_of_catalogues_from_memory() -> None:
    lib = rt.MaterialLibrary()
    for bad in ("", "A:B", "A B", "CONST"):
        with pytest.raises(ValueError):
            lib.add_catalog_text(ALL, bad)
    assert lib.catalogs() == []
    first = lib.glasses
    lib.add_catalog_text(ALL.encode("utf-8"), "BYTES")
    g = first("BYTES")[0]
    assert not g.transmission.flags.writeable
