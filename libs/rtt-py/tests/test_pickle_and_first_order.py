"""#336: System pickles (for multiprocessing) through its canonical file form, and
paraxial.first_order / first_order_section take a System like the other analyses."""

from __future__ import annotations

import pickle

import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR, REPO_ROOT

import raytatouille as rt
from raytatouille import paraxial


def singlet() -> rt.System:
    return rt.load(REFERENCE_DIR / "m1" / "singlet_const.rtt.json")


def test_a_system_pickles_to_an_equal_system() -> None:
    system = rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")
    system.name = "changed in Python"  # a change made from Python travels too
    back = pickle.loads(pickle.dumps(system))
    assert isinstance(back, rt.System)
    assert back.to_json() == system.to_json()
    assert back.name == "changed in Python"


def test_a_pickled_system_compiles_to_the_same_results() -> None:
    # What a worker process of multiprocessing does with it: unpickle and compile.
    system = singlet()
    back = pickle.loads(pickle.dumps(system, protocol=pickle.HIGHEST_PROTOCOL))
    efl = paraxial.first_order(back).efl
    assert efl is not None and efl == paraxial.first_order(system).efl


def test_first_order_takes_a_system_and_a_compiled_system() -> None:
    system = singlet()
    cs = rt.compile(system)
    a = paraxial.first_order(system)
    b = paraxial.first_order(cs, path="main")
    assert a.efl is not None and a.efl == b.efl  # bit for bit: the same compile
    assert a.bfl == b.bfl
    x = paraxial.first_order_section(system, paraxial.Section.Y)
    assert x.efl == a.efl  # rotationally symmetric: the section equals first_order


def test_the_wrappers_pass_every_argument_on() -> None:
    # Review of #394: each argument must reach the compiled call.
    system = singlet()
    cs = rt.compile(system)
    # wavelength: index 0 is not the reference of the singlet, so the EFL differs.
    assert paraxial.first_order(system, wavelength=0).efl == paraxial.first_order(cs, 0, 0).efl
    assert paraxial.first_order(system, wavelength=0).efl != paraxial.first_order(system).efl
    # path: an unknown index is an error of the compiled call.
    with pytest.raises(rt.ParaxialError):
        paraxial.first_order(system, path=5)
    # section: an anamorphic system, X differs from Y.
    biconic = rt.load(REPO_ROOT / "examples" / "anamorphic" / "biconic_singlet.rtt.json")
    x = paraxial.first_order_section(biconic, paraxial.Section.X)
    y = paraxial.first_order_section(biconic, paraxial.Section.Y)
    assert x.efl is not None and y.efl is not None and x.efl != y.efl
    assert x.efl == paraxial.first_order_section(rt.compile(biconic), paraxial.Section.X).efl
    # materials: a catalogue glass needs the library; without it the System does not compile.
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    bk7 = rt.load(REFERENCE_DIR / "m7" / "planoconvex_bk7.rtt.json")
    with_lib = paraxial.first_order(bk7, materials=lib)
    assert with_lib.efl is not None
    assert with_lib.efl == paraxial.first_order(rt.compile(bk7, lib)).efl
    with pytest.raises(rt.CompileError):
        paraxial.first_order(bk7)
    # a CompiledSystem takes no materials (as prescription and seidel).
    with pytest.raises(ValueError):
        paraxial.first_order(cs, materials=lib)
