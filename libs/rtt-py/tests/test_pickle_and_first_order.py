"""#336: System pickles (for multiprocessing) through its canonical file form, and
paraxial.first_order / first_order_section take a System like the other analyses."""

from __future__ import annotations

import pickle

from conftest import REFERENCE_DIR

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
