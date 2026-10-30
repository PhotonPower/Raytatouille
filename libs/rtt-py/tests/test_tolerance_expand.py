"""Tolerance classes from Python (#261, part 1; ADR 0035, point 5; ADR 0036, points 9 and 10):
rt.tolerance.expand, class_limit and abbe_share.

The limits of the class table are decided in ADR 0036, addendum #260. The abbe limit of a
generated item is abbe_share(class) * nu_d of its glass; nu_d is not available in Python (the
line wavelengths are not bound), so these tests check it for plausibility against V(d) of the
catalogue and the bitwise comparison with C++ checks it exactly (test_bitwise_tolerance.py)."""

from __future__ import annotations

import pytest
from conftest import CATALOG_DIR, REFERENCE_DIR

import raytatouille as rt
from raytatouille import model, tolerance

FILE = REFERENCE_DIR / "m7" / "achromat_tolerances.rtt.json"

C = model.ToleranceClass
K = model.ToleranceKind
L = tolerance.ToleranceLevel

# The items the classes generate in the reference file, in tree order (ADR 0035, point 5):
# class normal; objective precise (and so L1 and its surfaces); image none.
GENERATED: list[tuple[str, model.ToleranceClass, tolerance.ToleranceLevel]] = [
    ("class:node:stop:thickness", C.NORMAL, L.ELEMENT),
    ("class:node:stop:decenter_x", C.NORMAL, L.ELEMENT),
    ("class:node:stop:decenter_y", C.NORMAL, L.ELEMENT),
    ("class:node:stop:tilt_x_deg", C.NORMAL, L.ELEMENT),
    ("class:node:stop:tilt_y_deg", C.NORMAL, L.ELEMENT),
    ("class:node:objective:decenter_x", C.PRECISE, L.ASSEMBLY),
    ("class:node:objective:decenter_y", C.PRECISE, L.ASSEMBLY),
    ("class:node:objective:shift_z", C.PRECISE, L.ASSEMBLY),
    ("class:node:objective:tilt_x_deg", C.PRECISE, L.ASSEMBLY),
    ("class:node:objective:tilt_y_deg", C.PRECISE, L.ASSEMBLY),
    ("class:node:L1:thickness", C.PRECISE, L.ELEMENT),
    ("class:node:L1:decenter_x", C.PRECISE, L.ELEMENT),
    ("class:node:L1:decenter_y", C.PRECISE, L.ELEMENT),
    ("class:node:L1:tilt_x_deg", C.PRECISE, L.ELEMENT),
    ("class:node:L1:tilt_y_deg", C.PRECISE, L.ELEMENT),
    ("class:node:L1:index:0", C.PRECISE, L.ELEMENT),
    ("class:node:L1:index:1", C.PRECISE, L.ELEMENT),
    ("class:node:L1:abbe:0", C.PRECISE, L.ELEMENT),
    ("class:node:L1:abbe:1", C.PRECISE, L.ELEMENT),
    ("class:surface:L1.S2:curvature_relative", C.PRECISE, L.SURFACE),
    ("class:surface:L1.S2:thickness", C.PRECISE, L.SURFACE),
    ("class:surface:L1.S3:curvature_relative", C.PRECISE, L.SURFACE),
    ("class:surface:L1.S3:thickness", C.PRECISE, L.SURFACE),
]
EXPLICIT = ["S3-wedge", "S1-conic", "R1"]

# V(d) of the NM lines in tests/catalogs/schott.agf: N-BK7 (segment 0), F2 (segment 1).
CATALOG_VD = {0: 64.17, 1: 36.37}


def schott() -> rt.MaterialLibrary:
    lib = rt.MaterialLibrary()
    lib.add_catalog(CATALOG_DIR / "schott.agf")
    return lib


def find(r: tolerance.ToleranceItems, id: str) -> tolerance.ClassTolerance:
    found = [c for c in r.classes.items if c.id == id]
    assert len(found) == 1, id
    return found[0]


def test_items_of_the_reference_file() -> None:
    system = rt.load(FILE)
    r = tolerance.expand(system, materials=schott())
    assert isinstance(r, tolerance.ToleranceItems)
    assert [i.id for i in r.items] == [g[0] for g in GENERATED] + EXPLICIT
    assert len(r.origins) == len(r.items)
    for item, origin, (gid, cls, level) in zip(r.items, r.origins, GENERATED):
        assert origin.source == tolerance.ItemOriginSource.CLASS, gid
        c = r.classes.items[origin.index]
        assert (c.id, c.generated, c.tolerance_class, c.level) == (gid, True, cls, level)
        assert (item.surface, item.node, item.parameter, item.param, item.segment, item.kind) == (
            c.surface, c.node, None, None, c.segment, c.kind)
        if gid == "class:surface:L1.S2:curvature_relative":  # single value of the file
            assert (item.limit, item.distribution) == (0.005, model.ToleranceDistribution.NORMAL)
        elif gid == "class:node:L1:index:1":  # single value of the file
            assert (item.limit, item.distribution) == (0.001, model.ToleranceDistribution.UNIFORM)
        elif item.kind == K.ABBE:
            assert item.segment is not None
            nu_d = item.limit / tolerance.abbe_share(cls)
            assert abs(nu_d - CATALOG_VD[item.segment]) < 0.1, gid
            assert item.distribution == model.ToleranceDistribution.NORMAL
        else:
            assert item.limit == tolerance.class_limit(cls, item.kind, level), gid
            assert item.distribution == model.ToleranceDistribution.NORMAL
    for i, (item, origin) in enumerate(zip(r.items[len(GENERATED):], r.origins[len(GENERATED):])):
        assert item == system.tolerances.items[i]
        assert (origin.source, origin.index) == (tolerance.ItemOriginSource.EXPLICIT, i)

    # Skipped and replaced kinds stay in `classes`, with the reason.
    S = tolerance.ClassSkip
    assert find(r, "class:node:stop:index").skip == S.NO_MATERIAL_IN_BEAM
    assert find(r, "class:surface:STO:curvature_relative").skip == S.PLANE
    assert find(r, "class:surface:L1.S1:curvature_relative").skip == S.BOUND
    conic = find(r, "class:surface:L1.S1:conic")
    assert (conic.skip, conic.replaced_by, conic.generated) == (S.SPHERE, 1, False)
    assert [(v.index, v.effect, v.skip) for v in r.classes.value_overrides] == [
        (1, tolerance.OverrideEffect.APPLIED, S.NONE),
        (2, tolerance.OverrideEffect.APPLIED, S.NONE)]


def test_class_limits_and_the_abbe_share() -> None:
    # ADR 0036, addendum #260 (decided by the maintainer): a few rows, all classes.
    rows = [
        (K.THICKNESS, L.SURFACE, (0.25, 0.1, 0.025)),
        (K.TILT_X_DEG, L.ASSEMBLY, (0.2, 0.1, 0.05)),
        (K.CURVATURE_RELATIVE, L.SURFACE, (5e-3, 2e-3, 5e-4)),
        (K.CONIC, L.SURFACE, (0.05, 0.02, 0.005)),
        (K.INDEX, L.ELEMENT, (1e-3, 5e-4, 2e-4)),
    ]
    for kind, level, limits in rows:
        for cls, limit in zip((C.COARSE, C.NORMAL, C.PRECISE), limits):
            assert tolerance.class_limit(cls, kind, level) == limit
    # abbe is a share of nu_d (SCHOTT TIE-29, table 5-1, steps 4, 3, 1).
    assert [tolerance.abbe_share(c) for c in (C.COARSE, C.NORMAL, C.PRECISE)] == [
        0.008, 0.005, 0.002]
    with pytest.raises(ValueError, match="abbe_share"):
        tolerance.class_limit(C.NORMAL, K.ABBE, L.ELEMENT)
    with pytest.raises(ValueError):
        tolerance.abbe_share(C.NONE)
    with pytest.raises(ValueError):
        tolerance.class_limit(C.NONE, K.THICKNESS, L.ELEMENT)
    with pytest.raises(ValueError):
        tolerance.class_limit(C.NORMAL, K.RADIUS, L.SURFACE)  # no class generates it


def test_the_glass_of_a_generated_abbe_item_must_resolve() -> None:
    # Without the catalogue: material.unknown at the material of L1 (child 0 of the objective).
    with pytest.raises(rt.CompileError) as e:
        tolerance.expand(rt.load(FILE))
    assert [(d.code, d.location) for d in e.value.diagnostics] == [
        ("material.unknown", "/root/children/1/children/0/material/0"),
        ("material.unknown", "/root/children/1/children/0/material/1")]


def test_an_invalid_set_raises_compile_error() -> None:
    # thickness at surface 0 (no class generates it), ADR 0035, addendum #253.
    bad = rt.apply_patch(rt.load(FILE), [{
        "op": "add", "path": "/tolerances/overrides/-",
        "value": {"surface": "L1.S1", "kind": "thickness", "limit": 0.01}}],
        check="structure_only")
    with pytest.raises(rt.CompileError) as e:
        tolerance.expand(bad, materials=schott())
    assert ("tolerances.target_invalid", "/tolerances/overrides/5/kind") in [
        (d.code, d.location) for d in e.value.diagnostics]


def test_without_a_class_only_the_explicit_items() -> None:
    system = rt.apply_patch(rt.load(FILE), [
        {"op": "replace", "path": "/tolerances/class", "value": "none"},
        {"op": "remove", "path": "/tolerances/overrides/4"},  # image none
        {"op": "remove", "path": "/tolerances/overrides/0"},  # objective precise
    ])
    r = tolerance.expand(system)  # no glass is needed
    assert r.classes.items == []
    assert [i.id for i in r.items] == EXPLICIT
    assert [(o.source, o.index) for o in r.origins] == [
        (tolerance.ItemOriginSource.EXPLICIT, i) for i in range(3)]


def test_an_edit_changes_the_items_and_undo_restores_them() -> None:
    # The tolerance set is edited by patch (#253); expand follows the edited system.
    lib = schott()
    editor = rt.Editor(rt.load(FILE))
    before = tolerance.expand(editor.system, materials=lib)
    editor.apply([{"op": "add", "path": "/tolerances/overrides/-",
                   "value": {"node": "L1", "class": "coarse"}}])
    after = tolerance.expand(editor.system, materials=lib)
    old = {i.id: i for i in before.items}
    new = {i.id: i for i in after.items}
    assert old.keys() == new.keys()
    assert new["class:node:L1:decenter_x"].limit == 0.25  # coarse instead of precise
    assert old["class:node:L1:decenter_x"].limit == 0.05
    assert new["class:node:stop:decenter_x"] == old["class:node:stop:decenter_x"]
    assert new["class:node:L1:index:1"].limit == 0.001  # the single value goes first
    editor.undo()
    again = tolerance.expand(editor.system, materials=lib)
    assert again.to_json() == before.to_json()


def test_the_types_are_read_only() -> None:
    r = tolerance.expand(rt.load(FILE), materials=schott())
    with pytest.raises(AttributeError):
        setattr(r, "items", [])
    with pytest.raises(AttributeError):
        setattr(r.classes.items[0], "skip", tolerance.ClassSkip.NONE)
    with pytest.raises(AttributeError):
        setattr(r.origins[0], "index", 3)
