"""Read access to the whole model tree (#82, ADR 0024 point 1): immutable typed copies.

Oracle: the JSON text of every reference system. Each key of the file must appear with its
value in the typed tree. A key that is missing in the file has the value the C++ parser
(rtt-io, json_io.cpp) gives a missing key. The test does not repeat these defaults: it reads
them from DEFAULTS_TEXT, a system in which every optional key is missing, through the same
parser. The mapping of the file's type strings to classes and enum members is the file format
(schema/raytatouille.schema.json) and is written out below.
"""

from __future__ import annotations

import gc
import inspect
import json
import sys
from enum import Enum
from pathlib import Path
from typing import Any

import pytest
from conftest import REFERENCE_DIR, reference_files

import raytatouille as rt
from raytatouille import model

#: Every optional key is missing, so each first instance of a class in this system carries the
#: parser's defaults. Required keys have arbitrary values; the stop-size aperture has no value.
DEFAULTS_TEXT = """{
  "schema_version": "0.2.0",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [{"um": 0.5}],
  "aperture": {"type": "stop_size"},
  "fields": {"points": [{}]},
  "root": {
    "type": "assembly",
    "name": "root",
    "children": [
      {
        "type": "lens",
        "name": "L",
        "surfaces": [
          {"id": "S0"},
          {"id": "S1", "shape": {"base": {"type": "conic", "radius": 1.0}}},
          {"id": "S2", "shape": {"base": {"type": "even_asphere", "radius": 1.0},
                                 "terms": [{"type": "zernike_sag",
                                            "normalization_radius": 1.0}]}},
          {"id": "S3", "aperture": {"type": "circular", "radius": 1.0},
           "phases": [{"type": "linear_grating", "lines_per_mm": 1.0},
                      {"type": "radial_phase", "normalization_radius": 1.0}]},
          {"id": "S4", "interaction": {"type": "ideal_beam_splitter"}},
          {"id": "S5", "interaction": {"type": "ideal_polarizer"}},
          {"id": "S6", "interaction": {"type": "ideal_retarder"}}
        ]
      }
    ]
  },
  "paths": [{"name": "p", "events": [{"surface": "S0"}]}]
}
"""

#: Type strings of the file format and the enum members they stand for.
ENUM_NAMES: dict[type, dict[str, str]] = {
    model.ElementKind: {"lens": "LENS", "mirror": "MIRROR", "plate": "PLATE",
                        "thin_element": "THIN_ELEMENT", "stop": "STOP", "detector": "DETECTOR"},
    model.EventKind: {"refract": "REFRACT", "reflect": "REFLECT", "transmit": "TRANSMIT",
                      "ordinary": "ORDINARY", "extraordinary": "EXTRAORDINARY",
                      "diffract": "DIFFRACT"},
    model.FieldType: {"angle_deg": "ANGLE_DEG", "object_height": "OBJECT_HEIGHT",
                      "paraxial_image_height": "PARAXIAL_IMAGE_HEIGHT"},
    model.SystemApertureType: {"epd": "ENTRANCE_PUPIL_DIAMETER",
                               "image_fnumber": "IMAGE_SPACE_F_NUMBER",
                               "object_na": "OBJECT_SPACE_NA", "stop_size": "STOP_SIZE"},
}

#: Type strings of the variants (and of an assembly) and their classes.
VARIANT_CLASSES: dict[str, type] = {
    "assembly": model.Assembly,
    "plane": model.Plane, "conic": model.Conic, "even_asphere": model.EvenAsphere,
    "zernike_sag": model.ZernikeSag,
    "circular": model.CircularAperture, "rectangular": model.RectangularAperture,
    "elliptical": model.EllipticalAperture,
    "linear_grating": model.LinearGrating, "radial_phase": model.RadialPhase,
    "fresnel": model.Fresnel, "ideal_mirror": model.IdealMirror,
    "ideal_anti_reflection": model.IdealAntiReflection, "absorber": model.Absorber,
    "ideal_beam_splitter": model.IdealBeamSplitter, "coating": model.CoatingRef,
    "ideal_polarizer": model.IdealPolarizer, "ideal_retarder": model.IdealRetarder,
}

#: The new read-only classes of raytatouille.model (ADR 0024 point 1).
READ_ONLY_CLASSES: list[type] = [
    model.Param, model.Pose, model.ShapeStack, model.Surface, model.Element, model.Assembly,
    model.Event, model.Path, model.FieldSet, model.SystemAperture, model.ObjectSpace,
    *(c for c in VARIANT_CLASSES.values() if c is not model.Assembly),
]

#: File keys whose attribute has another name (the JSON key "object" would shadow the builtin
#: in the stubs, ADR 0024 addendum).
ATTRIBUTE_OF_KEY = {"object": "object_space"}


def attributes(obj: Any) -> list[str]:
    """Names of the data attributes (properties) of a bound object."""
    return sorted(name for name, value in inspect.getmembers(type(obj))
                  if not name.startswith("_") and inspect.isdatadescriptor(value))


def is_bound(obj: Any) -> bool:
    return type(obj).__module__.startswith("raytatouille")


def walk(obj: Any) -> list[Any]:
    """All bound objects in the tree below obj (obj included), depth first in attribute order."""
    found: list[Any] = []
    stack = [obj]
    while stack:
        o = stack.pop()
        if isinstance(o, (list, tuple)):
            stack.extend(reversed(o))
        elif is_bound(o) and not isinstance(o, Enum):
            found.append(o)
            stack.extend(reversed([getattr(o, a) for a in attributes(o)]))
    return found


def first_instances(system: Any) -> dict[type, Any]:
    out: dict[type, Any] = {}
    for o in walk(system):
        out.setdefault(type(o), o)
    return out


DEFAULTS = first_instances(rt.System.from_json(DEFAULTS_TEXT))


def check(expected: Any, actual: Any, where: str) -> None:
    """Compares a value of the file (JSON) with the typed copy; see the module docstring."""
    if isinstance(actual, Enum):
        assert actual.name == ENUM_NAMES[type(actual)][expected], where
        return
    if isinstance(actual, list):
        assert isinstance(expected, list) and len(expected) == len(actual), where
        for i, (e, a) in enumerate(zip(expected, actual)):
            check(e, a, f"{where}/{i}")
        return
    if not is_bound(actual):
        assert type(actual) is type(expected) or isinstance(actual, float), where
        assert actual == expected, where
        return
    if isinstance(actual, model.Param) and not isinstance(expected, dict):
        expected = {"value": expected}  # a plain number is a Param with the defaults
    assert isinstance(expected, dict), where
    covered: set[str] = set()
    for key, value in expected.items():
        at = f"{where}/{key}"
        if key == "type":
            if hasattr(actual, "kind"):  # Element
                check(value, actual.kind, at)
                covered.add("kind")
            elif hasattr(actual, "type"):  # FieldSet, SystemAperture
                check(value, actual.type, at)
                covered.add("type")
            else:  # a variant or an assembly
                assert type(actual) is VARIANT_CLASSES[value], at
        elif key == "units" and isinstance(actual, rt.System):
            assert value == {"length": "mm", "wavelength": "um"}, at
        elif key == "material" and isinstance(actual, model.Element):
            if isinstance(value, str):
                assert (actual.material, actual.segment_materials) == (value, []), at
            else:
                assert (actual.material, actual.segment_materials) == (None, value), at
            covered |= {"material", "segment_materials"}
        elif key == "events" and isinstance(actual, model.Path):
            assert actual.automatic == (value == "auto"), at
            check([] if value == "auto" else value, actual.events, at)
            covered |= {"automatic", "events"}
        else:
            name = ATTRIBUTE_OF_KEY.get(key, key)
            assert name in attributes(actual), f"{at}: no attribute {name}"
            check(value, getattr(actual, name), at)
            covered.add(name)
    missing = [name for name in attributes(actual) if name not in covered]
    if missing:  # classes without optional keys need no entry in DEFAULTS
        default = DEFAULTS[type(actual)]
        for name in missing:
            assert getattr(actual, name) == getattr(default, name), f"{where}/{name} (default)"


@pytest.mark.parametrize("file", reference_files(), ids=lambda f: f.name)
def test_typed_tree_matches_the_file(file: Path) -> None:
    check(json.loads(file.read_text(encoding="utf-8")), rt.load(file), "")


def test_the_comparison_detects_differences() -> None:
    # Counter-check of check(): a changed value, a removed key with a non-default value, an
    # extra key, another variant, another enum value and a shorter list must each fail.
    file = REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json"
    system = rt.load(file)
    text = file.read_text(encoding="utf-8")
    check(json.loads(text), system, "")
    changes = ("wavelength", "name", "pose", "radius", "extra key", "variant", "interaction",
               "element type", "event kind", "fewer events")
    for change in changes:
        data = json.loads(text)
        lens = data["root"]["children"][1]  # L2, placed by a pose
        surface = lens["surfaces"][0]
        events = data["paths"][0]["events"]
        if change == "wavelength":
            data["wavelengths"][0]["um"] = 0.5
        elif change == "name":
            lens["name"] = "L9"
        elif change == "pose":
            del lens["pose"]
        elif change == "radius":
            surface["shape"]["base"]["radius"] = 50.0
        elif change == "extra key":
            surface["extra"] = 1
        elif change == "variant":
            surface["shape"]["base"]["type"] = "even_asphere"
        elif change == "interaction":
            surface["interaction"] = {"type": "absorber"}
        elif change == "element type":
            lens["type"] = "plate"
        elif change == "event kind":
            events[0]["kind"] = "reflect"
        else:
            del events[-1]
        with pytest.raises(AssertionError):
            check(data, system, change)


def test_keys_missing_in_the_reference_files() -> None:
    # No reference file sets fast_axis or retardance_waves; check them on a changed copy of
    # the defaults system.
    data = json.loads(DEFAULTS_TEXT)
    retarder = data["root"]["children"][0]["surfaces"][6]["interaction"]
    retarder.update(fast_axis=[0.0, 1.0, 0.0], retardance_waves=0.5)
    system = rt.System.from_json(json.dumps(data))
    check(data, system, "")
    lens = system.root.children[0]
    assert isinstance(lens, model.Element)
    interaction = lens.surfaces[6].interaction
    assert isinstance(interaction, model.IdealRetarder)
    assert interaction.fast_axis == [0.0, 1.0, 0.0] and interaction.retardance_waves == 0.5


def json_objects(value: Any, pointer: str = "") -> list[tuple[str, dict[str, Any]]]:
    """All JSON objects below value with their JSON pointers."""
    out: list[tuple[str, dict[str, Any]]] = []
    if isinstance(value, dict):
        out.append((pointer, value))
        for k, v in value.items():
            out += json_objects(v, f"{pointer}/{k}")
    elif isinstance(value, list):
        for i, v in enumerate(value):
            out += json_objects(v, f"{pointer}/{i}")
    return out


def test_defaults_system_has_no_optional_values() -> None:
    # The source of the defaults is only valid if every key the parser does not require is a
    # carrier of a later instance, never a value of a first instance. A key is optional if the
    # system still parses without it; the parser decides, not a list in this test.
    data = json.loads(DEFAULTS_TEXT)
    optional = []
    for pointer, obj in json_objects(data):
        for key in list(obj):
            value = obj.pop(key)
            try:
                rt.System.from_json(json.dumps(data))
                optional.append(f"{pointer}/{key}")
            except rt.ParseError:
                pass
            obj[key] = value
    s = "/root/children/0/surfaces"
    assert sorted(optional) == sorted([
        f"{s}/1/shape", f"{s}/1/shape/base", f"{s}/2/shape", f"{s}/2/shape/base",
        f"{s}/2/shape/terms", f"{s}/3/aperture", f"{s}/3/phases", f"{s}/4/interaction",
        f"{s}/5/interaction", f"{s}/6/interaction",
    ])
    # The carriers come after the first instances of their classes.
    assert DEFAULTS[model.Surface].id == "S0"
    assert DEFAULTS[model.ShapeStack].base == DEFAULTS[model.Plane]
    assert isinstance(DEFAULTS[model.Surface].interaction, model.Fresnel)
    assert DEFAULTS[model.Surface].aperture is None and DEFAULTS[model.Surface].phases == []


def test_every_variant_class_appears_in_the_reference_systems() -> None:
    seen = {type(o) for f in reference_files() for o in walk(rt.load(f))}
    assert set(VARIANT_CLASSES.values()) <= seen
    tour = {type(o) for o in walk(rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json"))}
    assert set(VARIANT_CLASSES.values()) - {model.IdealMirror} <= tour


def test_copies_are_read_only_and_not_constructible() -> None:
    instances = {type(o): o for f in reference_files() for o in walk(rt.load(f))}
    for cls in READ_ONLY_CLASSES:
        obj = instances[cls]
        for name in attributes(obj):
            with pytest.raises(AttributeError):
                setattr(obj, name, getattr(obj, name))
        with pytest.raises(TypeError):
            cls()
        with pytest.raises(TypeError):
            hash(obj)


def test_copies_hold_no_reference_to_the_system() -> None:
    system = rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json")
    count = sys.getrefcount(system)
    kept = [system.root, system.paths, system.fields, system.aperture, system.object_space,
            system.wavelengths]
    assert sys.getrefcount(system) == count  # copies, not views kept alive by the system
    assert len(kept) == 6


def test_copies_outlive_the_system() -> None:
    system = rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json")
    root, paths = system.root, system.paths
    pose = system.root.children[1].pose  # nested: its parents are temporaries
    expected = [c.name for c in root.children]
    del system
    gc.collect()
    assert [c.name for c in root.children] == expected
    assert paths[0].name == "main"
    assert pose.position[2].value == 10.691


def test_field_points_are_copies() -> None:
    # Field is mutable (since before ADR 0024); FieldSet.points must hand out copies so that
    # the FieldSet and the System stay unchanged.
    system = rt.load(REFERENCE_DIR / "m2" / "cooke_triplet.rtt.json")
    fields = system.fields
    point = fields.points[0]
    point.x = 99.0
    fields.points.append(rt.Field(1.0, 2.0))
    assert fields.points[0].x == 0.0 and len(fields.points) == 3
    assert system.fields == fields


def test_each_access_is_an_independent_equal_copy() -> None:
    system = rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")
    a, b = system.root, system.root
    assert a is not b
    assert a == b
    assert a.children[0] != a.children[1]
    system.name = "renamed"  # changing the system does not touch a copy taken before
    assert a == b
    assert system.fields == rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json").fields


def element(file: str, name: str) -> model.Element:
    found = [n for n in rt.load(REFERENCE_DIR / file).root.children if n.name == name]
    assert len(found) == 1 and isinstance(found[0], model.Element)
    return found[0]


def test_element_segment_material() -> None:
    achromat = element("m2/achromat.rtt.json", "L1")
    assert achromat.segment_materials == ["SCHOTT:N-BK7", "SCHOTT:F2"]
    assert achromat.segment_material(0) == "SCHOTT:N-BK7"
    assert achromat.segment_material(1) == "SCHOTT:F2"
    assert achromat.segment_material(2) is None  # no third segment
    singlet = element("m0/singlet.rtt.json", "L1")
    assert singlet.segment_material(0) == singlet.material == "SCHOTT:N-BK7"
    assert singlet.segment_material(1) is None


def test_wavelength_equality_and_hash() -> None:
    first = rt.load(REFERENCE_DIR / "m0" / "singlet.rtt.json").wavelengths
    second = rt.load(REFERENCE_DIR / "m0" / "singlet.rtt.json").wavelengths
    assert first == second and first[0] is not second[0]
    assert first[0] != first[1]
    assert {hash(w) for w in first} == {hash(w) for w in second}
    assert len(set(first) | set(second)) == len(first)


def test_object_space_attribute() -> None:
    # ADR 0024 addendum: "object" would shadow the builtin in the stubs (mypy valid-type).
    tour = rt.load(REFERENCE_DIR / "m0" / "feature_tour.rtt.json")
    assert not tour.object_space.at_infinity
    assert tour.object_space.distance.value == 250.0
    assert tour.object_space.distance.variable
    assert not hasattr(tour, "object")
