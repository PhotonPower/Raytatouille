"""Read access to the model tree (ADR 0024): immutable typed copies.

``System.root``, ``paths``, ``fields``, ``aperture`` and ``object_space`` return copies of the
model. The classes below have no constructor and only read-only attributes; a copy holds no
reference to its System and stays valid when the System changes or is freed. Units: lengths in
mm, angles in degree where the name ends in ``_deg``; a Param has the unit of the attribute
that holds it. Variants (shape, aperture, phase, interaction) appear as the object of the active
alternative, e.g. ``Conic`` or ``Plane``.

Example::

    system = rt.load("singlet.rtt.json")
    for node in system.root.children:
        if isinstance(node, rt.model.Element):
            print(node.name, [s.id for s in node.surfaces])
"""

from __future__ import annotations

from ._core import (
    Absorber,
    Assembly,
    CircularAperture,
    CoatingRef,
    Conic,
    CrystalMaterial,
    DiffractionEfficiency,
    Element,
    ElementKind,
    EllipticalAperture,
    Event,
    EventKind,
    EvenAsphere,
    FieldSet,
    FieldType,
    Fresnel,
    IdealAntiReflection,
    IdealBeamSplitter,
    IdealMirror,
    IdealPolarizer,
    IdealRetarder,
    LinearGrating,
    ObjectSpace,
    Param,
    Path,
    Plane,
    Pose,
    RadialPhase,
    RectangularAperture,
    ShapeStack,
    Surface,
    SystemAperture,
    SystemApertureType,
    ZernikeSag,
)

__all__ = [
    "Absorber",
    "Assembly",
    "CircularAperture",
    "CoatingRef",
    "Conic",
    "CrystalMaterial",
    "DiffractionEfficiency",
    "Element",
    "ElementKind",
    "EllipticalAperture",
    "EvenAsphere",
    "Event",
    "EventKind",
    "FieldSet",
    "FieldType",
    "Fresnel",
    "IdealAntiReflection",
    "IdealBeamSplitter",
    "IdealMirror",
    "IdealPolarizer",
    "IdealRetarder",
    "LinearGrating",
    "ObjectSpace",
    "Param",
    "Path",
    "Plane",
    "Pose",
    "RadialPhase",
    "RectangularAperture",
    "ShapeStack",
    "Surface",
    "SystemAperture",
    "SystemApertureType",
    "ZernikeSag",
]
