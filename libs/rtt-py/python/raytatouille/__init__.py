"""raytatouille: Python API of the Raytatouille raytracing engine for optical design.

A thin mirror of the C++ library (ADR 0002). Units: lengths in mm, wavelengths in um
(vacuum), temperature in degC, pressure in atm; right-handed global coordinates with the
optical axis along +z. Results come as NumPy arrays.

Example::

    import raytatouille as rt

    system = rt.load("singlet.rtt.json")
    compiled = rt.compile(system)  # coatings=rt.CoatingLibrary() for coated surfaces
    fo = rt.paraxial.first_order(compiled, path="main")
    spot = rt.analysis.spot(compiled, path="main", field=1, rays="hexapolar:12")
    rays = rt.trace.make_rays(compiled, rt.trace.HexapolarPupil(rings=6))
    stats = rt.trace.trace(compiled, rays)
"""

from . import analysis, diagnostics, errors, layout, materials, paraxial, plot, polar, results, trace
from ._core import (
    CoatingLibrary,
    CompiledMedium,
    CompiledSystem,
    Diagnostic,
    Environment,
    Field,
    Severity,
    System,
    Wavelength,
    compile,
    load,
    save,
    validate,
)
from .materials import MaterialLibrary
from .errors import (
    AgfError,
    AnalysisError,
    CoatingCatalogError,
    CompileError,
    NoStopError,
    ParaxialError,
    ParseError,
    RaytatouilleError,
    RaytatouilleWarning,
    UnknownMaterial,
)

__version__ = "0.4.0"

__all__ = [
    "AgfError",
    "AnalysisError",
    "CoatingCatalogError",
    "CoatingLibrary",
    "CompileError",
    "CompiledMedium",
    "CompiledSystem",
    "Diagnostic",
    "Environment",
    "Field",
    "MaterialLibrary",
    "NoStopError",
    "ParaxialError",
    "ParseError",
    "RaytatouilleError",
    "RaytatouilleWarning",
    "Severity",
    "System",
    "UnknownMaterial",
    "Wavelength",
    "analysis",
    "compile",
    "diagnostics",
    "errors",
    "layout",
    "load",
    "materials",
    "paraxial",
    "plot",
    "polar",
    "results",
    "save",
    "trace",
    "validate",
]
