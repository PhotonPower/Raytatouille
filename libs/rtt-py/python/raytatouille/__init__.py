"""raytatouille: Python API of the Raytatouille raytracing engine for optical design.

A thin mirror of the C++ library (ADR 0002). Units: lengths in mm, wavelengths in um
(vacuum), temperature in degC, pressure in atm; right-handed global coordinates with the
optical axis along +z. Results come as NumPy arrays.

Example::

    import raytatouille as rt

    system = rt.load("singlet.rtt.json")
    compiled = rt.compile(system)
    fo = rt.paraxial.first_order(compiled, path="main")
    spot = rt.analysis.spot(compiled, path="main", field=1, rays="hexapolar:12")
    rays = rt.trace.make_rays(compiled, rt.trace.HexapolarPupil(rings=6))
    stats = rt.trace.trace(compiled, rays)
"""

from . import analysis, errors, paraxial, plot, trace
from ._core import (
    CompiledMedium,
    CompiledSystem,
    Diagnostic,
    Environment,
    Field,
    MaterialLibrary,
    Severity,
    System,
    Wavelength,
    compile,
    load,
    save,
    validate,
)
from .errors import (
    AgfError,
    AnalysisError,
    CompileError,
    ParaxialError,
    ParseError,
    RaytatouilleError,
    UnknownMaterial,
)

__version__ = "0.2.0"

__all__ = [
    "AgfError",
    "AnalysisError",
    "CompileError",
    "CompiledMedium",
    "CompiledSystem",
    "Diagnostic",
    "Environment",
    "Field",
    "MaterialLibrary",
    "ParaxialError",
    "ParseError",
    "RaytatouilleError",
    "Severity",
    "System",
    "UnknownMaterial",
    "Wavelength",
    "analysis",
    "compile",
    "errors",
    "load",
    "paraxial",
    "plot",
    "save",
    "trace",
    "validate",
]
