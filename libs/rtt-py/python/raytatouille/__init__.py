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

from . import (
    analysis,
    diagnostics,
    edit,
    errors,
    layout,
    materials,
    model,
    paraxial,
    plot,
    polar,
    results,
    trace,
)
from ._core import (
    CoatingLibrary,
    CompiledMedium,
    CompiledSystem,
    Diagnostic,
    Environment,
    Field,
    GhostPaths,
    GhostSystem,
    LoadWarning,
    Severity,
    System,
    Wavelength,
    compile,
    compile_with_ghosts,
    ghost_paths,
    load,
    save,
    validate,
)
from .edit import Editor, PatchResult, apply_patch, apply_patch_with_inverse
from .materials import MaterialLibrary
from .trace import CancelToken
from .errors import (
    AgfError,
    AnalysisError,
    Cancelled,
    CoatingCatalogError,
    CompileError,
    EditError,
    NoStopError,
    ParaxialError,
    ParseError,
    RaytatouilleError,
    RaytatouilleWarning,
    UnknownMaterial,
)

__version__ = "0.5.0"

__all__ = [
    "AgfError",
    "AnalysisError",
    "CancelToken",
    "Cancelled",
    "CoatingCatalogError",
    "CoatingLibrary",
    "CompileError",
    "CompiledMedium",
    "CompiledSystem",
    "Diagnostic",
    "EditError",
    "Editor",
    "Environment",
    "Field",
    "GhostPaths",
    "GhostSystem",
    "LoadWarning",
    "MaterialLibrary",
    "NoStopError",
    "ParaxialError",
    "ParseError",
    "PatchResult",
    "RaytatouilleError",
    "RaytatouilleWarning",
    "Severity",
    "System",
    "UnknownMaterial",
    "Wavelength",
    "analysis",
    "apply_patch",
    "apply_patch_with_inverse",
    "compile",
    "compile_with_ghosts",
    "diagnostics",
    "edit",
    "errors",
    "ghost_paths",
    "layout",
    "load",
    "materials",
    "model",
    "paraxial",
    "plot",
    "polar",
    "results",
    "save",
    "trace",
    "validate",
]
