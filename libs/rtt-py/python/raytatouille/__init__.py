"""raytatouille: raytracing engine for optical design (agents: raytatouille.llm_reference()).

Agents and language models: start with ``raytatouille.llm_reference()``, a compact reference
(workflow, conventions, file format, recipes, API index) shipped with the package as
``raytatouille/llm_reference.md``.

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
    design,
    diagnostics,
    edit,
    errors,
    layout,
    materials,
    model,
    optim,
    paraxial,
    plot,
    polar,
    reports,
    results,
    schemas,
    tolerance,
    trace,
)
from ._core import (
    SCHEMA_VERSION,
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
    ApertureError,
    Cancelled,
    CoatingCatalogError,
    CompileError,
    DesignError,
    EditError,
    LensCatalogError,
    NoStopError,
    OptimError,
    ParaxialError,
    ParseError,
    RaytatouilleError,
    RaytatouilleWarning,
    UnknownMaterial,
)

from .schemas import schema


def llm_reference() -> str:
    """The compact reference for language models and agents (Markdown, English, #408).

    Workflow, conventions, the essentials of the file format, runnable recipes and a generated
    index of all public functions; the file ``raytatouille/llm_reference.md`` in the package.
    """
    from importlib import resources

    return (resources.files(__package__) / "llm_reference.md").read_text(encoding="utf-8")


__version__ = "0.10.0"

__all__ = [
    "AgfError",
    "AnalysisError",
    "ApertureError",
    "CancelToken",
    "Cancelled",
    "CoatingCatalogError",
    "CoatingLibrary",
    "CompileError",
    "CompiledMedium",
    "CompiledSystem",
    "Diagnostic",
    "DesignError",
    "EditError",
    "Editor",
    "Environment",
    "Field",
    "GhostPaths",
    "GhostSystem",
    "LensCatalogError",
    "LoadWarning",
    "MaterialLibrary",
    "NoStopError",
    "OptimError",
    "ParaxialError",
    "ParseError",
    "PatchResult",
    "RaytatouilleError",
    "RaytatouilleWarning",
    "SCHEMA_VERSION",
    "Severity",
    "System",
    "UnknownMaterial",
    "Wavelength",
    "analysis",
    "apply_patch",
    "apply_patch_with_inverse",
    "compile",
    "compile_with_ghosts",
    "design",
    "diagnostics",
    "edit",
    "errors",
    "ghost_paths",
    "layout",
    "llm_reference",
    "load",
    "materials",
    "model",
    "optim",
    "paraxial",
    "plot",
    "polar",
    "reports",
    "results",
    "save",
    "schema",
    "schemas",
    "tolerance",
    "trace",
    "validate",
]
