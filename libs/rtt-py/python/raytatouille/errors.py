"""Exceptions of raytatouille.

Every error of the C++ library arrives as one of these classes (or as the built-in exception
nanobind maps a standard C++ exception to: ValueError, IndexError, ...). The message is the
C++ message. Ray problems during tracing are not exceptions but ray status flags (RayStatus).
Where the place in the system file is known, the errors carry it as a JSON pointer
(``location``) and, for a surface, its id (``surface``); see ADR 0022. Warnings arrive as
RaytatouilleWarning with a stable code (docs/diagnostics.md).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._core import Diagnostic, RayStatus

__all__ = [
    "AgfError",
    "AnalysisError",
    "Cancelled",
    "CoatingCatalogError",
    "CompileError",
    "EditError",
    "NoStopError",
    "ParaxialError",
    "ParseError",
    "RaytatouilleError",
    "RaytatouilleWarning",
    "UnknownMaterial",
]


class RaytatouilleError(Exception):
    """Base class of all raytatouille errors."""


class ParseError(RaytatouilleError, ValueError):
    """Structural error in a system file; ``pointer`` is a JSON pointer to the value."""

    def __init__(self, message: str, pointer: str) -> None:
        super().__init__(message)
        self.pointer = pointer


class CompileError(RaytatouilleError, ValueError):
    """Invalid model, unknown material or unsupported feature found by ``compile``.

    ``diagnostics`` lists the errors with JSON pointers into the system file.
    """

    def __init__(self, message: str, diagnostics: list[Diagnostic]) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics

    @property
    def codes(self) -> list[str]:
        """The stable codes of ``diagnostics``, in the same order."""
        return [d.code for d in self.diagnostics]


class EditError(RaytatouilleError, ValueError):
    """A JSON Patch (raytatouille.apply_patch, Editor) that cannot be applied (ADR 0024).

    ``code`` is a registered code: edit.patch_invalid, edit.path_not_found, edit.read_only,
    edit.test_failed, edit.invalid_value, edit.base_not_representable, or the code of the
    first error the patch would add to ``validate`` (e.g. surface.id_duplicate).
    ``location`` is the JSON pointer into the edit form (empty for the whole patch or system),
    ``op_index`` the index of the failing operation (None for errors found after applying),
    ``diagnostics`` the added errors of ``validate``. The system or Editor is unchanged.
    """

    def __init__(self, message: str, code: str, location: str, op_index: int | None,
                 diagnostics: list[Diagnostic]) -> None:
        super().__init__(message)
        self.code = code
        self.location = location
        self.op_index = op_index
        self.diagnostics = diagnostics


class ParaxialError(RaytatouilleError):
    """Paraxial data requested for a path that is not rotationally symmetric, or for an
    unknown path or wavelength index.

    ``surface`` is the id of the surface the error is about and ``location`` the JSON pointer
    of the place in the system file (a surface, or e.g. "/fields/points/3"); None if unknown.
    """

    def __init__(self, message: str, surface: str | None = None,
                 location: str | None = None) -> None:
        super().__init__(message)
        self.surface = surface
        self.location = location


class UnknownMaterial(RaytatouilleError, KeyError):
    """A material reference that the MaterialLibrary cannot resolve."""

    def __str__(self) -> str:
        # KeyError would show the repr of the message.
        return str(self.args[0]) if self.args else ""


class AnalysisError(RaytatouilleError):
    """An analysis has no defined result, e.g. the chief ray does not reach the image surface
    or no ray arrives.

    For a single lost ray (chief or zone ray): ``surface`` is the last surface it reached (for
    VIGNETTED, ABSORBED, TIR, EVENT_IMPOSSIBLE and EVANESCENT the surface where it stopped;
    for MISSED and NO_CONVERGENCE the surface before the one it did not reach), ``location``
    that surface's JSON pointer, ``ray_status`` its RayStatus, ``field`` the field index (None for a field
    given as value) and ``wavelength`` the wavelength index. All None otherwise.
    """

    def __init__(self, message: str, surface: str | None = None, location: str | None = None,
                 ray_status: RayStatus | None = None, field: int | None = None,
                 wavelength: int | None = None) -> None:
        super().__init__(message)
        self.surface = surface
        self.location = location
        self.ray_status = ray_status
        self.field = field
        self.wavelength = wavelength


class NoStopError(ParaxialError, AnalysisError, ValueError):
    """The path has no stop, but aiming, the pupils, Seidel sums or OPD need one (ADR 0022).

    One class for every analysis; it is also a ParaxialError, an AnalysisError and a
    ValueError, so existing handlers keep catching it. ``path_name`` names the path and
    ``location`` is its JSON pointer, e.g. "/paths/0".
    """

    def __init__(self, message: str, path_name: str, location: str) -> None:
        Exception.__init__(self, message)
        self.path_name = path_name
        self.surface = None
        self.location = location
        self.ray_status = None
        self.field = None
        self.wavelength = None


class RaytatouilleWarning(UserWarning):
    """A warning of the library with its stable ``code`` (docs/diagnostics.md) and the JSON
    pointer ``location`` into the system file. compile() issues one per warning in
    CompiledSystem.diagnostics; collect them with warnings.catch_warnings(record=True) or read
    the diagnostics directly.
    """

    def __init__(self, message: str, code: str, location: str) -> None:
        super().__init__(message)
        self.message = message
        self.code = code
        self.location = location

    def __str__(self) -> str:
        return f"warning [{self.code}] {self.location}: {self.message}"


class CoatingCatalogError(RaytatouilleError, ValueError):
    """Malformed coating catalogue (ADR 0019); ``file`` and ``pointer`` (JSON pointer to the
    offending value, empty for the whole file)."""

    def __init__(self, message: str, file: str, pointer: str) -> None:
        super().__init__(message)
        self.file = file
        self.pointer = pointer


class AgfError(RaytatouilleError, ValueError):
    """Malformed AGF glass catalogue; ``file`` and ``line`` (1-based, 0 = whole file)."""

    def __init__(self, message: str, file: str, line: int) -> None:
        super().__init__(message)
        self.file = file
        self.line = line


class Cancelled(RaytatouilleError):
    """A run was cancelled through its CancelToken (#83); the run's outputs (e.g. a partly
    traced RayBatch) are undefined."""
