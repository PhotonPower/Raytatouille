"""Exceptions of raytatouille.

Every error of the C++ library arrives as one of these classes (or as the built-in exception
nanobind maps a standard C++ exception to: ValueError, IndexError, ...). The message is the
C++ message. Ray problems during tracing are not exceptions but ray status flags (RayStatus).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ._core import Diagnostic

__all__ = [
    "AgfError",
    "AnalysisError",
    "CoatingCatalogError",
    "CompileError",
    "ParaxialError",
    "ParseError",
    "RaytatouilleError",
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


class ParaxialError(RaytatouilleError):
    """Paraxial data requested for a path that is not rotationally symmetric, or for an
    unknown path or wavelength index."""


class UnknownMaterial(RaytatouilleError, KeyError):
    """A material reference that the MaterialLibrary cannot resolve."""

    def __str__(self) -> str:
        # KeyError would show the repr of the message.
        return str(self.args[0]) if self.args else ""


class AnalysisError(RaytatouilleError):
    """An analysis has no defined result, e.g. the chief ray does not reach the image surface
    or no ray arrives."""


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
