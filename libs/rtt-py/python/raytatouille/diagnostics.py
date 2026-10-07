"""Registry of the stable diagnostic codes (ADR 0022).

Every Diagnostic, CompileError and RaytatouilleWarning carries a code such as
"material.unknown". Codes stay the same across versions: they are never renamed, reused or
given another meaning, only added. ``CODES`` maps each code to its severity and a one-line
summary; docs/diagnostics.md documents the same list.

Example::

    import raytatouille as rt

    for d in rt.validate(system):
        info = rt.diagnostics.CODES[d.code]
        print(d.code, info.severity, info.summary, d.location)
"""

from __future__ import annotations

from typing import NamedTuple

from . import _core
from ._core import Severity

__all__ = ["CODES", "CodeInfo"]


class CodeInfo(NamedTuple):
    """One registered code."""

    code: str
    severity: Severity
    summary: str


#: All registered codes, sorted by code.
CODES: dict[str, CodeInfo] = {
    code: CodeInfo(code, severity, summary) for code, severity, summary in _core.diagnostic_codes()
}
