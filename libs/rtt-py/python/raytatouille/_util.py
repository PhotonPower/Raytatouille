"""Helpers shared by the Python modules of raytatouille (private)."""

from __future__ import annotations

from typing import Union

from ._core import ChromaticPair, CompiledSystem, MaterialLibrary, System, compile

#: A System (compiled on each call) or a CompiledSystem (used as is).
SystemLike = Union[System, CompiledSystem]


def compiled(system: SystemLike, materials: MaterialLibrary | None) -> CompiledSystem:
    """CompiledSystem for an analysis: a System is compiled with ``materials`` (None: only
    VACUUM, AIR and CONST: resolve). For several analyses, compile once with rt.compile()."""
    if isinstance(system, CompiledSystem):
        if materials is not None:
            raise ValueError("materials only applies to a System; the CompiledSystem is used as is")
        return system
    return compile(system, materials)


def chromatic_pair(
    pair: ChromaticPair | tuple[int, int] | None,
) -> ChromaticPair | None:
    """Wavelength pair (first, second) as ChromaticPair; None stays None."""
    if pair is None or isinstance(pair, ChromaticPair):
        return pair
    first, second = pair
    return ChromaticPair(first, second)
