"""The JSON schemas of Raytatouille's file formats, shipped with the package (#386).

The files are those under ``schema/`` in the repository; the build copies them into the package
(``raytatouille/schema/``). Tools that edit system files use them for choices, descriptions and
limits (e.g. the node workflow of Raytatouille Fusion).
"""

from __future__ import annotations

import json
from importlib import resources
from typing import Any

_FILES = {
    "system": "raytatouille.schema.json",
    "result": "raytatouille-result.schema.json",
    "coatings": "raytatouille-coatings.schema.json",
}


def schema(name: str = "system") -> dict[str, Any]:
    """The JSON schema (draft 2020-12) of a file format as a dict.

    ``name``: ``"system"`` (system file ``.rtt.json``), ``"result"`` (result file) or
    ``"coatings"`` (coating catalogue). Raises ``ValueError`` for another name. A new dict on
    every call; changing it does not change the package's files.
    """
    file = _FILES.get(name)
    if file is None:
        raise ValueError(f"unknown schema {name!r}; one of: {', '.join(_FILES)}")
    text = resources.files(__package__).joinpath("schema", file).read_text(encoding="utf-8")
    loaded: dict[str, Any] = json.loads(text)
    return loaded


__all__ = ["schema"]
