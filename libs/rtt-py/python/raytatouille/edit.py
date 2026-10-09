"""Changing a system with JSON Patch, with undo and redo (ADR 0024).

A command is an RFC 6902 patch document: a list of operations ``add``, ``remove``, ``replace``,
``move``, ``copy`` and ``test`` whose JSON pointers (RFC 6901) address the edit form of the
system (``System.to_dict()``: every value written, every Param an object). A patch is applied
completely or not at all and is one undo step. Errors raise ``EditError`` with a code, a JSON
pointer and the index of the operation; the system or Editor stays unchanged.

Example::

    ed = rt.Editor(rt.load("singlet.rtt.json"))
    ed.set(ed.system.locate_surface("L1.S1") + "/shape/base/radius", 50.0)  # keeps "variable"
    ed.insert("/wavelengths/-", {"um": 0.55})
    ed.undo()
    rt.save(ed.system, "changed.rtt.json")
    replay = rt.Editor.from_history(ed.export_history())
"""

from __future__ import annotations

import copy
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, Union

from . import _core
from ._core import Diagnostic, System, validate

__all__ = ["Editor", "Patch", "PatchResult", "apply_patch", "apply_patch_with_inverse"]

#: A patch document: a list of operations, or the same as JSON text.
Patch = Union[Sequence[Mapping[str, Any]], str]
Check = Literal["no_new_errors", "structure_only"]


@dataclass(frozen=True)
class PatchResult:
    """Result of apply_patch_with_inverse."""

    system: System  #: the changed system
    #: RFC 6902 patch that turns the changed system back into the input; apply it with
    #: check="structure_only".
    inverse: list[dict[str, Any]]
    diagnostics: list[Diagnostic]  #: all diagnostics of validate(system)


def _text(patch: Patch) -> str:
    if isinstance(patch, str):
        return patch
    if isinstance(patch, Mapping) or not isinstance(patch, Sequence):
        raise TypeError("patch must be a list of operations (or its JSON text), not "
                        f"{type(patch).__name__}; wrap a single operation in a list")
    ops: list[dict[str, Any]] = []
    for op in patch:
        if not isinstance(op, Mapping):
            raise TypeError(f"each operation of a patch must be a dict, not {type(op).__name__}")
        ops.append(dict(op))
    # allow_nan=True: NaN and infinity become invalid JSON that the C++ parser reports as
    # edit.patch_invalid, the only error path for a bad patch (ADR 0024 addendum).
    return json.dumps(ops, allow_nan=True)


def apply_patch_with_inverse(system: System, patch: Patch, *,
                             check: Check = "no_new_errors") -> PatchResult:
    """Applies an RFC 6902 patch to the edit form of ``system`` and returns the changed system,
    the inverse patch and the diagnostics (ADR 0024). ``system`` is not changed.

    ``check="no_new_errors"`` (commands) rejects a patch that adds errors of ``validate``;
    systems with errors stay editable. ``check="structure_only"`` only reads the result strictly
    and is meant for the inverse of an accepted patch (undo, redo): it skips the check of
    ``validate``.

    ``patch`` is a LIST of operations (dicts) whose values are JSON values: int, float, str,
    bool, None, list, dict. Convert NumPy scalars with ``.item()`` first. It may also be the
    JSON text of such a list.

    Raises EditError if the patch cannot be applied, TypeError if ``patch`` is not a list (e.g.
    a single operation), ValueError for an unknown ``check``.
    """
    if check not in ("no_new_errors", "structure_only"):
        raise ValueError(f"check must be 'no_new_errors' or 'structure_only', got {check!r}")
    new, inverse, diagnostics = _core.apply_patch_json(system, _text(patch),
                                                       check == "structure_only")
    return PatchResult(new, json.loads(inverse), diagnostics)


def apply_patch(system: System, patch: Patch, *, check: Check = "no_new_errors") -> System:
    """Applies an RFC 6902 patch to the edit form of ``system`` and returns the changed system
    (``system`` is not changed). A script command: replaying the patches of an Editor's history
    on its base system gives the same model. See apply_patch_with_inverse for ``check``.

    Raises EditError if the patch cannot be applied.
    """
    return apply_patch_with_inverse(system, patch, check=check).system


def _is_param(value: Any) -> bool:
    return isinstance(value, dict) and "value" in value and "variable" in value


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _major_minor(version: str) -> str:
    return version.rsplit(".", 1)[0]


class Editor:
    """A system with undo and redo (ADR 0024). Each command is one patch and one undo step.

    The Editor works on its own copy; ``system`` returns a copy as well (changing it does not
    change the Editor). A system that cannot be written (e.g. a NaN set through
    ``System.environment``) raises ValueError already in the constructor. Commands that would add errors of ``validate`` raise EditError; systems
    with errors stay editable. No RaytatouilleWarning is issued; ``diagnostics`` has them as
    data.
    """

    def __init__(self, system: System) -> None:
        self._base = system.to_json()
        self._system = System.from_json(self._base)
        self._diagnostics = validate(self._system)
        #: (patch, inverse) per step; history() is the patches of _undo.
        self._undo: list[tuple[list[dict[str, Any]], list[dict[str, Any]]]] = []
        self._redo: list[tuple[list[dict[str, Any]], list[dict[str, Any]]]] = []

    @property
    def system(self) -> System:
        """A copy of the current system. The copy goes through the canonical JSON text, so a
        value of -0.0 that equals its default (0.0) comes back as +0.0 (physically the same)."""
        return System.from_json(self._system.to_json())

    @property
    def diagnostics(self) -> list[Diagnostic]:
        """All diagnostics of validate() for the current system (errors and warnings)."""
        return list(self._diagnostics)

    @property
    def can_undo(self) -> bool:
        """True if there is a step to take back."""
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        """True if there is an undone step to apply again."""
        return bool(self._redo)

    @property
    def history(self) -> list[list[dict[str, Any]]]:
        """The patches applied now, in order (the undo stack without the redo entries).
        Replayed on the base system they give the current system."""
        return copy.deepcopy([patch for patch, _ in self._undo])

    def apply(self, patch: Patch) -> None:
        """Applies an RFC 6902 patch as one undo step; drops the redo entries. ``patch`` is a
        list of operations with JSON values (see apply_patch_with_inverse) or its JSON text.

        Raises EditError, or TypeError if ``patch`` is not a list (the Editor is unchanged)."""
        # The text goes to the strict C++ parser unchanged (duplicate keys, invalid JSON and
        # NaN are edit.patch_invalid); the history keeps what was applied, as plain JSON.
        text = _text(patch)
        r = apply_patch_with_inverse(self._system, text)
        ops: list[dict[str, Any]] = json.loads(text)
        self._system, self._diagnostics = r.system, r.diagnostics
        self._undo.append((ops, r.inverse))
        self._redo.clear()

    def set(self, pointer: str, value: Any) -> None:
        """Sets the value at ``pointer``: ``replace``, or ``add`` for a missing member of an
        object (an optional member such as a surface aperture). A number set on an unbound Param
        object goes to ``pointer + "/value"`` so that ``variable`` and the bounds stay; on a Param
        bound to a parameter row ({"param": ...}, ADR 0029) it replaces the whole object, which
        unbinds it. A missing array element is not added (that is ``insert``): it raises
        EditError edit.path_not_found.

        Raises ValueError for an invalid pointer, EditError if the patch fails (the Editor is
        unchanged)."""
        try:
            current = self._system.json_at(pointer)
        except KeyError:
            parent_pointer = pointer.rsplit("/", 1)[0]
            try:
                parent = self._system.json_at(parent_pointer)
            except KeyError:
                parent = None
            op = "add" if isinstance(parent, dict) else "replace"
            self.apply([{"op": op, "path": pointer, "value": value}])
            return
        if _is_param(current) and _is_number(value):
            pointer += "/value"
        self.apply([{"op": "replace", "path": pointer, "value": value}])

    def insert(self, pointer: str, value: Any) -> None:
        """Inserts ``value`` at an array index or at the end ("-"): ``add``."""
        self.apply([{"op": "add", "path": pointer, "value": value}])

    def remove(self, pointer: str) -> None:
        """Removes the value at ``pointer``."""
        self.apply([{"op": "remove", "path": pointer}])

    def move(self, from_: str, to: str) -> None:
        """Moves the value at ``from_`` to ``to`` (RFC 6902: remove, then add)."""
        self.apply([{"op": "move", "from": from_, "path": to}])

    def undo(self) -> None:
        """Takes back the last step. Raises IndexError if there is none; if it fails
        (EditError), state and stacks stay unchanged."""
        if not self._undo:
            raise IndexError("nothing to undo")
        _, inverse = self._undo[-1]
        r = apply_patch_with_inverse(self._system, inverse, check="structure_only")
        self._system, self._diagnostics = r.system, r.diagnostics
        self._redo.append(self._undo.pop())

    def redo(self) -> None:
        """Applies the last undone step again. Raises IndexError if there is none; if it fails
        (EditError), state and stacks stay unchanged."""
        if not self._redo:
            raise IndexError("nothing to redo")
        patch, _ = self._redo[-1]
        r = apply_patch_with_inverse(self._system, patch, check="structure_only")
        self._system, self._diagnostics = r.system, r.diagnostics
        self._undo.append(self._redo.pop())

    def export_history(self) -> dict[str, Any]:
        """{"base": canonical JSON text of the starting system, "patches": history} as plain
        JSON data. Replaying needs the same base system (ADR 0024 point 5)."""
        return {"base": self._base, "patches": self.history}

    @classmethod
    def from_history(cls, data: Mapping[str, Any]) -> Editor:
        """An Editor with the base system of ``data`` and its patches applied in order, each
        checked like a new command (a changed history cannot add errors).

        Raises ValueError if the base has another schema major.minor than this build (the
        pointers of the patches belong to that version), EditError if a patch fails."""
        try:
            base = str(data["base"])
            version = str(json.loads(base)["schema_version"])
            patches = list(data["patches"])
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError("not an exported Editor history: needs \"base\" (system JSON text "
                             f"with a schema_version) and \"patches\" ({e})") from e
        ours = System().schema_version
        if _major_minor(version) != _major_minor(ours):
            raise ValueError(f"history of schema version {version}; this build edits {ours} "
                             "(the patches address the edit form of their version)")
        editor = cls(System.from_json(base))
        for patch in patches:
            editor.apply(patch)
        return editor
