"""The reference for language models shipped with the package (#408).

``raytatouille/llm_reference.md`` must stay true: its code blocks run in order as one session,
and its generated index names exactly the public functions (tools/llm_reference.py writes it).
"""

from __future__ import annotations

import re
import types
import warnings
from pathlib import Path

import pytest

import raytatouille as rt

_INDEX = re.compile(r"<!-- BEGIN GENERATED.*?-->\n(.*?)<!-- END GENERATED -->", re.S)

# Subpackages of rt.__all__ without public functions of their own in the index.
_NOT_INDEXED = {"errors", "edit", "materials", "model", "schemas"}

# The warnings the recipes cause on purpose; the text explains each of them.
_EXPECTED_WARNINGS = {"psf.aliasing", "psf.pupil_undersampled", "pupil.mapping_nonlinear",
                      "optim.parameter_at_bound"}


def _public_functions() -> set[str]:
    modules = [getattr(rt, n) for n in rt.__all__ if isinstance(getattr(rt, n), types.ModuleType)]
    names: set[str] = set()
    for module in [rt, *(m for m in modules if m.__name__.rpartition(".")[2] not in _NOT_INDEXED)]:
        for n in getattr(module, "__all__", []):
            obj = getattr(module, n, None)
            if callable(obj) and not isinstance(obj, type) and not n[0].isupper():
                names.add(n)
    for cls in (rt.System, rt.Editor, rt.MaterialLibrary):
        names.update(n for n in vars(cls) if not n.startswith("_") and callable(getattr(cls, n))
                     and n not in ("to_dict", "to_json"))
    return names


def test_llm_reference_is_the_packaged_file() -> None:
    text = rt.llm_reference()
    assert text.startswith("# raytatouille")
    assert f"Version {rt.__version__}." in text
    path = Path(rt.__file__).with_name("llm_reference.md")
    assert path.read_text(encoding="utf-8") == text
    assert "llm_reference()" in (rt.__doc__ or "").splitlines()[0]


def test_index_names_exactly_the_public_functions() -> None:
    match = _INDEX.search(rt.llm_reference())
    assert match is not None
    listed = set(re.findall(r"^([a-z_][a-z_0-9]*)\(", match.group(1), re.M))
    public = _public_functions()
    assert len(public) > 100  # content guard: the modules were found
    assert public - listed == set(), "missing in the index: run tools/llm_reference.py --write"
    assert listed - public == set(), "not public any more: run tools/llm_reference.py --write"


def test_code_blocks_run_as_one_session(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    blocks = re.findall(r"```python\n(.*?)```", rt.llm_reference(), re.S)
    assert len(blocks) >= 8
    monkeypatch.chdir(tmp_path)
    namespace: dict[str, object] = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", rt.RaytatouilleWarning)
        for k, block in enumerate(blocks):
            exec(compile(block, f"<llm_reference block {k}>", "exec"), namespace)  # noqa: S102
    codes = {getattr(w.message, "code", None) for w in caught
             if issubclass(w.category, rt.RaytatouilleWarning)}
    assert codes <= _EXPECTED_WARNINGS, f"unexpected warnings: {codes - _EXPECTED_WARNINGS}"
    assert (tmp_path / "singlet.rtt.json").exists()
