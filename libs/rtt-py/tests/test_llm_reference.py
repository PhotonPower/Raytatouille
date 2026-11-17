"""The reference for language models shipped with the package (#408).

``raytatouille/llm_reference.md`` must stay true: its code blocks run in order as one session,
and its generated index names exactly the public functions (tools/llm_reference.py writes it).
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pytest

import raytatouille as rt

_INDEX = re.compile(r"<!-- BEGIN GENERATED.*?-->\n(.*?)<!-- END GENERATED -->", re.S)

# The modules and classes of the generated index, as tools/llm_reference.py lists them.
_MODULES = ["paraxial", "trace", "analysis", "polar", "optim", "tolerance", "layout", "reports",
            "results", "plot", "diagnostics"]


def _public_functions() -> set[str]:
    names: set[str] = set()
    for module in [rt, *(getattr(rt, m) for m in _MODULES)]:
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
    assert "llm_reference()" in (rt.__doc__ or "").split("\n\n")[1]


def test_index_names_exactly_the_public_functions() -> None:
    match = _INDEX.search(rt.llm_reference())
    assert match is not None
    listed = set(re.findall(r"^([a-z_][a-z_0-9]*)\(", match.group(1), re.M))
    public = _public_functions()
    assert len(public) > 100  # content guard: the modules were found
    assert public - listed == set(), "missing in the index: run tools/llm_reference.py --write"
    assert listed - public == set(), "not public any more: run tools/llm_reference.py --write"


# The recipes run on a strongly aberrated singlet; the text explains the warnings that come
# with it (psf.*, pupil.mapping_nonlinear, optim.parameter_at_bound).
@pytest.mark.filterwarnings("ignore::raytatouille.errors.RaytatouilleWarning")
def test_code_blocks_run_as_one_session(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("numpy")
    blocks = re.findall(r"```python\n(.*?)```", rt.llm_reference(), re.S)
    assert len(blocks) >= 8
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("MPLBACKEND", os.environ.get("MPLBACKEND", "Agg"))
    namespace: dict[str, object] = {}
    for k, block in enumerate(blocks):
        exec(compile(block, f"<llm_reference block {k}>", "exec"), namespace)  # noqa: S102
    assert (tmp_path / "singlet.rtt.json").exists()
