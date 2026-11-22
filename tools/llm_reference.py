"""Generates the API index of raytatouille/llm_reference.md from the installed package (#408).

The reference for language models ends with one line per public function (and method of the
main classes), written between the markers ``<!-- BEGIN GENERATED`` and ``<!-- END GENERATED``.
This script imports ``raytatouille`` with the interpreter it runs in and rewrites that part:

    python tools/llm_reference.py --write    # update libs/rtt-py/python/raytatouille/llm_reference.md
    python tools/llm_reference.py            # print the index

Run it with the build of the version the reference describes (e.g. PYTHONPATH=build/<preset>/
python). libs/rtt-py/tests/test_llm_reference.py checks that the index names exactly the public
functions; signatures are informative and may differ in typing details between Python versions.
"""

from __future__ import annotations

import enum
import inspect
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "libs" / "rtt-py" / "python" / "raytatouille" / "llm_reference.md"
BEGIN = "<!-- BEGIN GENERATED: python tools/llm_reference.py --write -->"
END = "<!-- END GENERATED -->"

# (heading, module attribute path below raytatouille or "", class name or "") in output order.
SECTIONS: list[tuple[str, str, str]] = [
    ("raytatouille", "", ""),
    ("rt.paraxial", "paraxial", ""),
    ("rt.trace", "trace", ""),
    ("rt.analysis", "analysis", ""),
    ("rt.polar", "polar", ""),
    ("rt.optim", "optim", ""),
    ("rt.tolerance", "tolerance", ""),
    ("rt.layout", "layout", ""),
    ("rt.design", "design", ""),
    ("rt.reports", "reports", ""),
    ("rt.results", "results", ""),
    ("rt.plot (matplotlib)", "plot", ""),
    ("rt.diagnostics", "diagnostics", ""),
    ("rt.System methods", "", "System"),
    ("rt.Editor methods", "", "Editor"),
    ("rt.MaterialLibrary methods", "", "MaterialLibrary"),
]


# Keywords that many functions share; the index writes them as one "…" (explained in the text).
COMMON = {"materials", "coatings", "threads", "cancel", "progress", "lost_warning_fraction",
          "undersampled_waves", "mapping_warning", "mapping_limit", "max_recorded_rays"}


def _clean(t: str) -> str:
    t = re.sub(r"raytatouille\._core\.|raytatouille\.[a-z_]+\.|collections\.abc\.|typing\.", "", t)
    t = re.sub(r"numpy\.ndarray\[[^\]]*\]", "ndarray", t)
    t = t.replace("SinglePupilPoint | HexapolarPupil | GridPupil | FanXPupil | FanYPupil | "
                  "RandomPupil | GaussPupil | PointListPupil", "PupilSampling")
    return t.replace("NoneType", "None")


def _split_params(inner: str) -> list[str]:
    out, depth, cur = [], 0, ""
    for ch in inner:
        depth += ch in "([{"
        depth -= ch in ")]}"
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def _ann(a: object) -> str:
    return a if isinstance(a, str) else inspect.formatannotation(a)


def _params(obj: object) -> tuple[list[str], str]:
    """Parameters as text and the return annotation (" -> T" or "")."""
    try:
        sig = inspect.signature(obj)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        sig = None
    if sig is not None and not any(p.kind == p.VAR_POSITIONAL and p.name == "args"
                                   for p in sig.parameters.values()):
        ps: list[str] = []
        star = False
        for p in sig.parameters.values():
            if p.name == "self":
                continue
            if p.kind == p.KEYWORD_ONLY and not star:
                ps.append("*")
                star = True
            if p.kind == p.VAR_POSITIONAL:
                ps.append("*" + p.name)
                star = True
                continue
            if p.kind == p.VAR_KEYWORD:
                ps.append("**" + p.name)
                continue
            ann = "" if p.annotation is p.empty else ": " + _ann(p.annotation)
            default = "" if p.default is p.empty else " = " + repr(p.default)
            ps.append(p.name + ann + default)
        ret = ("" if sig.return_annotation is sig.empty
               else " -> " + _ann(sig.return_annotation))
        return ps, ret
    # nanobind functions: the first docstring line is "name(params) -> ret"
    doc = (getattr(obj, "__doc__", None) or "").strip().splitlines()
    m = re.match(r"\w+\((.*)\)(\s*->\s*(.*))?$", doc[0]) if doc else None
    if m is None:
        return ["..."], ""
    ps = [x for x in _split_params(m.group(1)) if x != "self" and not x.startswith("self:")]
    return ps, (" -> " + m.group(3) if m.group(3) else "")


def _signature(name: str, obj: object) -> str:
    ps, ret = _params(obj)
    kept: list[str] = []
    keyword_only = dropped = False
    for p in ps:
        keyword_only = keyword_only or p == "*" or p.startswith("*")
        if keyword_only and p.split(":")[0].split("=")[0].strip() in COMMON:
            dropped = True
            continue
        kept.append(p)
    if dropped:
        kept.append("…")
    if kept and kept[-1] == "*":
        kept.pop()
    return _clean(f"{name}({', '.join(kept)}){ret}")


def public_functions(module: str, cls: str) -> list[tuple[str, object]]:
    import raytatouille as rt

    if cls:
        owner = getattr(rt, cls)
        names = [n for n in vars(owner) if not n.startswith("_") and callable(getattr(owner, n))]
        return [(n, getattr(owner, n)) for n in names if n not in ("to_dict", "to_json")]
    mod = getattr(rt, module) if module else rt
    names = getattr(mod, "__all__", [n for n in dir(mod) if not n.startswith("_")])
    out = []
    for n in names:
        obj = getattr(mod, n, None)
        if obj is None or isinstance(obj, type) or not callable(obj) or n[0].isupper():
            continue
        out.append((n, obj))
    return out


def classes(module: str) -> list[str]:
    import raytatouille as rt

    mod = getattr(rt, module) if module else rt
    names = getattr(mod, "__all__", [])
    return [n for n in names if isinstance(getattr(mod, n, None), type)
            and not issubclass(getattr(mod, n), (BaseException, enum.Enum))]


def index() -> str:
    lines: list[str] = []
    for heading, module, cls in SECTIONS:
        funcs = public_functions(module, cls)
        if not funcs:
            continue
        lines.append(f"### {heading}")
        lines.append("```text")
        lines.extend(_signature(n, o) for n, o in funcs)
        lines.append("```")
        if not cls:
            names = classes(module)
            if names:
                lines.append("Classes: " + ", ".join(f"`{n}`" for n in names) + ".")
        lines.append("")
    import raytatouille as rt

    lines.append("Model classes (`rt.model`, read-only views of `System`): "
                 + ", ".join(f"`{n}`" for n in rt.model.__all__ if n[0].isupper()) + ".")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    text = index()
    if "--write" not in argv:
        print(text)
        return 0
    ref = REFERENCE.read_text(encoding="utf-8")
    a = ref.index(BEGIN) + len(BEGIN)
    b = ref.index(END)
    REFERENCE.write_text(ref[:a] + "\n" + text + "\n" + ref[b:], encoding="utf-8")
    print(f"{REFERENCE.relative_to(ROOT)}: index written ({text.count(chr(10)) + 1} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
