"""Changelog fragments: one file per change in changelog.d/, merged into CHANGELOG.md at release.

Every pull request adds a file changelog.d/<issue>.<kind>.md instead of editing CHANGELOG.md,
so parallel pull requests no longer conflict. <issue> is the issue number, <kind> one of
added, changed, fixed, removed (Keep a Changelog). The file holds one or more Markdown list
items ("- ..."), written as in CHANGELOG.md.

    python tools/changelog.py check                  validate all fragments
    python tools/changelog.py preview                print the [Unreleased] section with fragments
    python tools/changelog.py release 0.5.0 "Title"  move [Unreleased] and all fragments into a
                                                     new version section, delete the fragments

Only the standard library is used.
"""

from __future__ import annotations

import datetime
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FRAGMENTS = ROOT / "changelog.d"
CHANGELOG = ROOT / "CHANGELOG.md"

# Fragment kind -> section heading, in the order of the sections in CHANGELOG.md.
KINDS = {
    "added": "Hinzugefügt",
    "changed": "Geändert",
    "fixed": "Behoben",
    "removed": "Entfernt",
}
NAME = re.compile(r"^(\d+)\.(" + "|".join(KINDS) + r")(?:\.(\d+))?\.md$")


def fragments() -> list[tuple[int, int, str, Path]]:
    """All fragments as (issue, sequence, kind, path), sorted by kind order, issue, sequence."""
    result = []
    for path in sorted(FRAGMENTS.glob("*.md")):
        if path.name == "README.md":
            continue
        m = NAME.match(path.name)
        if not m:
            raise SystemExit(
                f"{path.relative_to(ROOT)}: name must be <issue>.<kind>.md or "
                f"<issue>.<kind>.<n>.md with kind in {', '.join(KINDS)}"
            )
        result.append((int(m.group(1)), int(m.group(3) or 0), m.group(2), path))
    order = list(KINDS)
    return sorted(result, key=lambda f: (order.index(f[2]), f[0], f[1]))


def text_of(path: Path) -> str:
    text = path.read_text(encoding="utf-8").strip("\n")
    if not text.startswith("- "):
        raise SystemExit(f"{path.relative_to(ROOT)}: must start with a list item '- '")
    for line in text.splitlines():
        if line and not (line.startswith("- ") or line.startswith("  ")):
            raise SystemExit(
                f"{path.relative_to(ROOT)}: every line is a list item '- ' or indented by two "
                f"spaces: {line!r}"
            )
    return text


def split_unreleased(changelog: str) -> tuple[str, dict[str, list[str]], str]:
    """Header up to and including '## [Unreleased]', its sections, and the rest."""
    start = changelog.index("## [Unreleased]")
    head_end = changelog.index("\n", start) + 1
    nxt = changelog.find("\n## [", head_end)
    body = changelog[head_end:] if nxt < 0 else changelog[head_end : nxt + 1]
    rest = "" if nxt < 0 else changelog[nxt + 1 :]
    sections: dict[str, list[str]] = {}
    current = None
    for line in body.splitlines():
        if line.startswith("### "):
            current = line[4:].strip()
            sections.setdefault(current, [])
        elif current is not None and line.strip():
            sections[current].append(line)
        elif line.strip():
            raise SystemExit(f"CHANGELOG.md: text outside a section in [Unreleased]: {line!r}")
    return changelog[:head_end], sections, rest


def merged_sections() -> dict[str, list[str]]:
    _, sections, _ = split_unreleased(CHANGELOG.read_text(encoding="utf-8"))
    for _, _, kind, path in fragments():
        sections.setdefault(KINDS[kind], []).extend(text_of(path).splitlines())
    return sections


def render(sections: dict[str, list[str]]) -> str:
    order = list(KINDS.values()) + [s for s in sections if s not in KINDS.values()]
    parts = []
    for name in order:
        lines = sections.get(name)
        if lines:
            parts.append(f"### {name}\n" + "\n".join(lines) + "\n")
    return "\n".join(parts)


def check() -> None:
    for _, _, _, path in fragments():
        text_of(path)
    split_unreleased(CHANGELOG.read_text(encoding="utf-8"))
    print(f"changelog: {len(fragments())} fragment(s) ok")


def preview() -> None:
    print("## [Unreleased]\n")
    print(render(merged_sections()))


def release(version: str, title: str) -> None:
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise SystemExit("version must be MAJOR.MINOR.PATCH")
    changelog = CHANGELOG.read_text(encoding="utf-8")
    head, _, rest = split_unreleased(changelog)
    body = render(merged_sections())
    if not body:
        raise SystemExit("nothing to release")
    date = datetime.date.today().isoformat()
    heading = f"## [{version}] – {title}" if title else f"## [{version}] – {date}"
    CHANGELOG.write_text(f"{head}\n{heading}\n\n{body}\n{rest}", encoding="utf-8", newline="\n")
    for _, _, _, path in fragments():
        path.unlink()
    print(f"CHANGELOG.md: [{version}] written, fragments removed")


def main(argv: list[str]) -> None:
    if len(argv) >= 1 and argv[0] == "check":
        check()
    elif len(argv) >= 1 and argv[0] == "preview":
        preview()
    elif len(argv) in (2, 3) and argv[0] == "release":
        release(argv[1], argv[2] if len(argv) == 3 else "")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
