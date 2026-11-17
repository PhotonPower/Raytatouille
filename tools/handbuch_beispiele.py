"""Runs the examples of the user manual and compares their output with the printed one.

docs/handbuch/raytatouille-handbuch.tex holds every example as a titled ``pycode`` block
(``\\begin{pycode}[title={\\code{ex01_schnellstart.py}}]``), followed by its output in an
``ausgabe`` block. This script writes all titled blocks into a temporary directory (the helper
modules rtt_bausteine.py, glas_katalog.py and achromat.py included), runs every ``ex*.py`` there
with the given interpreter and compares stdout with the printed output.

    python tools/handbuch_beispiele.py --python /pfad/zur/venv/bin/python
    python tools/handbuch_beispiele.py --python ... --update   # write the new outputs back

The interpreter needs the raytatouille version the manual describes (and matplotlib for the
plot examples). Numbers in the last digits can differ between platforms (libm); the printed
outputs come from the build named in the header of the manual. Exit status 0 if every output
matches, 1 otherwise. Only the standard library is used.
"""

from __future__ import annotations

import argparse
import difflib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEX = ROOT / "docs" / "handbuch" / "raytatouille-handbuch.tex"
BLOCK = re.compile(
    r"\\begin\{pycode\}\[title=\{\\code\{(?P<name>[^}]+)\}\}\]\n(?P<code>.*?)\\end\{pycode\}\n"
    r"(?P<out>\\begin\{ausgabe\}(?:\[[^\]]*\])?\n(?P<text>.*?)\\end\{ausgabe\})?",
    re.S,
)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--python", default=sys.executable, help="interpreter with raytatouille")
    ap.add_argument("--update", action="store_true", help="write the new outputs into the manual")
    ap.add_argument("--keep", type=Path, help="run in this directory and keep the files (figures)")
    args = ap.parse_args(argv)

    tex = TEX.read_text(encoding="utf-8")
    blocks = list(BLOCK.finditer(tex))
    work = args.keep or Path(tempfile.mkdtemp(prefix="rtt-handbuch-"))
    work.mkdir(parents=True, exist_ok=True)
    for m in blocks:
        (work / m["name"]).write_text(m["code"], encoding="utf-8")

    env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    failed = 0  # examples whose output differs
    crashed = 0  # examples that raised or exited with a non-zero status
    new_tex = tex
    for m in blocks:
        name = m["name"]
        if not name.startswith("ex"):
            continue
        run = subprocess.run([args.python, name], cwd=work, env=env, capture_output=True,
                             text=True, encoding="utf-8")
        if run.returncode != 0:
            print(f"ERROR  {name}: exit {run.returncode}\n{run.stderr.strip()[-2000:]}")
            crashed += 1
            continue
        got = run.stdout.rstrip("\n")
        if m["out"] is None:
            print(f"ok     {name} (no printed output)")
            continue
        want = m["text"].rstrip("\n")
        if got == want:
            print(f"ok     {name}")
            continue
        failed += 1
        print(f"DIFF   {name}:")
        for line in difflib.unified_diff(want.splitlines(), got.splitlines(), "manual", "run",
                                         lineterm="", n=0):
            print("   " + line)
        if args.update:
            old_block = m["out"]
            new_block = old_block.replace(m["text"], got + "\n", 1)
            new_tex = new_tex.replace(m.group(0), m.group(0).replace(old_block, new_block, 1), 1)

    if args.update and new_tex != tex:
        TEX.write_text(new_tex, encoding="utf-8")
        print(f"{TEX.relative_to(ROOT)}: new outputs written")
    total = len([m for m in blocks if m["name"].startswith("ex")])
    print(f"{total} examples, {failed} differing, {crashed} failed")
    if total == 0:  # content guard: the pattern found nothing (changed markup)
        print(f"no titled ex*.py blocks found in {TEX.relative_to(ROOT)}")
        crashed = 1
    if args.keep is None:
        shutil.rmtree(work, ignore_errors=True)
    else:
        print(f"files kept in {work}")
    return 1 if crashed or (failed and not args.update) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
