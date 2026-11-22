#!/usr/bin/env python3
"""Relative benchmark check against a baseline release (#323).

Runs the rtt_bench of the baseline (e.g. 0.9.0) and of the current build alternately in several
rounds on the same machine, each with ``--quick --runs N`` (one thread, axis field). Compared is
the spot of the Cooke triplet: per round the minimum over the runs, then the median over the
rounds. Alternating the order (base, new / new, base) and taking minima damps the noise of shared
CI runners. The check fails (exit 1) if the current build is slower than ``--threshold`` times
the baseline (default 1.25: more than 25 % slower).

The baseline's JSON has no "system" key (before #323 rtt_bench measured only the Cooke triplet);
an entry without it counts as the Cooke triplet.

Usage:
    python tools/bench_compare.py --base OLD/rtt_bench --base-arg=--reference-dir=DIR ... \\
        --new build/ci-linux-gcc/apps/rtt-bench/rtt_bench [--rounds 7] [--runs 21] \\
        [--threshold 1.25] [--summary FILE] [--github]
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Sample:
    """The spot of the Cooke triplet (one thread, axis field) of one rtt_bench run."""

    spot_min_ms: float
    rays: int


def cooke_spot(output: str) -> Sample:
    """The sample of the Cooke triplet from the JSON line of rtt_bench (any version)."""
    data = json.loads(output.strip().splitlines()[-1])
    for r in data["results"]:
        if r.get("system", "cooke_triplet") == "cooke_triplet" and r["threads"] == 1 and \
                r["field"] == 0:
            if "arrived" in r and r["arrived"] != r["rays"]:
                raise ValueError(f"only {r['arrived']} of {r['rays']} rays arrived")
            return Sample(float(r["spot_min_ms"]), int(r["rays"]))
    raise ValueError("no result for the Cooke triplet with one thread on axis")


def run(exe: str, extra: Sequence[str], runs: int) -> Sample:
    done = subprocess.run([exe, "--quick", "--runs", str(runs), *extra], capture_output=True,
                          text=True, check=False)
    if done.returncode != 0:
        raise RuntimeError(f"{exe} failed ({done.returncode}): {done.stderr.strip()}")
    return cooke_spot(done.stdout)


@dataclass(frozen=True)
class Comparison:
    base_ms: float
    new_ms: float
    base_rounds: list[float]
    new_rounds: list[float]

    @property
    def ratio(self) -> float:
        return self.new_ms / self.base_ms


def compare(base: list[Sample], new: list[Sample]) -> Comparison:
    """Median over the rounds of the per-round minima; the ray counts must agree."""
    rays = {s.rays for s in base + new}
    if len(rays) != 1:
        raise ValueError(f"the two builds traced different ray counts: {sorted(rays)}")
    b = [s.spot_min_ms for s in base]
    n = [s.spot_min_ms for s in new]
    return Comparison(statistics.median(b), statistics.median(n), b, n)


def summary_text(c: Comparison, base_name: str, threshold: float, rays: int) -> str:
    verdict = "ok" if c.ratio <= threshold else "**slower than allowed**"
    rounds = len(c.base_rounds)
    return (
        f"### Benchmark against {base_name}: spot of the Cooke triplet ({rays} rays, 1 thread)\n\n"
        f"Median over {rounds} alternating rounds of the minimum per round; allowed: "
        f"{threshold:.2f} x the baseline.\n\n"
        "| build | spot (ms) | rounds (ms) |\n| --- | --- | --- |\n"
        f"| {base_name} | {c.base_ms:.3f} | {', '.join(f'{x:.3f}' for x in c.base_rounds)} |\n"
        f"| current | {c.new_ms:.3f} | {', '.join(f'{x:.3f}' for x in c.new_rounds)} |\n\n"
        f"Ratio current / {base_name}: {c.ratio:.3f} ({verdict})\n\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", required=True, help="rtt_bench of the baseline")
    parser.add_argument("--base-arg", action="append", default=[],
                        help="extra argument for the baseline (repeatable)")
    parser.add_argument("--base-name", default="0.9.0")
    parser.add_argument("--new", required=True, help="rtt_bench of the current build")
    parser.add_argument("--new-arg", action="append", default=[],
                        help="extra argument for the current build (repeatable)")
    parser.add_argument("--rounds", type=int, default=7)
    parser.add_argument("--runs", type=int, default=21)
    parser.add_argument("--threshold", type=float, default=1.25)
    parser.add_argument("--summary", help="append a Markdown table to this file")
    parser.add_argument("--github", action="store_true", help="GitHub error annotation")
    args = parser.parse_args(argv)
    if args.rounds < 1 or args.runs < 1 or args.threshold <= 0:
        parser.error("--rounds and --runs must be >= 1, --threshold > 0")

    base: list[Sample] = []
    new: list[Sample] = []
    for k in range(args.rounds):
        order = [("base", args.base, args.base_arg), ("new", args.new, args.new_arg)]
        if k % 2 == 1:
            order.reverse()
        for which, exe, extra in order:
            (base if which == "base" else new).append(run(exe, extra, args.runs))
    c = compare(base, new)
    text = summary_text(c, args.base_name, args.threshold, base[0].rays)
    print(text, end="")
    if args.summary:
        with open(args.summary, "a", encoding="utf-8") as f:
            f.write(text)
    if c.ratio > args.threshold:
        if args.github:
            print(f"::error title=Benchmark::spot of the Cooke triplet {c.new_ms:.2f} ms is "
                  f"{c.ratio:.2f} x {args.base_name} ({c.base_ms:.2f} ms), more than "
                  f"{args.threshold:.2f} x")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
