"""Show the progress of a long run and cancel it, as a GUI would (cancel= and progress=, #83).

First computes a polychromatic spot diagram of the singlet with a progress callback that
prints every report (stages "aim" and "trace" per wavelength, at most every 50 ms, each stage
ending with done == total). Then starts the same spot in a worker thread and cancels it from
the main thread, as a "Cancel" button would: the run raises raytatouille.errors.Cancelled.

A GUI callback would only update a progress bar. Here the first report of the second run waits
until the main thread has cancelled, so that the example cancels reproducibly.

Run from the repository root after building the package (pip install .):

    python examples/python/cancel_progress.py
"""

from __future__ import annotations

import sys
import threading
from pathlib import Path

import raytatouille as rt

REPO_ROOT = Path(__file__).resolve().parents[2]
FILE = REPO_ROOT / "tests" / "reference" / "m1" / "singlet_const.rtt.json"


def main(file: Path = FILE) -> int:
    compiled = rt.compile(rt.load(file))

    finals: list[tuple[int, int, str]] = []

    def show(done: int, total: int, stage: str) -> None:
        print(f"  {stage}: {done} of {total}")
        if done == total:
            finals.append((done, total, stage))

    spot = rt.analysis.spot(compiled, field=2, rays="hexapolar:30", progress=show)
    print(f"spot complete: {spot.rays_arrived} of {spot.rays_launched} rays, "
          f"rms {spot.stats.rms_centroid * 1e3:.2f} um")

    token = rt.CancelToken()
    reported = threading.Event()
    cancelled = threading.Event()
    outcome: list[str] = []

    def wait_for_cancel(done: int, total: int, stage: str) -> None:
        if not reported.is_set():
            reported.set()
            cancelled.wait(timeout=30.0)

    def run() -> None:
        try:
            rt.analysis.spot(compiled, field=2, rays="hexapolar:30", cancel=token,
                             progress=wait_for_cancel)
            outcome.append("finished")
        except rt.errors.Cancelled:
            outcome.append("cancelled")

    worker = threading.Thread(target=run)
    worker.start()
    reported.wait(timeout=30.0)
    token.cancel()  # the "Cancel" button
    cancelled.set()
    worker.join()
    print(f"second run: {outcome[0]}")
    # Every stage of the complete run ended with done == total: aim and trace per wavelength.
    complete = {stage for _, _, stage in finals} == {"aim", "trace"}
    return 0 if complete and outcome == ["cancelled"] else 1


if __name__ == "__main__":
    sys.exit(main())
