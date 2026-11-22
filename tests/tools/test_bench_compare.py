"""Tests of tools/bench_compare.py (#323): parsing both JSON forms, the comparison and the exit
code, with fake rtt_bench programs."""

from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
import bench_compare  # noqa: E402

OLD = ('{"benchmark": "cooke_triplet_spot", "runs": 3, "cpus": 4, "target_ms": 10.0, '
       '"warn_factor": 3.0, "results": [{"threads": 1, "field": 0, "rays": 1027, '
       '"aim_ms": 1.0, "trace_ms": 1.0, "spot_ms": 5.0, "spot_min_ms": 4.0}, {"threads": 4, '
       '"field": 0, "rays": 1027, "aim_ms": 1.0, "trace_ms": 1.0, "spot_ms": 2.0, '
       '"spot_min_ms": 1.0}]}')


def new_json(spot: float, arrived: int = 1027) -> str:
    results = [
        {"system": "cooke_triplet_asphere", "threads": 1, "field": 0, "rays": 1027,
         "arrived": 1027, "spot_ms": 9.0, "spot_min_ms": 9.0},
        {"system": "cooke_triplet", "threads": 1, "field": 0, "rays": 1027,
         "arrived": arrived, "spot_ms": spot + 1.0, "spot_min_ms": spot},
    ]
    return json.dumps({"benchmark": "cooke_triplet_spot", "results": results})


def test_both_json_forms() -> None:
    assert bench_compare.cooke_spot(OLD) == bench_compare.Sample(4.0, 1027)
    assert bench_compare.cooke_spot(new_json(4.5)) == bench_compare.Sample(4.5, 1027)


def test_lost_rays_and_missing_result_are_errors() -> None:
    with pytest.raises(ValueError, match="1000 of 1027"):
        bench_compare.cooke_spot(new_json(4.5, arrived=1000))
    with pytest.raises(ValueError, match="no result"):
        bench_compare.cooke_spot('{"results": []}')


def test_median_of_the_round_minima() -> None:
    s = bench_compare.Sample
    c = bench_compare.compare([s(4.0, 1027), s(9.0, 1027), s(4.2, 1027)],
                              [s(5.0, 1027), s(4.8, 1027), s(20.0, 1027)])
    assert (c.base_ms, c.new_ms) == (4.2, 5.0)
    assert c.ratio == pytest.approx(5.0 / 4.2)
    with pytest.raises(ValueError, match="ray counts"):
        bench_compare.compare([s(4.0, 1027)], [s(4.0, 1000)])


def fake(tmp_path: Path, name: str, output: str) -> str:
    script = tmp_path / name
    script.write_text(f"#!{sys.executable}\nprint({output!r})\n", encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(script)


@pytest.mark.skipif(sys.platform == "win32", reason="fake programs are scripts with a shebang")
@pytest.mark.parametrize(("spot", "code"), [(4.9, 0), (5.0, 0), (5.1, 1)])
def test_exit_code_and_summary(tmp_path: Path, spot: float, code: int,
                               capsys: pytest.CaptureFixture[str]) -> None:
    summary = tmp_path / "summary.md"
    argv = ["--base", fake(tmp_path, "old", OLD), "--new", fake(tmp_path, "new", new_json(spot)),
            "--rounds", "3", "--runs", "1", "--summary", str(summary), "--github"]
    assert bench_compare.main(argv) == code  # threshold 1.25 x 4.0 = 5.0
    text = summary.read_text(encoding="utf-8")
    assert "Benchmark against 0.9.0" in text and "| current |" in text
    assert ("::error" in capsys.readouterr().out) == (code == 1)
