#!/usr/bin/env python3
"""Replay the fast core-cutoff pressure points and check the recorded trend."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from measure_core_cutoff_trend import QUICK, core_words, measure

INPUT = Path("evidence/core-cutoff-pressure-trend.json")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--full", action="store_true",
                        help="Replay every recorded cutoff instead of the quick set.")
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved["completed"] is True
    assert saved["threshold"] == 2 ** -0.5
    assert abs(saved["height38_floor"] - 29 / 41) < 1e-15
    rows = saved["rows"]
    assert [row["H"] for row in rows] == sorted(row["H"] for row in rows)
    assert len({row["H"] for row in rows}) == len(rows)

    for row in rows:
        assert row["q"] == len(core_words(row["H"]))
        assert row["pressure"]["completed"] is True
        assert row["pressure"]["rho"] > saved["height38_floor"]
        assert row["pressure"]["rho"] > saved["threshold"]
    for left, right in zip(rows, rows[1:]):
        assert right["pressure"]["rho"] < left["pressure"]["rho"]

    replay = rows if args.full else [row for row in rows if row["H"] in QUICK]
    replayed = []
    for row in replay:
        fresh = measure(row["H"])
        assert fresh["q"] == row["q"]
        assert abs(fresh["pressure"]["rho"] - row["pressure"]["rho"]) < 1e-12
        assert abs(fresh["pressure"]["scale_t"] - row["pressure"]["scale_t"]) < 1e-12
        replayed.append(fresh)
    receipt = {
        "passed": True,
        "replayed_heights": [row["H"] for row in replayed],
        "recorded_heights": [row["H"] for row in rows],
        "all_above_floor": True,
        "all_above_inverse_sqrt_two": True,
        "monotone_decreasing": True,
        "input_sha256": sha256(args.input.read_bytes()).hexdigest(),
        "scope": (
            "Exact replay of the fast identity-core profiles and structural "
            "checks of the recorded full trend. Pressures are numerical; the "
            "complete exponent remains 1/2."
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"passed": True, "replayed": [row["H"] for row in replayed]}))


if __name__ == "__main__":
    main()
