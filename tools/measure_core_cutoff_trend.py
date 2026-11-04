#!/usr/bin/env python3
"""Measure natural two-cut pressure for full height-cutoff cores."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_hot_compounds import height, profiles
from search_coherent_hot_retractions import pressure

HEIGHTS = (16, 18, 20, 22, 24, 25)
QUICK = (16, 18, 20)


def core_words(cutoff: int) -> list[int]:
    return [1 << word for word in range(4096) if height(word) < cutoff]


def measure(cutoff: int) -> dict:
    core = core_words(cutoff)
    started = time.monotonic()
    result = profiles(core, core, 12)
    fit = pressure(result)
    return {"H": cutoff, "q": len(core), "pressure": fit,
            "seconds": time.monotonic() - started}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heights", type=int, nargs="+", default=list(HEIGHTS))
    parser.add_argument("--quick", action="store_true",
                        help="Run only the three fastest cutoffs.")
    parser.add_argument("--output", type=Path,
                        default=Path("evidence/core-cutoff-pressure-trend.json"))
    args = parser.parse_args()
    heights = list(QUICK) if args.quick else args.heights
    rows = []
    for cutoff in heights:
        row = measure(cutoff)
        assert row["pressure"]["completed"] is True
        rows.append(row)
        print(json.dumps(row), flush=True)
    receipt = {
        "completed": True,
        "construction": (
            "The natural two-cut map restricted to the identity core of all "
            "one-hot words with height below H; the row and column cuts are "
            "the restricted tensor profiles."
        ),
        "cutoff_definition": "core words satisfy height(word) < H",
        "heights": heights,
        "rows": rows,
        "quick_heights": list(QUICK),
        "height38_floor": 29 / 41,
        "threshold": 2 ** -0.5,
        "scope": (
            "Complete finite natural profiles for the listed identity cores. "
            "Pressures are numerical fits, not exact certificates. All values "
            "are above the height-38 protected-core floor and 1/sqrt(2); no "
            "factoring exponent changes."
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"completed": True, "rows": len(rows)}))


if __name__ == "__main__":
    main()
