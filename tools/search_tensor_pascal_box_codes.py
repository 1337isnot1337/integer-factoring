#!/usr/bin/env python3
"""Exact finite homogeneous retraction screen for the smallest Pascal box.

Thirty cases exhaust all binary rank-one splits of one of the five repeated
occupation-and-energy blocks. Another 120 cases are a deterministic sample
with all five blocks restricted. No completeness is claimed for that sample.
Every retained case is rejected at every scale by an integer polynomial.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random

from tensor_pascal_box_crossings import (
    box_kernel, exterior_tensor, physical_profiles, small_box_code,
    three_quarter_certificate,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    seed = box_kernel(2, 2)
    tensor = exterior_tensor(seed)
    rng = random.Random(20261101)
    choices = [{index: choice} for index in range(5) for choice in range(6)]
    choices += [{index: rng.randrange(6) for index in range(5)} for _ in range(120)]
    cases = []
    for number, code in enumerate(choices):
        projected = small_box_code(tensor, code)
        profiles = physical_profiles(seed, projected)
        certificate = three_quarter_certificate(profiles["profiles"], profiles["total_heights"])
        coefficients = certificate.pop("translated_coefficients")
        certificate.update({"minimum_translated_coefficient": min(coefficients),
                            "positive_translated_coefficients": sum(value > 0 for value in coefficients),
                            "translated_coefficients_sha256": sha256(json.dumps(coefficients).encode()).hexdigest()})
        cases.append({"choices": code, "physical_letters": len(projected["alphabet"]),
                      "profiles": [sorted(Counter(profile).items()) for profile in profiles["profiles"]],
                      "total_heights": profiles["total_heights"],
                      "absolute_precisions": profiles["absolute_precisions"],
                      "three_quarter_certificate": certificate})
        if (number + 1) % 30 == 0:
            print("Completed exact all-scale certificates:", number + 1, flush=True)
    sources = [Path(__file__), Path(__file__).with_name("tensor_pascal_box_crossings.py"),
               Path(__file__).with_name("cube_symmetric_exterior_blocks.py")]
    report = {"status": "pass", "cases": cases, "single_block_cases": 30,
              "sampled_all_block_cases": 120, "sample_random_seed": 20261101,
              "strict_pressure_floor": [3, 4], "factoring_score": "1/2",
              "source_sha256": {path.name: sha256(path.read_bytes()).hexdigest() for path in sources},
              "scope": "Exact actual tensor profiles and all-scale integer certificates for 30 exhaustive single-block maps and 120 sampled all-block maps, at inherited centered gauges."}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
