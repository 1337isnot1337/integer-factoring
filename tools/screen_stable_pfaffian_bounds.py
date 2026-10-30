#!/usr/bin/env python3
"""Reproduce complete centered-jet upper costs near the dense-grade optimum."""

import argparse
import hashlib
import json
from pathlib import Path

from stable_pfaffian_bounds import centered_score_case


RATIOS = ((7, 20), (2, 5), (17, 40), (9, 20), (19, 40),
          (1, 2), (11, 20), (3, 5), (13, 20))


def lengths(pairs):
    return sorted({max(2, (pairs * numerator + denominator // 2) // denominator)
                   for numerator, denominator in RATIOS})


def compact_case(report):
    result = dict(report)
    profile = result.pop("complete_centered_jet_bounds")
    encoded = json.dumps(profile, separators=(",", ":")).encode()
    result["complete_coefficient_count"] = len(profile)
    result["complete_profile_sha256"] = hashlib.sha256(encoded).hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, nargs="+", default=(16, 32, 64, 128))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "running",
        "scope": "complete six-ancilla centered-jet upper-cost grid; no improved factoring exponent",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "module_sha256": hashlib.sha256(Path(__file__).with_name("stable_pfaffian_bounds.py").read_bytes()).hexdigest(),
        "length_ratios_to_complement_pairs": [list(ratio) for ratio in RATIOS],
        "cases": [], "minimum_by_pair_count": [],
    }
    for pairs in args.pairs:
        local = []
        for length in lengths(pairs):
            case = compact_case(centered_score_case(pairs, length))
            local.append(case)
            report["cases"].append(case)
            print(json.dumps({"pairs": pairs, "length": length,
                              "rho": case["rho_of_upper_bound"]}), flush=True)
            if args.output:
                args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        report["minimum_by_pair_count"].append(min(local, key=lambda row: row["rho_of_upper_bound"]))
    report["status"] = "pass"
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["minimum_by_pair_count"]), flush=True)


if __name__ == "__main__":
    main()
