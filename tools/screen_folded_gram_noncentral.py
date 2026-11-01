#!/usr/bin/env python3
"""Exact integer upper-bound screen over selected noncentral shell energies."""

import argparse
import hashlib
import json
from math import isqrt, log2, sqrt
from pathlib import Path

from folded_gram_aggregate import complete_gram_channels, shell_dimension


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", default="256:32,256:64")
    parser.add_argument("--numerators", default="1,2,3,4,5,6,8,10,12,14,16,18,20")
    parser.add_argument("--denominator", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    numerators = [int(value) for value in args.numerators.split(',')]
    if args.denominator <= 0 or any(not 0 <= value <= args.denominator for value in numerators):
        raise ValueError("energy fractions must lie between zero and one")
    report = {
        "status": "pass",
        "scope": "finite upper-bound screen, not an exclusion of lower actual ranks",
        "score": "sum of exact integer ceilings of square roots of capped rank-bound products",
        "criterion": "S < d^sqrt(2), equivalently rho < 1/sqrt(2)",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "aggregate_module_sha256": hashlib.sha256(Path(__file__).with_name(
            "folded_gram_aggregate.py").read_bytes()).hexdigest(),
        "rank_module_sha256": hashlib.sha256(Path(__file__).with_name(
            "folded_gram_ranks.py").read_bytes()).hexdigest(),
        "numerators": numerators, "denominator": args.denominator, "cases": [],
    }
    for description in args.cases.split(','):
        modulus, length = map(int, description.split(':'))
        lower, upper = length ** 2, modulus * length // 2
        tested = set()
        for numerator in numerators:
            weight = lower + (upper - lower) * numerator // args.denominator
            weight -= (weight - length) % 2
            if weight in tested:
                continue
            tested.add(weight)
            dimension = shell_dimension(modulus, length, weight)
            if dimension < 2:
                continue
            cap = dimension ** 2
            cost = 0
            profile = complete_gram_channels(modulus, length, weight)
            for _, r, c in profile:
                product = min(cap, r) * min(cap, c)
                root = isqrt(product)
                cost += root + int(root * root < product)
            row = {
                "mode_cutoff": modulus, "length": length, "operator_weight": weight,
                "partner_weight": length * modulus - weight, "labels": dimension,
                "integer_rank_sum_upper_bound": cost,
                "rho_of_upper_bound": log2(cost) / (2 * log2(dimension)),
                "log2_seed_margin": sqrt(2) * log2(dimension) - log2(cost),
                "complete_uncapped_channels": [list(piece) for piece in profile],
            }
            report["cases"].append(row)
            if args.output:
                args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        best = min((row for row in report["cases"] if row["mode_cutoff"] == modulus
                    and row["length"] == length), key=lambda row: row["rho_of_upper_bound"])
        print(json.dumps({key: best[key] for key in (
            "mode_cutoff", "length", "operator_weight", "labels",
            "integer_rank_sum_upper_bound", "rho_of_upper_bound", "log2_seed_margin")}), flush=True)


if __name__ == "__main__":
    main()
