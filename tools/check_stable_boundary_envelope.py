#!/usr/bin/env python3
"""Check a limitation of the stabilized centered bound, not of actual ranks."""

from functools import lru_cache
import argparse
import hashlib
import json
from math import comb, log2
from pathlib import Path

from folded_gram_ranks import tensor_power_blocks
from stable_pfaffian_bounds import centered_block_profile


@lru_cache(None)
def tensor_multiplicity(power, highest):
    return next((count for m, _, count in tensor_power_blocks(power) if m == highest), 0)


def odd_path_lower(power, highest):
    assert power & 1 and highest & 1
    steps, level = (power - 1) // 2, (highest - 1) // 2
    if not 1 <= level <= steps:
        return 0
    return 2 ** (steps - level) * comb(steps - 1, level - 1)


def exact_boundary_multiplicity(pairs, length, jet):
    modes = 2 * pairs
    first = 2 * (length - 2 * jet) + 1
    second = modes + 1 - 2 * length + 2 * jet
    total = 0
    for chosen in range(0, modes + 1, 2):
        total += comb(modes, chosen) * tensor_multiplicity(chosen + 1, first) * tensor_multiplicity(
            modes - chosen + 1, second)
    return first, second, total


def verify_path_counts():
    checked = 0
    for steps in range(1, 65):
        for level in range(1, steps + 1):
            assert odd_path_lower(2 * steps + 1, 2 * level + 1) <= tensor_multiplicity(
                2 * steps + 1, 2 * level + 1)
            checked += 1
    return checked


def verify_boundary_blocks():
    checked = 0
    ancillas = ((1, 0), (1, 1), (1, 2), (1, 3), (0, 0), (0, 2))
    for pairs, length in ((2, 2), (4, 2), (8, 4), (16, 7), (32, 14)):
        for jet in range(length // 2 + 1):
            first = 2 * (length - 2 * jet) + 1
            second = 2 * pairs + 1 - 2 * length + 2 * jet
            for row, col in ancillas:
                profile, dimension = centered_block_profile(first, second, 0,
                                                            2 * length - 2 * pairs, row, col)
                assert dimension > 0
                assert {(a, b): rank for a, b, rank in profile}.get((0, jet)) == 1
                checked += 1
    return checked


def endpoint_certificates():
    first_left, first_right = 3 ** 87 * 5 ** 15, 2 ** 16 * 7 ** 28 * 11 ** 22
    second_left, second_right = 45 ** 15, 2 ** 29 * 14 ** 14
    assert first_left > first_right and second_left > second_right
    assert 29 ** 2 > 2 * 20 ** 2
    return {"first_entropy_endpoint_integer_pair": [first_left, first_right],
            "second_entropy_endpoint_integer_pair": [second_left, second_right],
            "certified_computation_exponent_lower_bound": [29, 20],
            "certified_computation_rho_lower_bound": [29, 40]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-pairs", type=int, default=128)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "lower bounds only on the stabilized six-ancilla centered upper-bound computation; no lower bound on actual tensor ranks",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fusion_path_multiplicity_checks": verify_path_counts(),
        "all_six_ancillary_boundary_block_checks": verify_boundary_blocks(),
        "exact_integer_entropy_certificates": endpoint_certificates(),
        "finite_boundary_cases": [],
    }
    for pairs in (16, 32, 64, 128, 256, 512):
        if pairs > args.maximum_pairs:
            continue
        length, jet = round(.22550832 * 2 * pairs), round(.05657121 * 2 * pairs)
        first, second, count = exact_boundary_multiplicity(pairs, length, jet)
        labels = comb(2 * pairs, length)
        record = {"complement_pairs": pairs, "length": length, "symmetric_jet": jet,
                  "first_highest_weight": first, "second_highest_weight": second,
                  "determinant_twist": jet, "exact_stable_multiplicity": count,
                  "labels": labels,
                  "rho_lower_bound_on_stabilized_computation": log2(min(count, labels ** 2)) / (2 * log2(labels))}
        report["finite_boundary_cases"].append(record)
        print(json.dumps({"pairs": pairs, "length": length, "jet": jet,
                          "rho_lower_on_computation": record["rho_lower_bound_on_stabilized_computation"]}), flush=True)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
