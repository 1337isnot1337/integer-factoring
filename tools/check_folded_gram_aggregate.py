#!/usr/bin/env python3
"""Independent exact checks for weighted exterior-block aggregation."""

from collections import defaultdict
from itertools import combinations
import argparse
import hashlib
import json
from math import comb, isqrt, log2
from pathlib import Path

from folded_gram_aggregate import (
    complete_gram_channels, exterior_blocks, gaussian_distributions,
    shell_dimension, unweighted_complete_gram_channels, weighted_exterior_blocks,
)
from folded_gram_ranks import (
    complete_gram_channels as reference_weighted,
    residue_occupancy_counts, tensor_power_blocks,
    unweighted_complete_gram_channels as reference_unweighted,
)


def verify_gaussian():
    checked = 0
    for modes in range(1, 11):
        distributions = gaussian_distributions(modes)
        for length in range(modes + 1):
            counts = defaultdict(int)
            for subset in combinations(range(1, 2 * modes, 2), length):
                counts[(sum(subset) - length ** 2) // 2] += 1
            assert list(distributions[length]) == [
                counts[degree] for degree in range(length * (modes - length) + 1)]
            checked += 1
    return checked


def verify_unweighted_blocks():
    checked = 0
    for modes in range(1, 25):
        for k in range(min(modes, 8) + 1):
            expected = defaultdict(int)
            for i in range(k + 1):
                p = k - i
                if i + 2 * p > modes:
                    continue
                count = comb(modes, i) * comb(modes - i, 2 * p)
                for highest, twist, multiplicity in tensor_power_blocks(2 * p):
                    expected[highest, twist + i] += count * multiplicity
            actual = {(highest, twist): multiplicity
                      for highest, twist, multiplicity in exterior_blocks(modes, k)}
            assert actual == dict(expected)
            checked += 1
    return checked


def verify_weighted_blocks():
    checked = 0
    for modulus, length in ((4, 2), (8, 3), (16, 4)):
        expected = defaultdict(lambda: defaultdict(int))
        for i, p, energy, count in residue_occupancy_counts(modulus, length):
            k = i + p
            for highest, twist, multiplicity in tensor_power_blocks(2 * p):
                expected[k, energy][highest, twist + i] += count * multiplicity
        for k in range(length + 1):
            for energy in range(2 * k * (modulus - k) + 1):
                actual = {(highest, twist): multiplicity
                          for highest, twist, multiplicity in weighted_exterior_blocks(
                              modulus // 2, k, energy)}
                assert actual == dict(expected[k, energy]), (modulus, k, energy)
                checked += 1
    return checked


def score_case(modulus, length, weight):
    dimension = shell_dimension(modulus, length, weight)
    profile = complete_gram_channels(modulus, length, weight)
    cap = dimension ** 2
    total = 0
    for _, r, c in profile:
        product = min(cap, r) * min(cap, c)
        root = isqrt(product)
        total += root + int(root * root < product)
    return {
        "mode_cutoff": modulus, "length": length, "operator_weight": weight,
        "labels": dimension, "complete_Gram_channels": [list(row) for row in profile],
        "integer_rank_sum_upper_bound": total,
        "rho_of_upper_bound": log2(total) / (2 * log2(dimension)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "exact weighted aggregation of proved complete Gram channel bounds; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "aggregate_module_sha256": hashlib.sha256(Path(__file__).with_name(
            "folded_gram_aggregate.py").read_bytes()).hexdigest(),
        "rank_module_sha256": hashlib.sha256(Path(__file__).with_name(
            "folded_gram_ranks.py").read_bytes()).hexdigest(),
        "explicit_subset_distribution_checks": verify_gaussian(),
        "unweighted_multiplicity_checks": verify_unweighted_blocks(),
        "weighted_multiplicity_checks": verify_weighted_blocks(),
    }
    unweighted = []
    for modes, length in ((4, 2), (6, 3), (8, 4), (16, 8),
                          (32, 8), (32, 16), (64, 16), (64, 32)):
        assert unweighted_complete_gram_channels(modes, length) == reference_unweighted(modes, length)
        unweighted.append([modes, length])
    report["complete_unweighted_profile_comparisons"] = unweighted
    source = Path(__file__).resolve().parents[1] / "evidence/complement-shell-ranks.json"
    previous = json.loads(source.read_text())["cases"]
    weighted = []
    for case in previous:
        args_ = case["mode_cutoff"], case["length"], case["operator_weight"]
        assert shell_dimension(*args_) == case["labels"]
        profile = complete_gram_channels(*args_)
        assert profile == reference_weighted(*args_)
        bounds = {shift: (r, c) for shift, r, c in profile}
        for shift, r, c in case["nonzero_profile"]:
            upper_r, upper_c = bounds[shift // case["mode_cutoff"]]
            assert r <= upper_r and c <= upper_c
        weighted.append(list(args_))
    # A far-noncentral case checks delta independently of central symmetry.
    for args_ in ((16, 3, 15), (16, 4, 24), (32, 4, 32), (32, 8, 128), (64, 16, 512)):
        assert complete_gram_channels(*args_) == reference_weighted(*args_)
        weighted.append(list(args_))
    report["complete_weighted_profile_comparisons"] = weighted
    report["large_exact_examples"] = [
        score_case(*args_) for args_ in ((128, 32, 2048), (256, 32, 4096), (256, 64, 8192))]
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: report[key] for key in (
        "explicit_subset_distribution_checks", "unweighted_multiplicity_checks",
        "weighted_multiplicity_checks", "complete_unweighted_profile_comparisons",
        "complete_weighted_profile_comparisons")}), flush=True)
    print(json.dumps([{"M": row["mode_cutoff"], "L": row["length"], "j": row["operator_weight"],
                       "d": row["labels"], "rho": row["rho_of_upper_bound"]}
                      for row in report["large_exact_examples"]]), flush=True)


if __name__ == "__main__":
    main()
