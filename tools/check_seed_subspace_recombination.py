#!/usr/bin/env python3
"""Check dual label restrictions and exact split-seed rank-sum searches."""

import argparse
from itertools import permutations, product
import json
from pathlib import Path
import random

from crossing_rank_sum import (
    RankOracle, certifies_lower_power, certifies_upper_power,
    crossing_tensor, dual_pair, independent_generators, optimize_span,
    restrict_tensor, xor_shortest_path,
)
from check_permutation_homogeneous_labels import compose
from split_permutation_seed import coefficient_tensors


def random_dual_pair(old_labels, labels, rng):
    first = [1 << i for i in range(old_labels)]
    second = first[:]
    for _ in range(6 * old_labels):
        i, j = rng.sample(range(old_labels), 2)
        first[j] ^= first[i]
        second[i] ^= second[j]
    first, second = tuple(first[:labels]), tuple(second[:labels])
    assert dual_pair(first, second)
    return first, second


def pair_representatives(tensors):
    """All S4 x S4 orbits of two subset-degree indices, with duplicates removed."""
    seen = set()
    for s, t in product(range(5), repeat=2):
        first = ((1 << s) - 1) + (((1 << t) - 1) << 4)
        for si, so, ti, to in product(range(s + 1), range(5 - s),
                                      range(t + 1), range(5 - t)):
            second_s = ((1 << si) - 1) + (((1 << so) - 1) << s)
            second_t = ((1 << ti) - 1) + (((1 << to) - 1) << t)
            second = second_s + (second_t << 4)
            key = tuple(sorted((tensors[first], tensors[second])))
            if key not in seen:
                seen.add(key)
                yield first, second


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extended", action="store_true",
                        help="Also search the uniform span plus every pair of coefficient pieces.")
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    rng = random.Random(739812)

    shortest_path_checks = 0
    for _ in range(24):
        costs = [0] + [rng.randrange(1, 30) for _ in range(7)]
        brute = [sum(costs) + 1] * 8
        for subset in range(1 << 7):
            target = cost = 0
            for vector in range(1, 8):
                if subset >> (vector - 1) & 1:
                    target ^= vector
                    cost += costs[vector]
            brute[target] = min(brute[target], cost)
        for target in range(8):
            actual, witness = xor_shortest_path(costs, target)
            assert actual == brute[target]
            combined = 0
            for vector in witness:
                combined ^= vector
            assert combined == target and sum(costs[v] for v in witness) == actual
            shortest_path_checks += 1

    restriction_checks = 0
    for old_labels in range(2, 8):
        for labels in range(2, old_labels + 1):
            for _ in range(4):
                a, b = random_dual_pair(old_labels, labels, rng)
                u, y = random_dual_pair(old_labels, labels, rng)
                assert restrict_tensor(crossing_tensor(old_labels), old_labels, a, b, u, y) == (
                    crossing_tensor(labels))
                restriction_checks += 1

    # A tiny complete coordinate decomposition checks both the optimizer
    # and the strict direction of the integer certificate inequalities.
    d = 2
    target = crossing_tensor(d)
    coordinate_pieces = [1 << i for i in range(d**4) if target >> i & 1]
    tiny = optimize_span(coordinate_pieces, RankOracle(d))
    assert tiny["lower_scaled"] == tiny["upper_scaled"] == 4 * tiny["scale"]
    assert not certifies_upper_power(tiny, 2, 1)
    assert certifies_lower_power(tiny, 2, 1)
    assert certifies_upper_power(tiny, 3, 1)

    labels, tensors = coefficient_tensors(4)
    old_labels = len(labels)
    curve = [0] * 9
    for subset, tensor in enumerate(tensors):
        curve[subset.bit_count()] ^= tensor
    full_generators, _ = independent_generators(tensors)
    assert len(full_generators) == 141
    oracle = RankOracle(old_labels, coordinate_generators=tensors)
    uniform = optimize_span(curve, oracle)
    assert uniform["span_dimension"] == 6
    assert uniform["lower_scaled"] == uniform["upper_scaled"] == 325 * oracle.scale
    assert certifies_lower_power(uniform, 3, 2)

    restrictions = []
    for d in (7, 8, 10, 12, 16, 20):
        embedding = tuple(1 << i for i in range(d))
        compressed = [restrict_tensor(tensor, old_labels, embedding, embedding,
                                      embedding, embedding) for tensor in curve]
        result = optimize_span(compressed, RankOracle(d))
        baseline = (d + 1) * (d + 2) // 2
        assert result["lower_scaled"] == result["upper_scaled"] == baseline * result["scale"]
        assert certifies_lower_power(result, 3, 2)
        result["restriction"] = "coordinate prefix"
        restrictions.append(result)

    index = {p: i for i, p in enumerate(labels)}
    cycle = (1, 2, 0, 3)
    subgroup = (tuple(range(4)), cycle, compose(cycle, cycle))
    for side in ("left", "right"):
        remaining = set(range(old_labels))
        embedding = []
        while remaining:
            representative = min(remaining)
            group = tuple(index[compose(h, labels[representative]) if side == "left"
                                else compose(labels[representative], h)] for h in subgroup)
            remaining.difference_update(group)
            embedding.append(sum(1 << i for i in group))
        assert dual_pair(embedding, embedding)
        compressed = [restrict_tensor(tensor, old_labels, embedding, embedding,
                                      embedding, embedding) for tensor in curve]
        result = optimize_span(compressed, RankOracle(len(embedding)))
        assert result["lower_scaled"] == result["upper_scaled"] == 45 * result["scale"]
        result["restriction"] = "C3 " + side + " coset sums"
        restrictions.append(result)

    crossing = crossing_tensor(old_labels)
    pair_results = []
    for first, second in pair_representatives(tensors):
        result = optimize_span((crossing, tensors[first], tensors[second]), oracle)
        assert result["lower_scaled"] >= 325 * oracle.scale
        result["first_piece"], result["second_piece"] = first, second
        pair_results.append(result)
    assert len(pair_results) == 766
    assert min(result["upper_scaled"] for result in pair_results) == 325 * oracle.scale

    single_extensions = []
    for s, t in product(range(5), repeat=2):
        piece = ((1 << s) - 1) + (((1 << t) - 1) << 4)
        result = optimize_span(curve + [tensors[piece]], oracle)
        assert result["lower_scaled"] == result["upper_scaled"] == 325 * oracle.scale
        result["subset_sizes"] = [s, t]
        single_extensions.append(result)

    extended = []
    if options.extended:
        for first, second in pair_representatives(tensors):
            result = optimize_span(curve + [tensors[first], tensors[second]], oracle)
            result["first_piece"], result["second_piece"] = first, second
            extended.append(result)

    result = {
        "shortest_path_exhaustive_checks": shortest_path_checks,
        "dual_restriction_checks": restriction_checks,
        "strict_certificate_negative_control": "passed",
        "full_coefficient_span_dimension": len(full_generators),
        "uniform_span": uniform,
        "label_restrictions": restrictions,
        "crossing_plus_two_piece_spans": pair_results,
        "uniform_plus_one_piece_spans": single_extensions,
        "uniform_plus_two_piece_spans": extended,
        "scope": "exact optimization only within the explicitly listed coefficient spans",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(json.dumps({
        "shortest_path_exhaustive_checks": shortest_path_checks,
        "dual_restriction_checks": restriction_checks,
        "full_coefficient_span_dimension": len(full_generators),
        "uniform_span_dimension": uniform["span_dimension"],
        "uniform_optimal_rank_sum": uniform["upper_scaled"] // oracle.scale,
        "label_restrictions": len(restrictions),
        "crossing_plus_two_piece_spans": len(pair_results),
        "uniform_plus_one_piece_spans": len(single_extensions),
        "uniform_plus_two_piece_spans": len(extended),
        "scope": result["scope"],
    }, indent=2))


if __name__ == "__main__":
    main()
