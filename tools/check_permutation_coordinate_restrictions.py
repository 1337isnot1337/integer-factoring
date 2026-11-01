#!/usr/bin/env python3
"""Exact cycle-block and cost certificates for coordinate restrictions of S4."""

import argparse
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from math import isqrt
from pathlib import Path
import random

from split_permutation_seed import split_polynomial


def full_cycle(left, right):
    inverse = [0]*len(right)
    for index, value in enumerate(right):
        inverse[value] = index
    current = 0
    visited = set()
    while current not in visited:
        visited.add(current)
        current = inverse[left[current]]
    return current == 0 and len(visited) == len(left)


def graph_controls():
    labels = tuple(permutations(range(4)))
    edges = tuple((i, j) for i, j in combinations(range(24), 2)
                  if full_cycle(labels[i], labels[j]))
    assert len(edges) == 72
    degrees = [sum(vertex in edge for edge in edges) for vertex in range(24)]
    assert degrees == [6]*24
    unions = [frozenset(enumerate(labels[i])) | frozenset(enumerate(labels[j]))
              for i, j in edges]
    assert len(set(unions)) == len(edges)
    p_polynomial = 1 << 15
    q_polynomial = sum(1 << (subset << 4) for subset in range(16))
    pq_polynomial = q_polynomial << 15
    entries = 0
    for ai, yi in edges:
        orientations = ((ai, yi), (yi, ai))
        for i, j in product(range(2), repeat=2):
            aa, yy = orientations[i]
            bb, uu = orientations[j]
            assert split_polynomial(labels[aa], labels[bb], labels[uu], labels[yy]) == (
                p_polynomial if i == j else pq_polynomial)
            aa, uu = orientations[i]
            bb, yy = orientations[j]
            assert split_polynomial(labels[aa], labels[bb], labels[uu], labels[yy]) == (
                p_polynomial if i == j else q_polynomial)
            entries += 2

    def check_deleted(deleted):
        deleted = set(deleted)
        retained = sum(i not in deleted and j not in deleted for i, j in edges)
        assert retained >= 72-6*len(deleted)
        return retained

    exhaustive_deletions = 0
    for count in range(4):
        for deleted in combinations(range(24), count):
            check_deleted(deleted)
            exhaustive_deletions += 1
    rng = random.Random(424301)
    sampled_deletions = 0
    sharp_cases = []
    even = [i for i, label in enumerate(labels)
            if sum(label[a] > label[b] for a in range(4)
                   for b in range(a+1, 4)) % 2 == 0]
    assert len(even) == 12
    for d in range(15, 25):
        retained = check_deleted(even[:24-d])
        assert retained == 6*d-72
        sharp_cases.append({"labels": d, "retained_edges": retained})
        for _ in range(128):
            check_deleted(rng.sample(range(24), 24-d))
            sampled_deletions += 1
    return {"vertices": len(labels), "edges": len(edges),
            "degree_sequence": degrees,
            "distinct_union_blocks": len(set(unions)),
            "polynomial_block_entries": entries,
            "exhaustive_deletion_sets_through_size_three": exhaustive_deletions,
            "seeded_larger_deletion_sets": sampled_deletions,
            "sharp_degree_bound_controls": sharp_cases}


def cost_certificates():
    scale = 10**7
    assert 577**2-2*408**2 == 1

    def lower_root(value):
        numerator = isqrt(value*scale*scale)
        assert numerator*numerator <= value*scale*scale < (numerator+1)**2
        return Fraction(numerator, scale)

    records = []
    for d in range(15, 25):
        n, m = d*d, 6*d-72
        alternatives = [Fraction(n),
                        lower_root(m*(n-1))+lower_root(n-m),
                        m+2*lower_root(n-m-1)]
        lower_bound = min(alternatives)
        assert lower_bound**408 > d**577
        records.append({"labels": d,
                        "crossing_rank": n,
                        "forced_rank_in_each_cut": m,
                        "rational_alternative_lower_bounds": [str(x) for x in alternatives],
                        "rational_lower_bound": str(lower_bound),
                        "decimal_summary": float(lower_bound),
                        "exact_power_408_exceeds_d_power_577": True})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"scope": "all common coordinate subsets of the S4 identity-kernel coefficient span",
              "graph_and_polynomial_controls": graph_controls(),
              "rank_sum_certificates": cost_certificates(),
              "small_label_dependency": "minimum-cut-rank-budget.md, universal exclusion through at least d=14",
              "linear_label_restrictions_not_excluded": True,
              "tensor_power_recombinations_not_excluded": True}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass",
                      "graph_edges": result["graph_and_polynomial_controls"]["edges"],
                      "block_entries": result["graph_and_polynomial_controls"]["polynomial_block_entries"],
                      "strict_cost_certificates": len(result["rank_sum_certificates"])}))


if __name__ == "__main__":
    main()
