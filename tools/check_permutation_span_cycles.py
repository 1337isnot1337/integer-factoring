#!/usr/bin/env python3
"""Verify the single-cycle obstruction for the full permutation span."""

import argparse
from itertools import product
import json
from math import factorial
from pathlib import Path
import random

from crossing_rank_sum import RankOracle
from check_split_permutation_seed import direct_pfaffian, embed_even
from fermionic_piece_polynomial import polynomial_product
from split_permutation_seed import coefficient_tensors, split_polynomial


def single_cycle_graphs(labels):
    size = len(labels[0])
    result = {}
    for ai, a in enumerate(labels):
        for yi, y in enumerate(labels):
            if any(a[i] == y[i] for i in range(size)):
                continue
            inverse_y = {value: i for i, value in enumerate(y)}
            visited = set()
            current = 0
            while current not in visited:
                visited.add(current)
                current = inverse_y[a[current]]
            if len(visited) != size or current != 0:
                continue
            graph = frozenset(enumerate(a)) | frozenset(enumerate(y))
            result.setdefault(graph, (ai, yi))
    return tuple(result.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    rng = random.Random(198537)
    block_entries = source_entries = rank_checks = 0
    cases = []
    for size in (2, 3, 4):
        labels, tensors = coefficient_tensors(size)
        d = len(labels)
        expected_blocks = factorial(size) * factorial(size - 1) // 2
        graphs = single_cycle_graphs(labels)
        assert len(graphs) == expected_blocks
        p_mask = (1 << size) - 1
        p_polynomial = 1 << p_mask
        q_polynomial = sum(1 << (subset << size) for subset in range(1 << size))
        pq_polynomial = q_polynomial << p_mask
        for graph_index, (ai, yi) in enumerate(graphs):
            orientations = ((ai, yi), (yi, ai))
            for i, j in product(range(2), repeat=2):
                aa, yy = orientations[i]
                bb, uu = orientations[j]
                row_value = split_polynomial(labels[aa], labels[bb], labels[uu], labels[yy])
                assert row_value == (p_polynomial if i == j else pq_polynomial)
                aa, uu = orientations[i]
                bb, yy = orientations[j]
                column_value = split_polynomial(labels[aa], labels[bb], labels[uu], labels[yy])
                assert column_value == (p_polynomial if i == j else q_polynomial)
                block_entries += 2
                if graph_index < 2:
                    assert embed_even(column_value, 2 * size) == direct_pfaffian(
                        labels[aa], labels[bb], labels[uu], labels[yy])
                    source_entries += 1
                    polynomial_product.cache_clear()
        oracle = RankOracle(d, coordinate_generators=tensors)
        profiles = []
        for sample in range(24):
            functional = rng.getrandbits(1 << (2 * size))
            if sample % 2 == 0:
                functional |= 1 << p_mask
            tensor = 0
            for index, coefficient_tensor in enumerate(tensors):
                if functional >> index & 1:
                    tensor ^= coefficient_tensor
            a = (functional >> p_mask) & 1
            ai, yi = graphs[0]
            diagonal = (ai * d + yi) * d * d + ai * d + yi
            assert (tensor >> diagonal) & 1 == a
            ranks = oracle(tensor)
            if a:
                assert ranks.row_rank >= expected_blocks
                assert ranks.column_rank >= expected_blocks
                rank_checks += 1
            profiles.append([a, ranks.row_rank, ranks.column_rank])
        if size >= 3:
            # The nonzero functional condition is necessary.
            ranks = oracle(tensors[0])
            assert ranks.row_rank == 1 < expected_blocks
        cases.append({
            "permutation_size": size,
            "labels": d,
            "full_single_cycle_blocks": len(graphs),
            "recombination_rank_profiles": profiles,
        })

    d, total_rank, block_rank, threshold = 24, 576, 72, 96
    maximum_distinguished_rank = (threshold**2 - 1) // block_rank
    remaining_rank_budget = total_rank - maximum_distinguished_rank
    assert maximum_distinguished_rank == 127
    assert remaining_rank_budget == 449
    assert 4 * (remaining_rank_budget - 1) > (threshold - block_rank)**2
    assert 10**2 > 2 * 7**2
    assert threshold**7 > d**10
    assert factorial(5) > 4 * 5**2
    for size in range(5, 101):
        assert size**2 >= size + 1
        assert factorial(size) > 4 * size**2

    result = {
        "single_cycle_block_entries_checked": block_entries,
        "independent_Pfaffian_entries_checked": source_entries,
        "nonzero_functional_rank_bounds_checked": rank_checks,
        "cases": cases,
        "twenty_four_label_rank_sum_lower_bound": threshold,
        "twenty_four_label_residual_rank_budget": remaining_rank_budget,
        "strict_threshold_comparison": "96^7 > 24^10 and 10/7 > sqrt(2)",
        "factorial_induction_base": "5! > 4*5^2",
        "scope": "entire original full permutation coefficient span; dual label restrictions excluded",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(json.dumps({key: value for key, value in result.items() if key != "cases"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
