#!/usr/bin/env python3
"""Uniform binary-cost certificate for every cancelling P6 boundary section."""

import argparse
import hashlib
import heapq
from itertools import product
import json
import math
from pathlib import Path

from check_modular_crossing import rank_two
from check_noncommutative_crossing_pfaffian_isotropic import PAIRS, construction


def binary_coefficients(pairs=6):
    dimension = 1 << pairs
    rows = [[0] * (dimension * dimension) for _ in range(pairs + 1)]
    columns = [[0] * (dimension * dimension) for _ in range(pairs + 1)]
    for x, y in product(range(dimension), repeat=2):
        for subset in range(dimension):
            degree = subset.bit_count()
            change = (x ^ y) & subset
            v, u = x ^ change, y ^ change
            rows[degree][dimension * x + y] ^= 1 << (dimension * v + u)
            columns[degree][dimension * x + u] ^= 1 << (dimension * v + y)
    return rows, columns


def correction_tensors(gamma, even):
    dimension = 64
    rows0, columns0 = [0] * 4096, [0] * 4096
    rows1, columns1 = [0] * 4096, [0] * 4096
    for x, y in product(range(dimension), repeat=2):
        for (i, j), block in zip(PAIRS, gamma):
            xi, xj = (x >> i) & 1, (x >> j) & 1
            for b in range(4):
                if not block[b] >> y & 1:
                    continue
                bi, bj = b & 1, (b >> 1) & 1
                u = (x & ~((1 << i) | (1 << j))) | (bi << i) | (bj << j)
                remaining = block[xi | (xj << 1)]
                while remaining:
                    bit = remaining & -remaining
                    v = bit.bit_length() - 1
                    rows0[dimension * x + y] ^= 1 << (dimension * v + u)
                    columns0[dimension * x + u] ^= 1 << (dimension * v + y)
                    remaining ^= bit
                terms = (
                    ((x & ~(1 << i)) | (bi << i), xi | (bj << 1)),
                    ((x & ~(1 << j)) | (bj << j), bi | (xj << 1)),
                )
                for u, a in terms:
                    remaining = block[a]
                    while remaining:
                        bit = remaining & -remaining
                        v = bit.bit_length() - 1
                        rows1[dimension * x + y] ^= 1 << (dimension * v + u)
                        columns1[dimension * x + u] ^= 1 << (dimension * v + y)
                        remaining ^= bit
        remaining = even[y]
        while remaining:
            bit = remaining & -remaining
            v = bit.bit_length() - 1
            rows1[dimension * x + y] ^= 1 << (dimension * v + x)
            columns1[dimension * x + x] ^= 1 << (dimension * v + y)
            remaining ^= bit
    return rows0, columns0, rows1, columns1


def sum_rows(first, second):
    return [a ^ b for a, b in zip(first, second)]


def integer_shortest_path(weights, target=1):
    size = len(weights)
    infinity = size * max(weights) + 1
    distances = [infinity] * size
    predecessors = [None] * size
    distances[0] = 0
    queue = [(0, 0)]
    edges = [(weight, mask) for mask, weight in enumerate(weights) if mask]
    while queue:
        cost, state = heapq.heappop(queue)
        if cost != distances[state]:
            continue
        if state == target:
            masks = []
            current = target
            while current:
                current, mask = predecessors[current]
                masks.append(mask)
            return cost, masks
        for weight, mask in edges:
            successor = state ^ mask
            candidate = cost + weight
            if candidate < distances[successor]:
                distances[successor] = candidate
                predecessors[successor] = state, mask
                heapq.heappush(queue, (candidate, successor))
    raise AssertionError("the crossing must be reachable")


def audit():
    _, _, _, gamma, even = construction()
    rows, columns = binary_coefficients()
    r0, c0, r1, c1 = correction_tensors(gamma, even)
    correction_rows = [[0] * 4096, r0, r1, sum_rows(r0, r1)]
    correction_columns = [[0] * 4096, c0, c1, sum_rows(c0, c1)]
    correction_ranks = [[rank_two(r), rank_two(c)]
                        for r, c in zip(correction_rows, correction_columns)]
    shifted_ranks = [[rank_two(sum_rows(r, rows[6])),
                      rank_two(sum_rows(c, columns[6]))]
                     for r, c in zip(correction_rows, correction_columns)]
    assert correction_ranks == [[0, 0], [3840, 154], [1696, 18], [4096, 154]]
    assert shifted_ranks == [[4096, 1], [1696, 154], [3840, 18], [240, 154]]

    baseline_ranks = [[0, 0] for _ in range(128)]
    working_row, working_column = [0] * 4096, [0] * 4096
    previous = 0
    for index in range(1, 128):
        gray = index ^ (index >> 1)
        changed = (gray ^ previous).bit_length() - 1
        working_row = sum_rows(working_row, rows[changed])
        working_column = sum_rows(working_column, columns[changed])
        baseline_ranks[gray] = [rank_two(working_row), rank_two(working_column)]
        previous = gray

    # Ranks are preserved under every orthogonal change of the input basis.
    # No nonzero correction combination can equal a binary baseline tensor.
    assert all(pair not in baseline_ranks for pair in correction_ranks[1:])
    lower_pairs = []
    weights = []
    old_weights = []
    for delta in range(4):
        for baseline in range(128):
            if delta == 0:
                lower = baseline_ranks[baseline]
            else:
                lower = [
                    max(1,
                        abs(baseline_ranks[baseline][cut] - correction_ranks[delta][cut]),
                        abs(baseline_ranks[baseline ^ 64][cut] - shifted_ranks[delta][cut]))
                    for cut in range(2)
                ]
            lower_pairs.append(lower)
            weights.append(math.isqrt(lower[0] * lower[1]))
            old_lower = [abs(baseline_ranks[baseline][cut] - correction_ranks[delta][cut])
                         for cut in range(2)]
            if delta:
                old_lower = [max(1, value) for value in old_lower]
            old_weights.append(math.isqrt(old_lower[0] * old_lower[1]))
    lower, witness = integer_shortest_path(weights)
    old_lower, old_witness = integer_shortest_path(old_weights)
    assert lower == 834 and old_lower == 133
    xor_sum = 0
    for mask in witness:
        xor_sum ^= mask
    assert xor_sum == 1
    assert sum(weights[mask] for mask in witness) == lower
    assert lower ** 70 > 64 ** 99
    return {
        "labels": 64,
        "scope": "Every p^6-cancelling P6 boundary section, with arbitrary binary coefficient recombinations",
        "coefficient_order": "Q0 through Q6, Z0, Z1; Q6 is C and crossing is Q0",
        "baseline_ranks": baseline_ranks,
        "correction_ranks": correction_ranks,
        "C_shifted_correction_ranks": shifted_ranks,
        "uniform_binary_independence_by_rank_pairs": True,
        "uniform_lower_rank_pairs": lower_pairs,
        "integer_piece_lower_costs": weights,
        "integer_shortest_lower_bound": lower,
        "lower_weight_shortest_path_masks": witness,
        "old_unshifted_integer_lower_bound": old_lower,
        "old_lower_weight_shortest_path_masks": old_witness,
        "finite_seed_threshold_float": 64 ** math.sqrt(2),
        "integer_exclusion_comparison": "834^70 > 64^99",
        "all_cancelling_P6_boundary_sections_fail_binary_seed_criterion": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit()
    result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dependency = Path(__file__).with_name("check_noncommutative_crossing_pfaffian_isotropic.py")
    result["isotropic_constructor_sha256"] = hashlib.sha256(dependency.read_bytes()).hexdigest()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(json.dumps({key: result[key] for key in (
        "labels", "integer_shortest_lower_bound", "old_unshifted_integer_lower_bound",
        "all_cancelling_P6_boundary_sections_fail_binary_seed_criterion"
    )}))


if __name__ == "__main__":
    main()
