#!/usr/bin/env python3
"""Exact generic-block and coefficient-span controls for typed Gaussians.

The general bound concerns this one Gaussian seed and its coefficient span.
It does not exclude proportional occupation types or other crossings.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import heapq
from itertools import product
import json
from math import comb, factorial, isqrt
from pathlib import Path

from check_schur_gaussian_reduction import polynomial_multiply, rank_binary, seed_data


ROOT = Path(__file__).resolve().parent.parent
SEED = {tuple(entry["indices"]): entry["matrix"][0][0] for entry in seed_data()["local"]}


def typed_words(counts):
    if not any(counts):
        return [()]
    result = []
    for letter, count in enumerate(counts):
        if count:
            remaining = list(counts)
            remaining[letter] -= 1
            result.extend((letter,) + word for word in typed_words(tuple(remaining)))
    return result


def tensor_entry(x, y, u, v):
    value = 1
    for letters in zip(x, y, u, v):
        value = polynomial_multiply(value, SEED[letters])
        if not value:
            break
    return value


def disjoint(first, second):
    return not any(a and b for a, b in zip(first, second))


def signature(first, second, cut):
    plus = minus = both = 0
    for i, (a, b) in enumerate(zip(first, second)):
        assert not (a and b)
        if a == 3 or b == 3:
            both |= 1 << i
        elif (a == 1 or b == 2) if cut == 0 else (a == 1 or b == 1):
            plus |= 1 << i
        elif a or b:
            minus |= 1 << i
    return plus, minus, both


def check_generic_blocks(counts):
    a, b, c, d = counts
    m, occupied = sum(counts), b + c + d
    assert m >= 2 * occupied and occupied
    words = typed_words(tuple(counts))
    dimension = len(words)
    expected_dimension = factorial(m) // (factorial(a) * factorial(b) * factorial(c) * factorial(d))
    assert dimension == expected_dimension
    pairs = [(x, y) for x in words for y in words if disjoint(x, y)]
    central = 2 * (b + c) + 4 * d
    falling = factorial(m) // factorial(m - 2 * occupied)
    expected_blocks = (falling // (factorial(b + c) ** 2 * factorial(2 * d)),
                       falling // (factorial(2 * b) * factorial(2 * c) * factorial(2 * d)))
    expected_sizes = (comb(b + c, b) ** 2 * comb(2 * d, d),
                      comb(2 * b, b) * comb(2 * c, c) * comb(2 * d, d))
    for x, y in pairs:
        assert tensor_entry(x, y, y, x) == 1 << central
    cuts = []
    for cut in range(2):
        row_signatures = [signature(x, y, cut) for x, y in pairs]
        column_signatures = [signature(v, u, cut) if cut == 0 else signature(u, v, cut)
                             for u, v in pairs]
        counts_by_signature = Counter(row_signatures)
        assert len(counts_by_signature) == expected_blocks[cut]
        assert set(counts_by_signature.values()) == {expected_sizes[cut]}
        assert Counter(column_signatures) == counts_by_signature
        rows = [0] * len(pairs)
        zero_checks = nonzero_blocks = 0
        block_entries = Counter()
        for i, first in enumerate(pairs):
            for j, second in enumerate(pairs):
                if cut == 0:
                    x, y = first
                    u, v = second
                else:
                    x, u = first
                    y, v = second
                value = tensor_entry(x, y, u, v)
                if row_signatures[i] != column_signatures[j]:
                    assert value == 0
                    zero_checks += 1
                if value >> central & 1:
                    rows[i] |= 1 << j
                    block_entries[row_signatures[i]] += 1
        nonzero_blocks = len(block_entries)
        actual_rank = rank_binary(rows)
        assert nonzero_blocks == expected_blocks[cut]
        assert actual_rank >= expected_blocks[cut]
        cuts.append({"cut": cut, "blocks": expected_blocks[cut],
                     "block_size": expected_sizes[cut], "off_signature_zero_entries": zero_checks,
                     "generic_central_rank": actual_rank})
    return {"counts": counts, "m": m, "dimension": dimension,
            "generic_pairs": len(pairs), "central_degree": central,
            "pure_monomial_entries": len(pairs), "cuts": cuts}


def coefficient_matrices(counts):
    """Build actual typed coefficients, using simultaneous site permutations."""
    words = typed_words(tuple(counts))
    dimension, m = len(words), sum(counts)
    lookup = {word: i for i, word in enumerate(words)}
    canonical = words[0]
    actions = []
    for word in words:
        permutation = tuple(i for letter in range(4) for i, value in enumerate(word) if value == letter)
        assert tuple(word[i] for i in permutation) == canonical
        actions.append([lookup[tuple(other[i] for i in permutation)] for other in words])
    cube = [0] * dimension ** 3
    for yi, ui, vi in product(range(dimension), repeat=3):
        cube[(yi * dimension + ui) * dimension + vi] = tensor_entry(canonical, words[yi], words[ui], words[vi])
    degree = max(value.bit_length() - 1 for value in cube)
    row = [[0] * dimension ** 2 for _ in range(degree + 1)]
    column = [[0] * dimension ** 2 for _ in range(degree + 1)]
    for xi in range(dimension):
        action = actions[xi]
        for yi, ui, vi in product(range(dimension), repeat=3):
            value = cube[(action[yi] * dimension + action[ui]) * dimension + action[vi]]
            while value:
                bit = value & -value
                k = bit.bit_length() - 1
                value ^= bit
                row[k][dimension * xi + yi] ^= 1 << (dimension * ui + vi)
                column[k][dimension * xi + ui] ^= 1 << (dimension * yi + vi)
    endpoint = [0] * dimension ** 2
    for matrix in row:
        endpoint = [a ^ b for a, b in zip(endpoint, matrix)]
    assert endpoint == [1 << (dimension * y + x) for x in range(dimension) for y in range(dimension)]
    return dimension, row, column


def small_one_plus_minus(m):
    dimension, row, column = coefficient_matrices((m - 2, 1, 1, 0))
    degrees = (0, 2, 4, 6, 8)
    assert len(row) == 9
    assert all(not any(row[k]) for k in (1, 3, 5, 7))
    flattened = [sum(value << (i * dimension ** 2) for i, value in enumerate(row[k])) for k in degrees]
    assert rank_binary(flattened) == 5
    products = []
    all_ranks = []
    for code in range(32):
        r = [0] * dimension ** 2
        c = [0] * dimension ** 2
        for i, degree in enumerate(degrees):
            if code >> i & 1:
                r = [a ^ b for a, b in zip(r, row[degree])]
                c = [a ^ b for a, b in zip(c, column[degree])]
        pair = rank_binary(r), rank_binary(c)
        all_ranks.append({"code": code, "row_rank": pair[0], "column_rank": pair[1]})
        if code & 4:
            products.append(pair[0] * pair[1])
    minimum_product = min(products)
    # Every real rank cost dominates the corresponding integer floor cost.
    # Exact shortest paths therefore give a rigorous lower bound on every
    # F_2 coefficient-span decomposition, even when one central piece alone
    # has too small a lower bound.
    weights = [isqrt(item["row_rank"] * item["column_rank"]) for item in all_ranks]
    distance = [None] * 32
    distance[0] = 0
    queue = [(0, 0)]
    while queue:
        cost, state = heapq.heappop(queue)
        if cost != distance[state]:
            continue
        for code in range(1, 32):
            target, candidate = state ^ code, cost + weights[code]
            if distance[target] is None or candidate < distance[target]:
                distance[target] = candidate
                heapq.heappush(queue, (candidate, target))
    lower_bound = distance[31]
    assert lower_bound is not None
    # lower_bound > D^(17/12) > D^sqrt(2), entirely over integers.
    assert lower_bound ** 12 > dimension ** 17
    assert 17 ** 2 > 2 * 12 ** 2
    return {"m": m, "dimension": dimension, "coefficient_span_dimension": 5,
            "central_coefficient_minimum_rank_product": minimum_product,
            "integer_floor_shortest_path_lower_bound": lower_bound,
            "finite_field_scope": "F_2 coefficient recombinations",
            "integer_post_half_comparison": True, "all_recombination_ranks": all_ranks}


def generic_entry_polynomials():
    # Four distinct positions, independently enumerating compatible columns.
    def word(plus, minus):
        return tuple(int(i == plus) + 2 * int(i == minus) for i in range(4))
    from itertools import permutations
    x, y = word(0, 1), word(2, 3)
    row = []
    for first, second in product(permutations((0, 3)), permutations((1, 2))):
        v, u = word(first[0], second[0]), word(second[1], first[1])
        row.append(tensor_entry(x, y, u, v))
    x, u = word(0, 1), word(2, 3)
    column = []
    for first, second in product(permutations((0, 2)), permutations((1, 3))):
        v, y = word(first[0], second[0]), word(first[1], second[1])
        column.append(tensor_entry(x, y, u, v))
    assert row == [16, 80, 80, 272]
    assert column == [16, 20, 20, 17]
    return {"row_polynomial_bits": row, "column_polynomial_bits": column}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert 5 ** 7 > 42 ** 3 and 10 ** 2 > 2 * 7 ** 2
    assert (7 * 6 * 5 * 4 // 4) ** 7 > (7 * 6) ** 10
    result = {
        "status": "exact generic-block obstruction for this typed Gaussian coefficient span; no new factoring exponent",
        "generic_blocks": [check_generic_blocks(counts) for counts in (
            (2, 1, 1, 0), (2, 2, 0, 0), (2, 1, 0, 1), (2, 0, 0, 1), (3, 1, 1, 1),
        )],
        "generic_polynomials": generic_entry_polynomials(),
        "small_one_plus_minus": [small_one_plus_minus(m) for m in range(3, 7)],
        "all_m_ge_7_threshold": {"integer_comparison_5pow7_gt_42pow3": True,
                                  "rational_exponent_10_over_7_gt_sqrt2": True,
                                  "first_case": {"m": 7, "dimension": 42, "blocks": 210}},
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("tools/check_typed_gaussian_generic_blocks.py",
                                       "tools/check_schur_gaussian_reduction.py")},
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
