#!/usr/bin/env python3
"""Exact source-level checks for the isotropic three-label laser obstruction."""

import argparse
from collections import Counter
from itertools import combinations, product
import json
from math import comb, factorial
from pathlib import Path

from check_finite_seed_amplification import PIECES, add, column, crossing, multiply, rank
from check_modular_crossing import inverse_two


def transpose(matrix):
    return tuple(sum(((row >> j) & 1) << i for i, row in enumerate(matrix))
                 for j in range(len(matrix)))


def kronecker(left, right):
    size = len(right)
    return tuple(sum(((x >> j) & 1) * y << (size * j) for j in range(len(left)))
                 for x in left for y in right)


def transformed_seed():
    # Columns are e0, e2, e0+e1.  Matched legs use the dual basis change.
    basis = (5, 4, 2)
    inverse = tuple(inverse_two(list(basis), 3))
    left = kronecker(transpose(basis), transpose(basis))
    right = kronecker(transpose(inverse), transpose(inverse))
    return tuple(multiply(left, multiply(piece, right)) for piece in PIECES)


def support(matrix):
    return tuple((x, y, u, v) for x, y, u, v in product(range(3), repeat=4)
                 if matrix[3 * x + y] >> (3 * u + v) & 1)


def balanced_words(m):
    size = 3 * m
    for ones in combinations(range(size), m):
        remaining = [i for i in range(size) if i not in ones]
        for twos in combinations(remaining, m):
            yield tuple(1 if i in ones else 2 if i in twos else 0 for i in range(size))


def literal_type_matrix(m, seed):
    """Expand the ORIGINAL transformed16/17/19 entries, then restrict.

    This does not use the asserted permutation/all-ones block formulas.
    Prefix composition pruning only removes words that cannot finish with
    m copies of each symbol.  Repeated allocations are added modulo two.
    """
    size = 3 * m
    words = tuple(balanced_words(m))
    labels = {word: index for index, word in enumerate(words)}
    dimension = len(words)
    rows = [0] * (dimension * dimension)
    local = {i: support(seed[i]) for i in (16, 17, 19)}
    assignments = leaves = additions = 0
    prefixes = [[], [], [], []]
    counts = [[0, 0, 0] for _ in range(4)]

    def walk(types, position):
        nonlocal leaves, additions
        if position == size:
            leaves += 1
            x, y, u, v = (labels[tuple(word)] for word in prefixes)
            rows[x * dimension + y] ^= 1 << (u * dimension + v)
            additions += 1
            return
        for entry in local[types[position]]:
            if any(counts[leg][symbol] == m for leg, symbol in enumerate(entry)):
                continue
            for leg, symbol in enumerate(entry):
                prefixes[leg].append(symbol)
                counts[leg][symbol] += 1
            walk(types, position + 1)
            for leg, symbol in enumerate(entry):
                prefixes[leg].pop()
                counts[leg][symbol] -= 1

    for first in combinations(range(size), m):
        remaining = [i for i in range(size) if i not in first]
        for second in combinations(remaining, m):
            types = tuple(16 if i in first else 17 if i in second else 19 for i in range(size))
            walk(types, 0)
            assignments += 1
    return tuple(rows), dimension, {"type_allocations": assignments,
                                   "balanced_literal_leaves": leaves,
                                   "mod_two_entry_additions": additions}


def verify():
    from functools import reduce
    seed = transformed_seed()
    assert reduce(add, seed, (0,) * 9) == crossing(3)
    all_entries = first_filtered = second_filtered = 0
    first_survivors = []
    second = []
    for i, piece in enumerate(seed):
        first_piece = [0] * 9
        second_piece = [0] * 9
        for x, y, u, v in support(piece):
            all_entries += 1
            charge = (x == 2) + (y == 2) - (u == 2) - (v == 2)
            assert charge <= 0
            if charge:
                continue
            first_filtered += 1
            first_piece[3 * x + y] |= 1 << (3 * u + v)
            charge_one = (x == 1) + (y == 1) - (u == 1) - (v == 1)
            assert charge_one <= 0
            if charge_one:
                continue
            second_filtered += 1
            second_piece[3 * x + y] |= 1 << (3 * u + v)
        if any(first_piece):
            first_survivors.append(i)
        second.append(tuple(second_piece))
    assert first_survivors == [0, 5, 10, 12, 15, 16, 17, 19]
    assert [i for i, piece in enumerate(second) if any(piece)] == [12, 15, 16, 17, 19]
    assert support(second[12]) == ((1, 1, 1, 1),)
    assert set(support(second[15])) == set(product((0, 2), repeat=4)) - {
        q for q in product((0, 2), repeat=4) if q[0] == q[1] or q[2] == q[3]}
    assert set(support(second[16])) == {(a, b, a, b) for a, b in product((0, 2), repeat=2)}
    assert set(support(second[17])) == {(a, 1, 1, a) for a in (0, 2)}
    assert set(support(second[19])) == {(1, a, a, 1) for a in (0, 2)}
    assert reduce(add, second, (0,) * 9) == crossing(3)

    profiles = []
    for m in (1, 2):
        rows, dimension, checks = literal_type_matrix(m, seed)
        expected_dimension = factorial(3 * m) // factorial(m) ** 3
        expected_row = dimension * sum(comb(m, k) ** 4 for k in range(m + 1))
        expected_column = dimension * comb(2 * m, m)
        assert dimension == expected_dimension
        row_rank, column_rank = rank(rows), rank(column(rows, dimension))
        assert (row_rank, column_rank) == (expected_row, expected_column)
        assert all(row.bit_count() <= 1 for row in rows)
        assert len({row for row in rows if row}) == row_rank
        profiles.append({"m": m, "length": 3 * m, "dimension": dimension,
                         "row_rank": row_rank, "column_rank": column_rank,
                         "rank_product": row_rank * column_rank, **checks})

    entropy_checks = 0
    for m in range(1, 257):
        largest = max(comb(m, k) for k in range(m + 1))
        fourth = sum(comb(m, k) ** 4 for k in range(m + 1))
        square = comb(2 * m, m)
        assert largest * (m + 1) >= 2 ** m
        assert fourth * (m + 1) ** 4 >= 16 ** m
        assert fourth <= (m + 1) * 16 ** m
        assert square * (2 * m + 1) >= 4 ** m
        assert square <= 4 ** m
        entropy_checks += 1
    # log(6)/(2log(3)) > 4/5 > 1/sqrt(2), with integer checks.
    assert 6 ** 5 > 3 ** 8 and 2 * 4 ** 2 > 5 ** 2
    return {"scope": "exact isotropic-basis type ranks; no improved factoring exponent",
            "source_seed_commit": "e96ee50560284301428d72421330d8773155d25b",
            "basis_rows": [5, 4, 2],
            "transformed_nonzero_entries": all_entries,
            "first_charge_survivor_entries": first_filtered,
            "second_charge_survivor_entries": second_filtered,
            "first_charge_surviving_pieces": first_survivors,
            "second_charge_surviving_pieces": [12, 15, 16, 17, 19],
            "literal_rank_profiles": profiles,
            "exact_entropy_bound_cases": entropy_checks,
            "asymptotic_rho_lower_bound": "log(6)/(2log(3))",
            "limitation": "the proved restriction uses equal counts of all three symbols"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    report = verify()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(rendered)
    print(rendered, end="")
