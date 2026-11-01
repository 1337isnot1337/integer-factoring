#!/usr/bin/env python3
"""Check the compact permutation seed, exact row-rank formula, and limits."""

import argparse
from functools import lru_cache
from itertools import permutations, product
import json
from math import factorial
from pathlib import Path
import random

from check_permutation_homogeneous_labels import rank_two
from fermionic_piece_polynomial import polynomial_product
from split_permutation_seed import (
    coefficient_tensors, column_matrix, curve_matrices, split_polynomial,
    unpack_matrix,
)


def direct_pfaffian(a, b, u, y):
    size = len(a)
    particles = [
        (group, (i, label[i]))
        for group, label in (("U", u), ("B", b), ("A", a), ("Y", y))
        for i in range(size)
    ]
    count = len(particles)
    matrix = [[0] * count for _ in range(count)]
    for i, (g, e) in enumerate(particles):
        z = 1 << (5 ** e[0])
        t = 1 ^ (1 << (5 ** (size + e[1])))
        for j in range(i + 1, count):
            h, f = particles[j]
            if e != f:
                continue
            pair = g, h
            if pair in (("U", "B"), ("A", "Y")):
                matrix[i][j] = t
            elif pair in (("U", "A"), ("B", "Y")):
                matrix[i][j] = polynomial_product(z, t)
            elif pair in (("U", "Y"), ("B", "A")):
                matrix[i][j] = z

    @lru_cache(None)
    def recurse(remaining):
        if not remaining:
            return 1
        bit = remaining & -remaining
        first = bit.bit_length() - 1
        rest = remaining ^ bit
        result = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            second = bit.bit_length() - 1
            if matrix[first][second]:
                result ^= polynomial_product(matrix[first][second],
                                             recurse(rest ^ bit))
        return result
    return recurse((1 << count) - 1)


def embed_even(polynomial, variables):
    result = 0
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        subset = bit.bit_length() - 1
        code = sum(2 * 5**j for j in range(variables) if subset >> j & 1)
        result ^= 1 << code
    return result


def cycle_graphs(labels):
    size = len(labels[0])
    graphs = {}
    for a, y in product(labels, repeat=2):
        difference = frozenset(enumerate(a)) ^ frozenset(enumerate(y))
        if difference in graphs:
            continue
        rows = sum(1 << i for i in range(size) if a[i] != y[i])
        columns = sum(1 << a[i] for i in range(size) if a[i] != y[i])
        inverse_y = {value: i for i, value in enumerate(y)}
        remaining = {i for i in range(size) if rows >> i & 1}
        cycles = []
        while remaining:
            start = next(iter(remaining))
            current = start
            cycle_columns = 0
            while current in remaining:
                remaining.remove(current)
                cycle_columns |= 1 << a[current]
                current = inverse_y[a[current]]
            assert current == start
            cycles.append(cycle_columns)
        graphs[difference] = (rows, columns, tuple(cycles))
    return tuple(graphs.values())


def falling(total, chosen):
    if chosen > total:
        return 0
    return factorial(total) // factorial(total - chosen)


def row_rank_formula(s, t, graphs):
    result = 0
    for rows, columns, cycles in graphs:
        if rows & ~s:
            continue
        available = (t & ~columns).bit_count()
        selected = (s & ~rows).bit_count()
        result += (1 << sum(bool(t & cycle) for cycle in cycles)) * falling(available, selected)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    rng = random.Random(915820)
    source_checks = 0
    for size in (2, 3, 4):
        labels = tuple(permutations(range(size)))
        cases = product(labels, repeat=4) if size <= 3 else (
            tuple(rng.choice(labels) for _ in range(4)) for _ in range(256))
        for a, b, u, y in cases:
            compact = split_polynomial(a, b, u, y)
            assert embed_even(compact, 2 * size) == direct_pfaffian(a, b, u, y)
            assert compact.bit_count() % 2 == int(a == b and u == y)
            source_checks += 1
            polynomial_product.cache_clear()

    cases = []
    rank_formula_checks = 0
    curve_checks = 0
    for size in (2, 3, 4):
        labels, tensors = coefficient_tensors(size)
        count = len(labels)
        dimension = count * count
        summed = 0
        for tensor in tensors:
            summed ^= tensor
        assert summed == sum(1 << (dimension * row + row) for row in range(dimension))
        graphs = cycle_graphs(labels)
        profiles = {}
        maximum = 0
        complete_profiles = []
        for subset, tensor in enumerate(tensors):
            s = subset & ((1 << size) - 1)
            t = subset >> size
            rows = unpack_matrix(tensor, dimension)
            r = rank_two(rows)
            c = rank_two(column_matrix(rows, count))
            assert r == row_rank_formula(s, t, graphs)
            rank_formula_checks += 1
            key = s.bit_count(), t.bit_count()
            if key in profiles:
                assert profiles[key] == (r, c)
            profiles[key] = (r, c)
            maximum = max(maximum, r * c)
            complete_profiles.append([subset, r, c])
        expected_maximum = {2: 4, 3: 150, 4: 34560}[size]
        assert maximum == expected_maximum
        full_support = sum(rows == (1 << size) - 1 for rows, _, _ in graphs)
        no_fixed = sum(all(p[i] != i for i in range(size)) for p in labels)
        bad_profile = profiles[(size, 0)]
        assert bad_profile[0] == full_support
        assert bad_profile[1] >= count * no_fixed
        for _ in range(16):
            coordinates = [rng.choice((1, 2, 4, 7)) for _ in range(2 * size)]
            row_curve, _ = curve_matrices(tensors, coordinates, count)
            crossing = [0] * dimension
            for matrix in row_curve:
                for i, row in enumerate(matrix):
                    crossing[i] ^= row
            assert crossing == [1 << i for i in range(dimension)]
            curve_checks += 1
        cases.append({
            "particles": size, "labels": count,
            "nonzero_pieces": sum(bool(t) for t in tensors),
            "nonzero_degree_orbits": sum(bool(r) for r, _ in profiles.values()),
            "maximum_rank_product": maximum,
            "worst_orbits": [[s, t, r, c] for (s, t), (r, c) in sorted(profiles.items())
                            if r * c == maximum],
            "full_rows_empty_columns_ranks": list(bad_profile),
            "full_support_union_graphs": full_support,
            "ordered_disjoint_pairs": count * no_fixed,
            "orbit_profiles": [[s, t, r, c] for (s, t), (r, c) in sorted(profiles.items())],
            "complete_profiles": complete_profiles,
        })
    result = {
        "independent_pfaffian_checks": source_checks,
        "exact_row_rank_formula_checks": rank_formula_checks,
        "nonlinear_curve_crossing_checks": curve_checks,
        "cases": cases,
        "scope": "finite six-label seed; full permutation-domain growth is obstructed",
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(encoded)
    print(json.dumps({**result, "cases": [
        {key: value for key, value in case.items()
         if key not in ("orbit_profiles", "complete_profiles")} for case in cases
    ]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
