#!/usr/bin/env python3
"""Arbitrary-mode typed endpoints and independent-species pressure controls.

This checks exact constructions and counting identities. It does not
establish a stronger integer-factoring exponent.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import random

from check_schur_downset_crossing import cross_blocks, exact_orders, infinity_gauge, sparse_kernel
from check_schur_gaussian_reduction import (
    evaluate, flip, partition_from_frobenius, polynomial_determinant,
)
from schur_fermion import BinaryField, crossing_value_binary


ROOT = Path(__file__).resolve().parent.parent


def binary_determinant(matrix):
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        return 0
    rows = [sum(value << j for j, value in enumerate(row)) for row in matrix]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i] >> j & 1), None)
        if pivot is None:
            return 0
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(j + 1, n):
            if rows[i] >> j & 1:
                rows[i] ^= rows[j]
    return 1


def bits(mask, n):
    return tuple(i for i in range(n) if mask >> i & 1)


def endpoint_value(kernel, n, x, y, u, v):
    # A label is (arm-mode mask, complemented-leg-mode mask).
    plus = (u[1], x[1], v[0], y[0])
    minus = (u[0], x[0], v[1], y[1])
    rows = [block * n + i for block in range(4) for i in bits(plus[block], n)]
    columns = [block * n + i for block in range(4) for i in bits(minus[block], n)]
    if len(rows) != len(columns):
        return 0
    return binary_determinant([[kernel[i][j] for j in columns] for i in rows])


def minor(p, row_mask, column_mask):
    rows, columns = bits(row_mask, len(p)), bits(column_mask, len(p))
    if len(rows) != len(columns):
        return 0
    return binary_determinant([[p[i][j] for j in columns] for i in rows])


def mode_type_words(length, count_a, count_b):
    modes = len(count_a)
    options = [tuple(combinations(range(length), k)) for k in (*count_a, *count_b)]
    result = []
    for choices in product(*options):
        word = [[0, 0] for _ in range(length)]
        for mode, occupied in enumerate(choices):
            for position in occupied:
                word[position][mode // modes] |= 1 << (mode % modes)
        result.append(tuple(map(tuple, word)))
    expected = 1
    for k in (*count_a, *count_b):
        expected *= comb(length, k)
    assert len(result) == expected and len(set(result)) == expected
    return tuple(result)


def check_mode_type(indices, m, length, count_a, count_b):
    n = len(indices)
    polynomial = sparse_kernel(indices, m)
    endpoint = [[value.bit_count() % 2 for value in row] for row in polynomial]
    p = [[int(b & ~a == 0) for b in indices] for a in indices]
    labels = mode_type_words(length, count_a, count_b)
    pair_checks = 0
    for x, v in product(labels, repeat=2):
        value = 1
        for xp, vp in zip(x, v):
            value &= minor(p, xp[1], vp[1]) & minor(p, vp[0], xp[0])
        assert value == int(x == v)
        pair_checks += 1
    rng = random.Random(90000 + m * 100 + n * 10 + length)
    if len(labels) <= 8:
        queries = list(product(range(len(labels)), repeat=4))
    else:
        queries = [tuple(rng.randrange(len(labels)) for _ in range(4)) for _ in range(512)]
        queries += [(i, j, j, i) for i in range(len(labels)) for j in range(len(labels))]
    for xi, yi, ui, vi in queries:
        x, y, u, v = [labels[j] for j in (xi, yi, ui, vi)]
        value = 1
        for positions in zip(x, y, u, v):
            value &= endpoint_value(endpoint, n, *positions)
        assert value == int(xi == vi and yi == ui)

    # Arbitrary per-mode heights, with no affine-in-index restriction.
    restored = []
    for _ in range(2):
        heights = [[[rng.randrange(-5, 6) for _ in range(n)] for _ in range(2)] for _ in range(4)]
        expected_a = [heights[0][1][a] + heights[1][1][a] + heights[2][0][a] + heights[3][0][a]
                      for a in range(n)]
        expected_b = [heights[0][0][a] + heights[1][0][a] + heights[2][1][a] + heights[3][1][a]
                      for a in range(n)]
        for boundary in range(4):
            species_counts = (count_b, count_a) if boundary < 2 else (count_a, count_b)
            expected = sum(heights[boundary][species][a] * species_counts[species][a]
                           for species in range(2) for a in range(n))
            for word in labels:
                actual = 0
                for a_mask, b_mask in word:
                    masks = (b_mask, a_mask) if boundary < 2 else (a_mask, b_mask)
                    actual += sum(heights[boundary][species][a]
                                  for species in range(2) for a in bits(masks[species], n))
                assert actual == expected
        restored.append(sum(k * g for k, g in zip(count_a, expected_a))
                        + sum(k * g for k, g in zip(count_b, expected_b)))
    return {"indices": indices, "M": m, "copies": length, "A_mode_counts": count_a,
            "B_mode_counts": count_b, "label_dimension": len(labels),
            "unequal_total_species_occupations": sum(count_a) != sum(count_b),
            "inclusion_pairings": pair_checks, "direct_gaussian_endpoint_entries": len(queries),
            "arbitrary_height_label_checks": 8 * len(labels), "sample_restored_shifts": restored}


def negative_control():
    indices, m = (0, 1), 2
    endpoint = [[value.bit_count() % 2 for value in row] for row in sparse_kernel(indices, m)]
    x, v, empty = (1, 0), (2, 0), (0, 0)
    assert x != v and x[0].bit_count() == v[0].bit_count()
    assert endpoint_value(endpoint, 2, x, empty, empty, v) == 1
    return {"indices": indices, "x": x, "v": v, "y_and_u": empty,
            "same_species_cardinalities": True, "different_A_index_sums": [0, 1],
            "off_diagonal_endpoint_value": 1}


def check_variable_energy(indices, m, r, field):
    n = len(indices)
    subsets = tuple(combinations(range(n), r))
    groups = defaultdict(list)
    for aa, bb in product(subsets, repeat=2):
        a, b = tuple(indices[i] for i in aa), tuple(indices[i] for i in bb)
        delta = sum(a) - sum(b)
        partition = partition_from_frobenius(a, [m - 1 - value for value in b])
        assert sum(partition) == m * r + delta
        groups[delta].append((aa, bb, partition))
    w = 2
    z = field.power(w, 2)
    matrix = [[evaluate(value, w, field) for value in row] for row in sparse_kernel(indices, m)]
    records = []
    for delta, labels in sorted(groups.items()):
        checks = 0
        for x, y, u, v in product(labels, repeat=4):
            plus = (u[1], x[1], v[0], y[0])
            minus = (u[0], x[0], v[1], y[1])
            rows = [block * n + i for block in range(4) for i in plus[block]]
            columns = [block * n + i for block in range(4) for i in minus[block]]
            value = field.determinant([[matrix[i][j] for j in columns] for i in rows])
            restored = field.multiply(field.power(w, 4 * delta), value)
            expected = crossing_value_binary(x[2], y[2], u[2], v[2], z, field)
            assert restored == expected
            checks += 1
        records.append({"delta": delta, "h": m * r + delta, "labels": len(labels), "entries": checks})
    return {"indices": indices, "M": m, "r": r, "energy_slices": records}


def check_gibbs_counts():
    binary_counts = 0
    for length in range(1, 13):
        for k in range(length + 1):
            entropy_power = Fraction(length ** length, (k ** k) * ((length - k) ** (length - k)))
            dimension = comb(length, k)
            assert dimension <= entropy_power <= (length + 1) * dimension
            binary_counts += 1
    product_checks = 0
    root, numerator, denominator = Fraction(3, 2), 3, 2
    s = root ** numerator
    for weights in ((0, 1), (-2, 1, 3, 0), (4, -1, 2, -3, 0, 1)):
        local = tuple(root ** (-denominator * weight) for weight in weights)
        partition = 1
        for value in local:
            partition *= 1 + value
        explicit = sum(root ** (-denominator * sum(weight * bit for weight, bit in zip(weights, word)))
                       for word in product(range(2), repeat=len(weights)))
        assert partition == explicit
        probabilities = tuple(value / (1 + value) for value in local)
        assert all(0 < p < 1 for p in probabilities)
        assert s == root ** numerator
        product_checks += 1

    balanced_checks = 0
    for n in range(1, 7):
        subsets = [tuple(a for a in range(n) if mask >> a & 1) for mask in range(1 << n)]
        sectors = Counter((len(a), sum(a) - sum(b)) for a, b in product(subsets, repeat=2) if len(a) == len(b))
        x = Fraction(2, 3)
        for constant, slope in product((-2, 0, 4), (-2, 0, 2)):
            centered = [Fraction(constant, 2) + slope * (Fraction(a) - Fraction(n - 1, 2)) for a in range(n)]
            assert all(g.denominator == 1 for g in centered)
            partition = 1
            for g in centered:
                partition *= (1 + x ** int(g)) ** 2
            weighted = {key: count * x ** (constant * key[0] + slope * key[1]) for key, count in sectors.items()}
            equal_count = sum(weighted.values())
            assert partition / (n + 1) <= equal_count <= partition
            assert max(weighted.values()) >= partition / ((n + 1) * len(sectors))
            balanced_checks += 1
    return {"exact_binary_type_bounds": binary_counts,
            "exact_product_partition_expansions": product_checks,
            "balanced_equal_count_and_energy_pigeonhole_checks": balanced_checks}


def check_exterior_profiles():
    records = []
    for m, index in ((2, 0), (4, 1), (8, 3)):
        original = sparse_kernel((index,), m)
        infinity = infinity_gauge(original, (index,), m)
        for end, matrix in (("zero", original), ("gauged_infinity", infinity)):
            tensor = {}
            for x, y, u, v in product(range(4), repeat=4):
                states = (flip(u), flip(x), v, y)
                rows = [i for i, state in enumerate(states) if state & 1]
                columns = [i for i, state in enumerate(states) if state & 2]
                tensor[x, y, u, v] = polynomial_determinant([[matrix[i][j] for j in columns] for i in rows]) if len(rows) == len(columns) else 0
            for cut in range(2):
                small = sorted(sum((exact_orders(block) for block in cross_blocks(matrix, 1, cut)), []))
                expected = sorted(sum(small[i] for i in range(4) if mask >> i & 1) for mask in range(16))
                flattened = [[0] * 16 for _ in range(16)]
                for (x, y, u, v), value in tensor.items():
                    row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
                    flattened[row][column] = value
                actual = exact_orders(flattened)
                assert actual == expected
                records.append({"M": m, "index": index, "end": end, "cut": cut,
                                "one_particle_orders": small, "full_tensor_orders": actual})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {
        "status": "exact arbitrary-mode typed crossing and pressure reduction; no improved exponent",
        "typed_endpoints": [check_mode_type(*case) for case in (
            ((0, 1), 4, 3, (1, 2), (0, 1)),
            ((0, 3), 4, 2, (1, 1), (1, 0)),
            ((1, 2), 4, 2, (1, 0), (0, 1)),
            ((0, 1, 2), 4, 2, (1, 1, 1), (0, 1, 0)),
            ((0, 3, 5), 8, 2, (1, 1, 0), (0, 1, 1)),
        )],
        "omitted_mode_weight_negative_control": negative_control(),
        "variable_energy_schur": [check_variable_energy(indices, m, r, BinaryField(7))
                                   for indices, m, r in (((0, 1, 2), 4, 1), ((0, 3), 4, 1),
                                                         ((0, 2, 5), 8, 1), ((0, 1, 2, 3), 8, 2))],
        "counting": check_gibbs_counts(), "exterior_profiles": check_exterior_profiles(),
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in ("tools/check_fixed_mode_crossings.py",
                                       "tools/check_schur_downset_crossing.py",
                                       "tools/check_schur_gaussian_reduction.py", "tools/schur_fermion.py")},
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
