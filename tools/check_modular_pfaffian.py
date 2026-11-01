"""Exact checks of the folded characteristic-two Pfaffian formula.

The reference construction uses the original normally ordered coefficient
recurrence and intermediate partition sums. The Pfaffian path uses neither.
Finite rank profiles do not assert an asymptotic width improvement.
"""

import argparse
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import random

from check_half_piece_ranks import matmul_mod, partitions, zeros
from check_modular_pieces import (
    quotient_afam, quotient_opcoef, regular_partitions,
)
from modular_pfaffian import (
    final_piece_profile, folded_matrix, one_plus_z_power,
    polynomial_exact_quotient, polynomial_pfaffian, polynomial_product,
)


def matching_pfaffian(matrix):
    """Independent perfect-matching expansion, only for tiny test matrices."""
    if not matrix:
        return 1
    answer = 0
    for partner in range(1, len(matrix)):
        keep = [i for i in range(1, len(matrix)) if i != partner]
        answer ^= polynomial_product(
            matrix[0][partner],
            matching_pfaffian([[matrix[i][j] for j in keep] for i in keep]))
    return answer


def binary_rank(rows):
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def reference_matrices(h, labels):
    """Source formula (54), evaluated independently of the Pfaffian path."""
    d = len(labels)
    result = []
    for shift in range(-2 * h, 2 * h + 1):
        channels = [(i, tau) for i in range(2 * h + 1)
                    for tau in regular_partitions(shift + i, 2)]
        left = zeros(d * d, len(channels))
        right = zeros(len(channels), d * d)
        for column, (i, tau) in enumerate(channels):
            for x, first in enumerate(labels):
                terms = [(mu, value * comb(h + sum(mu), i) % 2)
                         for mu, value in quotient_afam(first, 2)
                         if h + sum(mu) >= i]
                terms = [(mu, value) for mu, value in terms if value]
                for y, state in enumerate(labels):
                    left[x * d + y][column] = sum(
                        value * quotient_opcoef(mu, tau, state, 2)
                        for mu, value in terms) % 2
            for u, output in enumerate(labels):
                for v, second in enumerate(labels):
                    right[column][u * d + v] = quotient_opcoef(
                        second, output, tau, 2)
        result.append(matmul_mod(left, right, 2))
    return result


def source_case(h):
    labels = tuple(sorted(regular_partitions(h, 2)))
    d = len(labels)
    reference = reference_matrices(h, labels)
    entries = 0
    for x, y, u, v in product(range(d), repeat=4):
        polynomial = final_piece_profile(
            labels[x], labels[v], labels[u], labels[y])
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        for offset, matrix in enumerate(reference):
            assert ((polynomial >> offset) & 1) == matrix[x*d+y][u*d+v]
            entries += 1
    return {"weight": h, "labels": d, "source_piece_entries": entries,
            "crossing_entries": d ** 4}


def sea_labels(cardinality, excitation):
    """Partitions near the filled distinct-odd-part ground state."""
    answer = []
    for excess in partitions(excitation):
        if len(excess) > cardinality:
            continue
        excess += (0,) * (cardinality - len(excess))
        answer.append(tuple(2 * (cardinality - i) - 1 + 2 * e
                            for i, e in enumerate(excess)))
    return tuple(sorted(answer))


def sea_case(cardinality, excitation):
    labels = sea_labels(cardinality, excitation)
    d = len(labels)
    h = cardinality ** 2 + 2 * excitation
    assert all(sum(label) == h for label in labels)
    matrices = {}
    for x, y, u, v in product(range(d), repeat=4):
        polynomial = final_piece_profile(
            labels[x], labels[v], labels[u], labels[y])
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        while polynomial:
            bit = polynomial & -polynomial
            shift = bit.bit_length() - 1 - 2 * h
            if shift not in matrices:
                matrices[shift] = [0] * (d * d)
            matrices[shift][x*d+y] |= 1 << (u*d+v)
            polynomial ^= bit
    profiles = []
    for shift, matrix in sorted(matrices.items()):
        column = [sum(((matrix[x*d+y] >> (u*d+v)) & 1) << (y*d+v)
                      for y in range(d) for v in range(d))
                  for x in range(d) for u in range(d)]
        profiles.append({"s": shift, "row_rank": binary_rank(matrix),
                         "column_rank": binary_rank(column)})
    answer = {
        "cardinality": cardinality, "excitation": excitation,
        "weight": h, "labels": [list(label) for label in labels],
        "dimension": d, "crossing_entries": d ** 4,
        "nonzero_pieces": len(profiles), "ranks": profiles,
        "max_rank_product": max(p["row_rank"] * p["column_rank"]
                                for p in profiles),
    }
    polynomial_product.cache_clear()
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()
    random.seed(1729)
    small_pfaffians = 0
    for size in (0, 2, 4, 6, 8):
        for _ in range(16):
            matrix = zeros(size, size)
            for i in range(size):
                for j in range(i + 1, size):
                    matrix[i][j] = matrix[j][i] = random.randrange(16)
            assert polynomial_pfaffian(matrix) == matching_pfaffian(matrix)
            small_pfaffians += 1
    for left in range(32):
        for right in range(1, 32):
            assert polynomial_exact_quotient(
                polynomial_product(left, right), right) == left

    folded_blocks = 0
    for n in range(1, 32, 2):
        for t in range(1, 32, 2):
            ua = ay = ab = 0
            for j in range(1, n + 1, 2):
                if comb(n, j) % 2:
                    if j <= t and comb(t - 1, j - 1) % 2:
                        ua ^= one_plus_z_power(j) << (n + t - j)
                    if comb(t + j - 1, j - 1) % 2:
                        ay ^= one_plus_z_power(j) << (n - j)
                        ab ^= 1 << n
            assert ua == 3 << (n | t)
            assert ay == 3 << ((n & t) - 1)
            assert ab == ((1 << n) if (n - 1) & ~(t - 1) == 0 else 0)
            folded_blocks += 3

    # The output profile is not the Pfaffian after dropping the folded factor.
    label = (5, 3)
    matrix = folded_matrix(label, label, label, label)
    good = polynomial_pfaffian(matrix)
    matrix[0][2] = matrix[2][0] = 1 << (label[0] | label[0])
    assert polynomial_pfaffian(matrix) != good

    source_weights = (8, 12, 16) if args.full else (8, 12)
    source_cases = [source_case(h) for h in source_weights]
    seas = [(2, 2), (4, 2), (8, 2), (3, 3), (6, 3)]
    if args.full:
        seas += [(16, 2), (12, 3), (4, 4), (8, 4)]
    sea_cases = [sea_case(k, e) for k, e in seas]
    print(json.dumps({
        "characteristic": 2,
        "tiny_pfaffians_checked_by_matchings": small_pfaffians,
        "polynomial_product_division_cases": 32 * 31,
        "folded_block_binomial_checks": folded_blocks,
        "missing_fold_factor_negative_control": True,
        "source_comparisons": source_cases,
        "source_piece_entries": sum(c["source_piece_entries"] for c in source_cases),
        "fermi_sea_cases": sea_cases,
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "implementation_sha256": sha256(Path(__file__).with_name(
            "modular_pfaffian.py").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "polynomial-time final-piece values and finite ranks; no improved asymptotic width",
        "status": "all exact checks passed",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
