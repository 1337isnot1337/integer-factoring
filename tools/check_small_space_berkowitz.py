#!/usr/bin/env python3
"""Independent finite checks for scalar Berkowitz queries and pivot masks."""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path
import random


class ScalarBerkowitz:
    def __init__(self, entry, order, modulus):
        self.entry = entry
        self.n = order
        self.q = modulus
        self.entry_queries = 0

    def read(self, i, j):
        self.entry_queries += 1
        return self.entry(i, j) % self.q

    def power_entry(self, start, exponent, i, j):
        if exponent == 0:
            return int(i == j)
        if exponent == 1:
            return self.read(start + i, start + j)
        left = exponent // 2
        value = 0
        for k in range(self.n - start):
            first = self.power_entry(start, left, i, k)
            second = self.power_entry(start, exponent - left, k, j)
            value = (value + first * second) % self.q
        return value

    def toeplitz_entry(self, start, i, j):
        delta = i - j
        if delta < 0:
            return 0
        if delta == 0:
            return 1
        if delta == 1:
            return -self.read(start, start) % self.q
        value = 0
        size = self.n - start - 1
        for a in range(size):
            for b in range(size):
                value = (value + self.read(start, start + 1 + a)
                         * self.power_entry(start + 1, delta - 2, a, b)
                         * self.read(start + 1 + b, start)) % self.q
        return -value % self.q

    def product_entry(self, left, right, i, j):
        if left == right:
            return int(i == j)
        if right - left == 1:
            return self.toeplitz_entry(left, i, j)
        middle = (left + right) // 2
        value = 0
        for k in range(self.n - middle + 1):
            first = self.product_entry(left, middle, i, k)
            second = self.product_entry(middle, right, k, j)
            value = (value + first * second) % self.q
        return value

    def coefficient(self, degree_index):
        return self.product_entry(0, self.n, degree_index, 0)

    def determinant(self):
        return (-1)**self.n * self.coefficient(self.n) % self.q


def det_leibniz(matrix, modulus):
    n = len(matrix)
    value = 0
    for permutation in permutations(range(n)):
        sign = (-1)**sum(permutation[i] > permutation[j] for i in range(n) for j in range(i + 1, n))
        term = sign
        for i, j in enumerate(permutation):
            term *= matrix[i][j]
        value += term
    return value % modulus


def gaussian(matrix, modulus):
    a = [row[:] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    determinant = 1
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col] % modulus), None)
        if pivot is None:
            continue
        if pivot != rank:
            a[rank], a[pivot] = a[pivot], a[rank]
            determinant = -determinant
        diagonal = a[rank][col] % modulus
        determinant = determinant * diagonal % modulus
        inverse = pow(diagonal, -1, modulus)
        for i in range(rank + 1, rows):
            multiplier = a[i][col] * inverse % modulus
            for j in range(col, cols):
                a[i][j] = (a[i][j] - multiplier * a[rank][j]) % modulus
        rank += 1
        if rank == rows:
            break
    return rank, determinant % modulus if rows == cols == rank else 0


def select(mask, ordinal):
    while ordinal:
        mask &= mask - 1
        ordinal -= 1
    return (mask & -mask).bit_length() - 1


def masked_determinant(entry, row_mask, col_mask, q):
    order = row_mask.bit_count()
    assert order == col_mask.bit_count()
    return ScalarBerkowitz(lambda i, j: entry(select(row_mask, i), select(col_mask, j)), order, q).determinant()


def greedy_rank(entry, rows, cols, q):
    row_mask = col_mask = 0
    tests = 0
    while True:
        changed = False
        for i in range(rows):
            if row_mask >> i & 1:
                continue
            for j in range(cols):
                if col_mask >> j & 1:
                    continue
                tests += 1
                if masked_determinant(entry, row_mask | (1 << i), col_mask | (1 << j), q):
                    row_mask |= 1 << i
                    col_mask |= 1 << j
                    changed = True
                    break
            if changed:
                break
        if not changed:
            return row_mask.bit_count(), row_mask, col_mask, tests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(20261030)
    coefficients = determinants = rings = ranks = recoveries = consistency = 0
    oracle_queries = 0

    # Principal-minor coefficients are independent of the Toeplitz product.
    for q in (2, 3, 5, 7, 17):
        for n in range(7):
            for _ in range(3):
                a = [[rng.randrange(q) for _ in range(n)] for _ in range(n)]
                oracle = ScalarBerkowitz(lambda i, j: a[i][j], n, q)
                for k in range(n + 1):
                    expected = (-1)**k * sum(
                        det_leibniz([[a[i][j] for j in subset] for i in subset], q)
                        for subset in combinations(range(n), k)) % q
                    assert oracle.coefficient(k) == expected
                    coefficients += 1
                assert oracle.determinant() == gaussian(a, q)[1]
                determinants += 1
                oracle_queries += oracle.entry_queries

    # Division-free determinant identities also survive nonfield rings.
    for q in (4, 8, 9, 12):
        for n in range(1, 6):
            for _ in range(3):
                a = [[rng.randrange(q) for _ in range(n)] for _ in range(n)]
                oracle = ScalarBerkowitz(lambda i, j: a[i][j], n, q)
                assert oracle.determinant() == det_leibniz(a, q)
                rings += 1

    for q in (2, 3, 5, 7):
        for rows in range(1, 7):
            for cols in range(1, 7):
                inner = rng.randrange(min(rows, cols) + 1)
                u = [[rng.randrange(q) for _ in range(inner)] for _ in range(rows)]
                v = [[rng.randrange(q) for _ in range(cols)] for _ in range(inner)]
                a = [[sum(u[i][k] * v[k][j] for k in range(inner)) % q for j in range(cols)] for i in range(rows)]
                result, imask, jmask, tests = greedy_rank(lambda i, j: a[i][j], rows, cols, q)
                assert result == gaussian(a, q)[0]
                assert imask.bit_count() == jmask.bit_count() == result
                assert masked_determinant(lambda i, j: a[i][j], imask, jmask, q)
                assert tests <= sum((rows - k) * (cols - k) for k in range(min(rows, cols)))
                ranks += 1

    delayed = [[0] * 6 for _ in range(5)]
    delayed[4][5] = 1
    result = greedy_rank(lambda i, j: delayed[i][j], 5, 6, 2)
    assert result[:3] == (1, 1 << 4, 1 << 5) and result[3] >= 30
    ranks += 1

    for q in (2, 3, 5):
        for _ in range(20):
            rows, cols = rng.randrange(1, 6), rng.randrange(1, 5)
            a = [[rng.randrange(q) for _ in range(cols)] for _ in range(rows)]
            x = [rng.randrange(q) for _ in range(cols)]
            b = [sum(v * w for v, w in zip(row, x)) % q for row in a]
            r, imask, jmask, _ = greedy_rank(lambda i, j: a[i][j], rows, cols, q)
            augmented = [row + [target] for row, target in zip(a, b)]
            assert greedy_rank(lambda i, j: augmented[i][j], rows, cols + 1, q)[0] == r
            consistency += 1
            denominator = masked_determinant(lambda i, j: a[i][j], imask, jmask, q)
            for pivot in range(r):
                def cramer(i, j):
                    row = select(imask, i)
                    if j != pivot:
                        return a[row][select(jmask, j)]
                    return (b[row] - sum(a[row][col] * x[col] for col in range(cols) if not (jmask >> col & 1))) % q
                numerator = ScalarBerkowitz(cramer, r, q).determinant()
                assert numerator * pow(denominator, -1, q) % q == x[select(jmask, pivot)]
                recoveries += 1
            arbitrary = [row + [rng.randrange(q)] for row in a]
            greedy_consistent = greedy_rank(lambda i, j: arbitrary[i][j], rows, cols + 1, q)[0] == r
            assert greedy_consistent == (gaussian(arbitrary, q)[0] == r)
            consistency += 1

    result = {
        "status": "passed",
        "scope": "Finite identities and oracle/mask algorithms; packed workspace follows from the written analysis, not Python memory measurements.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "characteristic_coefficients_vs_principal_minors": coefficients,
        "field_determinants_vs_elimination": determinants,
        "nonfield_ring_determinants_vs_leibniz": rings,
        "greedy_mask_ranks_vs_elimination": ranks,
        "affine_consistency_comparisons": consistency,
        "cramer_coordinate_recoveries": recoveries,
        "uncached_input_entry_queries_in_coefficient_checks": oracle_queries,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
