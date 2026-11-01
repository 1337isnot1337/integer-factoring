#!/usr/bin/env python3
"""Finite checks for weighted-Gram rank queries and virtual prefix bases.

The rank and coordinate algorithms recompute their entry queries and
prefix ranks without memoization. Dense reference matrices and minor
lists are used only by the independent finite checks. Python memory is
not evidence of the bit-space theorem in the accompanying proof.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import random
import sys
import time


HERE = Path(__file__).resolve().parent
if not (HERE / "check_seed_field_selection.py").exists():
    sys.path.insert(0, str(HERE.parent / "integer-factoring" / "tools"))

from check_seed_field_selection import Field  # noqa: E402
from check_small_space_berkowitz import det_leibniz, gaussian  # noqa: E402


@dataclass
class Counters:
    rank_calls: int = 0
    coefficient_calls: int = 0
    base_entry_queries: int = 0
    field_presentations: int = 0
    extension_presentations: int = 0


class Arithmetic:
    """The existing Field presentation, with uncached binary arithmetic."""

    def __init__(self, field: Field):
        self.field = field
        self.p, self.e, self.q = field.p, field.e, field.q
        self.binary_modulus = sum(c << i for i, c in enumerate(field.f)) if self.p == 2 else 0

    def add(self, a, b):
        if self.p == 2:
            return a ^ b
        if self.e == 1:
            return (a + b) % self.p
        return self.field.add(a, b)

    def neg(self, a):
        if self.p == 2:
            return a
        if self.e == 1:
            return -a % self.p
        return self.field.neg(a)

    def mul(self, a, b):
        if self.e == 1:
            return a * b % self.p
        if self.p != 2:
            return self.field.mul(a, b)
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a >> self.e:
                a ^= self.binary_modulus
        return result

    def power(self, a, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            exponent >>= 1
        return result

    def inv(self, a):
        if not a:
            raise ZeroDivisionError
        return self.power(a, self.q - 2)


class ScalarBerkowitz:
    """Generic-field version of the existing scalar-entry implementation."""

    def __init__(self, entry, order, arithmetic):
        self.entry, self.n, self.f = entry, order, arithmetic

    def power_entry(self, start, exponent, i, j):
        if exponent == 0:
            return int(i == j)
        if exponent == 1:
            return self.entry(start + i, start + j)
        left, result = exponent // 2, 0
        for k in range(self.n - start):
            first = self.power_entry(start, left, i, k)
            second = self.power_entry(start, exponent - left, k, j)
            result = self.f.add(result, self.f.mul(first, second))
        return result

    def toeplitz_entry(self, start, i, j):
        delta = i - j
        if delta < 0:
            return 0
        if delta == 0:
            return 1
        if delta == 1:
            return self.f.neg(self.entry(start, start))
        result = 0
        size = self.n - start - 1
        for a in range(size):
            for b in range(size):
                term = self.f.mul(self.entry(start, start + 1 + a),
                                  self.power_entry(start + 1, delta - 2, a, b))
                term = self.f.mul(term, self.entry(start + 1 + b, start))
                result = self.f.add(result, term)
        return self.f.neg(result)

    def product_entry(self, left, right, i, j):
        if left == right:
            return int(i == j)
        if right - left == 1:
            return self.toeplitz_entry(left, i, j)
        middle, result = (left + right) // 2, 0
        for k in range(self.n - middle + 1):
            first = self.product_entry(left, middle, i, k)
            second = self.product_entry(middle, right, k, j)
            result = self.f.add(result, self.f.mul(first, second))
        return result

    def coefficient(self, index):
        return self.product_entry(0, self.n, index, 0)

    def determinant(self):
        answer = self.coefficient(self.n)
        return self.f.neg(answer) if self.n & 1 else answer


def degree_bound(rows, cols):
    return rows * (rows + 1) // 2 + cols * (cols + 1) // 2


def choose_field(prime, degree, counters):
    """Rebuild the presentation for each query; no cached field search."""
    counters.field_presentations += 1
    if prime > degree:
        return Arithmetic(Field(prime, (0, 1)))
    extension_degree = 1
    while prime**extension_degree <= degree:
        extension_degree += 1
    for coefficients in product(range(prime), repeat=extension_degree):
        try:
            field = Field(prime, coefficients + (1,))
        except ValueError:
            continue
        counters.extension_presentations += 1
        return Arithmetic(field)
    raise AssertionError("no irreducible polynomial found")


def weighted_entry(entry, rows, cols, arithmetic, point, i, j, counters):
    result = 0
    for column in range(cols):
        first = entry(i, column)
        second = entry(j, column)
        counters.base_entry_queries += 2
        term = arithmetic.mul(first, second)
        term = arithmetic.mul(arithmetic.power(point, column + 1), term)
        result = arithmetic.add(result, term)
    return arithmetic.mul(arithmetic.power(point, i + 1), result)


def weighted_coefficient(entry, rows, cols, arithmetic, point, index, counters):
    counters.coefficient_calls += 1
    oracle = lambda i, j: weighted_entry(entry, rows, cols, arithmetic, point, i, j, counters)
    return ScalarBerkowitz(oracle, rows, arithmetic).coefficient(index)


def rank_query(entry, rows, cols, prime, counters):
    counters.rank_calls += 1
    if not rows or not cols:
        return 0
    degree = degree_bound(rows, cols)
    arithmetic = choose_field(prime, degree, counters)
    # Testing larger coefficients first is valid and avoids redundant work
    # when a full-rank matrix is detected at an early evaluation point.
    for index in range(min(rows, cols), 0, -1):
        for point in range(degree + 1):
            if weighted_coefficient(entry, rows, cols, arithmetic, point, index, counters):
                return index
    return 0


def prefix_rank(entry, rows, cols, prime, axis, length, counters):
    if axis == "row":
        return rank_query(entry, length, cols, prime, counters)
    return rank_query(entry, rows, length, prime, counters)


def basis_index(entry, rows, cols, prime, axis, ordinal, counters):
    """First prefix with rank >= ordinal+1, found without a pivot list."""
    low, high = 1, rows if axis == "row" else cols
    while low < high:
        middle = (low + high) // 2
        if prefix_rank(entry, rows, cols, prime, axis, middle, counters) >= ordinal + 1:
            high = middle
        else:
            low = middle + 1
    return low - 1


def is_pivot_column(entry, rows, cols, prime, column, counters):
    before = prefix_rank(entry, rows, cols, prime, "col", column, counters)
    after = prefix_rank(entry, rows, cols, prime, "col", column + 1, counters)
    return after > before


def requested_coordinate(entry, rhs, free_value, rows, cols, prime, column, counters):
    order = rank_query(entry, rows, cols, prime, counters)
    augmented = lambda i, j: rhs(i) if j == cols else entry(i, j)
    if rank_query(augmented, rows, cols + 1, prime, counters) != order:
        return None
    if not is_pivot_column(entry, rows, cols, prime, column, counters):
        return free_value(column)
    ordinal = prefix_rank(entry, rows, cols, prime, "col", column, counters)
    arithmetic = Arithmetic(Field(prime, (0, 1)))

    def pivot_entry(i, j):
        actual_row = basis_index(entry, rows, cols, prime, "row", i, counters)
        actual_col = basis_index(entry, rows, cols, prime, "col", j, counters)
        return entry(actual_row, actual_col)

    denominator = ScalarBerkowitz(pivot_entry, order, arithmetic).determinant()
    assert denominator

    def numerator_entry(i, j):
        if j != ordinal:
            return pivot_entry(i, j)
        actual_row = basis_index(entry, rows, cols, prime, "row", i, counters)
        target = rhs(actual_row)
        for free_column in range(cols):
            if not is_pivot_column(entry, rows, cols, prime, free_column, counters):
                target = (target - entry(actual_row, free_column) * free_value(free_column)) % prime
        return target

    numerator = ScalarBerkowitz(numerator_entry, order, arithmetic).determinant()
    return arithmetic.mul(numerator, arithmetic.inv(denominator))


def minor_polynomials(matrix, prime):
    rows, cols = len(matrix), len(matrix[0])
    degree = degree_bound(rows, cols)
    output = [[0] * (degree + 1) for _ in range(rows + 1)]
    output[0][0] = 1
    for order in range(1, min(rows, cols) + 1):
        for row_set in combinations(range(rows), order):
            for col_set in combinations(range(cols), order):
                determinant = det_leibniz([[matrix[i][j] for j in col_set] for i in row_set], prime)
                exponent = sum(i + 1 for i in row_set) + sum(j + 1 for j in col_set)
                output[order][exponent] = (output[order][exponent] + (-1)**order * determinant**2) % prime
    return output


def evaluate(coefficients, point, arithmetic):
    answer = 0
    for coefficient in reversed(coefficients):
        answer = arithmetic.add(arithmetic.mul(answer, point), coefficient)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    counters = Counters()
    coefficient_checks = rank_checks = basis_checks = recoveries = consistency_checks = 0
    arithmetic_checks = 0

    for field in [Field(2, (1, 1, 1)), Field(2, (1, 1, 0, 1)), Field(3, (1, 0, 1))]:
        fast = Arithmetic(field)
        for a, b in product(range(field.q), repeat=2):
            assert fast.add(a, b) == field.add(a, b)
            assert fast.mul(a, b) == field.mul(a, b)
            arithmetic_checks += 2

    rng = random.Random(20261030)
    matrices = [(2, [list(values[:3]), list(values[3:])]) for values in product(range(2), repeat=6)]
    for rows in range(1, 5):
        for cols in range(1, 5):
            matrices.append((2, [[rng.randrange(2) for _ in range(cols)] for _ in range(rows)]))
    for prime, rows, cols in [(3, 2, 3), (3, 3, 2), (5, 3, 4), (23, 4, 3), (23, 4, 4)]:
        matrices.append((prime, [[rng.randrange(prime) for _ in range(cols)] for _ in range(rows)]))

    for case, (prime, matrix) in enumerate(matrices):
        rows, cols = len(matrix), len(matrix[0])
        entry = lambda i, j: matrix[i][j]
        expected = gaussian(matrix, prime)[0]
        assert rank_query(entry, rows, cols, prime, counters) == expected
        rank_checks += 1
        # All tiny exhaustive matrices get coefficient-by-coefficient checks.
        # Larger matrices use three distinct extension points plus the
        # independent full rank-query scan above.
        arithmetic = choose_field(prime, degree_bound(rows, cols), counters)
        polynomials = minor_polynomials(matrix, prime)
        points = range(degree_bound(rows, cols) + 1) if rows <= 2 else [0, 1, min(2, arithmetic.q - 1)]
        for point in points:
            for index in range(rows + 1):
                actual = weighted_coefficient(entry, rows, cols, arithmetic, point, index, counters)
                assert actual == evaluate(polynomials[index], point, arithmetic)
                coefficient_checks += 1
        if case >= 64:
            row_basis = [basis_index(entry, rows, cols, prime, "row", i, counters) for i in range(expected)]
            col_basis = [basis_index(entry, rows, cols, prime, "col", i, counters) for i in range(expected)]
            assert det_leibniz([[matrix[i][j] for j in col_basis] for i in row_basis], prime)
            # Compare to the lexicographically least full row/column bases.
            reference_rows = next(indices for indices in combinations(range(rows), expected)
                                  if gaussian([matrix[i] for i in indices], prime)[0] == expected)
            reference_cols = next(indices for indices in combinations(range(cols), expected)
                                  if gaussian([[row[j] for j in indices] for row in matrix], prime)[0] == expected)
            assert tuple(row_basis) == reference_rows and tuple(col_basis) == reference_cols
            basis_checks += 1

    # The r=1 coefficient t^2+t^4 vanishes at every point of F_2.
    all_ones = [[1, 1], [1, 1]]
    polynomials = minor_polynomials(all_ones, 2)
    assert [i for i, c in enumerate(polynomials[1]) if c] == [2, 4]
    base = Arithmetic(Field(2, (0, 1)))
    assert all(evaluate(polynomials[1], point, base) == 0 for point in range(2))
    assert rank_query(lambda i, j: all_ones[i][j], 2, 2, 2, counters) == 1
    extension = choose_field(2, degree_bound(2, 2), counters)
    witnesses = [point for point in range(7) if evaluate(polynomials[1], point, extension)]
    assert witnesses and all(point >= 2 for point in witnesses)

    # At t=1, A=[[1,0],[1,0]] gives nonzero B=J_2 with B^2=0.
    nilpotent_source = [[1, 0], [1, 0]]
    nilpotent = [[weighted_entry(lambda i, j: nilpotent_source[i][j], 2, 2, base, 1, i, j, counters)
                  for j in range(2)] for i in range(2)]
    assert nilpotent == [[1, 1], [1, 1]]
    assert all(sum(nilpotent[i][k] * nilpotent[k][j] for k in range(2)) % 2 == 0
               for i, j in product(range(2), repeat=2))
    assert all(ScalarBerkowitz(lambda i, j: nilpotent[i][j], 2, base).coefficient(k) == 0 for k in (1, 2))
    assert rank_query(lambda i, j: nilpotent_source[i][j], 2, 2, 2, counters) == 1

    affine_cases = [
        (2, [[0, 1, 0], [0, 0, 1]], [1, 1, 0]),
        (2, [[1, 1, 1], [1, 0, 1]], [1, 0, 1]),
        (3, [[1, 2, 0], [2, 1, 0], [0, 0, 0]], [2, 1, 2]),
    ]
    for prime, matrix, known in affine_cases:
        rows, cols = len(matrix), len(matrix[0])
        rhs = [sum(x*y for x, y in zip(row, known)) % prime for row in matrix]
        entry = lambda i, j: matrix[i][j]
        for column in range(cols):
            actual = requested_coordinate(entry, lambda i: rhs[i], lambda j: known[j],
                                          rows, cols, prime, column, counters)
            assert actual == known[column]
            recoveries += 1
        arbitrary_rhs = [rng.randrange(prime) for _ in range(rows)]
        actual_consistent = rank_query(entry, rows, cols, prime, counters) == rank_query(
            lambda i, j: arbitrary_rhs[i] if j == cols else matrix[i][j], rows, cols + 1, prime, counters)
        expected_consistent = gaussian(matrix, prime)[0] == gaussian(
            [row + [b] for row, b in zip(matrix, arbitrary_rhs)], prime)[0]
        assert actual_consistent == expected_consistent
        consistency_checks += 1

    result = {
        "status": "passed",
        "scope": "Uncached scalar entry/rank/basis recomputation with dense independent finite references; no Python memory claim.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "field_arithmetic_checks": arithmetic_checks,
        "characteristic_coefficients_vs_weighted_minor_sums": coefficient_checks,
        "rank_comparisons": rank_checks + 2,
        "exhaustive_f2_2_by_3_matrices": 64,
        "lexicographic_prefix_basis_comparisons": basis_checks,
        "requested_cramer_coordinates": recoveries,
        "affine_consistency_comparisons": consistency_checks,
        "extension_required_example": {"matrix": all_ones, "coefficient_exponents": [2, 4],
                                       "field_order": extension.q, "nonbase_witnesses": witnesses},
        "nonzero_nilpotent_evaluation": {"source": nilpotent_source, "point": 1, "weighted_matrix": nilpotent},
        "recomputation_counters": counters.__dict__,
        "seconds": time.monotonic() - started,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
