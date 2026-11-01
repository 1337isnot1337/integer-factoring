#!/usr/bin/env python3
"""Exact finite geometric-transpose crossing certificates over F2.

The default run checks q=2,4,8,16 using only the standard library.  The
optional polynomial screen exhausts all 2^14 normalized polynomial changes
for the q=8 geometric channel and the independent q=4 regular control.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_nilpotent_flag_crossings import (
    apply_columns, compose_columns, make_q4_channel, nilpotent_family,
    matmul, phi, rank, reverse_at_degree, smith_series,
)


def independent_basis(vectors):
    pivots = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return list(pivots.values())


def invert_columns(columns):
    n = len(columns)
    rows = [sum(((columns[j] >> i) & 1) << j for j in range(n))
            | (1 << (n + i)) for i in range(n)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if (rows[i] >> j) & 1), None)
        if pivot is None:
            raise ValueError("singular matrix")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(n):
            if i != j and ((rows[i] >> j) & 1):
                rows[i] ^= rows[j]
    return [sum(((rows[i] >> (n + j)) & 1) << i for i in range(n))
            for j in range(n)]


def transpose_mask(value, q):
    return sum(((value >> (i * q + j)) & 1) << (j * q + i)
               for i in range(q) for j in range(q))


def geometric_channel(q):
    """Columns of S=I+U for P=(I+N^T)^(-1), N the lower shift.

    U(Y) = sum_i N^i P Y^T P^(-1) N^i (I+N^2).
    Each packed operation is a literal binary matrix operation.
    """
    if q < 2 or q & (q - 1):
        raise ValueError("q must be a power of two at least two")
    n = q * q
    full = (1 << n) - 1
    row_mask = (1 << q) - 1
    no_first_column = sum(((row_mask ^ 1) << (i * q)) for i in range(q))
    no_last_column = sum(((row_mask >> 1) << (i * q)) for i in range(q))
    last_two_zero = sum(((row_mask >> 2) << (i * q)) for i in range(q))

    def operator(value):
        value = transpose_mask(value, q)
        # Multiplication on the left by the upper triangular all-ones P.
        suffix, left = 0, 0
        for i in range(q - 1, -1, -1):
            suffix ^= (value >> (q * i)) & row_mask
            left |= suffix << (q * i)
        # Multiplication on the right by P^{-1}=I+N^T.
        term = left ^ ((left << 1) & no_first_column)
        total = 0
        for _ in range(q):
            total ^= term
            term = (((term << q) & full) >> 1) & no_last_column
        # Multiplication on the right by I+N^2.
        return total ^ ((total >> 2) & last_two_zero)

    columns = [operator(1 << j) ^ (1 << j) for j in range(n)]
    identity = sum(1 << (i * q + i) for i in range(q))
    assert apply_columns(columns, identity) == 0
    assert all(sum((v >> (i * q + i)) & 1 for i in range(q)) % 2 == 0
               for v in columns)
    return columns


def geometric_channel_reference(q):
    """Independent literal matrix-product implementation for small controls."""
    n = q * q
    ident = [int(i == j) for i in range(q) for j in range(q)]
    shift = [int(i == j + 1) for i in range(q) for j in range(q)]
    upper_ones = [int(i <= j) for i in range(q) for j in range(q)]
    upper_bidiagonal = [int(i == j or j == i + 1)
                        for i in range(q) for j in range(q)]
    powers = [ident]
    for _ in range(1, q):
        powers.append(matmul(powers[-1], shift, q))
    shift_squared = matmul(shift, shift, q)
    normalization = [a ^ b for a, b in zip(ident, shift_squared)]
    result = []
    for column in range(n):
        unit = [int(j == column) for j in range(n)]
        transposed = [unit[j * q + i] for i in range(q) for j in range(q)]
        middle = matmul(matmul(upper_ones, transposed, q),
                        upper_bidiagonal, q)
        total = [0] * n
        for power in powers:
            term = matmul(matmul(power, middle, q), power, q)
            total = [a ^ b for a, b in zip(total, term)]
        total = matmul(total, normalization, q)
        result.append(sum((a ^ b) << j for j, (a, b) in
                          enumerate(zip(total, unit))))
    return result


def powers_of(columns, count):
    powers = [[1 << i for i in range(len(columns))]]
    for _ in range(count):
        powers.append(compose_columns(columns, powers[-1]))
    return powers


def coefficient_columns(polynomial, degree):
    n = len(polynomial)
    return [sum(((polynomial[i][j] >> degree) & 1) << i for i in range(n))
            for j in range(n)]


def joint_nilpotence_certificate(polynomial):
    degree = max(v.bit_length() - 1 for row in polynomial for v in row)
    n = len(polynomial)
    constant_inverse = invert_columns(coefficient_columns(polynomial, 0))
    generators = [compose_columns(constant_inverse,
                                  coefficient_columns(polynomial, k))
                  for k in range(1, degree + 1)]
    space = [1 << i for i in range(n)]
    dimensions = [n]
    while space:
        next_space = independent_basis(
            apply_columns(g, v) for g in generators for v in space)
        dimensions.append(len(next_space))
        assert len(next_space) < len(space), dimensions
        space = next_space
    # Products of this many generators are zero.  Thus I+sum z^i A_i is
    # invertible over F2[z], and its determinant is the constant one.
    return dimensions


def literal_endpoint(polynomial, q):
    n = q * q
    for i in range(n):
        for j in range(n):
            assert polynomial[i][j].bit_count() % 2 == int(
                j == (i % q) * q + i // q)


def phi_rows(columns, q):
    return [sum(v << j for j, v in enumerate(row)) for row in phi(columns, q)]


def polynomial_screen(columns, q, index):
    """Complete f(S)=S+sum_{j=2}^{index-1} c_j S^j rank search."""
    powers = powers_of(columns, index)
    assert not any(powers[-1])
    reshuffled = [phi_rows(p, q) for p in powers]
    current = [a ^ b for a, b in zip(reshuffled[0], reshuffled[1])]
    previous_gray = 0
    histogram, minimizers = {}, []
    minimum = q * q + 1
    for step in range(1 << (index - 2)):
        gray = step ^ (step >> 1)
        if step:
            changed_exponent = (gray ^ previous_gray).bit_length() + 1
            current = [a ^ b for a, b in zip(current,
                                             reshuffled[changed_exponent])]
        previous_gray = gray
        value = rank(current)
        histogram[value] = histogram.get(value, 0) + 1
        if value < minimum:
            minimum, minimizers = value, []
        if value == minimum:
            minimizers.append(gray)
    return {"candidate_count": 1 << (index - 2),
            "rank_histogram": histogram, "minimum_rank": minimum,
            "minimum_masks": minimizers}


def verify(polynomial_changes=False):
    records = []
    expected_geometric_orders = {
        2: [0, 0, 1, 11],
        4: [0] * 4 + [1] * 6 + [8] * 2 + [15] * 2 + [30] * 2,
        8: [0] * 8 + [1] * 28 + [v for v in
            [7, 8, 13, 18, 20, 28, 30, 31, 33, 34, 46, 56, 64, 78]
            for _ in range(2)],
    }
    for q in [2, 4, 8, 16]:
        n, m = q * q, 2 * q
        columns = geometric_channel(q)
        if q <= 4:
            assert columns == geometric_channel_reference(q)
        powers = powers_of(columns, m)
        ranks = [rank(p) for p in powers]
        assert ranks == [(q // 2) * (m - k) for k in range(m + 1)]
        polynomial = nilpotent_family(columns, q, m)
        literal_endpoint(polynomial, q)
        dimensions = joint_nilpotence_certificate(polynomial)
        # U^{-1}=I+S+...+S^{m-1}; the new nilpotent is its nonconstant part.
        inverse_nilpotent = [0] * n
        for power in powers[1:m]:
            inverse_nilpotent = [a ^ b for a, b in zip(inverse_nilpotent, power)]
        inverted = nilpotent_family(inverse_nilpotent, q, m)
        literal_endpoint(inverted, q)
        inverse_dimensions = joint_nilpotence_certificate(inverted)
        assert inverse_dimensions == dimensions
        inverse_u = [x ^ (1 << i) for i, x in enumerate(inverse_nilpotent)]
        u = [x ^ (1 << i) for i, x in enumerate(columns)]
        assert compose_columns(u, inverse_u) == powers[0]
        assert compose_columns(inverse_u, u) == powers[0]
        assert rank(phi_rows(u, q)) == q
        assert rank(phi_rows(inverse_u, q)) == 2
        leading = phi_rows(inverse_u, q)
        subleading = phi_rows(compose_columns(inverse_nilpotent, inverse_u), q)
        restricted_rank = rank(leading + subleading) - rank(leading)
        alternating_jordan_rank = n // 2
        protected_product = restricted_rank * alternating_jordan_rank
        record = {"q": q, "index": m, "jordan_blocks": [m] * (q // 2),
                  "joint_image_dimensions": dimensions,
                  "inverse_joint_image_dimensions": inverse_dimensions,
                  "leading_row_rank": q, "inverse_leading_row_rank": 2,
                  "inverse_subleading_row_rank": rank(subleading),
                  "inverse_subleading_rank_on_leading_kernel": restricted_rank,
                  "alternating_jordan_rank": alternating_jordan_rank,
                  "protected_subleading_product": protected_product,
                  "column_zero_orders": [0] * (q // 2)
                      + [m] * (n - q // 2)}
        if q >= 4:
            assert protected_product ** 35 > q ** 99
        if q <= 8:
            expected = expected_geometric_orders[q]
            actual = smith_series(reverse_at_degree(polynomial, m - 1),
                                  max(expected) + 1)
            assert actual == expected
            assert sum(actual) == n * (m - 1)
            record["row_infinity_orders"] = actual
        if q == 8:
            expected = [0] * 2 + [1] * 49 + [2] * 2 + [5, 9, 13] \
                + [96] * 4 + [120] * 2 + [128] * 2
            actual = smith_series(reverse_at_degree(inverted, m - 1), 129)
            assert actual == expected
            record["inverse_row_infinity_orders"] = actual
            if polynomial_changes:
                screen = polynomial_screen(columns, q, m)
                assert screen["minimum_rank"] == 2
                assert screen["minimum_masks"] == [(1 << (m - 2)) - 1]
                assert sum(count for r, count in screen["rank_histogram"].items()
                           if r <= 5) == 1
                record["polynomial_screen"] = screen
        records.append(record)
    result = {"geometric_instances": records,
              "claim_scope": "exact finite instances; no exponent improvement"}
    if polynomial_changes:
        _, regular, _ = make_q4_channel()
        screen = polynomial_screen(regular, 4, 16)
        assert screen["minimum_rank"] == 2 and screen["minimum_masks"] == [0]
        assert sum(count for r, count in screen["rank_histogram"].items()
                   if r <= 3) == 1
        regular_u = [x ^ (1 << i) for i, x in enumerate(regular)]
        leading = phi_rows(regular_u, 4)
        subleading = phi_rows(compose_columns(regular, regular_u), 4)
        restricted_rank = rank(leading + subleading) - rank(leading)
        assert (rank(leading), rank(subleading), restricted_rank) == (2, 14, 13)
        screen["protected_subleading_product"] = restricted_rank * 8
        assert (restricted_rank * 8) ** 35 > 4 ** 99
        # All other candidates are excluded by their leading coefficient.
        assert (16 * 4) ** 35 > 4 ** 99
        assert (64 * 6) ** 35 > 8 ** 99
        result["q4_regular_polynomial_screen"] = screen
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--polynomial-screen", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.polynomial_screen)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
