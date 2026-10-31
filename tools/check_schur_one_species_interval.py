#!/usr/bin/env python3
"""Exact odd Smith spectra and rational pressure witnesses for one species.

The displayed family approaches the existing half-rate boundary from above.
It does not establish a stronger factoring exponent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from check_block_smith_crossings import smith_profile, toeplitz_rank
from check_fixed_mode_crossings import check_mode_type
from check_schur_downset_crossing import sparse_kernel
from check_schur_gaussian_reduction import polynomial_multiply


ROOT = Path(__file__).resolve().parent.parent


def power(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_multiply(result, value)
        value = polynomial_multiply(value, value)
        exponent >>= 1
    return result


def multiply(left, right):
    result = [[0] * len(right[0]) for _ in left]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(right[k]):
                    if other:
                        result[i][j] ^= polynomial_multiply(value, other)
    return result


def normalized_matrix(n):
    """The infinity matrix after unit changes, in the variable t=xi^2."""
    matrix = [[0] * (2 * n) for _ in range(2 * n)]
    for a in range(n):
        for b in range(n):
            matrix[a][b] = (power(3, a) << b) if a & ~b == 0 else 0
            matrix[a][n + b] = (1 << a) if b & ~a == 0 else 0
            matrix[n + a][b] = (1 << a) if a == b else 0
            matrix[n + a][n + b] = 3 << (a | b)
    return matrix


def check_algebra(ambient, n):
    assert ambient >= 2 and ambient & (ambient - 1) == 0
    assert n in (ambient - 1, ambient)
    pascal = [[int(b & ~a == 0) for b in range(n)] for a in range(n)]
    transpose = list(map(list, zip(*pascal)))
    q = [[power(3, a) if a == b else 0 for b in range(n)] for a in range(n)]
    diagonal = [[(1 << a) if a == b else 0 for b in range(n)] for a in range(n)]
    union = [[1 << (a | b) for b in range(n)] for a in range(n)]
    dp = multiply(diagonal, pascal)
    lhs = multiply(multiply(q, transpose), union)
    rhs = [[polynomial_multiply(power(3, ambient - 1), dp[a][b])
            ^ ((power(3, a) << n) if n == ambient - 1 else 0)
            for b in range(n)] for a in range(n)]
    assert lhs == rhs
    assert multiply(pascal, pascal) == [[int(a == b) for b in range(n)] for a in range(n)]
    if n == ambient - 1:
        assert all(sum(pascal[a][b] for a in range(n)) % 2 == 1 for b in range(n))
    schur = [[dp[a][b] ^ polynomial_multiply(3, lhs[a][b])
              for b in range(n)] for a in range(n)]
    reduced = [[(value >> n) for value in row] for row in schur]
    assert all(value & ((1 << n) - 1) == 0 for row in schur for value in row)
    rank_one = multiply(reduced, pascal)
    if n == ambient - 1:
        assert rank_one == [[((1 << (a + 1)) if a == b else 0) ^ power(3, a + 1)
                             for b in range(n)] for a in range(n)]
    else:
        assert rank_one == diagonal
    matrix = normalized_matrix(n)
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    actual, precision = smith_profile(matrix, degree)
    assert actual == tuple(range(2 * n))
    independent = 0
    if n <= 3:
        for jets in range(1, 2 * n + 1):
            assert toeplitz_rank(matrix, jets) == sum(max(0, jets - order) for order in actual)
            independent += 1
    return {"M": ambient, "modes": n, "normalized_t_orders": actual,
            "infinity_w_orders": [2 * order + 1 - 2 * n for order in actual],
            "smith_precision": precision, "independent_toeplitz_ranks": independent,
            "pascal_identity": True,
            "schur_complement": "diagonal after a unit change" if n == ambient else "diagonal plus rank one"}


def reverse_gauge(value, height, common_shift):
    result = 0
    while value:
        bit = value & -value
        exponent = height - (bit.bit_length() - 1) + common_shift
        assert exponent >= 0
        result ^= 1 << exponent
        value ^= bit
    return result


def check_direct_gauges(ambient, n):
    kernel = sparse_kernel(tuple(range(n)), ambient)
    zero = ([0] * n, [-2 * ambient + 1 + 2 * a for a in range(n)],
            [0] * n, [1 + 2 * a for a in range(n)])
    infinity = ([2 * a - n + 1 for a in range(n)], [2 * ambient - n] * n,
                [2 * a - n + 1 for a in range(n)], [-n] * n)
    regular_entries = 0
    for row_group, column_group in ((2, 0), (3, 1)):
        for a in range(n):
            for b in range(n):
                value = kernel[row_group * n + a][column_group * n + b]
                if value:
                    valuation = (value & -value).bit_length() - 1
                    assert zero[row_group][a] + zero[column_group][b] + valuation >= 0
                    regular_entries += 1
    zero_orders = []
    for row_group, column_group in ((2, 1), (3, 0)):
        matrix = []
        for a in range(n):
            row = []
            for b in range(n):
                value = kernel[row_group * n + a][column_group * n + b]
                height = zero[row_group][a] + zero[column_group][b]
                if height < 0:
                    assert value & ((1 << -height) - 1) == 0
                    value >>= -height
                else:
                    value <<= height
                row.append(value)
            matrix.append(row)
        degree = max(value.bit_length() - 1 for row in matrix for value in row)
        actual, _ = smith_profile(matrix, degree)
        assert actual == tuple(range(1, 2 * n, 2))
        zero_orders.extend(actual)
    rows = [(group, a) for group in (2, 3) for a in range(n)]
    columns = [(group, a) for group in (0, 1) for a in range(n)]
    shift = 2 * n - 1
    matrix = [[reverse_gauge(kernel[g * n + a][h * n + b],
                             infinity[g][a] + infinity[h][b], shift)
               for h, b in columns] for g, a in rows]
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    actual, _ = smith_profile(matrix, degree)
    restored = tuple(order - shift for order in actual)
    assert restored == tuple(range(1 - 2 * n, 2 * n, 2))
    weights = [sum(zero[g][a] + infinity[g][a] for g in range(4)) for a in range(n)]
    assert weights == [8 * a - 4 * (n - 1) for a in range(n)]
    return {"M": ambient, "zero_R_orders": sorted(zero_orders),
            "infinity_C_orders": restored, "combined_A_heights": weights,
            "internal_regular_entries": regular_entries,
            "zero_heights_B_X_V_Y": zero, "infinity_heights_B_X_V_Y": infinity}


def balanced_product(values):
    values = list(values)
    if not values:
        return 1
    while len(values) > 1:
        values = [values[i] * values[i + 1] if i + 1 < len(values) else values[i]
                  for i in range(0, len(values), 2)]
    return values[0]


def rational_pressure(n, numerator, denominator, root_numerator, root_denominator):
    """Check B(s)^denominator < Z(s)^numerator after cancelling powers."""
    a, b, p, q = numerator, denominator, root_numerator, root_denominator
    k = (n - 1) // 2
    left = balanced_product((p ** (a * (2 * j - 1)) + q ** (a * (2 * j - 1))) ** (4 * b)
                            for j in range(1, n + 1))
    if n % 2:
        right = ((1 << a) * p ** (a * b * (2 * n * n + 1)) * q ** (a * b)
                 * balanced_product((p ** (8 * b * j) + q ** (8 * b * j)) ** (2 * a)
                                    for j in range(1, k + 1)))
    else:
        right = (p ** (2 * a * b * n * n)
                 * balanced_product((p ** (4 * b * (2 * j - 1)) + q ** (4 * b * (2 * j - 1))) ** (2 * a)
                                    for j in range(1, n // 2 + 1)))
    assert left < right
    difference = right - left
    probabilities = []
    root = Fraction(p, q)
    for mode in range(n):
        weight = 8 * mode - 4 * (n - 1)
        if n > 16:
            probabilities.append({"rational_root": [p, q], "power": -b * weight,
                                  "formula": "root^power / (1 + root^power)"})
            continue
        value = root ** (-b * weight)
        probability = value / (1 + value)
        probabilities.append({"numerator": probability.numerator, "denominator": probability.denominator})
    return {"modes": n, "rho_upper_numerator": a, "rho_upper_denominator": 4 * b,
            "R": [a, b], "s": {"root": [p, q], "power": a},
            "strict_integer_inequality": True,
            "left_bits": left.bit_length(), "right_bits": right.bit_length(),
            "positive_difference_sha256": sha256(difference.to_bytes((difference.bit_length() + 7) // 8, "big")).hexdigest(),
            "rational_A_probabilities": probabilities}


def even_two_mode_control():
    # D_2(log(u^3),3)>0 for u>=1. Expand the numerator at u=1+v.
    polynomial = [0] * 49
    for a in range(5):
        for b in range(5):
            polynomial[3 * a + 9 * b] += comb(4, a) * comb(4, b)
    for a in range(7):
        polynomial[24 + 4 * a] -= comb(6, a)
    while polynomial and polynomial[-1] == 0:
        polynomial.pop()
    translated = [sum(value * comb(i, j) for i, value in enumerate(polynomial) if i >= j)
                  for j in range(len(polynomial))]
    assert len(translated) == 46 and all(value > 0 for value in translated)
    assert translated[0] == 192
    return {"n": 2, "R": 3, "substitution": "s=u^3; u=1+v",
            "strictly_positive_integer_coefficients": translated,
            "rho_floor_strictly_above_three_quarters": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--large-certificate", action="store_true",
                        help="also verify the larger 63-mode rho<17/24 integer witness")
    args = parser.parse_args()
    certificates = [rational_pressure(7, 23, 8, 50, 49), rational_pressure(8, 43, 15, 100, 99)]
    if args.large_certificate:
        certificates.append(rational_pressure(63, 17, 6, 177, 176))
    result = {
        "status": "exact one-species interval crossing and signed-odd profiles; no exponent above one half",
        "all_size_proof_scope": "M a power of two and n=M-1 or n=M; explicit displayed gauges",
        "algebra_and_profiles": [check_algebra(m, n) for m in (2, 4, 8, 16) for n in (m - 1, m)],
        "direct_gauged_profiles": [check_direct_gauges(m, n) for m in (2, 4, 8) for n in (m - 1, m)],
        "typed_endpoint": check_mode_type((0, 1, 2), 4, 2, (1, 1, 1), (0, 0, 0)),
        "exact_pressure_certificates": certificates,
        "even_two_mode_control": even_two_mode_control(),
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in (
            "tools/check_schur_one_species_interval.py", "tools/check_fixed_mode_crossings.py",
            "tools/check_schur_downset_crossing.py", "tools/check_block_smith_crossings.py")},
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
