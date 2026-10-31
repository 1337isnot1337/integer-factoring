#!/usr/bin/env python3
"""Exact localized quartic identities and zero-place Smith profiles over F2.

The checks certify formulas and endpoint obstructions, not an improved
factoring exponent. No native library or numerical optimizer is required.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank, smith_profile, toeplitz_rank
from check_schur_downset_crossing import (
    binomial_polynomial, poly_matrix_multiply, sparse_kernel,
)
from check_schur_gaussian_reduction import polynomial_multiply as mul


ROOT = Path(__file__).resolve().parent.parent


def subset_sums(values):
    result = [0]
    for value in values:
        result += [entry + value for entry in result]
    return result


def zero_profile(m, a, b, c, valuation):
    """Orders for d=m-1, a<d, and ord_0(f)=valuation."""
    d = m - 1
    assert 0 <= a < d and 0 <= b < m and 0 <= c < m
    energy = {(leg, j): 2 * j + 1 for leg in (0, 1) for j in range(m)}
    first = {(0, a), (1, b)}
    second = {(0, d), (1, c)}
    active = first | second
    A = sum(energy[key] for key in first)
    B = sum(energy[key] for key in second)
    pivot = min(A, B, valuation + 2 * (a + b + 1 - m))
    result = Counter(subset_sums(energy.values()))
    for offset in subset_sums(value for key, value in energy.items() if key not in active):
        result[A + offset] -= 1
        result[B + offset] -= 1
        result[pivot + offset] += 1
        result[A + B - pivot + offset] += 1
    assert min(result.values()) >= 0
    return tuple(sorted(result.elements()))


def scalar_tensor(m):
    kernel = sparse_kernel(tuple(range(m)), m)
    one_particle = [row[:2 * m] for row in kernel[2 * m:]]

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        value = 0
        remaining = columns
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length() - 1
            value ^= mul(one_particle[i][j], minor(rows ^ bit, columns ^ cb))
        return value

    q = 1 << m
    return [minor(v | (y << m), u | (x << m))
            for v, u, y, x in product(range(q), repeat=4)]


def perturb(values, m, a, b, c, numerator, denominator):
    """Return denominator*G + numerator*QG as polynomial coefficient bits."""
    monomial = ((1 << a) | (1 << (m + b)) | (1 << (2 * m + c))
                | (1 << (4 * m - 1)))
    result = [mul(denominator, value) for value in values]
    for mask, value in enumerate(values):
        if not mask & monomial:
            result[mask | monomial] ^= mul(numerator, value)
    return result


def gauged_matrix(values, m, infinity):
    q = 1 << m
    if infinity:
        heights = [[2 * m - 1] * m, [1 - 2 * m] * m,
                   [2 * a for a in range(m)],
                   [2 * a + 2 - 2 * m for a in range(m)]]
    else:
        heights = [[2 * a for a in range(m)],
                   [2 * a + 2 - 2 * m for a in range(m)],
                   [2 * m - 1] * m, [1 - 2 * m] * m]
    letter_heights = [[sum(row[a] for a in range(m) if letter >> a & 1)
                       for letter in range(q)] for row in heights]
    entries = []
    for mask, polynomial in enumerate(values):
        if not polynomial:
            continue
        labels = tuple(mask >> (m * leg) & (q - 1) for leg in range(4))
        shift = sum(letter_heights[leg][labels[leg]] for leg in range(4))
        exponents = [shift + (-j if infinity else j)
                     for j in range(polynomial.bit_length()) if polynomial >> j & 1]
        entries.append((labels, exponents))
    minimum = min(min(exponents) for _, exponents in entries)
    matrix = [[0] * (q * q) for _ in range(q * q)]
    for (x, y, u, v), exponents in entries:
        row, column = (q * x + u, q * y + v) if infinity else (q * x + y, q * u + v)
        matrix[row][column] = sum(1 << (e - minimum) for e in exponents)
    return matrix, minimum


def raw_profile(values, m, infinity, denominator):
    matrix, minimum = gauged_matrix(values, m, infinity)
    degree = max(p.bit_length() - 1 for row in matrix for p in row)
    orders, precision = smith_profile(matrix, degree)
    denominator_order = (-(denominator.bit_length() - 1) if infinity
                         else (denominator & -denominator).bit_length() - 1)
    restored = tuple(value + minimum - denominator_order for value in orders)
    return restored, matrix, orders, precision


def check_small_raw_profiles():
    m = 2
    base = scalar_tensor(m)
    baseline_infinity = sum(raw_profile(base, m, True, 1)[0])
    # numerator, denominator, zero valuation, infinity determinant-ratio order.
    factors = [(3, 1 << 3, -3, 1), (5, 1 << 4, -4, 4),
               (1, 1 << 2, -2, 2), (1, 6, -1, 1), (1, 5, 0, None)]
    records = []
    toeplitz_checks = 0
    for b, c in product(range(m), repeat=2):
        for numerator, denominator, valuation, gain in factors:
            values = perturb(base, m, 0, b, c, numerator, denominator)
            actual, matrix, local_orders, precision = raw_profile(values, m, False, denominator)
            assert actual == zero_profile(m, 0, b, c, valuation)
            if b == c == 0 and numerator == 1:
                for length in (1, 2, max(local_orders) + 1):
                    assert toeplitz_rank(matrix, length) == sum(max(0, length - v) for v in local_orders)
                    toeplitz_checks += 1
            record = {"indices": [0, b, c, m - 1], "numerator_bits": numerator,
                      "denominator_bits": denominator, "zero_orders": actual,
                      "zero_precision": precision}
            if gain is not None:
                infinite = raw_profile(values, m, True, denominator)[0]
                expected_gain = gain * (1 << (2 * m - 2))
                assert sum(infinite) - baseline_infinity == expected_gain
                record["infinity_determinant_gain"] = expected_gain
            records.append(record)
    return records, toeplitz_checks


def check_four_mode_raw_profiles():
    """Replay the two motivating 256 by 256 raw matrices independently."""
    base = scalar_tensor(4)
    records = []
    for c in (0, 1):
        values = perturb(base, 4, 1, 2, c, 3, 8)
        actual, _, _, precision = raw_profile(values, 4, False, 8)
        assert actual == zero_profile(4, 1, 2, c, -3)
        records.append({"indices": [1, 2, c, 3], "orders": actual,
                        "zero_precision": precision})
    return records


def check_one_particle_inverse(m):
    kernel = sparse_kernel(tuple(range(m)), m)
    matrix = [row[:2 * m] for row in kernel[2 * m:]]
    pascal = [[int(b & ~a == 0) for b in range(m)] for a in range(m)]
    A = [[mul(3, binomial_polynomial(a, 2)) * pascal[b][a]
          for b in range(m)] for a in range(m)]
    B = [[3 << (2 * (m - 1 - (a | b))) for b in range(m)] for a in range(m)]
    BP = poly_matrix_multiply(B, pascal)
    PA = poly_matrix_multiply(pascal, A)
    inverse = [BP[a] + [(1 << (2 * m)) * int(a == b) for b in range(m)]
               for a in range(m)] + [pascal[a] + PA[a] for a in range(m)]
    expected = [[int(a == b) for b in range(2 * m)] for a in range(2 * m)]
    assert poly_matrix_multiply(matrix, inverse) == expected
    assert poly_matrix_multiply(inverse, matrix) == expected
    for a, b, c in product(range(m - 1), range(m), range(m)):
        d = m - 1
        contraction = mul(BP[c][d], PA[a][b])
        contraction ^= ((1 << (2 * m)) * int(c == b) * pascal[a][d])
        assert contraction == 5 << (2 * (a & b))
    # The top-index variable is fixed by the inverse Pascal change.
    assert [pascal[a][m - 1] for a in range(m)] == [0] * (m - 1) + [1]
    return {"width": m, "inverse_entries": 4 * m * m,
            "inverse_minor_checks": (m - 1) * m * m}


def matrix_product(left, right):
    result = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            row ^= bit
            value ^= right[bit.bit_length() - 1]
        result.append(value)
    return result


def check_creation_annihilation():
    modes = 4
    size = 1 << modes
    creations = [[0] * size for _ in range(modes)]
    annihilations = [[0] * size for _ in range(modes)]
    for j, mask in product(range(modes), range(size)):
        if not mask >> j & 1:
            creations[j][mask | (1 << j)] |= 1 << mask
        else:
            annihilations[j][mask ^ (1 << j)] |= 1 << mask
    count = 0
    for first, second in product(range(size), repeat=2):
        B = []
        for mask in (first, second):
            rows = [0] * size
            for j in range(modes):
                if mask >> j & 1:
                    rows = [a ^ b for a, b in zip(rows, annihilations[j])]
            B.append(rows)
        N = matrix_product(matrix_product(matrix_product(creations[0], creations[1]), B[0]), B[1])
        determinant = ((first & 1) * ((second >> 1) & 1)) ^ (((first >> 1) & 1) * (second & 1))
        assert matrix_product(N, N) == (N if determinant else [0] * size)
        if determinant:
            assert binary_rank(N) == size // 4
        count += 1
    return count


def exterior_product(left, right):
    result = set()
    for a, b in product(left, right):
        if a & b:
            continue
        mask = a | b
        if mask in result:
            result.remove(mask)
        else:
            result.add(mask)
    return result


def check_residual_endpoint():
    # Physical variable order x0,xT,y0,yT,u0,uT,v0,vT.
    factors = [(0, 6), (1, 7), (4, 2), (5, 3),
               (5, 6), (1, 2), (0, 3, 4, 7)]
    coefficients = {0}
    for factor in factors:
        coefficients = exterior_product(coefficients, {0, sum(1 << bit for bit in factor)})
    flattenings = []
    for cut in (0, 1):
        rows = [0] * 16
        for mask in coefficients:
            x, y, u, v = (mask >> (2 * leg) & 3 for leg in range(4))
            row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
            rows[row] |= 1 << column
        flattenings.append(binary_rank(rows))
    assert flattenings[1] == 12
    for x, y, u, v in product(range(2), repeat=4):
        mask = (1 << x) | ((1 << y) << 2) | ((1 << u) << 4) | ((1 << v) << 6)
        assert int(mask in coefficients) == int(x + y == u + v)
    words = ((1, 2), (2, 1), (1, 2), (2, 1))
    value = 1
    for site in range(2):
        mask = sum(words[leg][site] << (2 * leg) for leg in range(4))
        value *= int(mask in coefficients)
    assert value == 1 and not (words[0] == words[3] and words[1] == words[2])
    return {"local_flattening_ranks": flattenings,
            "singleton_alphabet": [1, 2], "singleton_sum_identity_checks": 16,
            "off_crossing_common_type_words": words,
            "off_crossing_value": value, "all_occupied_endpoint": int(255 in coefficients)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw, toeplitz_checks = check_small_raw_profiles()
    result = {"status": "exact localized quartic formulas; no improved factoring bound",
              "one_particle": [check_one_particle_inverse(m) for m in (2, 4, 8)],
              "creation_annihilation_checks": check_creation_annihilation(),
              "raw_two_mode_profiles": raw, "independent_toeplitz_checks": toeplitz_checks,
              "raw_four_mode_profiles": check_four_mode_raw_profiles(),
              "four_mode_predicted_zero_profiles": [
                  {"c": c, "valuation": valuation,
                   "orders": zero_profile(4, 1, 2, c, valuation)}
                  for c, valuation in product(range(4), (-3, -4, -2, -1, 0))],
              "rank_deficient_endpoint_control": check_residual_endpoint(),
              "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                  for name in ("tools/check_localized_quartic_profiles.py",
                               "tools/check_schur_downset_crossing.py",
                               "tools/check_block_smith_crossings.py")}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
