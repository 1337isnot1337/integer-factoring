#!/usr/bin/env python3
"""Conjugate F4 residual codes, polynomial endpoints, and exact F2 descent.

Polynomial coefficients use two adjacent bits for an F4 element a+b*omega,
with omega^2=omega+1. This checker contains no floating-point search.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank
from check_schur_downset_crossing import binomial_polynomial, sparse_kernel
from check_schur_gaussian_reduction import polynomial_multiply as multiply_binary
from check_localized_quartic_profiles import exterior_product


ROOT = Path(__file__).resolve().parent.parent


def scalar_multiply(polynomial, coefficient):
    if coefficient == 0:
        return 0
    if coefficient == 1:
        return polynomial
    mask = ((1 << (2 * ((polynomial.bit_length() + 1) // 2))) - 1) // 3
    low, high = polynomial & mask, (polynomial >> 1) & mask
    omega = high ^ ((low ^ high) << 1)
    return omega if coefficient == 2 else polynomial ^ omega


def multiply(left, right):
    answer = 0
    while right:
        degree = ((right & -right).bit_length() - 1) // 2
        coefficient = right >> (2 * degree) & 3
        right ^= coefficient << (2 * degree)
        answer ^= scalar_multiply(left, coefficient) << (2 * degree)
    return answer


def exact_divide(numerator, denominator):
    degree = (denominator.bit_length() - 1) // 2
    inverse = (0, 1, 3, 2)[denominator >> (2 * degree) & 3]
    answer = 0
    while numerator and (numerator.bit_length() - 1) // 2 >= degree:
        top = (numerator.bit_length() - 1) // 2
        coefficient = scalar_multiply(numerator >> (2 * top) & 3, inverse)
        shift = 2 * (top - degree)
        answer ^= coefficient << shift
        numerator ^= scalar_multiply(denominator, coefficient) << shift
    if numerator:
        raise ArithmeticError("polynomial division is not exact")
    return answer


def lift_binary(polynomial):
    return sum(1 << (2 * j) for j in range(polynomial.bit_length()) if polynomial >> j & 1)


def conjugate(polynomial):
    mask = ((1 << (2 * ((polynomial.bit_length() + 1) // 2))) - 1) // 3
    return polynomial ^ ((polynomial >> 1) & mask)


def evaluate_one(polynomial):
    result = 0
    while polynomial:
        result ^= polynomial & 3
        polynomial >>= 2
    return result


def check_residual_phases():
    residual = {0}
    for factor in ((0, 6), (1, 7), (4, 2), (5, 3),
                   (5, 6), (1, 2), (0, 3, 4, 7)):
        residual = exterior_product(residual, {0, sum(1 << a for a in factor)})
    checks = 0
    for p, q in product(range(4), repeat=2):
        input_maps = (((1, 1),), ((0, 1), (3, p)))
        output_maps = (((1, 1),), ((0, 1), (3, q)))
        phase = 1 ^ scalar_multiply(p, q)
        for x, y, u, v in product(range(2), repeat=4):
            value = 0
            for (xx, a), (yy, b), (uu, c), (vv, d) in product(
                    input_maps[x], output_maps[y], input_maps[u], output_maps[v]):
                mask = xx | (yy << 2) | (uu << 4) | (vv << 6)
                if mask in residual:
                    value ^= scalar_multiply(scalar_multiply(a, b), scalar_multiply(c, d))
            expected = (phase if x != u else 1) if x == v and y == u else 0
            assert value == expected
            checks += 1
    assert scalar_multiply(2, conjugate(2)) == 1
    assert scalar_multiply(3, conjugate(3)) == 1
    return checks


def one_copy_numerator(indices, ambient):
    """c^2 times the coded rational tensor; one bulk mode, two core states."""
    assert len(indices) == 3 and indices[0] == 0 and indices[-1] < ambient
    n = len(indices)
    kernel = sparse_kernel(indices, ambient)
    matrix = [row[:2 * n] for row in kernel[2 * n:]]

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        result, remaining = 0, columns
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length() - 1
            result ^= multiply_binary(matrix[i][j], minor(rows ^ bit, columns ^ cb))
        return result

    def code(label, output):
        bulk = (label & 1) << 1
        if label < 2:
            return ((bulk | 1, 1),)
        return ((bulk, 1), (bulk | 5, 2 if output else 1))

    quartic_rows, quartic_columns = (1 << 2) | (1 << 5), 1 | (1 << 3)
    answer = [0] * 256
    for x, y, u, v in product(range(4), repeat=4):
        total = 0
        for (xx, bx), (yy, by), (uu, bu), (vv, bv) in product(
                code(x, False), code(y, True), code(u, False), code(v, True)):
            rows, columns = vv | (yy << n), uu | (xx << n)
            polynomial = multiply_binary(5, minor(rows, columns))
            if rows & quartic_rows == quartic_rows and columns & quartic_columns == quartic_columns:
                polynomial ^= minor(rows ^ quartic_rows, columns ^ quartic_columns)
            coefficient = scalar_multiply(scalar_multiply(bx, by), scalar_multiply(bu, bv))
            total ^= scalar_multiply(lift_binary(polynomial), coefficient)
        answer[x + 4 * y + 16 * u + 64 * v] = total
    return answer


def grouped_seed(indices=(0, 1, 2), ambient=4, sign=1):
    """Bulk-copy bits occupy positions0,1; the repeated core bit is position2."""
    one = one_copy_numerator(indices, ambient)
    denominators = {j: lift_binary(binomial_polynomial(j)) for j in range(2, 7)}
    answer = [0] * 4096
    failures = 0
    for x, y, u, v in product(range(8), repeat=4):
        first = tuple((a & 1) + 2 * (a >> 2) for a in (x, y, u, v))
        second = tuple(((a >> 1) & 1) + 2 * (a >> 2) for a in (x, y, u, v))
        a = one[first[0] + 4 * first[1] + 16 * first[2] + 64 * first[3]]
        b = one[second[0] + 4 * second[1] + 16 * second[2] + 64 * second[3]]
        delta = sign * ((y >> 2) + (v >> 2) - (x >> 2) - (u >> 2))
        try:
            value = exact_divide(multiply(a, conjugate(b)), denominators[4 - delta])
        except ArithmeticError:
            failures += 1
            value = 0
        answer[x + 8 * y + 64 * u + 512 * v] = value
    return answer, failures


def swap_bulk(label):
    return (label & 4) | ((label & 1) << 1) | ((label & 2) >> 1)


def descent_matrices():
    B, C = [[0] * 8 for _ in range(8)], [[0] * 8 for _ in range(8)]
    for a in range(8):
        b = swap_bulk(a)
        if a == b:
            B[a][a] = C[a][a] = 1
        elif a < b:
            B[a][a], B[b][a], B[a][b], B[b][b] = 1, 1, 2, 3
            C[a][a], C[a][b], C[b][a], C[b][b] = 3, 2, 1, 1
    for a, b in product(range(8), repeat=2):
        value = 0
        for k in range(8):
            value ^= scalar_multiply(C[a][k], B[k][b])
        assert value == int(a == b)
        assert conjugate(B[a][b]) == B[swap_bulk(a)][b]
        assert conjugate(C[a][b]) == C[a][swap_bulk(b)]
    return B, C


def descend(values):
    B, C = descent_matrices()
    for x, y, u, v in product(range(8), repeat=4):
        old = x + 8 * y + 64 * u + 512 * v
        swapped = swap_bulk(x) + 8 * swap_bulk(y) + 64 * swap_bulk(u) + 512 * swap_bulk(v)
        assert conjugate(values[old]) == values[swapped]
    result = list(values)
    for leg in range(4):
        replacement = [0] * 4096
        for mask, polynomial in enumerate(result):
            if not polynomial:
                continue
            old = mask >> (3 * leg) & 7
            cleared = mask ^ (old << (3 * leg))
            for new in range(8):
                coefficient = B[old][new] if leg in (0, 2) else C[new][old]
                if coefficient:
                    replacement[cleared | (new << (3 * leg))] ^= scalar_multiply(polynomial, coefficient)
        result = replacement
    binary = []
    for polynomial in result:
        value = 0
        for j in range((polynomial.bit_length() + 1) // 2):
            coefficient = polynomial >> (2 * j) & 3
            assert coefficient in (0, 1)
            value |= coefficient << j
        binary.append(value)
    return binary


def coefficient_rank(values, degree, cut, binary=False):
    rows = [0] * 64
    for x, y, u, v in product(range(8), repeat=4):
        value = values[x + 8 * y + 64 * u + 512 * v]
        coefficient = (value >> degree & 1) if binary else (value >> (2 * degree) & 3)
        row, column = (8 * x + y, 8 * u + v) if cut == 0 else (8 * x + u, 8 * y + v)
        rows[row] |= coefficient << ((1 if binary else 2) * column)
    if binary:
        return binary_rank(rows)
    pivots = {}
    for row in rows:
        while row:
            column = ((row & -row).bit_length() - 1) // 2
            coefficient = row >> (2 * column) & 3
            row = scalar_multiply(row, (0, 1, 3, 2)[coefficient])
            if column not in pivots:
                pivots[column] = row
                break
            row ^= pivots[column]
    return len(pivots)


def endpoint_check(values, binary=False):
    for x, y, u, v in product(range(8), repeat=4):
        value = values[x + 8 * y + 64 * u + 512 * v]
        endpoint = (value.bit_count() & 1) if binary else evaluate_one(value)
        assert endpoint == int(x == v and y == u)


def digest(values):
    return sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed-directory", type=Path)
    args = parser.parse_args()
    records = []
    for indices in ((0, 1, 2), (0, 1, 3), (0, 2, 3)):
        values, failures = grouped_seed(indices)
        assert failures == 0
        endpoint_check(values)
        binary = descend(values)
        endpoint_check(binary, True)
        ranks = []
        for degree, cut in product((0, 1, 2, 12, 24, 48), (0, 1)):
            original = coefficient_rank(values, degree, cut)
            reduced = coefficient_rank(binary, degree, cut, True)
            assert original == reduced
            ranks.append({"degree": degree, "cut": cut, "rank": original})
        record = {"indices": indices, "ambient": 4, "alphabet": 8,
                  "exact_polynomial_divisions": 4096, "strict_endpoint_checks": 8192,
                  "degree_F4": max((v.bit_length() - 1) // 2 for v in values),
                  "degree_F2": max(v.bit_length() - 1 for v in binary),
                  "nonzero_F4_entries": sum(bool(v) for v in values),
                  "nonzero_F2_entries": sum(bool(v) for v in binary),
                  "F4_coefficients_sha256": digest(values),
                  "F2_coefficients_sha256": digest(binary),
                  "coefficient_rank_checks": ranks}
        records.append(record)
        if args.seed_directory:
            args.seed_directory.mkdir(parents=True, exist_ok=True)
            tag = "".join(map(str, indices))
            for field, coefficients in (("F4", values), ("F2", binary)):
                data = {"q": 8, "field": field, "values": coefficients,
                        "coefficient_encoding": "two-bit F4" if field == "F4" else "one-bit F2",
                        "indices": indices, "ambient": 4, "endpoint": 1,
                        "index_order": "x+8*y+64*u+512*v"}
                (args.seed_directory / f"conjugate-phase-{tag}-{field}.json").write_text(json.dumps(data) + "\n")
    _, wrong_sign = grouped_seed(sign=-1)
    assert wrong_sign == 29
    result = {"status": "exact polynomial crossings and fixed-field descent; no improved factoring bound",
              "residual_phase_checks": check_residual_phases(),
              "seeds": records, "wrong_gauge_sign_failed_divisions": wrong_sign,
              "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                  for name in ("tools/check_conjugate_residual_phase_codes.py",
                               "tools/check_schur_downset_crossing.py")}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
