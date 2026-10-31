#!/usr/bin/env python3
"""An exact one-parameter q4 crossing and two polynomial Wick obstructions.

No floating-point arithmetic or native library is used. Polynomial masks
represent coefficients in F2[w]. Optional F4 exports use two bit planes.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_schur_downset_crossing import binomial_polynomial, sparse_kernel
from weighted_smith_oracle import polynomial_multiply, polynomial_exact_divide


ROOT = Path(__file__).resolve().parent.parent


def parameter_basis():
    """Return A,B,C such that the even-core contraction is A+lambda B+lambda^2 C."""
    n = 4
    kernel = sparse_kernel(tuple(range(n)), n)
    matrix = [row[:2 * n] for row in kernel[2 * n:]]

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        answer = 0
        remaining = columns
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length() - 1
            answer ^= polynomial_multiply(matrix[i][j], minor(rows ^ bit, columns ^ cb))
        return answer

    basis = [[0] * 256 for _ in range(3)]
    quartic_rows = (1 << 3) | (1 << 7)
    quartic_columns = 1 | (1 << 4)
    for labels in product(range(4), repeat=4):
        target = sum(a << (2 * leg) for leg, a in enumerate(labels))
        for states in product((0, 1), repeat=4):
            old = [((a & 1) << 1) | ((a >> 1) << 2) | (9 if both else 0)
                   for a, both in zip(labels, states)]
            x, y, u, v = old
            rows, columns = v | (y << n), u | (x << n)
            numerator = polynomial_multiply(5, minor(rows, columns))
            if rows & quartic_rows == quartic_rows and columns & quartic_columns == quartic_columns:
                numerator ^= minor(rows ^ quartic_rows, columns ^ quartic_columns)
            basis[states[1] + states[3]][target] ^= numerator
    # The binary polynomial 5 is c^2 = 1+w^2.
    return [[polynomial_exact_divide(value, 5) for value in row] for row in basis]


def quartic_relation(basis, mask):
    bits = [1 << i for i in range(8) if mask >> i & 1]
    assert len(bits) == 4
    a, b, c, d = bits
    result = [0] * 5
    for left, right in ((0, mask), (a | b, c | d), (a | c, b | d), (a | d, b | c)):
        for i, j in product(range(3), repeat=2):
            result[i + j] ^= polynomial_multiply(basis[i][left], basis[j][right])
    return result


def gaussian_relation_failures(values):
    """Check every even Pfaffian identity beyond its two-variable coefficients."""
    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        bit = mask & -mask
        rest = mask ^ bit
        remaining = rest
        answer = 0
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            answer ^= polynomial_multiply(values[bit | cb], pfaffian(rest ^ cb))
        return answer

    failures = []
    checks = 0
    for mask in range(256):
        weight = mask.bit_count()
        if weight % 2 or weight < 4:
            continue
        value = values[mask]
        for _ in range(weight // 2 - 1):
            value = polynomial_multiply(value, values[0])
        if value != pfaffian(mask):
            failures.append(mask)
        checks += 1
    assert checks == 99
    return failures


def field_multiply(a, b):
    p = (a & 1) * (b & 1)
    q = (a >> 1) * (b >> 1)
    r = ((a & 1) ^ (a >> 1)) * ((b & 1) ^ (b >> 1))
    return (p ^ q) | ((r ^ p) << 1)


def specialize(basis, parameter):
    square = field_multiply(parameter, parameter)
    result = []
    for values in zip(*basis):
        pair = [0, 0]
        for polynomial, coefficient in zip(values, (1, parameter, square)):
            if coefficient & 1:
                pair[0] ^= polynomial
            if coefficient & 2:
                pair[1] ^= polynomial
        result.append(pair)
    return result


def digest(value):
    return sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def check(seed_directory=None):
    basis = parameter_basis()
    degrees = [max(value.bit_length() - 1 for value in row) for row in basis]
    assert degrees == [17, 32, 14]
    endpoint_checks = 0
    parity_checks = 0
    for power, row in enumerate(basis):
        for mask, polynomial in enumerate(row):
            x, y, u, v = (mask >> (2 * leg) & 3 for leg in range(4))
            assert polynomial.bit_count() % 2 == int(power == 0 and x == v and y == u)
            endpoint_checks += 1
            if mask.bit_count() % 2:
                assert polynomial == 0
                parity_checks += 1

    r23 = quartic_relation(basis, 23)
    r83 = quartic_relation(basis, 83)
    w2c4 = binomial_polynomial(4) << 2
    w2c7 = binomial_polynomial(7) << 2
    assert r23 == [0, 0, 0, w2c4 << 8, w2c4]
    assert r83 == [0, 0, 0, w2c7, w2c7]
    gaussian_control = gaussian_relation_failures(basis[0])
    binary_one = [a ^ b ^ c for a, b, c in zip(*basis)]
    non_gaussian_control = gaussian_relation_failures(binary_one)
    assert basis[0][0] == 1 and not gaussian_control
    assert binary_one[0] == 1 << 16
    assert len(non_gaussian_control) == 68
    assert r23[3] ^ r23[4] == binomial_polynomial(12) << 2

    specializations = []
    for parameter in range(4):
        values = specialize(basis, parameter)
        for mask, pair in enumerate(values):
            x, y, u, v = (mask >> (2 * leg) & 3 for leg in range(4))
            endpoint = (pair[0].bit_count() % 2) | ((pair[1].bit_count() % 2) << 1)
            assert endpoint == int(x == v and y == u)
        record = {
            "parameter": parameter,
            "degree": max(max(pair).bit_length() - 1 for pair in values),
            "coefficient_sha256": digest(values),
            "endpoint_entries_checked": 256,
        }
        specializations.append(record)
        if seed_directory:
            document = {"q": 4, "field": "F4", "coefficient_planes": "lo+omega*hi",
                        "parameter": parameter, "strict_endpoint": True,
                        "denominator_degree": 0, "values": values}
            (seed_directory / f"even-core-parameter-{parameter}-f4.json").write_text(
                json.dumps(document) + "\n")
            if parameter in (0, 1):
                assert all(not pair[1] for pair in values)
                document.update(field="F2", values=[pair[0] for pair in values])
                (seed_directory / f"even-core-parameter-{parameter}-binary.json").write_text(
                    json.dumps(document) + "\n")
    if seed_directory:
        (seed_directory / "even-core-parameter-basis.json").write_text(
            json.dumps({"q": 4, "parameter_coefficients": basis,
                        "parameter_degree": 2, "degrees": degrees}) + "\n")

    return {
        "claim": "A+lambda*B+lambda^2*C is a strict q4 polynomial crossing",
        "field": "every characteristic-two field",
        "parameter_coefficients_sha256": digest(basis),
        "parameter_coefficient_degrees": degrees,
        "exact_c_squared_divisions": 768,
        "coefficient_endpoint_checks": endpoint_checks,
        "odd_parity_zeros_checked": parity_checks,
        "wick_relations": {
            "mask23_lambda_coefficients": r23,
            "mask23_factor": "lambda^3*w^2*(1+w)^4*(w^8+lambda)",
            "mask83_lambda_coefficients": r83,
            "mask83_factor": "lambda^3*w^2*(1+w)^7*(1+lambda)",
            "lambda1_mask23_factor": "w^2*(1+w)^12",
        },
        "gaussian_control_lambda0_relations_checked": 99,
        "gaussian_control_lambda0_failures": len(gaussian_control),
        "lambda1_relations_checked": 99,
        "lambda1_relation_failures": len(non_gaussian_control),
        "constant_specializations": specializations,
        "scope": "No pressure or factoring-exponent improvement is asserted.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence/even-core-postselection.json")
    parser.add_argument("--seed-directory", type=Path)
    args = parser.parse_args()
    if args.seed_directory:
        args.seed_directory.mkdir(parents=True, exist_ok=True)
    report = check(args.seed_directory)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
