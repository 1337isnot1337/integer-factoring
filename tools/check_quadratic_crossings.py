#!/usr/bin/env python3
"""Independent truncated-series replay and exact all-scale catalogue certificate."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import comb
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def matrices(mask):
    infinity, zero = [[0] * 4 for _ in range(4)], [[0] * 4 for _ in range(4)]
    for x, y, u, v in product(range(2), repeat=4):
        identity, crossing = int(x == u and y == v), int(x == v and y == u)
        middle = mask >> (8*x + 4*y + 2*u + v) & 1
        leading = identity ^ crossing ^ middle
        infinity[2*x+y][2*u+v] = leading | middle << 1 | identity << 2
        zero[2*x+u][2*y+v] = identity | middle << 1 | leading << 2
    return infinity, zero


def binary_rank(rows):
    basis = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit not in basis:
                basis[bit] = row
                break
            row ^= basis[bit]
    return len(basis)


def orders(matrix):
    """Smith orders from F2 ranks of block Toeplitz maps modulo t^k."""
    previous_rank = previous_count = 0
    result = []
    for precision in range(1, 10):
        rows = []
        for output_degree in range(precision):
            for row in matrix:
                packed = 0
                for input_degree in range(max(0, output_degree - 2), output_degree + 1):
                    for j, entry in enumerate(row):
                        if entry >> (output_degree - input_degree) & 1:
                            packed |= 1 << (4*input_degree + j)
                rows.append(packed)
        rank = binary_rank(rows)
        count = rank - previous_rank
        result.extend([precision - 1] * (count - previous_count))
        if count == 4:
            return tuple(result)
        previous_count, previous_rank = count, rank
    raise AssertionError("endpoint implies determinant degree at most eight")


def convolution(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def profile_polynomial(key):
    polynomials = []
    for side in key:
        coefficients = [0] * (max(side) + 1)
        for order in side:
            coefficients[order] += 1
        polynomials.append(coefficients)
    return convolution(*polynomials)


def divide_exact(poly, factor):
    work = list(map(Fraction, poly))
    quotient = [Fraction(0)] * (len(work) - len(factor) + 1)
    for degree in range(len(work) - len(factor), -1, -1):
        coefficient = work[degree + len(factor) - 1] / factor[-1]
        quotient[degree] = coefficient
        for j, value in enumerate(factor):
            work[degree+j] -= coefficient * value
    assert not any(work)
    return quotient


def bernstein_certificate(poly):
    """Exact nonnegativity on [0,1] by rational Bernstein subdivisions."""
    degree = len(poly) - 1
    bernstein = [sum(Fraction(poly[j]) * Fraction(comb(k, j), comb(degree, j))
                     for j in range(k+1)) for k in range(degree+1)]
    pending, leaves, max_depth = [(bernstein, 0)], 0, 0
    while pending:
        coefficients, depth = pending.pop()
        if min(coefficients) >= 0:
            leaves += 1
            max_depth = max(max_depth, depth)
            continue
        assert depth < 32, "no certificate found"
        triangle = [coefficients]
        while len(triangle[-1]) > 1:
            last = triangle[-1]
            triangle.append([(a+b)/2 for a, b in zip(last, last[1:])])
        left = [row[0] for row in triangle]
        right = [row[-1] for row in reversed(triangle)]
        pending.extend(((left, depth+1), (right, depth+1)))
    return {"leaves": leaves, "max_depth": max_depth}


def certify(catalogue):
    keys = [tuple(tuple(row[side]) for side in ("infinity", "zero")) for row in catalogue]
    pareto = [key for key in keys if not any(
        other != key and all(a >= b for a, b in zip(sum(other, ()), sum(key, ())))
        for other in keys)]
    certificates = []
    for key in pareto:
        poly = [27 * value for value in profile_polynomial(key)]
        poly[2] -= 256
        # Exact equality cases are the elementary pencil and its z -> z^2 pullback.
        equality = None
        if key == ((1, 2, 2, 3), (0, 1, 1, 2)):
            poly = divide_exact(poly, [1, -6, 9])
            equality = "x=1/3"
        elif key == ((0, 2, 2, 4), (0, 2, 2, 4)):
            poly = divide_exact(poly, [1, 0, -6, 0, 9])
            equality = "x^2=1/3"
        certificates.append({"profile": key, "equality": equality,
                             **bernstein_certificate(poly)})
    return certificates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="independently enumerate all 65,536 tensors")
    args = parser.parse_args()
    data = json.loads((ROOT / "evidence/quadratic-crossing-catalogue.json").read_text())
    catalogue = data["profiles"]
    assert sum(row["count"] for row in catalogue) == 65536
    expected = {(tuple(row["infinity"]), tuple(row["zero"])): row for row in catalogue}
    assert len(expected) == len(catalogue) == 162
    rng = random.Random(6546458276)
    masks = range(65536) if args.full else sorted(
        {row["first_mask"] for row in catalogue} | {rng.randrange(65536) for _ in range(256)})
    counts, first = Counter(), {}
    for mask in masks:
        key = tuple(orders(matrix) for matrix in matrices(mask))
        assert key in expected, (mask, key)
        counts[key] += 1
        first.setdefault(key, mask)
        if mask == expected[key]["first_mask"]:
            assert key == (tuple(expected[key]["infinity"]), tuple(expected[key]["zero"]))
    if args.full:
        assert all(counts[key] == row["count"] and first[key] == row["first_mask"]
                   for key, row in expected.items())
    certificates = certify(catalogue)
    assert len(certificates) == 18
    print(json.dumps({"checked_tensors": len(masks), "profiles": len(catalogue),
                      "all_scale_lower_bound": "256/27", "certificates": certificates}, indent=2))


if __name__ == "__main__":
    main()
