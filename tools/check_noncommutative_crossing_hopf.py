#!/usr/bin/env python3
"""Symbolic audits of graded Hopf crossing factorizations and rank bounds."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from check_noncommutative_crossing_orbit import binary_rank, polynomial_mul


def add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            out[monomial] = out.get(monomial, 0) ^ coefficient
    return {m: c for m, c in out.items() if c}


def scale(polynomial, coefficient):
    return {m: polynomial_mul(c, coefficient) for m, c in polynomial.items()
            if polynomial_mul(c, coefficient)}


def multiply(left, right):
    out = {}
    for a, ca in left.items():
        for b, cb in right.items():
            if not a & b:
                out[a | b] = out.get(a | b, 0) ^ polynomial_mul(ca, cb)
    return {m: c for m, c in out.items() if c}


def variable(index):
    return {1 << index: 1}


def substitute(polynomial, images):
    result = {}
    for monomial, coefficient in polynomial.items():
        term = {0: coefficient}
        while monomial:
            bit = monomial & -monomial
            term = multiply(term, images[bit.bit_length() - 1])
            monomial ^= bit
        result = add(result, term)
    return result


def compose(left, right):
    return [substitute(image, left) for image in right]


def parameter_translate(polynomial):
    out = {}
    for monomial, coefficient in polynomial.items():
        result = 0
        while coefficient:
            bit = coefficient & -coefficient
            exponent = bit.bit_length() - 1
            sub = exponent
            while True:
                result ^= 1 << sub
                if not sub:
                    break
                sub = (sub - 1) & exponent
            coefficient ^= bit
        if result:
            out[monomial] = result
    return out


def evaluate_one(polynomial):
    return {m: 1 for m, c in polynomial.items() if c.bit_count() & 1}


def maps(m, t, columns):
    width, count = 2 * m + t, 4 * m + 2 * t
    identity = [variable(i) for i in range(count)]
    R, M, K = [], [], []
    for side in range(2):
        for j in range(width):
            own, other = variable(side * width + j), variable((1 - side) * width + j)
            if j < 2 * m:
                image = add(scale(own, 2), scale(other, 3))
                R.append(image)
                M.append(image)
                K.append(own)
                continue
            center = j - 2 * m
            linear = add(scale(own, 4), scale(other, 5))
            quadratic_R, quadratic_K = {}, {}
            for i, column in enumerate(columns):
                if not column >> center & 1:
                    continue
                xL, yL = variable(i), variable(m + i)
                xR, yR = variable(width + i), variable(width + m + i)
                terms = [multiply(xL, yL), multiply(xL, yR),
                         multiply(yL, xR), multiply(xR, yR)]
                rcoef = (0, 2, 4, 6) if side == 0 else (2, 0, 6, 4)
                kcoef = (10, 8, 14, 12) if side == 0 else (8, 10, 12, 14)
                quadratic_R = add(quadratic_R, *(scale(v, c) for v, c in zip(terms, rcoef)))
                quadratic_K = add(quadratic_K, *(scale(v, c) for v, c in zip(terms, kcoef)))
            R.append(add(linear, quadratic_R))
            M.append(linear)
            K.append(add(own, quadratic_K))
    assert compose(M, K) == compose(K, M) == R
    return identity, R, M, K


def normalized_triangular(m, t, columns, original_R):
    width, count = 2 * m + t, 4 * m + 2 * t
    # New linear coordinates are px,py,qx,qy,P,Q.
    change = []
    for component in range(2):
        for i in range(m):
            change.append(add(variable(component * m + i), variable(width + component * m + i)))
    for component in range(2):
        for i in range(m):
            change.append(variable(width + component * m + i))
    for j in range(t):
        change.append(add(variable(2 * m + j), variable(width + 2 * m + j)))
    for j in range(t):
        change.append(variable(width + 2 * m + j))
    inverse = []
    for side in range(2):
        for j in range(width):
            if j < 2 * m:
                image = variable(2 * m + j)
                if side == 0:
                    image = add(image, variable(j))
            else:
                image = variable(4 * m + t + j - 2 * m)
                if side == 0:
                    image = add(image, variable(4 * m + j - 2 * m))
            inverse.append(image)
    identity = [variable(i) for i in range(count)]
    assert [substitute(image, inverse) for image in change] == identity
    R = [substitute(substitute(image, original_R), inverse) for image in change]
    A, B, C = [], [], []
    for j in range(t):
        a = b = c = {}
        for i, column in enumerate(columns):
            if column >> j & 1:
                a = add(a, multiply(variable(i), variable(m + i)))
                b = add(b, multiply(variable(i), variable(3 * m + i)))
                c = add(c, multiply(variable(2 * m + i), variable(m + i)))
        A.append(a)
        B.append(b)
        C.append(c)
    expected = list(identity)
    for i in range(2 * m):
        expected[2 * m + i] = add(identity[2 * m + i], scale(identity[i], 3))
    for j in range(t):
        P, Q = identity[4 * m + j], identity[4 * m + t + j]
        expected[4 * m + j] = add(P, scale(A[j], 2))
        expected[4 * m + t + j] = add(Q, scale(P, 5), scale(A[j], 2),
                                           scale(B[j], 2), scale(C[j], 4))
    assert R == expected
    R_one = [evaluate_one(image) for image in R]
    assert compose(R_one, R_one) == identity
    F = compose([parameter_translate(image) for image in R], R_one)
    # Nonlinear coordinates Z=P+C and W=Q+P.
    nonlinear = list(identity)
    nonlinear_inverse = list(identity)
    for j in range(t):
        P, Q = identity[4 * m + j], identity[4 * m + t + j]
        nonlinear[4 * m + j] = add(P, C[j])
        nonlinear[4 * m + t + j] = add(Q, P)
        nonlinear_inverse[4 * m + j] = add(P, C[j])
        nonlinear_inverse[4 * m + t + j] = add(Q, P, C[j])
    assert compose(nonlinear_inverse, nonlinear) == identity
    triangular = [substitute(substitute(image, F), nonlinear_inverse) for image in nonlinear]
    expected = list(identity)
    for i in range(2 * m):
        expected[2 * m + i] = add(identity[2 * m + i], scale(identity[i], 2))
    for j in range(t):
        expected[4 * m + t + j] = add(identity[4 * m + t + j],
                                           scale(B[j], 2), scale(identity[4 * m + j], 4))
    assert triangular == expected
    return triangular


def coefficients(images):
    count = len(images)
    dimension = 1 << count
    degree = max((coefficient.bit_length() - 1 for image in images
                  for coefficient in image.values()), default=0) * count
    matrices = [[0] * dimension for _ in range(degree + 1)]
    values = [{0: 1}]
    for image in images:
        values += [multiply(value, image) for value in values]
    for source, polynomial in enumerate(values):
        for target, value in polynomial.items():
            while value:
                bit = value & -value
                matrices[bit.bit_length() - 1][source] ^= 1 << target
                value ^= bit
    while matrices and not any(matrices[-1]):
        matrices.pop()
    return matrices


def column_ranks(matrices, width):
    d = 1 << width
    ranks = []
    for matrix in matrices:
        rearranged = [0] * (d * d)
        for source, row in enumerate(matrix):
            left, right = source & (d - 1), source >> width
            while row:
                bit = row & -row
                target = bit.bit_length() - 1
                target_left, target_right = target & (d - 1), target >> width
                rearranged[left * d + target_left] ^= 1 << (right * d + target_right)
                row ^= bit
        ranks.append(binary_rank(rearranged))
    return ranks


def hilbert(m, t):
    n = 2 * m
    h = [0] * (n + 2 * t + 1)
    for a in range(n + 1):
        for b in range(t + 1):
            h[a + 2 * b] += math.comb(n, a) * math.comb(t, b)
    return h


def row_bounds(m, t):
    h = hilbert(m, t)
    D = len(h) - 1
    return [sum(min(h[a] * h[b], h[a + i] * h[b - i])
                for a in range(D - i + 1) for b in range(i, D + 1))
            for i in range(D + 1)]


def integer_poly_mul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def selection_counts(columns):
    independent = []
    for mask in range(1 << len(columns)):
        vectors = [v for i, v in enumerate(columns) if mask >> i & 1]
        if binary_rank(vectors) == len(vectors):
            independent.append(mask)
    counts = [0] * (len(columns) + 1)
    for left in independent:
        for right in independent:
            if not left & right:
                counts[(left | right).bit_count()] += 1
    return counts


def column_bound_polynomial(columns):
    counts = selection_counts(columns)
    result = [0] * (3 * len(columns) + 1)
    power = [1]
    for count in counts:
        for i, value in enumerate(power):
            result[i] += count * value
        power = integer_poly_mul(power, [2, 2, 2, 4])
    while result and not result[-1]:
        result.pop()
    return counts, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for m, t, columns in [(1, 1, (1,)), (2, 1, (1, 1)),
                          (3, 2, (1, 2, 3)), (3, 2, (1, 1, 2))]:
        identity, R, M, K = maps(m, t, columns)
        triangular = normalized_triangular(m, t, columns, R)
        counts, kappa = column_bound_polynomial(columns)
        record = {"horizontal_pairs": m, "center_bits": t, "columns": list(columns),
                  "symbolic_generator_factorizations": True,
                  "normalized_triangular_coordinates": True,
                  "row_upper_bounds": row_bounds(m, t),
                  "disjoint_independent_selection_counts": counts,
                  "translation_coefficient_Schmidt_upper_bounds": kappa}
        if 2 * m + t <= 5:
            original = coefficients([parameter_translate(image) for image in R])
            normal = coefficients(triangular)
            actual_rows = [binary_rank(matrix) for matrix in original]
            assert actual_rows == [binary_rank(matrix) for matrix in normal]
            assert all(a <= b for a, b in zip(actual_rows, record["row_upper_bounds"]))
            abelian = coefficients([parameter_translate(image) for image in M])
            abelian_rows = [binary_rank(matrix) for matrix in abelian]
            abelian_columns = column_ranks(abelian, 2 * m + t)
            assert abelian_columns == abelian_rows[::-1]
            actual_columns = column_ranks(original, 2 * m + t)
            D, dimension = 2 * m + 2 * t, 1 << (2 * m + t)
            bounds = [min(dimension * dimension, sum(
                kappa[a] * abelian_rows[D - i + a]
                for a in range(min(i, len(kappa) - 1) + 1)
                if 0 <= D - i + a <= D)) for i in range(D + 1)]
            assert all(a <= b for a, b in zip(actual_columns, bounds))
            record.update({"exact_row_ranks": actual_rows, "exact_column_ranks": actual_columns,
                           "abelian_row_ranks": abelian_rows,
                           "convolution_column_upper_bounds": bounds})
            translated_K = coefficients([parameter_translate(image) for image in K])
            exact_kappa = column_ranks(translated_K, 2 * m + t)
            assert all(a <= b for a, b in zip(exact_kappa, kappa))
            record["exact_translation_coefficient_Schmidt_ranks"] = exact_kappa
        records.append(record)
    result = {"status": "symbolic all-generator identities and finite exact rank-bound checks",
              "records": records,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
