#!/usr/bin/env python3
"""Local odd-field bases and the tensor-stable characteristic-three certificate."""

import argparse
from functools import lru_cache
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial
import json
from pathlib import Path
import random

from check_partition_channels import rank_mod
from check_word_symmetry_barrier import log_bounds


def add_vector(basis, vector, prime):
    v = [x % prime for x in vector]
    for j in range(len(v) - 1, -1, -1):
        if not v[j]:
            continue
        if j in basis:
            value = v[j]
            v = [(a - value * b) % prime for a, b in zip(v, basis[j])]
        else:
            scale = pow(v[j], -1, prime)
            basis[j] = [(x * scale) % prime for x in v]
            return True
    return False


def local_basis(d, prime):
    n = d * d
    w = [int(i // d == i % d) for i in range(n)]
    vectors, kinds, independent = [], [], {}

    def append(v):
        assert add_vector(independent, v, prime)
        vectors.append(v)

    append(w)
    nilpotent = d % prime == 0
    if nilpotent:
        v = [0] * n
        v[-1] = 1
        append(v)
        kinds.append(("nilpotent", 2))
    else:
        kinds.append(("special", 1))
    for i in range(d - 1):
        v = [0] * n
        v[i * d + i], v[-1] = 1, -1
        if add_vector(independent, v, prime):
            vectors.append(v)
            kinds.append(("scalar", 1))
    for i in range(d):
        for j in range(i + 1, d):
            v = [0] * n
            v[i * d + j] = v[j * d + i] = 1
            append(v)
            kinds.append(("scalar", 1))
    for i in range(d):
        for j in range(i + 1, d):
            v = [0] * n
            v[i * d + j], v[j * d + i] = 1, -1
            append(v)
            kinds.append(("anti", 1))
    assert len(vectors) == n
    assert rank_mod([[v[i] for v in vectors] for i in range(n)], prime) == n
    checked = 0
    for cut in (False, True):
        for letter in range(3):  # r,h,c
            index = 2 - letter if cut else letter
            expected, start = [[0] * n for _ in range(n)], 0
            for kind, length in kinds:
                if kind == "nilpotent":
                    block = (
                        [[0, 1], [0, 0]]
                        if index == 0
                        else [[0, -1], [0, 0]] if index == 1 else [[1, 0], [0, 1]]
                    )
                elif kind == "special":
                    block = [[(d, -d, 1)[index]]]
                elif kind == "scalar":
                    block = [[int(index == 2)]]
                else:
                    block = [[(0, 2, -1)[index]]]
                for i in range(length):
                    for j in range(length):
                        expected[start + i][start + j] = block[i][j] % prime
                start += length
            for column, vector in enumerate(vectors):
                trace = sum(vector[i * d + i] for i in range(d)) % prime
                r = [trace * x % prime for x in w]
                c = [vector[(i % d) * d + i // d] % prime for i in range(n)]
                h = [(a - b - e) % prime for a, b, e in zip(vector, r, c)]
                actual = (r, h, c)[index]
                nonzero = [
                    (j, expected[j][column]) for j in range(n) if expected[j][column]
                ]
                wanted = [
                    sum(value * vectors[j][i] for j, value in nonzero) % prime
                    for i in range(n)
                ]
                assert actual == wanted
                checked += n
    return checked


def multiply(left, right):
    result = {}
    for a, value in left.items():
        for b, other in right.items():
            key = tuple(x + y for x, y in zip(a, b))
            result[key] = (result.get(key, 0) + value * other) % 3
    return {key: value for key, value in result.items() if value}


def power(polynomial, exponent):
    result = {(0, 0, 0): 1}
    while exponent:
        if exponent & 1:
            result = multiply(result, polynomial)
        polynomial = multiply(polynomial, polynomial)
        exponent >>= 1
    return result


@lru_cache(maxsize=None)
def rank_floor(t, k, j):
    return (
        comb(t, j)
        * 3**j
        * sum(
            comb(t - j, a) * 2 ** (a - k) * 4 ** (t - j - a)
            for a in range(k, t - j + 1)
        )
    )


def certificate_checks():
    # Coordinates are c, alpha=r-h, beta=2h-c in characteristic three.
    c = {(1, 0, 0): 1}
    h = {(1, 0, 0): 2, (0, 0, 1): 2}
    r = {(0, 1, 0): 1, (1, 0, 0): 2, (0, 0, 1): 2}
    p = multiply(multiply(r, c), power(h, 2))
    expected = {(4, 0, 0): 2, (1, 0, 3): 2, (3, 1, 0): 1, (2, 1, 1): 2, (1, 1, 2): 1}
    assert p == expected
    # Check the full variable involution, not only point evaluations.
    swapped = (
        r,
        {(1, 0, 0): 2, (0, 0, 1): 1},
        {(1, 0, 0): 2, (0, 1, 0): 2, (0, 0, 1): 2},
    )
    image = {}
    for monomial, coefficient in p.items():
        term = {(0, 0, 0): coefficient}
        for variable, exponent in zip(swapped, monomial):
            term = multiply(term, power(variable, exponent))
        for key, value in term.items():
            image[key] = (image.get(key, 0) + value) % 3
    assert {key: value for key, value in image.items() if value} == p
    assert sum(value for (_, k, _), value in p.items() if k == 0) % 3 == 1
    floors = [rank_floor(4, k, j) for _, k, j in expected]
    assert floors == [1296, 648, 520, 912, 540]
    inequalities = 0
    for t1 in range(1, 6):
        for t2 in range(1, 6):
            for k1 in range(t1 + 1):
                for j1 in range(t1 - k1 + 1):
                    for k2 in range(t2 + 1):
                        for j2 in range(t2 - k2 + 1):
                            assert rank_floor(t1 + t2, k1 + k2, j1 + j2) >= (
                                rank_floor(t1, k1, j1) * rank_floor(t2, k2, j2)
                            )
                            inequalities += 1
    supported = 0
    for m in range(1, 17):
        for (_, k, j), value in power(p, m).items():
            assert value and rank_floor(4 * m, k, j) >= 520**m
            supported += 1
    return {
        "polynomial_terms_c_alpha_beta": [
            [*key, value] for key, value in expected.items()
        ],
        "five_rank_floors": floors,
        "supermultiplicative_count_checks": inequalities,
        "repeated_certificate_monomials_checked": supported,
    }


def determinant_checks():
    cases = 0
    for d in (2, 3, 4):
        n, anti = d * d, d * (d - 1) // 2
        sym = d * (d + 1) // 2 - 1
        for k in range(1, 7):
            t = 2 * n * k
            lower = factorial(t) // (factorial((n + 1) * k) * factorial(k) ** (n - 1))
            for a in range(k, t - (sym + anti) * k + 1):
                for b in range(sym * k, t - a - anti * k + 1):
                    c = t - a - b
                    multiplicity = (
                        factorial(t)
                        // (factorial(a) * factorial(b) * factorial(c))
                        * sym**b
                        * anti**c
                    )
                    assert multiplicity >= lower
                    cases += 1
    for d in range(6, 1001):
        a = d * (d - 1) // 2
        assert 4 * (a * a - d**3) == d * d * (d * d - 6 * d + 1) > 0
    return cases


def matrix_add(left, right):
    return [[(a + b) % 3 for a, b in zip(x, y)] for x, y in zip(left, right)]


def kron(left, right):
    return [[a * b % 3 for a in x for b in y] for x in left for y in right]


def type_matrices(left, right):
    return [
        (
            kron(left[i], right[j])
            if i == j
            else matrix_add(kron(left[i], right[j]), kron(left[j], right[i]))
        )
        for i in range(3)
        for j in range(i, 3)
    ]


def combine(matrices, coefficients):
    return [
        [
            sum(c * matrix[i][j] for c, matrix in zip(coefficients, matrices)) % 3
            for j in range(len(matrices[0]))
        ]
        for i in range(len(matrices[0]))
    ]


def direct_characteristic_three_checks():
    n, d = 9, 3
    w = [int(i // d == i % d) for i in range(n)]
    r = [[a * b for b in w] for a in w]
    c = [[int(j == (i % d) * d + i // d) for j in range(n)] for i in range(n)]
    h = [[(int(i == j) - r[i][j] - c[i][j]) % 3 for j in range(n)] for i in range(n)]
    nilpotent = [[[0, 1], [0, 0]], [[0, 2], [0, 0]], [[1, 0], [0, 1]]]
    scalar, anti = [[[0]], [[0]], [[1]]], [[[0]], [[2]], [[2]]]
    rng = random.Random(520)
    coefficients = [[0] * 6, [1] * 6] + [
        [int(i == j) for i in range(6)] for j in range(6)
    ]
    coefficients += [[rng.randrange(3) for _ in range(6)] for _ in range(24)]
    monomials = []
    for i in range(3):
        for j in range(i, 3):
            monomial = [0] * 3
            monomial[i] += 1
            monomial[j] += 1
            monomials.append(tuple(monomial))
    cases = 0
    for column in (False, True):
        generators = (r, h, c) if not column else (c, h, r)
        actual_types = type_matrices(generators, generators)
        kinds = [(nilpotent, 1), (scalar, 4), (anti, 3)]
        blocks = []
        for i, (left, count_left) in enumerate(kinds):
            for j in range(i, len(kinds)):
                right, count_right = kinds[j]
                if column:
                    a, b = [left[k] for k in (2, 1, 0)], [right[k] for k in (2, 1, 0)]
                else:
                    a, b = left, right
                blocks.append(
                    (
                        type_matrices(a, b),
                        count_left * count_right * (1 if i == j else 2),
                    )
                )
        base = {(0, 0, 1): 1} if not column else {(1, 0, 0): 1}
        alpha = (
            {(1, 0, 0): 1, (0, 1, 0): 2} if not column else {(0, 0, 1): 1, (0, 1, 0): 2}
        )
        beta = (
            {(0, 1, 0): 2, (0, 0, 1): 2} if not column else {(0, 1, 0): 2, (1, 0, 0): 2}
        )
        for coefficient in coefficients:
            actual = rank_mod(combine(actual_types, coefficient), 3)
            expected = sum(
                count * rank_mod(combine(matrices, coefficient), 3)
                for matrices, count in blocks
            )
            assert actual == expected
            for k in range(3):
                for j in range(3 - k):
                    polynomial = multiply(
                        multiply(power(base, 2 - k - j), power(alpha, k)),
                        power(beta, j),
                    )
                    response = (
                        sum(
                            value * polynomial.get(monomial, 0)
                            for value, monomial in zip(coefficient, monomials)
                        )
                        % 3
                    )
                    if response:
                        assert actual >= rank_floor(2, k, j)
            cases += 1
    return cases


def decimal_interval(lower, upper):
    scale = 10**12
    a = lower.numerator * scale // lower.denominator
    b = -((-upper.numerator * scale) // upper.denominator)
    return {
        "lower": f"{a//scale}.{a%scale:012d}",
        "upper": f"{b//scale}.{b%scale:012d}",
    }


def rate_checks():
    intervals, rates = {}, {}
    for d in (2, 3, 4):
        n = d * d
        lo, hi = log_bounds(2 * n)
        qlo, qhi = log_bounds(n + 1)
        dlo, dhi = log_bounds(n)
        coefficient = Fraction(n + 1, 2 * n)
        lower, upper = (lo - coefficient * qhi) / dhi, (hi - coefficient * qlo) / dlo
        rates[f"semisimple_seed_{d}"] = lower, upper
    lo, hi = log_bounds(520)
    dlo, dhi = log_bounds(3)
    rates["characteristic_three_seed_three"] = lo / (8 * dhi), hi / (8 * dlo)
    lo, hi = log_bounds(10)
    dlo, dhi = log_bounds(5)
    rates["large_seed_five"] = lo / (2 * dhi), hi / (2 * dlo)
    gamma = rates["semisimple_seed_4"]
    for name, (lower, upper) in rates.items():
        assert 2 * lower * lower > 1
        if name != "semisimple_seed_4":
            assert lower > gamma[1]
        intervals[name] = decimal_interval(lower, upper)
    l2, u2 = log_bounds(2)
    l14, u14 = log_bounds(14)
    assert (2 + l14 / u2) / 8 > gamma[1]
    return intervals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    local = sum(local_basis(d, p) for p in (3, 5, 7) for d in range(2, 9))
    result = {
        "status": "all-field full pairing-type obstruction; no factoring lower bound",
        "local_scalar_nilpotent_entries_checked": local,
        "characteristic_three_certificate": certificate_checks(),
        "integer_determinant_sector_checks": determinant_checks(),
        "independent_characteristic_three_global_rank_checks": direct_characteristic_three_checks(),
        "rho_intervals": rate_checks(),
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
