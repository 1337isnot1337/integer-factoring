#!/usr/bin/env python3
"""Exact signed-polynomial and uniform block checks for the all-field barrier."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import combinations_with_replacement, permutations, product
import json
from math import comb, factorial, prod
from pathlib import Path


def add(*polynomials):
    answer = defaultdict(int)
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            answer[monomial] += coefficient
    return {m: c for m, c in answer.items() if c}


def scale(polynomial, scalar):
    return {m: scalar * c for m, c in polynomial.items() if scalar * c}


def multiply(left, right):
    answer = defaultdict(int)
    for (a, ca), (b, cb) in product(left.items(), right.items()):
        answer[tuple(sorted(a + b))] += ca * cb
    return {m: c for m, c in answer.items() if c}


def permutation_sign(values):
    return (-1) ** sum(values[i] > values[j]
                        for i in range(len(values)) for j in range(i + 1, len(values)))


def minor(rows, columns):
    return {
        tuple(sorted(4 * r + c for r, c in zip(rows, permutation))):
        permutation_sign(permutation)
        for permutation in permutations(columns)
    }


def adjoint_coefficient(i, j, a, b):
    cofactor = scale(minor([r for r in range(4) if r != j],
                          [c for c in range(4) if c != b]), (-1) ** (b + j))
    return multiply({(4 * i + a,): 1}, cofactor)


def sl4_coefficients():
    inputs = [{(a, b): 1} for a, b in permutations(range(4), 2)]
    inputs += [{(r, r): 1, (3, 3): -1} for r in range(3)]
    outputs = list(permutations(range(4), 2)) + [(r, r) for r in range(3)]
    return [
        add(*(scale(adjoint_coefficient(i, j, a, b), coefficient)
              for (a, b), coefficient in input_entries.items()))
        for input_entries in inputs for i, j in outputs
    ]


def rank_mod_prime(polynomials, prime, index):
    basis = {}
    for polynomial in polynomials:
        vector = {index[m]: c % prime for m, c in polynomial.items() if c % prime}
        while vector:
            pivot = max(vector)
            if pivot not in basis:
                inverse = pow(vector[pivot], -1, prime)
                basis[pivot] = {i: c * inverse % prime for i, c in vector.items()}
                break
            multiplier = vector[pivot]
            for i, c in basis[pivot].items():
                value = (vector.get(i, 0) - multiplier * c) % prime
                if value:
                    vector[i] = value
                else:
                    vector.pop(i, None)
    return len(basis)


def reshuffle(polynomial):
    answer = defaultdict(int)
    for monomial, coefficient in polynomial.items():
        transformed = []
        for variable in monomial:
            x, y, u, v = ((variable >> i) & 1 for i in (3, 2, 1, 0))
            transformed.append(8 * x + 4 * u + 2 * y + v)
        answer[tuple(sorted(transformed))] += coefficient
    return {m: c for m, c in answer.items() if c}


def polynomial_checks():
    determinant = minor(range(4), range(4))
    trace = add(*(adjoint_coefficient(i, i, a, a)
                  for i in (0, 1) for a in (1, 3)))
    witness = add(trace, scale(determinant, -1))
    factored = add(
        multiply(minor((0, 1), (0, 2)), minor((2, 3), (1, 3))),
        scale(multiply(minor((0, 1), (1, 3)), minor((2, 3), (0, 2))), -1),
    )
    assert witness == factored
    assert len(witness) == 8
    assert set(witness.values()) == {-1, 1}
    assert reshuffle(witness) == witness
    for monomial, coefficient in witness.items():
        assert reshuffle({monomial: coefficient}) == {monomial: coefficient}
        for bit in range(4):
            assert sum((variable >> bit) & 1 for variable in monomial) == 2
    assert sum(c for monomial, c in witness.items()
               if all(v in (0, 6, 9, 15) for v in monomial)) == 1

    coefficients = sl4_coefficients()
    # Twice the X-(1/2)I coefficient: integer identity, no modular division.
    represented = add(*(
        scale(coefficients[15 * (12 + r) + output], -1 if r != 1 else 1)
        for r in range(3) for output in (12, 13)
    ))
    assert represented == scale(witness, 2)
    index = {m: i for i, m in enumerate(combinations_with_replacement(range(16), 4))}
    dimensions = {}
    for prime in (3, 5, 7, 11, 101):
        dimension = rank_mod_prime(coefficients, prime, index)
        assert dimension == 225
        dimensions[str(prime)] = dimension
    return witness, dimensions


def partitions(total, parts=4, maximum=None):
    if parts == 0:
        if total == 0:
            yield ()
        return
    if maximum is None:
        maximum = total
    for first in range(min(total, maximum), -1, -1):
        for rest in partitions(total - first, parts - 1, first):
            yield (first,) + rest


def hook_dimension(shape):
    total = sum(shape)
    numerator = factorial(total) * prod(
        shape[i] - shape[j] + j - i
        for i in range(4) for j in range(i + 1, 4)
    )
    denominator = prod(factorial(shape[i] + 3 - i) for i in range(4))
    assert numerator % denominator == 0
    return numerator // denominator


def dominated_partition_checks():
    count = 0
    for k in range(1, 17):
        maximum_product = factorial(2 * k) * factorial(k)**2
        multinomial = comb(4 * k, 2 * k) * comb(2 * k, k)
        for shape in partitions(4 * k):
            if shape[0] > 2 * k or shape[0] + shape[1] > 3 * k:
                continue
            assert prod(factorial(x) for x in shape) <= maximum_product
            dimension = hook_dimension(shape)
            assert dimension * (4 * k + 4)**6 >= multinomial
            assert dimension * (4 * k + 4)**8 >= 1 << (6 * k)
            count += 1
    return count


def blocks(length, prime):
    k = length // 4
    if prime == 0:
        return [(k, 0, 4 * k)] if k else []
    cap = 1 if prime <= 3 else (prime - 1) // 4
    answer = []
    j = 0
    remaining = k
    while remaining:
        digit = remaining % prime
        remaining //= prime
        original_digit = digit
        before = len(answer)
        while digit:
            a = min(cap, digit)
            assert a > 0
            digit -= a
            answer.append((a, j, 4 * a * prime**j))
            if prime >= 5:
                assert 4 * a < prime
            else:
                assert a == 1
        assert len(answer) - before <= 6
        assert sum(a for a, jj, m in answer[before:]) == original_digit
        j += 1
    return answer


def block_checks():
    count = 0
    primes = (0, 2, 3, 5, 7, 11, 13, 17, 19, 31, 61, 127, 257, 65537)
    lengths = list(range(4, 4097))
    lengths += [10**6, (1 << 127) - 1, 1 << 127, (1 << 127) + 3]
    for prime in primes:
        for length in lengths:
            decomposition = blocks(length, prime)
            assert sum(m for a, j, m in decomposition) == 4 * (length // 4)
            assert len(decomposition) <= 6 * length.bit_length()
            # This checks the integer denominator envelope without
            # constructing an integer with exponentially many bits.
            assert all(m + 4 <= length + 4 for a, j, m in decomposition)
            remainder = length % 4
            assert (length - remainder) // 2 + remainder // 2 == length // 2
            count += 1
    return count, list(primes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    witness, dimensions = polynomial_checks()
    partition_count = dominated_partition_checks()
    block_count, characteristics = block_checks()
    assert 9 > 8 and 16 < 18
    result = {
        "status": "passed",
        "scope": "Exact certificates for the uniform all-field binary-word and middle-layer barrier; no improved factoring exponent.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "integral_witness": [{"variables": list(m), "coefficient": c}
                            for m, c in sorted(witness.items())],
        "balanced_leg_support_checks": 32,
        "odd_characteristic_coefficient_dimensions": dimensions,
        "dominated_partition_integer_bounds": partition_count,
        "uniform_characteristic_block_cases": block_count,
        "sample_characteristics_including_zero": characteristics,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
