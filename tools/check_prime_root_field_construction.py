#!/usr/bin/env python3
"""Exact orbit/trace construction of a prime-order root over a binary subfield."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from check_carry_fourier_elimination import BinaryField, prime_factors
from check_frobenius_rational_trace import field_of_degree


COUNTS = Counter()


def trim(value):
    result = list(value)
    while result and result[-1] == 0:
        result.pop()
    return result


def add(left, right):
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] ^= value
    for index, value in enumerate(right):
        result[index] ^= value
    return trim(result)


def multiply(field, left, right):
    if not left or not right:
        return []
    result = [0] * (len(left) + len(right) - 1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i + j] ^= field.mul(value, other)
    return trim(result)


def remainder(field, left, right):
    value = trim(left)
    right = trim(right)
    if not right:
        raise ZeroDivisionError("zero polynomial modulus")
    inverse = field.inverse(right[-1])
    while len(value) >= len(right):
        scale = field.mul(value[-1], inverse)
        shift = len(value) - len(right)
        for index, coefficient in enumerate(right):
            value[shift + index] ^= field.mul(scale, coefficient)
        value = trim(value)
    return value


def gcd(field, left, right):
    left, right = trim(left), trim(right)
    while right:
        left, right = right, remainder(field, left, right)
    if left:
        inverse = field.inverse(left[-1])
        left = [field.mul(inverse, value) for value in left]
    return left


def power_mod(field, value, exponent, modulus):
    answer = [1]
    value = remainder(field, value, modulus)
    while exponent:
        if exponent & 1:
            answer = remainder(field, multiply(field, answer, value), modulus)
        value = remainder(field, multiply(field, value, value), modulus)
        exponent >>= 1
    return answer


def cyclic_square(field, vector):
    size = len(vector)
    result = [0] * size
    for index, value in enumerate(vector):
        result[(2 * index) % size] ^= field.mul(value, value)
    return result


def multiplicative_order(value, modulus):
    current = 1
    for order in range(1, modulus):
        current = current * value % modulus
        if current == 1:
            return order
    raise AssertionError("order did not divide the prime unit group")


def orbit_partition(prime, base):
    remaining = set(range(prime))
    orbits = []
    while remaining:
        start = min(remaining)
        orbit, value = [], start
        while value not in orbit:
            orbit.append(value)
            remaining.remove(value)
            value = base * value % prime
        assert value == start
        orbits.append(orbit)
    return orbits


def root_factor(field, prime):
    base, w = field.size, field.degree
    order = multiplicative_order(base % prime, prime)
    current = [1] * prime
    orbits = orbit_partition(prime, base)
    for orbit in orbits:
        original = [int(index in orbit) for index in range(prime)]
        raised = original[:]
        for _ in range(w):
            raised = cyclic_square(field, raised)
        assert raised == original
        COUNTS["exact_Frobenius_orbit_vectors"] += 1
        for bit in range(w):
            if len(current) - 1 == order:
                return current, order, orbits
            value = [(1 << bit) if index in orbit else 0 for index in range(prime)]
            trace = [0] * prime
            for _ in range(w):
                trace = [left ^ right for left, right in zip(trace, value)]
                value = cyclic_square(field, value)
            assert cyclic_square(field, trace) == trace
            COUNTS["exact_binary_trace_idempotents"] += 1
            divisor = gcd(field, current, trace)
            if 1 < len(divisor) < len(current):
                current = divisor
                COUNTS["proper_trace_gcd_splits"] += 1
    assert len(current) - 1 == order
    return current, order, orbits


def irreducibility_control(field, polynomial):
    degree = len(polynomial) - 1
    x = remainder(field, [0, 1], polynomial)
    powers = [x]
    for _ in range(degree):
        powers.append(power_mod(field, powers[-1], field.size, polynomial))
    assert powers[degree] == x
    for prime in prime_factors(degree):
        assert gcd(field, polynomial, add(powers[degree // prime], x)) == [1]
        COUNTS["independent_irreducibility_gcd_controls"] += 1
    COUNTS["independent_irreducibility_Frobenius_controls"] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    models = []
    for w in range(1, 6):
        field = BinaryField(0b11) if w == 1 else field_of_degree(w)
        for prime in [3, 5, 7, 11, 13, 17, 19, 31]:
            factor, order, orbits = root_factor(field, prime)
            assert len(factor) - 1 == order and (prime - 1) % order == 0
            assert all(len(orbit) == order for orbit in orbits if orbit != [0])
            assert remainder(field, [1] * prime, factor) == []
            irreducibility_control(field, factor)
            root = remainder(field, [0, 1], factor)
            assert root != [1]
            assert power_mod(field, root, prime, factor) == [1]
            COUNTS["exact_primitive_prime_root_controls"] += 1
            if order == 1:
                COUNTS["root_already_in_subfield_controls"] += 1
            # Constant elements retain their original multiplication in
            # the tower representation, including zero and all basis bits.
            for left in [0, 1] + [1 << bit for bit in range(w)]:
                for right in [0, 1] + [1 << bit for bit in range(w)]:
                    expected = trim([field.mul(left, right)])
                    assert remainder(field, multiply(field, [left], [right]), factor) == expected
                    COUNTS["exact_prescribed_subfield_products"] += 1
            models.append({"w": w, "Q": prime, "extension_degree": order,
                           "binary_degree": w * order,
                           "ground_polynomial": bin(field.polynomial),
                           "root_factor_coefficients": factor})
    result = {"status": "passed", "claim": "deterministic prescribed-subfield prime-root construction",
              "counts": dict(sorted(COUNTS.items())), "models": models,
              "scope": "Root-field preparation only; no improved whole-network contraction is claimed."}
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(json.dumps({key: value for key, value in result.items() if key != "models"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
