#!/usr/bin/env python3
"""Exact exploratory p-typical Witt pairing and vertex vacuum calculations."""

from fractions import Fraction
from functools import lru_cache
from math import factorial
import argparse
import json
import time

from weighted_energy_pencil import rank_two


def artin_hasse_coefficients(N):
    coefficients = [Fraction(1)]
    for n in range(1, N + 1):
        total = Fraction(0)
        power = 1
        while power <= n:
            total += coefficients[n - power]
            power *= 2
        value = total / n
        assert value.denominator % 2 == 1
        coefficients.append(value)
    return coefficients


@lru_cache(None)
def monomials(energy, variables):
    if variables == 0:
        return ((),) if energy == 0 else ()
    weight = 1 << (variables - 1)
    return tuple(low + (a,) for a in range(energy // weight + 1)
                 for low in monomials(energy - a * weight, variables - 1))


def energy(monomial):
    return sum(a << j for j, a in enumerate(monomial))


def multiply_two(A, B, cutoff):
    result = set()
    for x in A:
        for y in B:
            z = tuple(a + b for a, b in zip(x, y))
            if energy(z) <= cutoff:
                if z in result:
                    result.remove(z)
                else:
                    result.add(z)
    return result


def vacuum_map(h, neutral=False):
    variables = h.bit_length()
    parity = [value.numerator & 1 for value in artin_hasse_coefficients(h)]
    all_monomials = [m for e in range(h + 1) for m in monomials(e, variables)]
    inverse_coherent = None
    if neutral:
        inverse_parity = [1]
        for n in range(1, h + 1):
            inverse_parity.append(sum(parity[k] * inverse_parity[n-k]
                                      for k in range(1, n + 1)) & 1)
        inverse_coherent = {m for m in all_monomials
                            if all(inverse_parity[a] for a in m)}
    H = [{m for m in monomials(e, variables)
          if all(parity[a] for a in m)} for e in range(h + 1)]
    basis = list(monomials(h, variables))
    positions = {m: i for i, m in enumerate(basis)}
    B = []
    for j in range(variables):
        n = 1 << j
        B.append({m for k in range(1, h + 1) if ((-k) >> j) & 1 for m in H[k]})
    rows = []
    for label in basis:
        value = {(0,) * variables}
        for j, occupation in enumerate(label):
            bit = 0
            while occupation:
                if occupation & 1:
                    factor = {tuple(a << bit for a in m) for m in B[j]
                              if energy(m) << bit <= h}
                    value = multiply_two(value, factor, h)
                occupation >>= 1
                bit += 1
        if neutral:
            length, bit = sum(label), 0
            while length:
                if length & 1:
                    factor = {tuple(a << bit for a in m) for m in inverse_coherent
                              if energy(m) << bit <= h}
                    value = multiply_two(value, factor, h)
                length >>= 1
                bit += 1
        rows.append(sum(1 << positions[m] for m in value if energy(m) == h))
    return basis, rows


def multiply_rational(A, B):
    result = {}
    for x, a in A.items():
        for y, b in B.items():
            z = tuple(i + j for i, j in zip(x, y))
            result[z] = result.get(z, Fraction(0)) + a * b
    return {m: c for m, c in result.items() if c}


def power_rational(A, exponent, variables):
    result = {(0,) * variables: Fraction(1)}
    while exponent:
        if exponent & 1:
            result = multiply_rational(result, A)
        exponent >>= 1
        if exponent:
            A = multiply_rational(A, A)
    return result


def witt_ghost_generators(variables):
    result = []
    for j in range(variables):
        key = tuple(int(i == j) for i in range(variables))
        polynomial = {key: Fraction(1)}
        for i in range(j):
            for m, c in power_rational(result[i], 1 << (j - i), variables).items():
                polynomial[m] = polynomial.get(m, Fraction(0)) - (1 << i) * c
        result.append({m: c / (1 << j) for m, c in polynomial.items() if c})
    return result


def hall_gram(h):
    variables = h.bit_length()
    generators = witt_ghost_generators(variables)
    basis = list(monomials(h, variables))
    transformed = []
    for label in basis:
        polynomial = {(0,) * variables: Fraction(1)}
        for j, a in enumerate(label):
            if a:
                polynomial = multiply_rational(polynomial,
                                              power_rational(generators[j], a, variables))
        transformed.append(polynomial)
    norms = {m: __import__('math').prod(factorial(a) * (1 << j) ** a
                                      for j, a in enumerate(m)) for m in basis}
    gram, rows = [], []
    for A in transformed:
        row = []
        for B in transformed:
            value = sum((a * B.get(m, 0) * norms[m] for m, a in A.items()), Fraction(0))
            assert value.denominator == 1
            row.append(value.numerator)
        gram.append(row)
        rows.append(sum((value & 1) << j for j, value in enumerate(row)))
    return basis, gram, rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h", type=int, default=16)
    parser.add_argument("--gram", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    started = time.monotonic()
    basis, vacuum = vacuum_map(args.h)
    result = {"h": args.h, "dimension": len(basis),
              "normal_ordered_derivative_vacuum_rank": rank_two(vacuum),
              "vacuum_rows_hex": [hex(x) for x in vacuum]}
    assert result["normal_ordered_derivative_vacuum_rank"] == len(basis)
    if args.gram:
        other, gram, rows = hall_gram(args.h)
        assert other == basis
        result["hall_gram_rank_mod2"] = rank_two(rows)
        assert result["hall_gram_rank_mod2"] == len(basis)
        result["hall_gram_rows_hex"] = [hex(x) for x in rows]
        result["hall_gram_integral"] = True
    result["elapsed_seconds"] = round(time.monotonic() - started, 3)
    if args.output:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    print(json.dumps({k: v for k, v in result.items() if not k.endswith("hex")}, indent=2))


if __name__ == "__main__":
    main()
