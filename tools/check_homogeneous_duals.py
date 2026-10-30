#!/usr/bin/env python3
"""Exact checks for the homogeneous Appendix B duals.

The interpolation formula is checked against the independent integer
operator recursion in check_partition_channels.py. This does not test an
asymptotic factoring runtime or claim a stronger logarithmic exponent.
"""

from collections import Counter, defaultdict
from functools import lru_cache
import json
from math import comb, factorial, prod

from check_partition_channels import operator, partition_denominator, partitions


Polynomial = dict[tuple[int, ...], int]


def add_scaled(target: Polynomial, source: Polynomial, scale: int) -> None:
    for monomial, value in source.items():
        target[monomial] = target.get(monomial, 0) + scale * value
        if target[monomial] == 0:
            del target[monomial]


def shift(poly: Polynomial, direction: int) -> Polynomial:
    """L_1 lowers a mode; L_-1 raises a mode; each coefficient is its old mode."""
    if direction not in (-1, 1):
        raise ValueError("direction must be minus or plus one")
    result: Polynomial = defaultdict(int)
    for monomial, coefficient in poly.items():
        for mode, multiplicity in Counter(monomial).items():
            if direction == -1 and mode == 1:
                continue
            changed = list(monomial)
            changed.remove(mode)
            changed.append(mode + direction)
            result[tuple(sorted(changed))] += coefficient * mode * multiplicity
    return {key: value for key, value in result.items() if value}


def x_operator(poly: Polynomial) -> Polynomial:
    return shift(shift(poly, -1), 1)


def multiply_linear(poly: list[int], constant: int) -> list[int]:
    """Multiply an ordinary univariate polynomial by X+constant."""
    result = [0] * (len(poly) + 1)
    for degree, value in enumerate(poly):
        result[degree] += constant * value
        result[degree + 1] += value
    return result


@lru_cache(None)
def inverse_polynomial(h: int) -> tuple[int, tuple[int, ...]]:
    if h < 1:
        raise ValueError("positive weight required")
    eigen_x = [h * (h - 1) - j * (j - 1) for j in range(1, h + 1)]
    eigen_b = [comb(2 * h - 1, h - j) for j in range(1, h + 1)]
    scale = prod(eigen_b) * prod(eigen_x[j] - eigen_x[k]
                               for j in range(h) for k in range(j + 1, h))
    assert scale > 0
    coefficients = [0] * h
    for j in range(h):
        denominator = eigen_b[j] * prod(eigen_x[j] - eigen_x[k]
                                       for k in range(h) if k != j)
        assert scale % denominator == 0
        polynomial = [1]
        for k in range(h):
            if k != j:
                polynomial = multiply_linear(polynomial, -eigen_x[k])
        for degree, value in enumerate(polynomial):
            coefficients[degree] += scale // denominator * value
    for x, b in zip(eigen_x, eigen_b):
        assert sum(value * x**degree for degree, value in enumerate(coefficients)) * b == scale
    return scale, tuple(coefficients)


def apply_polynomial(coefficients: tuple[int, ...], monomial: tuple[int, ...]) -> Polynomial:
    result: Polynomial = {}
    for value in reversed(coefficients):
        result = x_operator(result)
        add_scaled(result, {monomial: 1}, value)
    return result


@lru_cache(None)
def homogeneous_dual(monomial: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], int], ...]:
    h = sum(monomial)
    _, coefficients = inverse_polynomial(h)
    sign = (-1) ** (h - len(monomial))
    return tuple(sorted((key, sign * value)
                        for key, value in apply_polynomial(coefficients, monomial).items()))


@lru_cache(None)
def dual_operator(monomial: tuple[int, ...], output: tuple[int, ...],
                  input_: tuple[int, ...]) -> int:
    return sum(coefficient * operator(term, output, input_)
               for term, coefficient in homogeneous_dual(monomial))


def piece_entry(h: int, s: int, x: tuple[int, ...], y: tuple[int, ...],
                u: tuple[int, ...], v: tuple[int, ...], exchanged: bool) -> int:
    scale = factorial(h) // partition_denominator(v)
    value = 0
    for i in range(2 * h + 1):
        coefficient = (-1) ** i * comb(2 * h, i) * scale
        weight = s + i if exchanged else 2 * h - s - i
        for intermediate in partitions(weight):
            if exchanged:
                term = dual_operator(x, intermediate, y) * operator(v, u, intermediate)
            else:
                term = dual_operator(x, u, intermediate) * operator(v, intermediate, y)
            value += coefficient * term
    return value


def verify() -> dict:
    pairing_cases = commutator_cases = 0
    weight_records = []
    for h in range(1, 9):
        basis = partitions(h)
        scale, polynomial = inverse_polynomial(h)
        height_bound = h * scale * (1 + h * h) ** (2 * h - 2)
        max_bits = 0
        for mu in basis:
            # Check the sl2 commutator independently on each basis vector.
            vector = {mu: 1}
            commutator = shift(shift(vector, 1), -1)
            add_scaled(commutator, shift(shift(vector, -1), 1), -1)
            assert commutator == {mu: 2 * h}
            commutator_cases += 1
            dual = homogeneous_dual(mu)
            assert all(sum(term) == h and len(term) == len(mu) for term, _ in dual)
            assert sum(abs(value) for _, value in dual) <= height_bound
            max_bits = max(max_bits, max((abs(value).bit_length() for _, value in dual), default=0))
            for nu in basis:
                b_scale = factorial(h) // partition_denominator(nu)
                actual = dual_operator(mu, (), nu) * b_scale
                assert actual == (scale * factorial(h) if mu == nu else 0)
                pairing_cases += 1
        assert scale.bit_length() <= h * (2 * h - 1) + h * (h - 1) * h.bit_length() + 1
        weight_records.append({
            "h": h,
            "dimension": len(basis),
            "scale": str(scale),
            "scale_bits": scale.bit_length(),
            "largest_dual_coefficient_bits": max_bits,
            "inverse_polynomial_degree": len(polynomial) - 1,
        })

    exchange_cases = crossing_cases = 0
    for h in range(1, 4):
        scale, _ = inverse_polynomial(h)
        basis = partitions(h)
        for x in basis:
            for y in basis:
                for u in basis:
                    for v in basis:
                        total = 0
                        for s in range(-2 * h, 2 * h + 1):
                            row = piece_entry(h, s, x, y, u, v, True)
                            column = piece_entry(h, s, x, y, u, v, False)
                            assert row == column
                            total += row
                            exchange_cases += 1
                        assert total == (scale * factorial(h) if x == v and y == u else 0)
                        crossing_cases += 1

    # Removing the grading sign breaks the prescribed positive dual pairing.
    h, mu = 2, (2,)
    scale, coefficients = inverse_polynomial(h)
    unsigned = sum(value * operator(term, (), mu)
                   for term, value in apply_polynomial(coefficients, mu).items())
    assert unsigned == -scale * factorial(h)
    # Omitting the Gram inverse gives a non-uniform diagonal, not the desired common scale.
    assert (-1) ** (2 - len(mu)) * operator(mu, (), mu) == 6
    assert operator((1, 1), (), (1, 1)) == 2

    return {
        "scope": "exact finite duality and crossing checks; all-weight and cost proofs are separate",
        "weights": weight_records,
        "duality_entries": pairing_cases,
        "sl2_commutator_basis_cases": commutator_cases,
        "order_exchange_entries": exchange_cases,
        "crossing_sum_entries": crossing_cases,
        "negative_controls": ["omitted grading sign", "omitted Gram inverse"],
        "new_factoring_exponent_claimed": False,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
