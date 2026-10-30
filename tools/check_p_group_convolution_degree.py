#!/usr/bin/env python3
"""Exact controls for the abelian p-group convolution degree lemma.

This checks the matrix-only extension.  It makes no four-leg crossing
identification for label groups other than elementary binary groups.
Only the Python standard library is required.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb, prod
from pathlib import Path
import random


def rank_mod(rows, p):
    basis = {}
    for source in rows:
        row = [x % p for x in source]
        while any(row):
            pivot = next(j for j, x in enumerate(row) if x)
            if pivot not in basis:
                inverse = pow(row[pivot], -1, p)
                basis[pivot] = [(x * inverse) % p for x in row]
                break
            scale = row[pivot]
            row = [(x - scale * y) % p for x, y in zip(row, basis[pivot])]
    return len(basis)


class PGroup:
    def __init__(self, p, exponents):
        self.p = p
        self.exponents = tuple(exponents)
        self.orders = tuple(p**e for e in exponents)
        self.points = tuple(product(*(range(q) for q in self.orders)))
        self.size = len(self.points)
        self.index = {x: i for i, x in enumerate(self.points)}
        self.weights = tuple(
            sum((Fraction(e * a, q - 1)
                 for e, a, q in zip(self.exponents, point, self.orders)), Fraction())
            for point in self.points
        )
        self.keys = tuple(zip(self.weights, self.points))
        # The rows are functions, and the columns are points of the group.
        self.binomial = tuple(
            tuple(prod(comb(x, a) for x, a in zip(point, exponent)) % p
                  for point in self.points)
            for exponent in self.points
        )
        self.addition = tuple(
            tuple(self.index[tuple((a + b) % q
                                   for a, b, q in zip(x, y, self.orders))]
                  for y in self.points)
            for x in self.points
        )

    def evaluate(self, coefficients):
        return tuple(sum(c * self.binomial[i][j]
                         for i, c in enumerate(coefficients)) % self.p
                     for j in range(self.size))

    def leading_index(self, coefficients):
        return max((i for i, c in enumerate(coefficients) if c % self.p),
                   key=self.keys.__getitem__)

    def convolution_rank(self, values):
        return rank_mod([[values[k] for k in row] for row in self.addition], self.p)

    def difference(self, values, coordinate):
        shift = [0] * len(self.orders)
        shift[coordinate] = 1
        row = self.addition[self.index[tuple(shift)]]
        return tuple((values[row[j]] - values[j]) % self.p
                     for j in range(self.size))


def digits(number, p, size):
    result = []
    for _ in range(size):
        number, digit = divmod(number, p)
        result.append(digit)
    assert not number
    return result


def check_group(p, exponents, rng):
    group = PGroup(p, exponents)
    d = group.size
    assert rank_mod(group.binomial, p) == d
    periodicity_checks = 0
    for q in group.orders:
        for a in range(q):
            for x in range(2 * q):
                assert (comb(x + q, a) - comb(x, a)) % p == 0
                periodicity_checks += 1

    difference_basis_checks = 0
    for i, exponent in enumerate(group.points):
        for coordinate in range(len(group.orders)):
            actual = group.difference(group.binomial[i], coordinate)
            if exponent[coordinate]:
                previous = list(exponent)
                previous[coordinate] -= 1
                expected = group.binomial[group.index[tuple(previous)]]
            else:
                expected = (0,) * d
            assert actual == expected
            difference_basis_checks += 1

    function_checks = 0
    derivative_checks = 0
    smallest_rank_margin = d
    # Enumerate every function in the declared groups, via its unique
    # binomial coefficient vector.  All groups here have at most nine points.
    for number in range(1, p**d):
        coefficients = digits(number, p, d)
        lead = group.leading_index(coefficients)
        exponent = group.points[lead]
        bound = prod(a + 1 for a in exponent)
        values = group.evaluate(coefficients)
        actual_rank = group.convolution_rank(values)
        assert actual_rank >= bound, (p, exponents, coefficients, exponent,
                                      actual_rank, bound)
        smallest_rank_margin = min(smallest_rank_margin, actual_rank - bound)
        function_checks += 1

        # The complete derivative family is directly checked on a reproducible
        # sample, including every individual binomial basis function.
        if sum(c != 0 for c in coefficients) == 1 or number % 251 == 1:
            derivatives = []
            for b in product(*(range(a + 1) for a in exponent)):
                actual = values
                for coordinate, count in enumerate(b):
                    for _ in range(count):
                        actual = group.difference(actual, coordinate)
                shifted = [0] * d
                for j, a in enumerate(group.points):
                    if all(x >= y for x, y in zip(a, b)):
                        smaller = tuple(x - y for x, y in zip(a, b))
                        shifted[group.index[smaller]] = coefficients[j]
                assert actual == group.evaluate(shifted)
                expected_lead = tuple(x - y for x, y in zip(exponent, b))
                assert group.points[group.leading_index(shifted)] == expected_lead
                derivatives.append(actual)
                derivative_checks += 1
            assert len(derivatives) == bound
            assert rank_mod(derivatives, p) == bound

    # The exact rational weights suppress their common log_2(p) factor.
    expected_degree_sum = Fraction(d * sum(exponents), 2)
    assert sum(group.weights) == expected_degree_sum
    basis_checks = 0
    for _ in range(100):
        coefficients = [[int(i == j) for j in range(d)] for i in range(d)]
        for _ in range(10 * d):
            a, b = rng.sample(range(d), 2)
            scalar = rng.randrange(1, p)
            coefficients[a] = [(x + scalar * y) % p
                               for x, y in zip(coefficients[a], coefficients[b])]
        assert rank_mod(coefficients, p) == d
        total_weight = sum(group.weights[group.leading_index(row)]
                           for row in coefficients)
        assert total_weight >= expected_degree_sum
        basis_checks += 1

    return {
        "characteristic": p,
        "cyclic_orders": list(group.orders),
        "group_size": d,
        "periodicity_checks": periodicity_checks,
        "difference_basis_checks": difference_basis_checks,
        "all_nonzero_functions_checked": function_checks,
        "direct_derivative_checks": derivative_checks,
        "random_function_bases_checked": basis_checks,
        "monomial_basis_total_weight": str(expected_degree_sum),
        "minimum_convolution_rank_minus_product_bound": smallest_rank_margin,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(20261030)
    groups = [(2, (2,)), (2, (3,)), (3, (1,)), (3, (2,)),
              (2, (1, 2)), (3, (1, 1))]
    results = [check_group(p, exponents, rng) for p, exponents in groups]
    report = {
        "status": "PASS",
        "scope": "Abelian p-group convolution-matrix lemma in matching characteristic; no nonbinary crossing claim",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "groups": results,
        "totals": {key: sum(x[key] for x in results) for key in (
            "periodicity_checks", "difference_basis_checks",
            "all_nonzero_functions_checked", "direct_derivative_checks",
            "random_function_bases_checked")},
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
