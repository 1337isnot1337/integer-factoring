#!/usr/bin/env python3
"""Exact degree-sensitive shared-minor pressure certificates."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def balanced_polynomial(count, total):
    quotient, remainder = divmod(total, count)
    coefficients = [0]*(quotient+2)
    coefficients[quotient] = count-remainder
    coefficients[quotient+1] = remainder
    while len(coefficients) > 1 and not coefficients[-1]:
        coefficients.pop()
    return coefficients


def add(left, right):
    return [(left[i] if i < len(left) else 0)+(right[i] if i < len(right) else 0)
            for i in range(max(len(left), len(right)))]


def multiply(left, right):
    output = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i+j] += a*b
    return output


def profile_polynomial(q, degree, budget):
    n = q*q
    return add([1], add(balanced_polynomial(q-1, budget),
                        balanced_polynomial(n-q, n*degree-budget)))


def evaluate(polynomial, point):
    answer = Fraction(0)
    for coefficient in reversed(polynomial):
        answer = answer*point+coefficient
    return answer


def threshold(q, denominator=10**6):
    """Choose the least p/denominator strictly above q^(577/204)."""
    target = q**577 * denominator**204
    low, high = 0, q**3*denominator
    while low+1 < high:
        middle = (low+high)//2
        if middle**204 > target:
            high = middle
        else:
            low = middle
    assert high**204 > target >= (high-1)**204
    return Fraction(high, denominator)


def bernstein_coefficients(polynomial):
    degree = len(polynomial)-1
    return [sum((Fraction(coefficient*comb(k, j), comb(degree, j))
                 for j, coefficient in enumerate(polynomial[:k+1]) if coefficient), Fraction(0))
            for k in range(degree+1)]


def subdivide(coefficients):
    row = list(coefficients)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(a+b)/2 for a, b in zip(row, row[1:])]
        left.append(row[0])
        right.append(row[-1])
    return left, list(reversed(right))


def positive_certificate(polynomial):
    leaves = []
    pending = [(bernstein_coefficients(polynomial), 0, 0)]
    while pending:
        coefficients, depth, index = pending.pop()
        if all(value > 0 for value in coefficients):
            leaves.append({"depth": depth, "interval_index": index,
                           "minimum_bernstein_coefficient": str(min(coefficients))})
            continue
        if depth >= 28:
            raise ArithmeticError("The declared Bernstein depth did not certify positivity.")
        left, right = subdivide(coefficients)
        pending.append((right, depth+1, 2*index+1))
        pending.append((left, depth+1, 2*index))
    assert sum((Fraction(1, 1 << leaf["depth"]) for leaf in leaves), Fraction(0)) == 1
    return leaves


def certify_case(q, degree, bound):
    records = []
    maximum_degree = 0
    for budget in range(q*degree//2+1):
        polynomial = multiply(profile_polynomial(q, degree, budget),
                              profile_polynomial(q, degree, q*degree-budget))
        polynomial = [bound.denominator*value for value in polynomial]
        polynomial[degree] -= bound.numerator
        maximum_degree = max(maximum_degree, len(polynomial)-1)
        leaves = positive_certificate(polynomial)
        records.append({
            "prefix_budget": budget,
            "polynomial_sha256": sha256(json.dumps(polynomial, separators=(",", ":")).encode()).hexdigest(),
            "degree": len(polynomial)-1,
            "leaves": leaves,
        })
    return {"alphabet": q, "degree": degree, "pressure_threshold": str(bound),
            "maximum_bernstein_degree": maximum_degree, "budgets": records}


def finite_controls():
    balancing = 0
    for count in range(1, 5):
        for total in range(9):
            expected = evaluate(balanced_polynomial(count, total), Fraction(2, 3))
            for orders in product(range(total+1), repeat=count):
                if sum(orders) <= total:
                    assert sum((Fraction(2, 3)**order for order in orders), Fraction(0)) >= expected
                    balancing += 1
    monotonicity = affine = 0
    for q in range(2, 25):
        for degree in range(1, 13):
            for budget in range((q*degree+1)//2):
                first = budget//(q-1)
                second = (q*q*degree-budget-1)//(q*q-q)
                assert second >= first
                monotonicity += 1
        for point in (Fraction(1), Fraction(1, 2), Fraction(2, 3), Fraction(9, 10)):
            products = [evaluate(profile_polynomial(q, 1, budget), point)
                        * evaluate(profile_polynomial(q, 1, q-budget), point)/point
                        for budget in range(q+1)]
            expected = q**4*point + q**3*(1-point)**2 + (q-1)*(1-point)**4/point
            assert min(products) == expected >= q**3
            affine += 1
    return {"balanced_integer_allocations_checked": balancing,
            "budget_monotonicity_cases": monotonicity,
            "affine_exact_minimum_cases": affine}


def run():
    assert 2*577**2 > 816**2
    ranges = {7: range(2, 12), 8: range(2, 6), 9: range(2, 4), 10: range(2, 4), 11: (3,)}
    bounds = {q: threshold(q) for q in ranges}
    records = [certify_case(q, degree, bounds[q])
               for q, degrees in ranges.items() for degree in degrees]
    leaves = [leaf for case in records for budget in case["budgets"] for leaf in budget["leaves"]]
    return {
        "scope": "bare ungauged polynomial local Smith pressures, all fields and fixed auxiliary dimension",
        "all_alphabet_degree_one_floor": "rho >= 3/4",
        "finite_cases_rho_floor": "rho > 577/816 > 1/sqrt(2)",
        "exact_thresholds": {str(q): str(value) for q, value in bounds.items()},
        "finite_controls": finite_controls(),
        "alphabet_degree_cases": len(records),
        "budget_polynomials": sum(len(case["budgets"]) for case in records),
        "bernstein_leaves": len(leaves),
        "maximum_subdivision_depth": max(leaf["depth"] for leaf in leaves),
        "certificates": records,
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence" / "integer-shared-minor-floor.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "certificates"}, indent=2))
