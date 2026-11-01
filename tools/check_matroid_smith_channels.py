#!/usr/bin/env python3
"""Exact basis-rotation gain over every scale of the previous word bound."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path

import check_cancellation_module_automata as fixed
import check_exact_module_orbits as exact
import check_minimal_boundary_modules as minimal
from matroid_smith_channels import MatroidSmithLattices, binary_plane_lattices, greedy_basis


def apply(matrix, vector):
    return tuple(sum(value * vector[target] for target, value in row.items())
                 for row in matrix)


def dense(matrix):
    return [[str(row.get(j, Fraction(0))) for j in range(len(matrix))] for row in matrix]


def poly_add(left, right):
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, Fraction(0)) + value
    return {key: value for key, value in result.items() if value}


def poly_mul(left, right):
    result = {}
    for (i, j), value in left.items():
        for (k, ell), other in right.items():
            key = i + k, j + ell
            result[key] = result.get(key, Fraction(0)) + value * other
    return {key: value for key, value in result.items() if value}


def characteristic_polynomial(quotient):
    """det(x I - T), with y=s^-1, using an independent permutation sum."""
    size = len(quotient)
    matrix = [[{} for _ in range(size)] for _ in range(size)]
    for i, entries in enumerate(quotient):
        matrix[i][i][1, 0] = Fraction(1)
        for (j, power), multiplicity in entries.items():
            assert power <= 0
            matrix[i][j] = poly_add(matrix[i][j], {(0, -power): Fraction(-multiplicity)})
    result = {}
    for permutation in permutations(range(size)):
        value = {(0, 0): Fraction(1)}
        for i, j in enumerate(permutation):
            value = poly_mul(value, matrix[i][j])
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(size) for j in range(i + 1, size))
        if inversions % 2:
            value = {key: -entry for key, entry in value.items()}
        result = poly_add(result, value)
    return result


def midpoint_split(coefficients):
    left, right = [coefficients[0]], [coefficients[-1]]
    while len(coefficients) > 1:
        coefficients = [(a + b) / 2 for a, b in zip(coefficients, coefficients[1:])]
        left.append(coefficients[0])
        right.append(coefficients[-1])
    return left, right[::-1]


def bernstein_positive(coefficients, left=Fraction(0), right=Fraction(1), depth=0):
    """An exact finite cover with strictly positive Bernstein coefficients."""
    assert depth <= 12, "positivity search did not terminate within its declared budget"
    if min(coefficients) > 0:
        return [{"interval": [str(left), str(right)],
                 "coefficients": list(map(str, coefficients)), "depth": depth}]
    assert max(coefficients) > 0, "polynomial is nonpositive on a subinterval"
    first, second = midpoint_split(coefficients)
    middle = (left + right) / 2
    return (bernstein_positive(first, left, middle, depth + 1)
            + bernstein_positive(second, middle, right, depth + 1))


def global_fixed_lower(local):
    """The previous exact word-pressure product exceeds 129/10 for all s>=1."""
    quotients, state_counts = [], []
    for cut, infinity, classes in ((0, True, (0, 1)),
                                  (1, False, (0, 1, 2, 2, 3, 1, 2))):
        factors = fixed.local_factors(local, cut, 0, (0, 0), infinity)
        representation = minimal.regularize(
            minimal.minimize(factors["channels"], (fixed.ONE, fixed.ZERO),
                             (fixed.ONE, fixed.ZERO)), factors["exponents"])
        orbit = exact.ExactModuleOrbit(representation["channels"])
        state_counts.append(len(orbit.states))
        quotient = exact.symbolic_quotient(orbit, representation["exponents"], classes)
        quotients.append(quotient)
    assert state_counts == [2, 7]
    assert quotients[0] == [{(1, 1): 4}, {(0, 0): 3, (1, 1): 1}]
    assert quotients[1] == [
        {(0, -1): 1, (1, 0): 1, (2, -1): 2},
        {(0, -1): 1, (2, -1): 3},
        {(1, 0): 1, (2, -1): 2, (3, -1): 1},
        {(1, 0): 1, (1, -1): 1, (2, -2): 2},
    ]
    polynomial = characteristic_polynomial(quotients[1])
    assert polynomial == {
        (4, 0): 1, (3, 1): -3, (2, 1): -4, (2, 2): 2,
        (2, 3): -2, (1, 3): -3, (1, 4): 2, (0, 3): 1, (0, 4): 3,
    }
    bound = Fraction(129, 10)
    # If h=rho(T_infinity), then h>=4 and s=h^2/(h+12).
    # Clear h's denominator from P(bound/h, (h+12)/h^2).
    denominator_power = max(i + 2 * j for i, j in polynomial)
    numerator = {}
    for (i, j), value in polynomial.items():
        shift = denominator_power - i - 2 * j
        for ell in range(j + 1):
            exponent = shift + ell
            numerator[exponent] = numerator.get(exponent, Fraction(0)) + (
                value * bound ** i * comb(j, ell) * 12 ** (j - ell))
    numerator = {key: value for key, value in numerator.items() if value}
    degree = max(numerator)
    assert denominator_power == 9 and degree == 6
    # Translate h=4+u, then u=t/(1-t).  Positivity of -P is enough:
    # det(xI-T)<0 forces a positive real eigenvalue larger than x.
    translated = [sum(value * comb(i, j) * 4 ** (i - j)
                      for i, value in numerator.items() if i >= j)
                  for j in range(degree + 1)]
    # Use the integer-scaled numerator displayed in the research note.
    coefficients = [-10000 * translated[i] / comb(degree, i) for i in range(degree + 1)]
    leaves = bernstein_positive(coefficients)
    assert len(leaves) == 5 and max(leaf["depth"] for leaf in leaves) == 4
    return {"strict_lower": str(bound), "scale_domain": "s >= 1",
            "states": state_counts,
            "symbolic_quotients": [[[j, power, multiplicity]
                                    for (j, power), multiplicity in sorted(row.items())]
                                   for quotient in quotients for row in quotient],
            "characteristic_polynomial": [[i, j, str(value)]
                                           for (i, j), value in sorted(polynomial.items())],
            "cleared_denominator_power": denominator_power,
            "displayed_numerator_scale": 10000,
            "numerator_coefficients": [str(numerator.get(i, 0)) for i in range(degree + 1)],
            "translated_coefficients": list(map(str, translated)),
            "projective_bernstein_coefficients": list(map(str, coefficients)),
            "positive_bernstein_cover": leaves}


def check(maximum_length=3):
    local, point = exact.companion_two(), Fraction(5)
    states = binary_plane_lattices(2)
    row = MatroidSmithLattices(local, 1, cut=0, infinity=True, seed_states=states)
    column = MatroidSmithLattices(local, 1, cut=1, infinity=False, seed_states=states)
    assert len(row.states) == len(column.states) == 10
    row_potential = (10000, 5320, 10000, 10000, 5320, 10000, 5320, 10000, 10000, 10000)
    col_potential = (10000, 9365, 6316, 6168, 9365, 6211, 9365, 6316, 5423, 6062)
    row_rate, col_rate = Fraction(266, 25), Fraction(181, 150)
    records, exhaustive_bases = [], 0
    for graph, potential, rate in ((row, row_potential, row_rate),
                                   (column, col_potential, col_rate)):
        certificate = graph.select(potential, point)
        image = apply(certificate["transfer"], potential)
        assert all(value <= rate * height for value, height in zip(image, potential))
        prefactor = potential[graph.start] * max(
            value / potential[i] for i, value in enumerate(certificate["terminal"]))
        for blocks in graph.blocks:
            for block in blocks:
                rank = len(block["lefts"])
                chosen = greedy_basis(block["candidates"], rank, potential)
                cost = sum(potential[target] for _, target in chosen)
                for alternative in combinations(block["candidates"], rank):
                    if len(fixed.binary_basis([mask for mask, _ in alternative])) == rank:
                        assert cost <= sum(potential[target] for _, target in alternative)
                        exhaustive_bases += 1
        records.append({"graph": graph.summary(), "certificate": certificate,
                        "potential": list(potential), "rate": str(rate),
                        "max_row_ratio": str(max(value / height for value, height in zip(image, potential))),
                        "prefactor": str(prefactor)})
    assert records[0]["prefactor"] == "1"
    assert records[1]["prefactor"] == "9365/6211"
    upper = row_rate * col_rate
    assert upper == Fraction(24073, 1875) < Fraction(129, 10)
    global_lower = global_fixed_lower(local)

    finite_records, literal_entries, rank_checks = [], 0, 0
    for length in range(1, maximum_length + 1):
        ranks, entries = fixed.coefficient_ranks(local, 1, length)
        literal_entries += entries
        words = []
        for graph, record in zip((row, column), records):
            certificate = record["certificate"]
            vector = tuple(certificate["terminal"])
            for _ in range(length):
                vector = apply(certificate["transfer"], vector)
            words.append(vector[graph.start])
            assert words[-1] <= Fraction(record["prefactor"]) * Fraction(record["rate"]) ** length
        for degree in range(length + 1):
            assert ranks[0][degree] <= (length - degree + 1) * point ** (-degree) * words[0]
            assert ranks[1][degree] <= (degree + 1) * point ** degree * words[1]
            rank_checks += 2
        finite_records.append({"length": length, "row_ranks": ranks[0],
                               "column_ranks": ranks[1], "word_sums": list(map(str, words))})
    for record in records:
        certificate = record["certificate"]
        certificate["transfer"] = dense(certificate["transfer"])
        certificate["terminal"] = list(map(str, certificate["terminal"]))
    return {"status": "pass", "point": str(point),
            "adaptive_product_upper": str(upper), "fixed_product_global_lower": global_lower,
            "row": records[0], "column": records[1],
            "exhaustive_competing_bases": exhaustive_bases,
            "literal_endpoint_entries": literal_entries, "literal_rank_inequalities": rank_checks,
            "finite_records": finite_records,
            "scope": ("Equal-order state-adaptive channel bases beat the previous fixed-factor "
                      "exact word-pressure product at every scale s>=1. This is not a new "
                      "factoring exponent or an optimum over all crossing representations.")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-length", type=int, default=3)
    args = parser.parse_args()
    assert 1 <= args.maximum_length <= 4
    report = check(args.maximum_length)
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    report["engine_sha256"] = sha256(Path(__file__).with_name("matroid_smith_channels.py").read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
