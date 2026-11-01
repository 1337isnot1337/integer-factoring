#!/usr/bin/env python3
"""Independent exact checks of Cauchy channels, minors, and mode counts."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import isqrt, log2, prod
from pathlib import Path
from time import monotonic

from cauchy_two_adic_smith import (
    cauchy_matrix, count_nonzero_channels, rational_ldu, residue,
    smith_exponent, valuation2, weighted_monomial_count,
)


def determinant(matrix):
    work = [row.copy() for row in matrix]
    value = Fraction(1)
    for k in range(len(work)):
        pivot = next((i for i in range(k, len(work)) if work[i][k]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k:
            work[k], work[pivot] = work[pivot], work[k]
            value = -value
        value *= work[k][k]
        for i in range(k + 1, len(work)):
            ratio = work[i][k] / work[k][k]
            for j in range(k + 1, len(work)):
                work[i][j] -= ratio * work[k][j]
    return value


def chain_ring_profile(matrix, precision):
    """Independent min-valuation row/column elimination modulo 2^precision."""
    modulus = 1 << precision
    work = [[residue(x, precision) for x in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    profile = []
    for k in range(min(rows, columns)):
        candidates = [(valuation2(work[i][j]), i, j)
                      for i in range(k, rows) for j in range(k, columns)
                      if work[i][j]]
        if not candidates:
            break
        exponent, i, j = min(candidates)
        work[k], work[i] = work[i], work[k]
        for row in work:
            row[k], row[j] = row[j], row[k]
        smaller_modulus = 1 << (precision - exponent)
        inverse = pow(work[k][k] >> exponent, -1, smaller_modulus)
        for i in range(k + 1, rows):
            assert work[i][k] % (1 << exponent) == 0
            multiplier = ((work[i][k] >> exponent) * inverse) % smaller_modulus
            for j in range(k, columns):
                work[i][j] = (work[i][j] - multiplier * work[k][j]) % modulus
            assert work[i][k] == 0
        for j in range(k + 1, columns):
            assert work[k][j] % (1 << exponent) == 0
            multiplier = ((work[k][j] >> exponent) * inverse) % smaller_modulus
            for i in range(k, rows):
                work[i][j] = (work[i][j] - multiplier * work[i][k]) % modulus
            assert work[k][j] == 0
        profile.append(exponent)
    assert profile == sorted(profile)
    return profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    started = monotonic()
    identity_entries = exact_minors = modular_profiles = modular_entries = 0

    for rows, columns in [(1, 1), (2, 5), (5, 2), (6, 6), (9, 7), (7, 9)]:
        for c in [-9, -3, -1, 1, 3, 7, 15]:
            matrix = cauchy_matrix(rows, columns, c)
            left, diagonal, right = rational_ldu(rows, columns, c)
            for k, d in enumerate(diagonal):
                assert d.denominator % 2 == 1
                assert valuation2(d.numerator) == smith_exponent(k)
                assert left[k][k] == right[k][k] == 1
            assert all(x.denominator % 2 for row in left + right for x in row)
            for a in range(rows):
                for b in range(columns):
                    assert sum((left[a][k] * diagonal[k] * right[k][b]
                                for k in range(len(diagonal))), Fraction(0)) == matrix[a][b]
                    identity_entries += 1

    minor_profiles = []
    for c in [-5, 1, 11]:
        matrix = cauchy_matrix(6, 6, c)
        profile = []
        for size in range(1, 7):
            minimum = None
            for selected_rows in combinations(range(6), size):
                for selected_columns in combinations(range(6), size):
                    minor = [[matrix[a][b] for b in selected_columns] for a in selected_rows]
                    value = determinant(minor)
                    assert value and value.denominator % 2
                    v = valuation2(value.numerator)
                    minimum = v if minimum is None else min(minimum, v)
                    exact_minors += 1
            assert minimum == sum(smith_exponent(k) for k in range(size))
            profile.append(minimum)
        minor_profiles.append({"odd_shift": c, "determinantal_valuations": profile})

    for dimension in [1, 2, 4, 7, 12, 24]:
        for c in [-7, 1, 9]:
            matrix = cauchy_matrix(dimension, dimension, c)
            left, diagonal, right = rational_ldu(dimension, dimension, c)
            for precision in [1, 2, 3, 7, 16, 40, 80]:
                expected = [smith_exponent(k) for k in range(dimension)
                            if smith_exponent(k) < precision]
                assert chain_ring_profile(matrix, precision) == expected
                assert count_nonzero_channels(dimension, precision) == len(expected)
                modular_profiles += 1
                modulus = 1 << precision
                dl = [[residue(left[a][k] * diagonal[k], precision)
                       for k in range(dimension)] for a in range(dimension)]
                ur = [[residue(right[k][b], precision) for b in range(dimension)]
                      for k in range(dimension)]
                for a in range(dimension):
                    for b in range(dimension):
                        value = sum(dl[a][k] * ur[k][b] for k in range(dimension)) % modulus
                        assert value == residue(matrix[a][b], precision)
                        modular_entries += 1

    # Independent explicit tuple enumeration for small weighted mode boxes.
    def literal_count(weights, budget, index=0):
        if index == len(weights):
            return 1
        weight = weights[index]
        return sum(literal_count(weights, budget - weight * count, index + 1)
                   for count in range(budget // weight + 1))

    tuple_counts = 0
    for dimension in range(1, 5):
        for copies in [1, 2]:
            for shift in [1, 3]:
                weights = [shift + 4 * k - 2 * k.bit_count()
                           for k in range(dimension) for _ in range(copies)]
                for budget in [0, 1, 3, 8, 17]:
                    assert weighted_monomial_count(dimension, copies, shift, budget) == literal_count(weights, budget)
                    tuple_counts += 1

    lower_box_checks = 0
    for u in range(1, 65):
        h = isqrt(u)
        h += h * h < u
        for factor in [16, 24, 32, 64]:
            dimension = factor * h
            for shift in [1, dimension // (4 * h)]:
                for copies in [1, 3, 7]:
                    K = dimension // (8 * h)
                    budget = copies * dimension * dimension // u
                    maximum_subset_weight = copies * sum(shift + smith_exponent(k) for k in range(K))
                    assert maximum_subset_weight <= budget
                    # Square to check K >= D/(32 sqrt(u)) without rounding.
                    assert (32 * K) ** 2 * u >= dimension * dimension
                    assert K <= dimension
                    lower_box_checks += 1

    entropy_rows = []
    for dimension, copies, u, shift in [(16, 1, 4, 4), (32, 1, 4, 5),
                                       (64, 1, 8, 6), (64, 2, 8, 6)]:
        budget = copies * dimension * dimension // u
        count = weighted_monomial_count(dimension, copies, shift, budget)
        entropy_rows.append({"dimension": dimension, "copies": copies, "u": u,
                             "base_valuation": shift, "budget": budget,
                             "weighted_monomials": count,
                             "log2_weighted_monomials": log2(count)})

    spectrum = [{"index": k, "valuation": smith_exponent(k)} for k in range(17)]
    report = {
        "status": "PASS", "exact_rational_identity_entries": identity_entries,
        "exact_minors": exact_minors, "minor_profiles": minor_profiles,
        "independent_chain_ring_profiles": modular_profiles,
        "modular_factorization_entries": modular_entries,
        "literal_tuple_count_comparisons": tuple_counts,
        "lower_entropy_box_checks": lower_box_checks,
        "initial_smith_exponents": spectrum, "entropy_examples": entropy_rows,
        "seconds": round(monotonic() - started, 3),
        "scope": "Exact finite checks support the accompanying all-size proof; no factoring exponent improvement is claimed.",
    }
    output = json.dumps(report, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
