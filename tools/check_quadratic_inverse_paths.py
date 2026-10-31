#!/usr/bin/env python3
"""Exact noncommuting quadratic-inverse paths and a complete first-jet screen."""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json

from check_quadratic_crossings import binary_rank, ROOT


def apply(columns, vector):
    result = 0
    while vector:
        bit = (vector & -vector).bit_length()-1
        result ^= columns[bit]
        vector &= vector-1
    return result


def compose(left, right):
    return [apply(left, column) for column in right]


def inverse(columns):
    n = len(columns)
    rows = [sum((columns[j] >> i & 1) << j for j in range(n)) | 1 << (n+i) for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i] >> j & 1)
        rows[j], rows[pivot] = rows[pivot], rows[j]
        for i in range(n):
            if i != j and rows[i] >> j & 1:
                rows[i] ^= rows[j]
    return [sum((rows[i] >> (n+j) & 1) << i for i in range(n)) for j in range(n)]


def kernel(columns):
    n = len(columns)
    rows = [sum((columns[j] >> i & 1) << j for j in range(n)) for i in range(n)]
    pivots, rank = [], 0
    for j in range(n):
        pivot = next((i for i in range(rank, n) if rows[i] >> j & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(n):
            if i != rank and rows[i] >> j & 1:
                rows[i] ^= rows[rank]
        pivots.append(j)
        rank += 1
    return [1 << free | sum((rows[i] >> free & 1) << pivot for i, pivot in enumerate(pivots))
            for free in range(n) if free not in pivots]


def kernel_flag(nilpotent):
    n = len(nilpotent)
    power, columns, levels, basis = [1 << i for i in range(n)], [], [], {}
    for level in range(1, n+1):
        power = compose(nilpotent, power)
        for vector in kernel(power):
            reduced = vector
            while reduced:
                bit = reduced.bit_length()-1
                if bit not in basis:
                    basis[bit] = reduced
                    columns.append(vector)
                    levels.append(level)
                    break
                reduced ^= basis[bit]
        if len(columns) == n:
            return columns, levels
    raise AssertionError("operator is not nilpotent")


def row_cut(columns, q):
    return [sum((columns[q*v+y] >> (q*x+u) & 1) << (q*u+v)
                for u in range(q) for v in range(q)) for x in range(q) for y in range(q)]


def quotient_jet(rows, constant):
    pivot_row = next(i for i, row in enumerate(constant) if row)
    pivot_column = (constant[pivot_row] & -constant[pivot_row]).bit_length()-1
    result = []
    for i, row in enumerate(rows):
        if i == pivot_row:
            continue
        if constant[i] >> pivot_column & 1:
            row ^= rows[pivot_row]
        if row >> pivot_column & 1:
            row ^= constant[pivot_row]
        result.append(row & ~(1 << pivot_column))
    return result


def construct(q, basis, levels, nilpotent_triangular, inverse_basis, one_levels):
    n = q*q
    # I-A(w) is strictly upper triangular off the diagonal, with edge weights 1 or 2.
    rows = [sum((nilpotent_triangular[j] >> i & 1) << j for j in range(n)) for i in range(n)]
    resolvent = [[int(i == j) for j in range(n)] for i in range(n)]
    for i in reversed(range(n)):
        for k in range(i+1, n):
            if rows[i] >> k & 1:
                weight = 1 if levels[k] in one_levels else 2
                for j in range(k, n):
                    resolvent[i][j] ^= resolvent[k][j] << weight
    unit = [column ^ 1 << j for j, column in enumerate(nilpotent_triangular)]
    # Multiply the resolvent by I-N, then change back to the physical matrix basis.
    current = [[0] * n for _ in range(n)]
    for j, column in enumerate(unit):
        for k in range(n):
            if column >> k & 1:
                for i in range(n):
                    current[i][j] ^= resolvent[i][k]
    left = [[0] * n for _ in range(n)]
    for k, column in enumerate(basis):
        for i in range(n):
            if column >> i & 1:
                for j in range(n):
                    left[i][j] ^= current[k][j]
    physical = [[0] * n for _ in range(n)]
    for j, column in enumerate(inverse_basis):
        for k in range(n):
            if column >> k & 1:
                for i in range(n):
                    physical[i][j] ^= left[i][k]
    return physical, resolvent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alphabet", type=int, default=7, choices=(3, 7, 9))
    args = parser.parse_args()
    q, n = args.alphabet, args.alphabet**2
    identity = [1 << i for i in range(n)]
    # U(X)=J_q(1) X^T is an invertible transpose channel, with row rank exactly one.
    unit = [(1 << (q*b+a)) ^ ((1 << (q*(b-1)+a)) if b else 0)
            for a in range(q) for b in range(q)]
    nilpotent = [a ^ b for a, b in zip(identity, unit)]
    basis, levels = kernel_flag(nilpotent)
    inverse_basis = inverse(basis)
    assert compose(basis, inverse_basis) == compose(inverse_basis, basis) == identity
    triangular = compose(compose(inverse_basis, nilpotent), basis)
    for j, column in enumerate(triangular):
        assert all(levels[i] < levels[j] for i in range(n) if column >> i & 1)
    constant = row_cut(unit, q)
    assert binary_rank(constant) == 1
    maximum_level = max(levels)
    possible = list(range(2, maximum_level+1))
    count_ones = {3: 3, 7: 6, 9: 1}[q]
    jets = {}
    for level in possible:
        selected = [column if levels[j] == level else 0 for j, column in enumerate(triangular)]
        physical = compose(compose(basis, selected), inverse_basis)
        jets[level] = quotient_jet(row_cut(compose(physical, unit), q), constant)
    histogram, best = Counter(), None
    for selected in combinations(possible, count_ones):
        rows = [0] * (n-1)
        for level in selected:
            rows = [a ^ b for a, b in zip(rows, jets[level])]
        multiplicity = binary_rank(rows)
        histogram[multiplicity] += 1
        if best is None or multiplicity < best[0]:
            best = multiplicity, selected
    multiplicity, selected = best
    physical, resolvent = construct(q, basis, levels, triangular, inverse_basis, selected)
    degree = max(value.bit_length()-1 for row in physical for value in row)
    assert degree == 2*(maximum_level-1)-count_ones
    assert all(physical[i][j].bit_count() % 2 == int(i == j) for i in range(n) for j in range(n))
    assert all(physical[i][j] & 1 == (unit[j] >> i & 1) for i in range(n) for j in range(n))
    # An independent Toeplitz rank recovers the best witness's order-one count.
    first_columns = [sum((physical[i][j] >> 1 & 1) << i for i in range(n)) for j in range(n)]
    first_rows = row_cut(first_columns, q)
    toeplitz = constant + [first_rows[i] | constant[i] << n for i in range(n)]
    assert binary_rank(toeplitz)-2 == multiplicity
    # Verify the entire original polynomial identity Q(w)(I-wS-w^2T)=I.
    for i in range(n):
        for j in range(n):
            value = resolvent[i][j]
            weight = 1 if levels[j] in selected else 2
            for k in range(n):
                if triangular[j] >> k & 1:
                    value ^= resolvent[i][k] << weight
            assert value == int(i == j)
    small_scale = Fraction(n, 1) / Fraction(3, 2)**2 * (1+Fraction(multiplicity, 1)/Fraction(3, 2))
    large_scale = Fraction(3, 2)**degree
    excluded = small_scale**7 > q**20 and large_scale**7 > q**20
    if q == 7:
        assert excluded, "a surviving first-jet candidate requires full profile evaluation"
    compatibility_controls = 0
    for name in ("quadratic-crossing-catalogue.json", "parity-cubic-crossing-catalogue.json",
                 "parity-quartic-crossing-catalogue.json", "cubic-crossing-million-sample.json",
                 "rank-one-constant-quadratic-catalogues.json"):
        data = json.loads((ROOT / "evidence" / name).read_text())
        for catalogue in data if isinstance(data, list) else [data]:
            regularization = catalogue["degree"]
            for row in catalogue["profiles"]:
                mu, nu = row["infinity"], row["zero"]
                if sum(nu) == 4*regularization and mu[0] == 0 < mu[1]:
                    assert mu[1] <= 2*(nu[-1]-regularization)
                    compatibility_controls += 1
    report = {"alphabet": q, "dimension": n, "kernel_level_counts": dict(Counter(levels)),
                      "paths_screened": sum(histogram.values()), "order_one_histogram": dict(sorted(histogram.items())),
                      "minimum_order_one_multiplicity": multiplicity, "best_one_weight_levels": selected,
                      "actual_degree": degree, "row_constant_rank": 1, "inverse_degree_at_most": 2,
                      "original_polynomial_resolvent_identity": "passed",
                      "physical_endpoint_and_constant": "passed",
                      "independent_best_first_jet_toeplitz": "passed",
                      "physical_polynomial_sha256": sha256(json.dumps(physical).encode()).hexdigest(),
                      "prior_catalogue_inverse_compatibility_controls": compatibility_controls,
                      "complete_first_jet_exclusion": excluded,
                      "small_scale_B_lower_bound": str(small_scale),
                      "large_scale_B_lower_bound": str(large_scale),
                      "scope": "This weighted kernel-flag family only; no factoring improvement"}
    if q == 9:
        from check_leading_rank_entropy_floor import entropy_dual
        alpha, _, _, denominator, entropy_power, _ = entropy_dual(9)
        penalty = alpha[0]-2*alpha[1]/7
        cutoff = Fraction(239, 200)
        assert (entropy_power/cutoff**int(denominator*penalty))**6 > 9**(17*denominator)
        assert cutoff**(6*35) > 9**17
        full_jet = [0] * (n-1)
        for jet in jets.values():
            full_jet = [a ^ b for a, b in zip(full_jet, jet)]
        full_rank = binary_rank(full_jet)
        top_three = sorted((rank for rank, count in histogram.items() for _ in range(count)), reverse=True)[:3]
        early_floor = full_rank-sum(top_three)
        assert full_rank == 80 and top_three == [22, 22, 18] and early_floor == 18
        split = Fraction(5, 4)
        small_case = Fraction(n)/split**2 * (1+Fraction(early_floor)/split)
        large_case = split**31
        assert small_case**7 > 9**20 and large_case**7 > 9**20
        report['all_binary_level_weights'] = {
            'assignments_covered': 2**31, 'fixed_kernel_flag': True,
            'degree': '31 + number of weight-two levels',
            'high_degree_case': 'at least four weight-two levels: degree >=35',
            'high_degree_rho_floor': '17/24', 'all_one_jet_rank': full_rank,
            'top_three_removed_jet_ranks': top_three,
            'few_twos_order_one_lower_bound': early_floor,
            'few_twos_scale_split': str(split),
            'few_twos_small_scale_pressure_floor': str(small_case),
            'few_twos_large_scale_pressure_floor': str(large_case),
            'few_twos_rho_floor': '5/7'}
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
