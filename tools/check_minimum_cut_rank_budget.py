#!/usr/bin/env python3
"""Exact minimum-cut rank budgets and the twenty-one-label rank-sum cutoff."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import random

from check_swap_product_rank_sum import SCALE, rank_budget, rank_mod
from check_distinguished_rank_budget import surds, compare


@lru_cache(None)
def root_numerator(value):
    result = isqrt(value * SCALE * SCALE)
    assert result * result <= value * SCALE * SCALE < (result + 1) ** 2
    return result


def determinant(first, second, third):
    return (first[0] * (second[1] * third[2] - second[2] * third[1])
            - first[1] * (second[0] * third[2] - second[2] * third[0])
            + first[2] * (second[0] * third[1] - second[1] * third[0]))


def feasible_vertices(d, a, b):
    """Return exact rational coordinates with one common integer denominator."""
    n = d * d
    constraints = [(n, 2, 2, root_numerator(2 * n)),
                   (2, n, 2, root_numerator(2 * n)),
                   (1, 1, 1, SCALE)]
    if a:
        constraints.append((a, 1, 1, root_numerator(a)))
    if b:
        constraints.append((1, b, 1, root_numerator(b)))
    planes = constraints + [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]
    for first, second, third in combinations(planes, 3):
        denominator = determinant(first, second, third)
        if not denominator:
            continue
        numerators = [
            determinant((first[3], first[1], first[2]),
                        (second[3], second[1], second[2]),
                        (third[3], third[1], third[2])),
            determinant((first[0], first[3], first[2]),
                        (second[0], second[3], second[2]),
                        (third[0], third[3], third[2])),
            determinant((first[0], first[1], first[3]),
                        (second[0], second[1], second[3]),
                        (third[0], third[1], third[3])),
        ]
        if denominator < 0:
            denominator = -denominator
            numerators = [-value for value in numerators]
        if min(numerators) < 0:
            continue
        if any(sum(line[i] * numerators[i] for i in range(3)) > line[3] * denominator
               for line in constraints):
            continue
        yield tuple(numerators), denominator


def check_unrounded_endpoints(d, a, b, numerators, determinant_denominator):
    n = d * d
    denominator = determinant_denominator * SCALE
    u, v, w = numerators
    assert min(numerators) >= 0
    endpoints = [(2, 2), (2, n), (n, 2), (n, n), (1, 1)]
    if a:
        endpoints.append((a, 1))
    if b:
        endpoints.append((1, b))
    for row, column in endpoints:
        value = row * u + column * v + min(row, column) * w
        assert value * value <= row * column * denominator * denominator
    return len(endpoints)


def profile_bound(d, a, b, use_nonlinear=True):
    row_budget, column_budget = rank_budget(d, a), rank_budget(d, b)
    best_numerator, best_denominator = 0, 1
    selected = None
    vertices = 0
    for numerators, denominator in feasible_vertices(d, a, b):
        u, v, w = numerators
        values = {"dual": row_budget * u + column_budget * v + d * w}
        if a:
            values["one_row_maximum"] = (root_numerator(a) * denominator
                + (row_budget - a) * u + (column_budget - 1) * v + (d - 1) * w)
        if b:
            values["one_column_maximum"] = (root_numerator(b) * denominator
                + (row_budget - 1) * u + (column_budget - b) * v + (d - 1) * w)
        if a and b and (a, b) != (1, 1):
            values["two_distinct_maxima"] = ((root_numerator(a) + root_numerator(b)) * denominator
                + (row_budget - a - 1) * u + (column_budget - b - 1) * v + (d - 2) * w)
        for method, value in values.items():
            if value * best_denominator > best_numerator * denominator:
                best_numerator, best_denominator = value, denominator
                selected = method, numerators, denominator
        vertices += 1
    assert selected is not None
    method, numerators, denominator = selected
    endpoints = check_unrounded_endpoints(d, a, b, numerators, denominator)
    point = [str(Fraction(value, SCALE * denominator)) for value in numerators]
    bound = Fraction(best_numerator, SCALE * best_denominator)
    if use_nonlinear and a == b and a in (0, 1):
        nonlinear = d - 4 + 2 * Fraction(root_numerator(2 * (d * d - d + 2)), SCALE)
        if nonlinear > bound:
            bound, method = nonlinear, "nonlinear_surplus"
    return bound, method, point, vertices, endpoints


def rank_sum_certificates():
    assert 577 ** 2 - 2 * 408 ** 2 == 1
    records = []
    profiles = vertices = endpoints = 0
    for d in range(2, 22):
        products = [0] + sorted({a * b for a in range(1, d + 1) for b in range(1, d + 1)})
        smallest = None
        worst = None
        methods = Counter()
        count = 0
        for index, a in enumerate(products):
            for b in products[index:]:
                if (a, b) == (0, 1):
                    continue
                bound, method, point, found_vertices, tested_endpoints = profile_bound(d, a, b)
                assert bound ** 408 > d ** 577
                if smallest is None or bound < smallest:
                    smallest = bound
                    worst = {"a": a, "b": b, "method": method,
                             "dual_u": point[0], "dual_v": point[1], "dual_w": point[2]}
                methods[method] += 1
                vertices += found_vertices
                endpoints += tested_endpoints
                count += 1
        profiles += count
        records.append({"dimension": d, "maximum_pairs": count,
                        "rational_lower_bound": str(smallest),
                        "decimal_summary": float(smallest), "worst_pair": worst,
                        "winning_bound_counts": dict(sorted(methods.items())),
                        "exact_power_408_exceeds_d_power_577": True})
    assert records[-1]["rational_lower_bound"] == "65335133917/878000000"
    # Record why the nonlinear addition matters, without calling the
    # earlier linear failure a tensor example.
    linear, _, _, _, _ = profile_bound(21, 0, 0, use_nonlinear=False)
    assert linear == Fraction(507589488, 6859375)
    assert not linear ** 408 > 21 ** 577
    failed, method, point, _, _ = profile_bound(22, 0, 0)
    assert failed == Fraction(98657731, 1250000)
    assert not failed ** 408 > 22 ** 577
    assert 239 ** 2 < 2 * 169 ** 2 and failed ** 169 < 22 ** 239
    numerical = ((18, 18), (2, 464), (464, 2))
    assert sum(pair[0] for pair in numerical) == sum(pair[1] for pair in numerical) == 484
    assert sum(min(pair) for pair in numerical) == 22
    root = root_numerator(928)
    ceiling = root if root * root == 928 * SCALE * SCALE else root + 1
    upper_cost = 18 + 2 * Fraction(ceiling, SCALE)
    assert upper_cost ** 169 < 22 ** 239
    required_complement = max(s * (22 - 2 * s) for s in range(23))
    assert required_complement == 60 > numerical[0][1] + numerical[2][1] == 20
    return {"maximum_pairs": profiles, "feasible_vertices": vertices,
            "unrounded_endpoint_checks": endpoints, "certificates": records,
            "nonlinear_improves_linear_frontier": {"dimension": 21,
                "old_linear_bound": str(linear), "exact_power_test_failed": True},
            "first_uncertified_profile": {"dimension": 22, "a": 0, "b": 0,
                "rational_lower_bound": str(failed), "method": method,
                "dual": point, "numerical_rank_pairs": numerical,
                "required_rank_two_complement": required_complement,
                "available_complement_rank_budget": 20,
                "is_tensor_decomposition": False}}


def nonlinear_controls():
    """Test the general scalar lemma, including its exact numerical extrema."""
    rng = random.Random(20261104)
    random_cases = attainment_cases = 0
    for d in range(4, 13):
        n = d * d
        radicand = 2 * (n - d + 2)
        lower_terms = (() if d == 4 else ((d - 4) ** 2,)) + (radicand, radicand)
        lower = surds(lower_terms)
        for _ in range(150):
            pairs = []
            for _ in range(rng.randrange(7)):
                if rng.randrange(5) == 0:
                    pairs.append((1, 1))
                    continue
                minimum = rng.randrange(2, n + 1)
                maximum = rng.randrange(minimum, n + 1)
                pair = (minimum, maximum)
                pairs.append(pair if rng.randrange(2) else tuple(reversed(pair)))
            deficit = max(0, d - sum(min(pair) for pair in pairs),
                          n - sum(pair[0] for pair in pairs),
                          n - sum(pair[1] for pair in pairs))
            if deficit:
                pairs.append((deficit, deficit))
            assert sum(pair[0] for pair in pairs) >= n and sum(pair[1] for pair in pairs) >= n
            assert sum(min(pair) for pair in pairs) >= d
            assert all(min(pair) >= 2 or pair == (1, 1) for pair in pairs)
            assert compare(surds(tuple(a * b for a, b in pairs)), lower) >= 0
            random_cases += 1
        # These are numerical rank-budget extrema, not tensor decompositions.
        unit_extremum = [(1, 1)] * (d - 4) + [(2, n - d + 2), (n - d + 2, 2)]
        assert sum(a for a, _ in unit_extremum) == sum(b for _, b in unit_extremum) == n
        assert sum(min(pair) for pair in unit_extremum) == d
        assert surds(tuple(a * b for a, b in unit_extremum)) == lower
        attainment_cases += 1
        if d >= 6:
            diagonal_extremum = [(d - 4, d - 4), (2, n - d + 2), (n - d + 2, 2)]
            assert surds(tuple(a * b for a, b in diagonal_extremum)) == lower
            attainment_cases += 1
    return {"scalar_budget_checks": random_cases,
            "exact_numerical_extrema": attainment_cases,
            "scope": "Rank-budget relaxation only; extrema need not be tensor-realizable."}


def column_matrix(row_matrix, d):
    result = [[0] * (d * d) for _ in range(d * d)]
    for x in range(d):
        for y in range(d):
            for u in range(d):
                for v in range(d):
                    result[x * d + u][y * d + v] = row_matrix[x * d + y][u * d + v]
    return result


def product_matrix(left, right, prime):
    d = len(left)
    return [[left[x][u] * right[y][v] % prime for u in range(d) for v in range(d)]
            for x in range(d) for y in range(d)]


def matrix_controls():
    rng = random.Random(20261103)
    quotients = decompositions = zero_product_cases = 0
    minima = []
    for prime in (2, 3, 5):
        for d in range(1, 7):
            n = d * d
            swap = [[int(x == v and y == u) for u in range(d) for v in range(d)]
                    for x in range(d) for y in range(d)]
            for summands in range(4):
                for _ in range(3):
                    pieces = []
                    total = [[0] * n for _ in range(n)]
                    for _ in range(summands):
                        left = [[rng.randrange(prime) for _ in range(d)] for _ in range(d)]
                        right = [[rng.randrange(prime) for _ in range(d)] for _ in range(d)]
                        piece = product_matrix(left, right, prime)
                        pieces.append(piece)
                        total = [[(x + y) % prime for x, y in zip(row, other)]
                                 for row, other in zip(total, piece)]
                    complement = [[(x - y) % prime for x, y in zip(row, other)]
                                  for row, other in zip(swap, total)]
                    actual_k = rank_mod(column_matrix(total, d), prime)
                    assert actual_k <= summands
                    lower = max(s * (d - actual_k * s) for s in range(d + 1))
                    assert rank_mod(complement, prime) >= lower
                    quotients += 1
                    if not actual_k:
                        zero_product_cases += 1
                    pieces.append(complement)
                    cost = sum(min(rank_mod(piece, prime),
                                   rank_mod(column_matrix(piece, d), prime)) for piece in pieces)
                    assert cost >= d
                    assert all(sum(piece[i][j] for piece in pieces) % prime == swap[i][j]
                               for i in range(n) for j in range(n))
                    decompositions += 1
            # A small sharp control for the universal constant in (1).
            if d == 2:
                identity = [[int(i == j) for j in range(n)] for i in range(n)]
                residual = [[(swap[i][j] - identity[i][j]) % prime for j in range(n)]
                            for i in range(n)]
                cost = sum(min(rank_mod(piece, prime), rank_mod(column_matrix(piece, d), prime))
                           for piece in (identity, residual))
                assert cost == d
                minima.append({"prime": prime, "dimension": d, "sum_minimum_ranks": cost})
    assert max(s * (13 - 2 * s) for s in range(14)) == 21
    return {"higher_rank_quotient_checks": quotients,
            "exact_crossing_decompositions": decompositions,
            "zero_product_cases": zero_product_cases,
            "two_label_sharp_controls": minima,
            "thirteen_label_rank_two_complement_bound": 21}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = rank_sum_certificates()
    result = {"status": "PASS", "scope": "Finite rank-sum amplification only; no factoring lower bound.",
              "minimum_cut_budget": "sum_i min(r_i,c_i)>=d",
              "minimum_seed_dimension_for_a_strict_rank_sum_win": 22,
              "rational_root_scale": SCALE, "rank_sum": certificate,
              "nonlinear_surplus_controls": nonlinear_controls(),
              "finite_field_controls": matrix_controls()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"],
                      "minimum_seed_dimension": 22,
                      "maximum_pairs": certificate["maximum_pairs"],
                      "feasible_vertices": certificate["feasible_vertices"],
                      "unrounded_endpoints": certificate["unrounded_endpoint_checks"],
                      "matrix_controls": result["finite_field_controls"]}))


if __name__ == "__main__":
    main()
