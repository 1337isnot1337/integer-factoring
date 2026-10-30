#!/usr/bin/env python3
"""Exact certificates for the swap/product lemma and rank-sum cutoff."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import random


SCALE = 10**7


def root_lower(value):
    numerator = isqrt(value*SCALE*SCALE)
    assert numerator*numerator <= value*SCALE*SCALE < (numerator+1)**2
    return Fraction(numerator, SCALE)


def rank_budget(d, maximum):
    n = d*d
    if not maximum:
        return n
    fitting_cap = max(min(left, maximum//left) for left in range(1, d+1)
                      if maximum % left == 0 and maximum//left <= d)
    complement = n-fitting_cap*(fitting_cap+1)//2
    return max(n, maximum+complement)


def dual_bound(d, a, b):
    n = d*d
    constraints = [(n, 2, root_lower(2*n)),
                   (2, n, root_lower(2*n))]
    if a or b:
        constraints.append((1, 1, Fraction(1)))
    if a:
        constraints.append((a, 1, root_lower(a)))
    if b:
        constraints.append((1, b, root_lower(b)))
    row_budget, col_budget = rank_budget(d, a), rank_budget(d, b)
    best, point = Fraction(0), (Fraction(0), Fraction(0))
    lines = constraints+[(1, 0, Fraction(0)), (0, 1, Fraction(0))]
    for (x1, y1, z1), (x2, y2, z2) in combinations(lines, 2):
        determinant = x1*y2-x2*y1
        if not determinant:
            continue
        u = (z1*y2-z2*y1)/determinant
        v = (x1*z2-x2*z1)/determinant
        if u < 0 or v < 0:
            continue
        if any(x*u+y*v > z for x, y, z in constraints):
            continue
        objective = row_budget*u+col_budget*v
        if objective > best:
            best, point = objective, (u, v)
    u, v = point
    endpoints = [(2, 2), (2, n), (n, 2), (n, n)]
    if a or b:
        endpoints.append((1, 1))
    if a:
        endpoints.append((a, 1))
    if b:
        endpoints.append((1, b))
    # An independent exact check against the unrounded endpoint costs.
    for row, col in endpoints:
        assert (row*u+col*v)**2 <= row*col
    assert best == row_budget*u+col_budget*v
    return best, point, len(endpoints)


def case_bound(d, a, b):
    n = d*d
    row_budget, col_budget = rank_budget(d, a), rank_budget(d, b)
    dual, point, endpoint_count = dual_bound(d, a, b)
    bounds = {"dual": dual}
    if a and row_budget-a >= 2 and col_budget-1 >= 2:
        bounds["one_row_maximum"] = (root_lower(a)
            + root_lower(row_budget-a-1)+root_lower(col_budget-2))
    if b and row_budget-1 >= 2 and col_budget-b >= 2:
        bounds["one_column_maximum"] = (root_lower(b)
            + root_lower(row_budget-2)+root_lower(col_budget-b-1))
    if a and b and (a, b) != (1, 1):
        if row_budget-a-1 >= 2 and col_budget-b-1 >= 2:
            bounds["two_distinct_maxima"] = (root_lower(a)+root_lower(b)
                + root_lower(row_budget-a-2)+root_lower(col_budget-b-2))
    if a == b == 0:
        bounds["both_cut_minima_two"] = 2*root_lower(2*(n-2))
    if a == b == 1:
        bounds["shared_rank_one_pieces"] = 1+2*root_lower(2*(n-3))
    if a and not b and row_budget-a >= 4:
        bounds["no_row_rank_one_remainder"] = (root_lower(a)
            + root_lower(row_budget-a-2)+root_lower(2*n-4))
    if b and not a and col_budget-b >= 4:
        bounds["no_column_rank_one_remainder"] = (root_lower(b)
            + root_lower(col_budget-b-2)+root_lower(2*n-4))
    method = max(bounds, key=bounds.get)
    return bounds[method], method, point, endpoint_count


def rank_sum_certificates():
    assert 577**2-2*408**2 == 1
    records = []
    total_cases = endpoint_checks = 0
    for d in range(2, 13):
        possible = [0]+sorted({a*b for a in range(1, d+1)
                              for b in range(1, d+1)})
        minimum = None
        worst = None
        methods = Counter()
        cases = 0
        for i, a in enumerate(possible):
            for b in possible[i:]:
                if (a, b) == (0, 1):
                    # A row-rank-one piece of column rank one also makes a>=1.
                    continue
                bound, method, point, endpoint_count = case_bound(d, a, b)
                assert bound**408 > d**577
                if minimum is None or bound < minimum:
                    minimum = bound
                    worst = {"a": a, "b": b, "method": method,
                             "dual_u": str(point[0]), "dual_v": str(point[1])}
                methods[method] += 1
                endpoint_checks += endpoint_count
                cases += 1
        total_cases += cases
        assert minimum**408 > d**577
        records.append({
            "dimension": d,
            "possible_maxima": len(possible),
            "unordered_maximum_pairs": cases,
            "rational_lower_bound": str(minimum),
            "decimal_summary": float(minimum),
            "worst_pair": worst,
            "winning_bound_counts": dict(sorted(methods.items())),
            "exact_lower_bound_power_408_exceeds_d_power_577": True,
        })
    return records, total_cases, endpoint_checks


def rank_mod(matrix, prime):
    rows = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(rank, len(rows))
                      if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(entry*inverse) % prime for entry in rows[rank]]
        for i in range(rank+1, len(rows)):
            factor = rows[i][column]
            if factor:
                rows[i] = [(x-factor*y) % prime
                           for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def multiply(left, right, prime):
    return [[sum(x*y for x, y in zip(row, column)) % prime
             for column in zip(*right)] for row in left]


def matrix_power(matrix, exponent, prime):
    size = len(matrix)
    result = [[int(i == j) for j in range(size)] for i in range(size)]
    while exponent:
        if exponent & 1:
            result = multiply(result, matrix, prime)
        matrix = multiply(matrix, matrix, prime)
        exponent >>= 1
    return result


def swap_difference(p_matrix, q_matrix, prime):
    d = len(p_matrix)
    return [[(int(x == v and y == u)-p_matrix[x][u]*q_matrix[y][v]) % prime
             for u in range(d) for v in range(d)]
            for x in range(d) for y in range(d)]


def verify_matrix_pair(p_matrix, q_matrix, prime):
    d = len(p_matrix)
    pq = multiply(p_matrix, q_matrix, prime)
    qp = multiply(q_matrix, p_matrix, prime)
    stable = rank_mod(matrix_power(qp, d, prime), prime)
    assert stable == rank_mod(matrix_power(pq, d, prime), prime)
    rank_p = rank_mod(p_matrix, prime)
    rank_q = rank_mod(q_matrix, prime)
    assert stable <= min(rank_p, rank_q)
    actual = rank_mod(swap_difference(p_matrix, q_matrix, prime), prime)
    assert actual >= d*d-stable*(stable+1)//2
    assert actual >= d*(d-1)//2
    if min(rank_p, rank_q) < d:
        assert actual >= d*(d+1)//2
    return actual, stable


def matrix_controls():
    exhaustive = random_cases = sharp_cases = nilpotent_cases = 0
    for prime in (2, 3):
        matrices = [[[a, b], [c, d]]
                    for a, b, c, d in product(range(prime), repeat=4)]
        for left, right in product(matrices, repeat=2):
            verify_matrix_pair(left, right, prime)
            exhaustive += 1
    rng = random.Random(710301)
    for prime, d in product((2, 3, 5, 7), range(1, 8)):
        identity = [[int(i == j) for j in range(d)] for i in range(d)]
        assert verify_matrix_pair(identity, identity, prime)[0] == d*(d-1)//2
        projection = [[int(i == j and i < d-1) for j in range(d)]
                      for i in range(d)]
        assert verify_matrix_pair(projection, projection, prime)[0] == d*(d+1)//2
        sharp_cases += 2
        nilpotent = [[int(j == i+1) for j in range(d)] for i in range(d)]
        assert verify_matrix_pair(identity, nilpotent, prime) == (d*d, 0)
        nilpotent_cases += 1
        for iteration in range(12):
            left = [[rng.randrange(prime) for _ in range(d)] for _ in range(d)]
            right = [[rng.randrange(prime) for _ in range(d)] for _ in range(d)]
            if iteration % 3 == 0:
                left[-1] = [0]*d
            if iteration % 3 == 1:
                for row in right:
                    row[-1] = 0
            verify_matrix_pair(left, right, prime)
            random_cases += 1
    # The singular strengthening must not be asserted for arbitrary invertible factors.
    assert rank_mod(swap_difference([[1, 0], [0, 1]],
                                    [[1, 0], [0, 1]], 2), 2) == 1 < 3
    return {"exhaustive_two_by_two_pairs": exhaustive,
            "seeded_larger_pairs": random_cases,
            "sharpness_controls": sharp_cases,
            "nilpotent_fitting_controls": nilpotent_cases,
            "negative_invertible_singular_bound_control": True}


def verify():
    certificates, cases, endpoints = rank_sum_certificates()
    return {
        "scope": "all-field finite rank-sum criterion; not a factoring lower bound",
        "all_field_lemma": "rank(Swap-P tensor Q)>=d^2-r(r+1)/2, r=rank((QP)^d)",
        "rational_root_scale": SCALE,
        "rank_maximum_pairs_checked": cases,
        "unrounded_dual_endpoints_checked": endpoints,
        "exact_certificates": certificates,
        "finite_field_controls": matrix_controls(),
        "minimum_seed_dimension_for_a_strict_rank_sum_win": 13,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass",
                      "rank_maximum_pairs": result["rank_maximum_pairs_checked"],
                      "dual_endpoints": result["unrounded_dual_endpoints_checked"],
                      "minimum_seed_dimension": 13,
                      "matrix_controls": result["finite_field_controls"]}))


if __name__ == "__main__":
    main()
