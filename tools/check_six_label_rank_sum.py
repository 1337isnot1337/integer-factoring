#!/usr/bin/env python3
"""Integer certificates for the universal six-label rank-sum obstruction."""

import argparse
import itertools
import json
from math import isqrt
from pathlib import Path
import random


def rank_mod(matrix, prime):
    rows = [row[:] for row in matrix]
    rank = 0
    if not rows:
        return 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows))
                      if rows[i][column] % prime), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][column], -1, prime)
        rows[rank] = [(x*inv) % prime for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            if rows[i][column]:
                factor = rows[i][column]
                rows[i] = [(x-factor*y) % prime
                           for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def swap_checks():
    rng = random.Random(630261)
    checks = 0
    for prime, dimension in itertools.product((2, 3, 5), range(2, 7)):
        for _ in range(12):
            matrices = []
            for _ in range(2):
                while True:
                    candidate = [[rng.randrange(prime) for _ in range(dimension)]
                                 for _ in range(dimension)]
                    if rank_mod(candidate, prime) == dimension:
                        matrices.append(candidate)
                        break
            p, q = matrices
            difference = [
                [(int(x == v and y == u)-p[x][u]*q[y][v]) % prime
                 for u in range(dimension) for v in range(dimension)]
                for x in range(dimension) for y in range(dimension)
            ]
            assert rank_mod(difference, prime) >= dimension-1
            checks += 1
    return checks


def verify():
    assert 99**2 > 2*70**2
    assert 1261**70 > 6**99*100**70
    assert 34*200**2 > 1161**2
    tested = 0
    surviving = []
    for r1 in range(1, 37):
        for r2 in range(max(1, 36-r1), 37):
            for c1 in range(1, 37):
                for c2 in range(max(1, 36-c1), 37):
                    tested += 1
                    lower = isqrt(10000*r1*c1)+isqrt(10000*r2*c2)
                    if lower < 1300:
                        surviving.append((r1, c1, r2, c2))
    expected = sorted(
        [(1, a, b, 1) for a in (35, 36) for b in (35, 36)]
        + [(a, 1, 1, b) for a in (35, 36) for b in (35, 36)]
    )
    assert sorted(surviving) == expected
    products = [(a, b, a*b) for a in range(1, 7) for b in range(1, 7)
                if a*b >= 35]
    assert products == [(6, 6, 36)]
    return {
        "scope": "universal six-label exclusion for the rank-sum criterion",
        "all_field_matrix_argument": "rank(Swap-P tensor Q)>=d-1 for invertible P,Q",
        "exact_threshold_certificates": {
            "sqrt2_less_than_99_over_70": True,
            "six_to_sqrt2_less_than_1261_over_100": True,
            "one_plus_two_sqrt34_greater_than_1261_over_100": True,
        },
        "integer_rank_quadruples_checked": tested,
        "rank_patterns_surviving_relaxed_threshold_13": surviving,
        "rank_products_at_least_35": products,
        "finite_field_swap_controls": swap_checks(),
        "minimum_possible_seed_dimension_for_this_criterion": 7,
        "does_not_exclude_other_amplification_methods": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass",
                      "rank_quadruples": result["integer_rank_quadruples_checked"],
                      "surviving_patterns": len(result["rank_patterns_surviving_relaxed_threshold_13"]),
                      "swap_controls": result["finite_field_swap_controls"]}))


if __name__ == "__main__":
    main()
