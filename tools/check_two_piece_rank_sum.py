#!/usr/bin/env python3
"""Exact controls for the all-dimensional two-piece rank-sum obstruction."""

import argparse
from itertools import product
import json
from pathlib import Path
import random

from check_swap_product_rank_sum import rank_mod


def root_sum_ge(a, b, c, d):
    """Compare sqrt(a)+sqrt(b) >= sqrt(c)+sqrt(d), using integers only."""
    assert min(a, b, c, d) >= 0
    delta = a+b-c-d
    left, right = a*b, c*d
    if delta >= 0:
        if left >= right:
            return True
        residual = 4*right-delta*delta-4*left
        return residual <= 0 or residual*residual <= 16*delta*delta*left
    if left <= right:
        return False
    delta = -delta
    residual = 4*left-delta*delta-4*right
    return residual >= 0 and residual*residual >= 16*delta*delta*right


def scalar_controls():
    assert root_sum_ge(1, 4, 4, 1)
    assert root_sum_ge(4, 9, 1, 4)
    assert not root_sum_ge(1, 4, 4, 9)
    assert root_sum_ge(0, 9, 1, 4)
    raw = normalized = sharp = 0
    for d in range(2, 5):
        n = d*d
        target = (n-d+1, (n-1)*(d-1))
        for r1, r2, c1, c2 in product(range(1, n+1), repeat=4):
            if r1+r2 < n or c1+c2 < n or r1+c2 < d or r2+c1 < d:
                continue
            assert root_sum_ge(r1*c1, r2*c2, *target)
            raw += 1
    for d in range(2, 13):
        n = d*d
        target = (n-d+1, (n-1)*(d-1))
        for x, y in product(range(1, n), repeat=2):
            if not d <= x+y <= 2*n-d:
                continue
            assert root_sum_ge(x*(n-y), y*(n-x), *target)
            normalized += 1
        assert root_sum_ge(n*n, 0, *target)
        # Exact algebra for the strict g(d)^2>d^3 comparison.
        assert (n-d+1)*(n-1) > d-1
        assert root_sum_ge(target[0], target[1], d**3, 0)
        assert not root_sum_ge(d**3, 0, *target)
        assert root_sum_ge(target[0], target[1], *target)
        sharp += 1
    # Omitting the cross constraints gives a spurious lower-cost rank tuple.
    d = 5
    n = d*d
    assert not root_sum_ge(n-1, n-1, n-d+1, (n-1)*(d-1))
    return {"raw_admissible_rank_quadruples": raw,
            "normalized_rank_pairs": normalized,
            "strict_algebraic_dimension_controls": sharp,
            "negative_missing_cross_constraint_control": True}


def other_cut(matrix, d):
    return [[matrix[x*d+y][u*d+v] for y in range(d) for v in range(d)]
            for x in range(d) for u in range(d)]


def matrix_controls():
    rng = random.Random(230130)
    cases = zero_piece_cases = 0
    for prime in (2, 3, 5):
        for d in range(2, 7):
            n = d*d
            swap = [[int(x == v and y == u) for u in range(d) for v in range(d)]
                    for x in range(d) for y in range(d)]
            target = (n-d+1, (n-1)*(d-1))
            for iteration in range(14):
                first = [[rng.randrange(prime) for _ in range(n)] for _ in range(n)]
                if iteration == 0:
                    first = [[0]*n for _ in range(n)]
                if iteration == 1:
                    first = [row[:] for row in swap]
                second = [[(x-y) % prime for x, y in zip(row, other)]
                          for row, other in zip(swap, first)]
                r1, r2 = rank_mod(first, prime), rank_mod(second, prime)
                c1 = rank_mod(other_cut(first, d), prime)
                c2 = rank_mod(other_cut(second, d), prime)
                assert r1+r2 >= n and c1+c2 >= n
                assert r1+c2 >= d and r2+c1 >= d
                assert root_sum_ge(r1*c1, r2*c2, *target)
                zero_piece_cases += int(not min(r1, r2))
                cases += 1
    return {"exact_two_piece_decompositions": cases,
            "zero_piece_controls": zero_piece_cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"scope": "all-field two-piece exclusion for finite rank-sum amplification",
              "universal_bound": "sqrt(d^2-d+1)+sqrt((d^2-1)(d-1)) > d^(3/2)",
              "scalar_controls": scalar_controls(),
              "matrix_controls": matrix_controls(),
              "minimum_nonzero_pieces_for_a_strict_rank_sum_win": 3}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass", **result["scalar_controls"],
                      **result["matrix_controls"]}))


if __name__ == "__main__":
    main()
