#!/usr/bin/env python3
"""Independent exact controls for the strict universal mixed-cut gate.

The all-field proof is in the companion note. The finite rank check exhausts
every one-term opposite-cut correction on Mat_3(F_2). Threshold comparisons
use only integer arithmetic.
"""
import argparse
import json
from math import isqrt
from pathlib import Path


def rank(rows):
    pivots = {}
    for v in rows:
        while v:
            k = v.bit_length()-1
            if k not in pivots:
                pivots[k] = v
                break
            v ^= pivots[k]
    return len(pivots)


def kron_rows(a, b, q):
    return [sum(b[j] << (q*k) for k in range(q) if a[i] >> k & 1)
            for i in range(q) for j in range(q)]


def exhaustive_one_term_q3():
    q = 3
    matrices = [[(mask >> (q*i)) & 7 for i in range(q)] for mask in range(1,512)]
    ranks = [rank(a) for a in matrices]
    swap = [1 << (q*j+i) for i in range(q) for j in range(q)]
    count = invertible_pairs = 0
    minimum = minimum_singular = minimum_invertible = 10
    for a, ra in zip(matrices, ranks):
        for b, rb in zip(matrices, ranks):
            correction = kron_rows(a, b, q)
            residual_rank = rank([x ^ y for x, y in zip(swap, correction)])
            assert residual_rank >= 3
            count += 1
            minimum = min(minimum, residual_rank)
            if ra == rb == 3:
                invertible_pairs += 1
                minimum_invertible = min(minimum_invertible, residual_rank)
            else:
                minimum_singular = min(minimum_singular, residual_rank)
    assert count == 511**2 and invertible_pairs == 168**2
    return dict(q=3, prime=2, nonzero_factor_pairs=count,
                invertible_factor_pairs=invertible_pairs,
                minimum_residual_row_rank=minimum,
                minimum_with_invertible_factors=minimum_invertible,
                minimum_with_a_singular_factor=minimum_singular)


def q2_exception():
    identity = [1,2]
    swap = [1 << (2*j+i) for i in range(2) for j in range(2)]
    product = kron_rows(identity, identity, 2)
    r = rank([x ^ y for x, y in zip(swap, product)])
    assert r == 1
    return dict(q=2, prime=2, residual_row_rank=r, correction_column_rank=1,
                mixed_budget=2, strict_q_plus_one_claim_does_not_apply=True)


def threshold_certificates():
    # sqrt(2) < 17/12, since 288 < 289.
    assert 2*12**2 < 17**2
    cases = []
    for q in range(2,12):
        floor = 3 if q == 2 else 3*q-3
        assert floor**12 > q**17
        cases.append(dict(q=q, rational_floor_numerator=floor,
                          rational_floor_denominator=1,
                          left_integer=floor**12, right_integer=q**17))
    # At q12, 11+2sqrt(132)>11+2(57/5)=169/5.
    assert 132*25 > 57**2
    assert 169**12 > 5**12*12**17
    cases.append(dict(q=12, rational_floor_numerator=169,
                      rational_floor_denominator=5,
                      left_integer=169**12, right_integer=5**12*12**17))
    return dict(exponent_upper_bound_numerator=17, exponent_upper_bound_denominator=12,
                excluded_physical_sizes=list(range(2,13)), comparisons=cases)


def relaxation_witness(q):
    n, h = q*q, q*q-q
    pairs = [(h,1), (1,h)]+[(1,1)]*(q-1)
    assert sum(r for r,c in pairs) == sum(c for r,c in pairs) == n
    assert sum(min(r,c) for r,c in pairs) == q+1
    assert h == q*(q-1)
    return dict(q=q, rank_pairs=[list(pair) for pair in pairs],
                row_budget=n, column_budget=n, mixed_budget=q+1,
                cost_expression=f'{q-1}+2*sqrt({h})',
                realization_as_one_crossing_decomposition_not_claimed=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='universal all-field obstruction through physical size12; no factoring exponent improvement',
                  all_field_claims=dict(strict_mixed_budget='sum min(row_rank,column_rank) >= q+1 for q>=3',
                    cost_floor='q-1+2*sqrt(q*q-q) for q>=3'),
                  exhaustive_one_term_control=exhaustive_one_term_q3(),
                  q2_exception=q2_exception(), threshold_certificates=threshold_certificates(),
                  optimal_budget_relaxation_controls=[relaxation_witness(q) for q in (3,12,13)])
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
