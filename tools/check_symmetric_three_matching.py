#!/usr/bin/env python3
"""Exact checks for the symmetric three-matching rank obstruction."""
from functools import lru_cache
from itertools import product
from math import comb
import json
from check_binary_coefficient_entropy import binary_rank, coefficient_rank, elementary_rank


def tail(v, degree):
    return sum(comb(v, j) for j in range((v-degree)//2+1)) if degree <= v else 0


@lru_cache(None)
def split_rank(t, degree, split):
    return sum(comb(t, a)*elementary_rank(a, split)*elementary_rank(t-a, degree-split)
               for a in range(t+1))


def walk_count(t, degree, split):
    return sum(comb(t, a)*tail(a, split)*tail(t-a, degree-split)
               for a in range(t+1))


def walk_count_rotated(t, degree, split):
    answer = 0
    for plus in range(t+1):
        first = 2*plus-t
        if first < degree:
            continue
        for other_plus in range(t+1):
            second = 2*other_plus-t
            if 2*split-first <= second <= first-2*degree+2*split:
                answer += comb(t, plus)*comb(t, other_plus)
    return answer


def literal_types(t):
    types = [(r, s, t-r-s) for r in range(t+1) for s in range(t-r+1)]
    indices = {(r, s): i for i, (r, s, _) in enumerate(types)}
    size = 1 << t
    rows = [[0]*(size*size) for _ in types]
    cols = [[0]*(size*size) for _ in types]
    for x, y, u, v in product(range(size), repeat=4):
        polynomial = {(0, 0): 1}
        for bit in range(t):
            a, b, c, d = ((value >> bit)&1 for value in (x, y, u, v))
            values = [int(a == b and c == d), int(a == c and b == d), int(a^b^c^d == 0)]
            updated = {}
            for (r, s), value in polynomial.items():
                if not value:
                    continue
                for j, enabled in enumerate(values):
                    if enabled:
                        key = r+int(j == 0), s+int(j == 1)
                        updated[key] = updated.get(key, 0)^1
            polynomial = updated
        for key, value in polynomial.items():
            if value:
                i = indices[key]
                rows[i][size*x+y] |= 1 << (size*v+u)
                cols[i][size*x+u] |= 1 << (size*v+y)
    return types, rows, cols


def verify():
    full_tensors = 0
    for t in range(1, 4):
        types, generators, other_generators = literal_types(t)
        size = 1 << t
        row, col = [0]*(size*size), [0]*(size*size)
        old = 0
        for index in range(1, 1 << len(types)):
            mask = index^(index >> 1)
            changed = (old^mask).bit_length()-1
            old = mask
            row = [a^b for a, b in zip(row, generators[changed])]
            col = [a^b for a, b in zip(col, other_generators[changed])]
            active = [typ for j, typ in enumerate(types) if (mask >> j)&1]
            max_s = max(s for r, s, f in active)
            min_r = min(r for r, s, f in active if s == max_s)
            max_r = max(r for r, s, f in active)
            min_s = min(s for r, s, f in active if r == max_r)
            assert binary_rank(row) >= split_rank(t, t-max_s, min_r)
            assert binary_rank(col) >= split_rank(t, t-max_r, min_s)
            full_tensors += 1

    walk_cases = filtration_tail_cases = rotated_cases = 0
    for t in range(1, 49):
        for degree in range(t+1):
            endpoint = walk_count(t, degree, 0)
            literal = coefficient_rank(t, t-degree)
            assert split_rank(t, degree, 0) == literal
            assert literal <= 4*endpoint
            for split in range(degree+1):
                count = walk_count(t, degree, split)
                assert count >= endpoint
                actual_rank = split_rank(t, degree, split)
                assert count <= (t+1)**2*actual_rank
                assert actual_rank <= 4*count
                assert literal <= 4*(t+1)**2*actual_rank
                walk_cases += 1
                filtration_tail_cases += 1
                if t <= 10:
                    assert count == walk_count_rotated(t, degree, split)
                    rotated_cases += 1

    # An ordered word need not satisfy the symmetric-space asymptotic rate.
    # R^(t/2) S^(t/2) has each rank 4^(t/2), product 4^t.
    t = 128
    non_symmetric_product = 1 << (2*t)
    a, b = t//2, t-t//2
    left_tail = sum(comb(2*t, j) for j in range(a+1))
    right_tail = sum(comb(2*t, j) for j in range(b+1))
    assert non_symmetric_product*16*(t+1)**6 < left_tail*right_tail
    return {'status': 'pass',
            'scope': 'all pieces in the permutation-symmetric three-matching span on the full cube',
            'literal_full_tensor_filtration_checks': full_tensors,
            'walk_unimodality_cases': walk_cases,
            'rank_tail_comparisons': filtration_tail_cases,
            'independent_rotated_walk_counts': rotated_cases,
            'nonsymmetric_word_negative_control_bits': t,
            'piece_count_restriction': None,
            'new_factoring_exponent_claimed': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
