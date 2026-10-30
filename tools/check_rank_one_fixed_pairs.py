#!/usr/bin/env python3
"""Exact finite controls for rank-one fixed pairs of rank-two channels.

The small exhaustive control uses all binary 4 by 4 linear operators, without
assuming a Kraus representation or using the theorem to generate candidates.
The larger controls are explicit constructions, not random searches.
"""
import argparse
import json
from pathlib import Path


def rank(vectors):
    pivots = {}
    for v in vectors:
        while v:
            k = v.bit_length() - 1
            if k not in pivots:
                pivots[k] = v
                break
            v ^= pivots[k]
    return len(pivots)


def apply(columns, vector):
    answer = 0
    while vector:
        bit = vector & -vector
        answer ^= columns[bit.bit_length()-1]
        vector ^= bit
    return answer


def compose(left, right):
    return [apply(left, v) for v in right]


def transpose(a, q):
    return sum(((a >> (q*i+j)) & 1) << (q*j+i)
               for i in range(q) for j in range(q))


def multiply(a, b, q):
    result = 0
    for i in range(q):
        for j in range(q):
            entry = 0
            for k in range(q):
                entry ^= ((a >> (q*i+k)) & 1) & ((b >> (q*k+j)) & 1)
            result |= entry << (q*i+j)
    return result


def physical_rank(a, q):
    return rank([(a >> (q*i)) & ((1 << q)-1) for i in range(q)])


def phi_rows(columns, q):
    return [sum(((columns[q*v+y] >> (q*x+u)) & 1) << (q*u+v)
                for u in range(q) for v in range(q))
            for x in range(q) for y in range(q)]


def fixed_pair_check(columns, a, b, q):
    bt = transpose(b, q)
    assert apply(columns, a) == a
    assert all(((bt & col).bit_count() & 1) == ((bt >> j) & 1)
               for j, col in enumerate(columns))
    assert not ((bt & a).bit_count() & 1)
    assert multiply(a, b, q) == multiply(b, a, q) == 0


def exhaustive_q2():
    rank_one = [a for a in range(16) if physical_rank(a, 2) == 1]
    invertible = low_rank = pairs = nonorthogonal = 0
    rank_histogram = {}
    for operator in range(1 << 16):
        columns = [(operator >> (4*j)) & 15 for j in range(4)]
        if rank(columns) != 4:
            continue
        invertible += 1
        r = rank(phi_rows(columns, 2))
        if r > 2:
            continue
        low_rank += 1
        rank_histogram[r] = rank_histogram.get(r, 0) + 1
        right = [a for a in rank_one if apply(columns, a) == a]
        left = []
        for b in rank_one:
            bt = transpose(b, 2)
            if all(((bt & col).bit_count() & 1) == ((bt >> j) & 1)
                   for j, col in enumerate(columns)):
                left.append(b)
        for a in right:
            for b in left:
                if (transpose(b, 2) & a).bit_count() & 1:
                    nonorthogonal += 1
                    continue
                fixed_pair_check(columns, a, b, 2)
                pairs += 1
    assert invertible == 20160
    return dict(all_operators=65536, invertible_operators=invertible,
                invertible_choi_rank_at_most_two=low_rank,
                choi_rank_histogram=rank_histogram,
                orthogonal_rank_one_fixed_pairs=pairs,
                nonorthogonal_fixed_pairs_outside_claim=nonorthogonal)


def anchored_q4():
    q = 4
    a, b, r = 522, 15319, 50673
    columns = []
    for j in range(q*q):
        xt = transpose(1 << j, q)
        columns.append(multiply(xt, r, q) ^ multiply(multiply(a, xt, q), b, q))
    s = [v ^ (1 << j) for j, v in enumerate(columns)]
    orbit = 1 << 5
    previous = 0
    for _ in range(16):
        previous, orbit = orbit, apply(s, orbit)
    assert previous == 1 and orbit == 0
    fixed_pair_check(columns, 1, 1 << 5, q)
    dyadic = []
    power = columns
    for _ in range(4):
        dyadic.append(rank(phi_rows(power, q)))
        power = compose(power, power)
    assert dyadic == [2, 12, 12, 16]
    return dict(q=q, formula='U(X)=X^T R+A X^T B',
                binary_row_major_masks=dict(A=a, B=b, R=r),
                nilpotence_index=16, right_fixed_mask=1,
                left_fixed_trace_matrix_mask=32, dyadic_choi_ranks=dyadic)


def cycle_join(q):
    assert q >= 2 and q & (q-1) == 0
    n = q*q
    permutation = list(range(q))
    swaps = []
    for k in range(q//2-1):
        a, b = q-1-k, k+1
        permutation[a], permutation[b] = permutation[b], permutation[a]
        swaps.append([a, b])
    state_map = []
    for i in range(q):
        for j in range(q):
            x, y = j, (i+1) % q
            if y == 0:
                x = permutation[x]
            state_map.append(q*x+y)
    visited = set()
    point = 0
    while point not in visited:
        visited.add(point)
        point = state_map[point]
    assert point == 0 and len(visited) == n
    columns = [1 << point for point in state_map]
    jmask = (1 << n)-1
    fixed_pair_check(columns, jmask, jmask, q)
    dyadic = []
    power = columns
    for _ in range(n.bit_length()-1):
        dyadic.append(rank(phi_rows(power, q)))
        power = compose(power, power)
    assert power == [1 << j for j in range(n)]
    assert dyadic[0] == (1 if q == 2 else 2)
    # A full binary n-cycle has one Jordan block: (I+U)^(n-1)e_0 is all ones.
    s = [v ^ (1 << j) for j, v in enumerate(columns)]
    vector = 1
    for _ in range(n-1):
        vector = apply(s, vector)
    assert vector == jmask and apply(s, vector) == 0
    return dict(q=q, controlled_row_swaps=swaps,
                full_cycle_length=n, nilpotence_index=n,
                right_and_left_fixed_physical_rank=1,
                dyadic_choi_ranks=dyadic)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = dict(status='exact structural controls; no factoring exponent improvement',
                  exhaustive_binary_q2=exhaustive_q2(),
                  nonpermutation_regular_q4=anchored_q4(),
                  controlled_cycle_family=[cycle_join(q) for q in (2,4,8,16)])
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
