#!/usr/bin/env python3
"""Fast independent controls for two completed binary q8 cycle searches.

This does not rerun the exhaustive native searches. It verifies their class
sizes, primitive feedbacks, and explicit best low-leading-rank controls.
"""
import argparse
import itertools
import json
from pathlib import Path


def rank(vectors):
    pivots = {}
    for v in vectors:
        while v:
            k = v.bit_length()-1
            if k not in pivots:
                pivots[k] = v
                break
            v ^= pivots[k]
    return len(pivots)


def xor(values):
    out = 0
    for v in values:
        out ^= v
    return out


def basis_map(row_basis, column_basis):
    assert rank(row_basis+column_basis) == 6
    embedding = [xor([row_basis[k] for k in range(3) if i >> k & 1]
                     +[column_basis[k] for k in range(3) if j >> k & 1])
                 for i in range(8) for j in range(8)]
    inverse = [0]*64
    for i, v in enumerate(embedding):
        inverse[v] = i
    return embedding, inverse


def choi(permutation):
    rows = [0]*64
    for j, output in enumerate(permutation):
        v, y = divmod(j, 8)
        x, u = divmod(output, 8)
        rows[8*x+y] |= 1 << (8*u+v)
    return rows


def check_control(state_map, row_basis, column_basis):
    embedding, inverse = basis_map(row_basis, column_basis)
    u = [inverse[state_map[x]] for x in embedding]
    visited = set()
    x = 0
    while x not in visited:
        visited.add(x)
        x = u[x]
    assert x == 0 and len(visited) == 64
    a = choi(u)
    c = choi([u[u[j]] for j in range(64)])
    r, r2 = rank(a), rank(c)
    first_jet = a+[(x ^ y) | (x << 64) for x, y in zip(a, c)]
    ones = rank(first_jet)-2*r
    return dict(full_cycle_length=64, leading_rank=r, square_rank=r2,
                order_one_count=ones, first_jet_cost_floor=64*r+63*ones)


def primitive_taps():
    output = []
    for taps in range(1,64,2):
        x, length = 1, 0
        while True:
            x = (x >> 1) | (((x & taps).bit_count() & 1) << 5)
            length += 1
            if x == 1:
                break
        if length == 63:
            output.append(taps)
    return output


def subspace_count():
    spaces = set()
    for a, b, c in itertools.combinations(range(1,64), 3):
        if c == a ^ b:
            continue
        spaces.add(frozenset([0,a,b,c,a^b,a^c,b^c,a^b^c]))
    assert len(spaces) == 1395
    return len(spaces)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    taps = primitive_taps()
    assert taps == [3,27,33,39,45,51]
    spaces = subspace_count()
    feedback, gray_index = 1, 103 ^ (103 >> 1)
    monomials = [m for m in range(32) if m.bit_count() <= 2]
    for j, m in enumerate(monomials):
        if gray_index >> j & 1:
            feedback ^= sum(1 << x for x in range(32) if x & m == m)
    assert feedback == 4244438269
    debruijn_map = [(x >> 1) | (((x & 1) ^ ((feedback >> (x >> 1)) & 1)) << 5)
                    for x in range(64)]
    debruijn_control = check_control(debruijn_map, [1,4,8], [2,16,32])
    assert debruijn_control == dict(full_cycle_length=64, leading_rank=5,
                                   square_rank=8, order_one_count=7,
                                   first_jet_cost_floor=761)
    linear = [(x >> 1) | (((x & 3).bit_count() & 1) << 5) for x in range(64)]
    linear[0], linear[1] = 32, 0
    lfsr_control = check_control(linear, [1,2,8], [4,16,32])
    assert lfsr_control == dict(full_cycle_length=64, leading_rank=5,
                               square_rank=9, order_one_count=8,
                               first_jet_cost_floor=824)
    result = dict(status='completed finite searches; no local Smith-profile survivor',
                  default_check_scope='class sizes and representative witnesses; full native reruns are optional',
                  quadratic_feedback=dict(feedback_choices=65536,
                      full_cycles=1952, balanced_partitions=10,
                      full_cycle_partition_cases=19520,
                      leading_rank_histogram={'3':80,'4':752,'5':1406,'at_least_6':17282},
                      first_jet_survivors=0, representative_feedback_mask=feedback,
                      representative_partition_mask=13,
                      representative=debruijn_control),
                  linear_gauges=dict(primitive_taps=taps,
                      row_subspaces=spaces, complements_per_subspace=512,
                      ordered_decompositions=spaces*512,
                      cases=spaces*512*len(taps),
                      conservative_linear_prefilter_passes=258552,
                      leading_rank_at_most_five=3360,
                      leading_rank_histogram={'3':120,'4':96,'5':3144},
                      minimum_first_jet_floor_by_leading_rank={'3':2145,'4':2272,'5':824},
                      first_jet_survivors=0, representative_taps=3,
                      representative_row_basis=[1,2,8], representative_column_basis=[4,16,32],
                      representative=lfsr_control))
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
