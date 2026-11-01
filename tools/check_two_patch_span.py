#!/usr/bin/env python3
"""Exact binary recombination optimum for small two-patch trace spans."""
import argparse
import hashlib
import heapq
import json
import math
import os
from pathlib import Path

from check_trace_bilinear_kernel import Field, irreducibility_by_division, rank


def coefficient_kernels(t, polynomial):
    field = Field(polynomial)
    assert field.degree == t
    irreducibility_by_division(polynomial)
    d = 1 << t
    alpha_powers = [field.power(field.alpha, j) for j in range(t)]
    coordinates = [sum(field.trace(field.mul(a, z)) << j
                       for j, a in enumerate(alpha_powers))
                   for z in range(d)]
    assert len(set(coordinates)) == d
    kernels = [0] * (2*t)
    weights_seen = set()
    for s in range(d):
        for w in range(d):
            error = s ^ w
            b0 = field.mul(s, error)
            b1 = field.mul(s ^ 1, error)
            polynomial_bits = 1
            for value in (coordinates[b0], coordinates[b1]):
                for j in range(t):
                    if value >> j & 1:
                        polynomial_bits ^= polynomial_bits << 1
            assert polynomial_bits.bit_length() <= 2*t
            assert polynomial_bits.bit_count() % 2 == int(s == w)
            weights_seen.add(coordinates[b0].bit_count()
                             + coordinates[b1].bit_count())
            for j in range(2*t):
                if polynomial_bits >> j & 1:
                    kernels[j] |= 1 << (s*d+w)
    assert weights_seen == set(range(2*t))
    assert rank(kernels) == 2*t
    tensors = [0] * (1 << (2*t))
    for mask in range(1, len(tensors)):
        bit = mask & -mask
        tensors[mask] = tensors[mask ^ bit] ^ kernels[bit.bit_length()-1]
    assert len(set(tensors)) == len(tensors)
    target = len(tensors)-1
    assert tensors[target] == sum(1 << (s*d+s) for s in range(d))
    return tensors, target


def cut_ranks(kernel, d, cache):
    def convolution_rank(bits):
        if bits not in cache:
            matrix = [sum(((bits >> (x ^ y)) & 1) << y
                          for y in range(d)) for x in range(d)]
            cache[bits] = rank(matrix)
        return cache[bits]
    row_mask = (1 << d)-1
    r = sum(convolution_rank((kernel >> (s*d)) & row_mask)
            for s in range(d))
    c = sum(convolution_rank(sum(((kernel >> (s*d+w)) & 1) << s
                                 for s in range(d))) for w in range(d))
    return [r, c]


def dense_cut_ranks(kernel, d):
    # Build literal four-leg cuts, without the block/circulant shortcut.
    rows, columns = [0]*(d*d), [0]*(d*d)
    for x in range(d):
        for y in range(d):
            for u in range(d):
                v = x ^ y ^ u
                if kernel >> ((x ^ y)*d+(x ^ u)) & 1:
                    rows[x*d+y] |= 1 << (u*d+v)
                    columns[x*d+u] |= 1 << (y*d+v)
    return [rank(rows), rank(columns)]


def shortest_paths(costs):
    n = len(costs)
    distances = [None]*n
    previous = [None]*n
    distances[0] = 0
    heap = [(0, 0)]
    while heap:
        distance, state = heapq.heappop(heap)
        if distance != distances[state]:
            continue
        for edge in range(1, n):
            successor = state ^ edge
            proposal = distance+costs[edge]
            if distances[successor] is None or proposal < distances[successor]:
                distances[successor] = proposal
                previous[successor] = (state, edge)
                heapq.heappush(heap, (proposal, successor))
    assert distances[0] == 0
    # A complete finite potential certificate independently verifies the
    # lower bound: every possible edge can increase potential by at most cost.
    for state in range(n):
        for edge in range(1, n):
            assert distances[state ^ edge] <= distances[state]+costs[edge]
    return distances, previous


def profile(t, polynomial):
    d = 1 << t
    tensors, target = coefficient_kernels(t, polynomial)
    cache = {0: 0}
    ranks = [cut_ranks(kernel, d, cache) for kernel in tensors]
    for kernel, expected in zip(tensors, ranks):
        assert dense_cut_ranks(kernel, d) == expected
    assert ranks[target] == [d*d, d*d]
    assert all(r > 0 and c > 0 for r, c in ranks[1:])
    # Begin with integer floors. Only refine if they fail to certify the
    # one-piece optimum. All comparisons below use exact integer arithmetic.
    denominator = 1
    while True:
        costs = [math.isqrt(r*c*denominator**2) for r, c in ranks]
        distances, previous = shortest_paths(costs)
        path = []
        state = target
        while state:
            parent, edge = previous[state]
            path.append(edge)
            state = parent
        assert len(path) <= 2*t
        mask = 0
        upper_numerator = 0
        for edge in path:
            mask ^= edge
            r, c = ranks[edge]
            square = r*c*denominator**2
            floor = math.isqrt(square)
            upper_numerator += floor + int(floor*floor != square)
        assert mask == target
        if distances[target] == d*d*denominator:
            status = 'exact_optimum_d_squared'
            assert upper_numerator == d*d*denominator
            break
        if upper_numerator < d*d*denominator:
            status = 'strict_improvement_over_one_piece'
            break
        if denominator >= 1 << 20:
            status = 'unresolved_interval'
            break
        denominator *= 2
    return {
        't': t,
        'field_polynomial': polynomial,
        'labels': d,
        'coefficient_span_dimension': 2*t,
        'span_tensors_checked': len(tensors),
        'dense_cut_comparisons': 2*len(tensors),
        'target_coefficient_mask': target,
        'rank_pairs_by_mask': ranks,
        'lower_cost_denominator': denominator,
        'lower_cost_numerators_by_mask': costs,
        'lower_potential_by_mask': distances,
        'potential_inequalities_checked': len(tensors)*(len(tensors)-1),
        'lower_optimum_numerator': distances[target],
        'upper_decomposition_masks': path,
        'upper_decomposition_rank_pairs': [ranks[edge] for edge in path],
        'upper_cost_numerator': upper_numerator,
        'one_piece_cost': d*d,
        'status': status,
        'cached_convolution_kernels': len(cache),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    profiles = [profile(3, 11), profile(4, 19)]
    here = Path(__file__).resolve().parent
    sources = [Path(__file__), here/'check_trace_bilinear_kernel.py']
    evidence = {
        'status': 'passed',
        'profiles': profiles,
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        'scope': 'Complete F2 coefficient spans at t=3,4 only. No extension-field or all-size exclusion and no factoring-score improvement.',
    }
    temporary = args.output.with_suffix(args.output.suffix+'.tmp')
    with temporary.open('w') as out:
        json.dump(evidence, out, indent=2, sort_keys=True)
        out.write('\n')
        out.flush()
        os.fsync(out.fileno())
    temporary.replace(args.output)
    assert len(json.loads(args.output.read_text())['profiles']) == 2
    print(json.dumps({'status': 'passed', 'profiles': [
        {k: p[k] for k in ('t', 'span_tensors_checked', 'status',
                           'lower_cost_denominator', 'lower_optimum_numerator',
                           'upper_cost_numerator', 'upper_decomposition_masks')}
        for p in profiles]}))


if __name__ == '__main__':
    main()
