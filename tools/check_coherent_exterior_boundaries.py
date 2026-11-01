#!/usr/bin/env python3
"""Exact checks for coherent exterior boundaries and their common projection.

The complete binary screen is finite evidence, not an asymptotic assertion.
Only the Python standard library and a C++17 compiler are required.
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile

STATES = (3, 5, 6, 9, 10, 12)
COORDINATES = tuple((i, i) for i in range(6)) + ((2, 3), (3, 2))


def rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def polynomial_product(a, b):
    answer = 0
    while b:
        bit = b & -b
        answer ^= a << (bit.bit_length() - 1)
        b ^= bit
    return answer


def local_operators():
    """Construct c_q a_p in the column convention from exterior subsets."""
    index = {s: i for i, s in enumerate(STATES)}
    result = {}
    for p, q in itertools.product(range(4), repeat=2):
        matrix = [[int(p == q and i == j) for j in range(6)]
                  for i in range(6)]
        for column, state in enumerate(STATES):
            if not (state >> p & 1):
                continue
            removed = state ^ (1 << p)
            if removed >> q & 1:
                continue
            row = index[removed | (1 << q)]
            matrix[row][column] ^= 3  # (1+z)c_q a_p.
        result[p, q] = matrix
    return result


def direct_product(t, p, q, operators):
    matrix = [[int(i == j) for j in range(6)] for i in range(6)]
    for position in range(t):
        local = operators[p >> (2*position) & 3, q >> (2*position) & 3]
        answer = [[0]*6 for _ in range(6)]
        for i, k, j in itertools.product(range(6), repeat=3):
            if matrix[i][k] and local[k][j]:
                answer[i][j] ^= polynomial_product(matrix[i][k], local[k][j])
        matrix = answer
    return matrix


def zip_bits(t, x, y):
    return sum(((x >> i & 1) + 2*(y >> i & 1)) << (2*i) for i in range(t))


def evaluate_boundary(matrix, mask):
    result = 0
    for k, (i, j) in enumerate(COORDINATES):
        if mask >> k & 1:
            result ^= matrix[i][j]
    return result


def tensor_profiles(t, boundary, matrices):
    d = 1 << t
    row = [[0]*(d*d) for _ in range(d*d)]
    column = [[0]*(d*d) for _ in range(d*d)]
    for x, y, u, v in itertools.product(range(d), repeat=4):
        value = evaluate_boundary(matrices[zip_bits(t, x, y), zip_bits(t, v, u)], boundary)
        row[x*d+y][u*d+v] = value
        column[x*d+u][y*d+v] = value
    profiles = []
    for cut in (row, column):
        profiles.append([rank([sum((((value & functional).bit_count() & 1) << j)
                                   for j, value in enumerate(values)) for values in cut])
                         for functional in range(1, 1 << (t+1))])
    return profiles


def occupation_marginals(state):
    return (sum(state >> p & 1 for p in (1, 3)),
            sum(state >> p & 1 for p in (2, 3)))


def check_projection(matrices, t):
    rng = random.Random(903)
    checks = 0
    changed = 0
    matching_survives = False
    for trial in range(16):
        boundary = [[rng.randrange(2) for _ in range(6)] for _ in range(6)]
        boundary[0][0] ^= (sum(boundary[i][i] for i in range(6)) & 1) ^ 1
        projected = [[boundary[i][j] if occupation_marginals(STATES[i]) ==
                      occupation_marginals(STATES[j]) else 0
                      for j in range(6)] for i in range(6)]
        assert sum(projected[i][i] for i in range(6)) % 2 == 1
        for p, q in matrices:
            matrix = matrices[p, q]
            raw = reduced = 0
            for i, j in itertools.product(range(6), repeat=2):
                if boundary[i][j]:
                    raw ^= matrix[i][j]
                if projected[i][j]:
                    reduced ^= matrix[i][j]
            x = sum((p >> (2*i) & 1) << i for i in range(t))
            y = sum((p >> (2*i+1) & 1) << i for i in range(t))
            v = sum((q >> (2*i) & 1) << i for i in range(t))
            u = sum((q >> (2*i+1) & 1) << i for i in range(t))
            selected = x.bit_count() == v.bit_count() and y.bit_count() == u.bit_count()
            assert (raw if selected else 0) == reduced
            changed += raw != reduced
            if selected and matrix[2][3]:
                matching_survives = True
            checks += 1
    assert changed and matching_survives
    return checks, changed


def incoming_edges(S):
    incoming = [0]*4
    for a, b in itertools.product(range(4), repeat=2):
        x, u, v, y = a % 2, a // 2, b % 2, b // 2
        if a != b and (S >> (v+2*u) & 1) and not (S >> (x+2*y) & 1):
            incoming[a] |= 1 << b
    return incoming


def acyclic(incoming, subset):
    left = subset
    while left:
        choices = [i for i in range(4) if left >> i & 1 and not (incoming[i] & left)]
        if not choices:
            return False
        left ^= 1 << choices[0]
    return True


def check_union_bounds(t, profiles):
    checks = 0
    for state_index, state in enumerate(STATES):
        row, col = profiles[1 << state_index]
        for subset in range(1, 16):
            if not acyclic(incoming_edges(state), subset):
                continue
            n, s = subset.bit_count(), (subset & state).bit_count()
            for functional in range(1, 1 << (t+1)):
                rr = cc = 0
                for k in range(t+1):
                    if functional >> k & 1:
                        rr += math.comb(t, k) * 2**t
                        cc += math.comb(t, k) * s**k * (n-s)**(t-k)
                assert row[functional-1] >= rr
                assert col[functional-1] >= cc
                checks += 1
    return checks


def shortest_integer_lower_bound(profile, precision=16):
    """Exact rational lower costs certify every binary coefficient decomposition."""
    rows, columns = profile['ranks']
    costs = [math.isqrt((r*c) << (2*precision)) for r, c in zip(rows, columns)]
    count = len(costs)+1
    distance = [None]*count
    distance[0] = 0
    queue = [(0, 0)]
    while queue:
        value, state = heapq.heappop(queue)
        if value != distance[state]:
            continue
        for functional, cost in enumerate(costs, 1):
            following = state ^ functional
            candidate = value + cost
            if distance[following] is None or candidate < distance[following]:
                distance[following] = candidate
                heapq.heappush(queue, (candidate, following))
    return distance[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output')
    parser.add_argument('--screen-output')
    args = parser.parse_args()
    source = Path(__file__).with_name('coherent_exterior.cpp')
    operators = local_operators()
    direct_comparisons = union_checks = projection_checks = projection_changes = 0
    with tempfile.TemporaryDirectory(prefix='coherent-exterior-') as temporary:
        temporary = Path(temporary)
        native = temporary/'profile'
        subprocess.run(['g++', '-O3', '-std=c++17', str(source), '-o', str(native)], check=True)
        for t in (1, 2, 3):
            destination = temporary/f't{t}.json'
            subprocess.run([str(native), str(t), 'all', 'span', str(destination)], check=True,
                           stdout=subprocess.DEVNULL)
            data = json.loads(destination.read_text())
            profiles = {p['boundary']: p['ranks'] for p in data['profiles']}
            matrices = {(p, q): direct_product(t, p, q, operators)
                        for p, q in itertools.product(range(1 << (2*t)), repeat=2)}
            boundaries = list(profiles) if t <= 2 else [1, 4, 7, 13, 65, 68, 129, 132, 193, 196]
            for boundary in boundaries:
                expected = tensor_profiles(t, boundary, matrices)
                assert expected == profiles[boundary]
                direct_comparisons += 2*((1 << (t+1))-1)
            union_checks += check_union_bounds(t, profiles)
            if t == 2:
                projection_checks, projection_changes = check_projection(matrices, t)
        destination = temporary/'screen.json'
        subprocess.run([str(native), '4', 'all', 'span', str(destination)], check=True,
                       stdout=subprocess.DEVNULL)
        screen = json.loads(destination.read_text())
        assert len(screen['profiles']) == 128
        target = 256 << 16
        for profile in screen['profiles']:
            assert profile['ranks'][0][-1] == profile['ranks'][1][-1] == 256
            assert shortest_integer_lower_bound(profile) == target
        if args.screen_output:
            Path(args.screen_output).write_text(json.dumps(screen, indent=2)+'\n')
    report = {
        'status': 'passed',
        'native_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'direct_cut_rank_comparisons': direct_comparisons,
        'common_projection_entry_checks': projection_checks,
        'projection_negative_control_changes': projection_changes,
        'triangular_union_rank_checks': union_checks,
        'screen': {
            'width': 4, 'labels': 16, 'normalized_boundaries': 128,
            'nonzero_coefficient_functionals_per_boundary': 31,
            'rank_pairs': 128*31, 'binary_span_optimum_for_every_boundary': 256,
            'lower_cost_denominator': 1 << 16,
            'scope': 'F2 coefficient span, finite width only; no extension-field or asymptotic exclusion',
        },
    }
    if args.output:
        Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
