#!/usr/bin/env python3
"""Independently replay constant and polynomial F4 Smith-flag records."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_block_smith_crossings import binary_rank
from check_f4_recentered_lattices import scalar_product
from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from f4_recentered_lattices import multiply_matrix
from search_f4_recentered_crossings import digest, root_for_profiles
from replay_source_provenance import resolve_source_hashes


def binary_field_rank(matrix):
    rows = []
    for row in matrix:
        low = high = 0
        for j, coefficient in enumerate(row):
            for basis in range(2):
                value = scalar_product(coefficient, 1 << basis)
                low |= (value & 1) << (2*j+basis)
                high |= (value >> 1) << (2*j+basis)
        rows.extend((low, high))
    return binary_rank(rows)


def scale_polynomial(value, coefficient):
    result = power = 0
    while value:
        result |= scalar_product(value & 3, coefficient) << (2*power)
        value >>= 2
        power += 1
    return result


def certify_module(matrix, weights):
    n = len(weights)
    assert len(matrix) == n and all(len(row) == n for row in matrix)
    constant = [[0]*n for _ in range(n)]
    for i, row in enumerate(matrix):
        for j, entry in enumerate(row):
            assert isinstance(entry, int) and entry >= 0
            difference = weights[i]-weights[j]
            if difference > 0:
                assert not (entry & ((1 << (2*difference))-1))
            if difference >= 0:
                constant[i][j] = (entry >> (2*difference)) & 3
    assert binary_field_rank(constant) == 2*n


def replay(path):
    started = time.monotonic()
    record = json.loads(path.read_text())
    kind = record['kind']
    assert kind in ('equal-order-F4-flag-search', 'polynomial-F4-flag-search')
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    q, degree = record['alphabet'], record['degree']
    n = q*q
    heights, energy = record['heights_infinity_then_zero'], record['paired_energies']
    floors, caps, orders = record['floor_targets'], record['ceiling_targets'], record['reference_targets']
    assert len(floors) == len(caps) == len(orders) == 2*n
    assert digest(record['reference_values']) == record['reference_digest']
    frames = []
    for side in range(2):
        matrix, scalar = literal_matrix(record['reference_values'], q, degree, heights, side)
        assert tuple(value+scalar for value in record['reference_profiles'][side]) == tuple(orders[side*n:(side+1)*n])
        frames.append(certify_profile(matrix, scalar, record['reference_profiles'][side], max(caps[side*n:(side+1)*n])))
    same_cost(root_for_profiles(record['reference_profiles'], *energy), record['reference_cost'])
    unit_certificates = 0
    for index, flag in enumerate(record['flags']):
        assert flag['flag'] == index
        assert 0 <= flag['completed_trials'] <= record['arguments']['trials_per_flag']
        matrices = flag['matrices']
        assert len(matrices) == 2
        if kind == 'polynomial-F4-flag-search':
            actual = [[[int(i == j) for j in range(n)] for i in range(n)] for _ in range(2)]
            for side, i, j, shift, coefficient in flag['operations']:
                assert side in (0, 1) and 0 <= i < n and 0 <= j < n and i != j
                assert coefficient in (1, 2, 3)
                a, b = side*n+i, side*n+j
                assert shift >= max(0, orders[a]-orders[b], floors[a]-floors[b], caps[a]-caps[b])
                mask = (1 << (2*max(caps[side*n:(side+1)*n])))-1
                actual[side][i] = [(first ^ (scale_polynomial(second, coefficient) << (2*shift))) & mask
                                   for first, second in zip(actual[side][i], actual[side][j])]
            assert actual == matrices
            assert flag['nonbinary'] == any(operation[4] in (2, 3) for operation in flag['operations'])
        else:
            for side, matrix in enumerate(matrices):
                assert all(0 <= entry < 4 for row in matrix for entry in row)
                assert all(not matrix[i][j] or orders[side*n+i] == orders[side*n+j]
                           for i in range(n) for j in range(n))
            assert flag['nonbinary'] == any(c in (2, 3) for matrix in matrices for row in matrix for c in row)
        for side, matrix in enumerate(matrices):
            for weights in ((0,)*n, floors[side*n:(side+1)*n], orders[side*n:(side+1)*n]):
                certify_module(matrix, weights)
                unit_certificates += 1
    targets_checked = profiles_checked = 0
    costs = [record['reference_cost']['rho']]
    for witness in record['retained_witnesses']:
        assert digest(witness['values']) == witness['digest']
        targets = witness['targets']
        assert len(targets) == 2*n and all(f <= t <= c for f, t, c in zip(floors, targets, caps))
        flag = record['flags'][witness['flag']]
        for side, (left, right) in enumerate(frames):
            matrix, scalar = literal_matrix(witness['values'], q, degree, heights, side)
            certify_profile(matrix, scalar, witness['profiles'][side])
            profiles_checked += 1
            precision = max(caps[side*n:(side+1)*n])
            changed_left = multiply_matrix(flag['matrices'][side], left, precision)
            transformed = multiply_matrix(multiply_matrix(changed_left, matrix), right)
            for row, target in zip(transformed, targets[side*n:(side+1)*n]):
                assert all(not (entry & ((1 << (2*target))-1)) for entry in row)
            targets_checked += 1
        cost = root_for_profiles(witness['profiles'], *energy)
        same_cost(cost, witness)
        costs.append(cost['rho'])
    assert sum(item['samples'] for item in record['profile_classes']) == record['complete_profile_queries']
    assert sum(flag['completed_trials'] for flag in record['flags']) == record['completed_target_walks']
    assert abs(min(costs)-record['best_recorded_rho']) < 2e-8
    assert record['strict_pressure_candidate_found'] == any(w['critical_pressure'] < -1e-9 for w in record['retained_witnesses'])
    if record['completed_requested_schedule']:
        assert len(record['flags']) == record['arguments']['flags']
        assert all(flag['completed_requested_trials'] and flag['completed_trials'] == record['arguments']['trials_per_flag']
                   for flag in record['flags'])
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(), 'status': 'passed',
            'source_provenance': provenance,
            'kind': kind, 'flags': len(record['flags']),
            'binary_rank_certificates_for_regular_modules': unit_certificates,
            'complete_profiles_certified_by_polynomial_frames': profiles_checked+2,
            'untruncated_target_products': targets_checked,
            'literal_endpoint_entries': (1+len(record['retained_witnesses']))*q**4,
            'retained_witnesses': len(record['retained_witnesses']),
            'completed_target_walks': record['completed_target_walks'],
            'complete_profile_queries': record['complete_profile_queries'],
            'best_recorded_rho': record['best_recorded_rho'],
            'strict_pressure_candidate_found': record['strict_pressure_candidate_found'],
            'elapsed_seconds': time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = []
    for path in args.records:
        item = replay(path)
        results.append(item)
        print(json.dumps(item), flush=True)
    if args.output:
        args.output.write_text(json.dumps({'scope': 'Exact literal row-frame and complete-profile replays; numerical cost minimizations are exploratory.',
                                         'records': results, 'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2)+'\n')


if __name__ == '__main__':
    main()
