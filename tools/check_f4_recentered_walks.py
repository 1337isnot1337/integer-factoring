#!/usr/bin/env python3
"""Replay saved nonbinary-frame searches by literal polynomial identities.

This checker does not solve the search's affine systems or use the native
Smith profile oracle.  It reconstructs exact polynomial regular frames
from each recorded complete profile, verifies full untruncated row products,
and separately repeats the exploratory numerical cost calculation.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from f4_recentered_lattices import multiply_matrix, regular_frame
from search_f4_recentered_crossings import digest, root_for_profiles
from replay_source_provenance import resolve_source_hashes


ROOT = Path(__file__).resolve().parent.parent


def literal_matrix(values, q, degree, heights, side):
    """Independently assemble the two packed F4 polynomial matrices."""
    h = heights[4*side:4*side+4]
    scalar = (degree if side == 0 else 0)-sum(min(row) for row in h)
    matrix = [[0]*(q*q) for _ in range(q*q)]
    assert len(values) == q**4
    for index, pair in enumerate(values):
        assert len(pair) == 2 and all(isinstance(a, int) and 0 <= a < 1 << (degree+1) for a in pair)
        x, y, u, v = bits = tuple((index//q**leg) % q for leg in range(4))
        assert (pair[0].bit_count() & 1) == int(x == v and y == u)
        assert pair[1].bit_count() % 2 == 0
        i, j = (q*x+u, q*y+v) if side == 0 else (q*x+y, q*u+v)
        height = sum(h[leg][letter] for leg, letter in enumerate(bits))
        for power in range(degree+1):
            coefficient = ((pair[0] >> power) & 1) | (((pair[1] >> power) & 1) << 1)
            if coefficient:
                at = scalar+height+(-power if side == 0 else power)
                assert at >= 0
                matrix[i][j] ^= coefficient << (2*at)
    return matrix, scalar


def certify_profile(matrix, scalar, claimed, precision=None):
    assert all(isinstance(order, int) for order in claimed)
    expected = tuple(order+scalar for order in claimed)
    assert len(expected) == len(matrix) and list(expected) == sorted(expected) and min(expected) >= 0
    precision = max(expected)+1 if precision is None else precision
    actual, left, right, diagonal = regular_frame(matrix, precision)
    assert actual == expected
    return left, right


def same_cost(actual, claimed):
    for name in ('rho', 'scale', 'critical_pressure'):
        assert abs(actual[name]-claimed[name]) < 2e-8*max(1., abs(actual[name]))


def replay(path, allow_source_mismatch=False):
    started = time.monotonic()
    record = json.loads(path.read_text())
    q, degree = record['alphabet'], record['degree']
    heights = record['heights_infinity_then_zero']
    energy_x, energy_y = record['paired_energies']
    assert len(heights) == 8 and all(len(row) == q for row in heights)
    assert len(energy_x) == len(energy_y) == q
    provenance, mismatches = resolve_source_hashes(record['source_sha256'])
    if not allow_source_mismatch:
        assert not mismatches, f'Recorded source hashes changed: {mismatches}'
    frames = {}
    certified_profiles = endpoints = target_products = 0
    queries = 0
    reference_costs = []
    for i, item in enumerate(record['rounds']):
        assert item['round'] == i
        values = item['reference_values']
        assert digest(values) == item['reference_digest']
        assert item['nonbinary_reference'] == any(pair[1] for pair in values)
        if i:
            assert record['rounds'][i-1]['next_reference_digest'] == item['reference_digest']
        assert 0 <= item['completed_trials'] <= record['arguments']['trials_per_round']
        queries += sum(row['samples'] for row in item['profile_classes'])
        per_side = []
        for side in range(2):
            matrix, scalar = literal_matrix(values, q, degree, heights, side)
            cap = max(item['ceiling_targets'][side*q*q:(side+1)*q*q])
            per_side.append(certify_profile(matrix, scalar, item['reference_profiles'][side], cap))
            certified_profiles += 1
        frames[i] = per_side
        endpoints += q**4
        cost = root_for_profiles(item['reference_profiles'], energy_x, energy_y)
        same_cost(cost, item['reference_cost'])
        reference_costs.append(cost['rho'])

    witness_costs = []
    for witness in record['retained_witnesses']:
        values, targets = witness['values'], witness['targets']
        assert digest(values) == witness['digest']
        source = record['rounds'][witness['round']]
        assert len(targets) == 2*q*q
        assert all(floor <= target <= cap for floor, target, cap in
                   zip(source['floor_targets'], targets, source['ceiling_targets']))
        for side in range(2):
            matrix, scalar = literal_matrix(values, q, degree, heights, side)
            certify_profile(matrix, scalar, witness['profiles'][side])
            certified_profiles += 1
            left, right = frames[witness['round']][side]
            # No jet truncation is used in this independent witness replay.
            transformed = multiply_matrix(multiply_matrix(left, matrix), right)
            for row, target in zip(transformed, targets[side*q*q:(side+1)*q*q]):
                mask = (1 << (2*target))-1
                assert all(not (entry & mask) for entry in row)
            target_products += 1
        endpoints += q**4
        cost = root_for_profiles(witness['profiles'], energy_x, energy_y)
        same_cost(cost, witness)
        witness_costs.append(cost['rho'])

    assert queries == record['complete_profile_queries']
    if reference_costs:
        assert abs(min(reference_costs+witness_costs)-record['best_recorded_rho']) < 2e-8
    assert record['strict_pressure_candidate_found'] == any(
        w['critical_pressure'] < -1e-9 for w in record['retained_witnesses'])
    if record['completed_requested_schedule']:
        assert len(record['rounds']) == record['arguments']['rounds']
        assert all(r['completed_requested_trials'] and r['completed_trials'] == record['arguments']['trials_per_round']
                   for r in record['rounds'])
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'status': 'passed', 'source_hash_mismatches': mismatches,
            'source_provenance': provenance,
            'references': len(record['rounds']),
            'nonbinary_references': sum(r['nonbinary_reference'] for r in record['rounds']),
            'retained_literal_witnesses': len(record['retained_witnesses']),
            'complete_profiles_certified_by_polynomial_frames': certified_profiles,
            'untruncated_target_products': target_products, 'literal_endpoint_entries': endpoints,
            'recorded_complete_profile_queries': queries,
            'completed_requested_schedule': record['completed_requested_schedule'],
            'best_recorded_rho': record['best_recorded_rho'],
            'strict_pressure_candidate_found': record['strict_pressure_candidate_found'],
            'elapsed_seconds': time.monotonic()-started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records', nargs='+', type=Path)
    parser.add_argument('--allow-source-mismatch', action='store_true',
                        help='Explicitly permit replay of historical records with newer source files.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = []
    for path in args.records:
        item = replay(path, args.allow_source_mismatch)
        results.append(item)
        print(json.dumps(item), flush=True)
    if args.output:
        report = {'scope': 'Exact literal endpoints, complete polynomial Smith frames, untruncated target products, and numerical cost replays; no new factoring bound.',
                  'records': results,
                  'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
        args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
