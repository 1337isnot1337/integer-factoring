#!/usr/bin/env python3
"""Independent four-map contraction and exact-profile replay for paired GL5 moves."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import time

from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from paired_physical_transvections import transform
from replay_source_provenance import resolve_source_hashes
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles
from search_saturated_dual_directions import GAUGES, complete_profiles


def field_multiply(first, second):
    result = 0
    while second:
        if second & 1:
            result ^= first
        second >>= 1
        first <<= 1
        if first & 4:
            first ^= 7
    return result


def physical_maps(operations):
    matrices = [[[int(i == j) for j in range(5)] for i in range(5)] for leg in range(4)]
    for pair, target, source, coefficient in operations:
        assert pair in (0, 1) and target != source and 0 <= target < 5 and 0 <= source < 5
        assert coefficient in (1, 2, 3)
        for leg, i, j in ((pair, target, source), (3-pair, source, target)):
            matrices[leg][i] = [a ^ field_multiply(coefficient, b)
                                for a, b in zip(matrices[leg][i], matrices[leg][j])]
    for pair in (0, 1):
        first, second = matrices[pair], matrices[3-pair]
        for i in range(5):
            for j in range(5):
                value = 0
                for k in range(5):
                    value ^= field_multiply(first[i][k], second[j][k])
                assert value == int(i == j)
    return matrices


def literal_contraction(base, operations):
    matrices = physical_maps(operations)
    supports = [[[(j, c) for j, c in enumerate(row) if c] for row in matrix] for matrix in matrices]
    result = []
    for index in range(625):
        letters = [(index//5**leg) % 5 for leg in range(4)]
        low = high = 0
        for terms in product(*(supports[leg][letters[leg]] for leg in range(4))):
            coefficient = 1
            old_index = 0
            for leg, (letter, value) in enumerate(terms):
                coefficient = field_multiply(coefficient, value)
                old_index += letter*5**leg
            old_low, old_high = base[old_index]
            # Use independent scalar multiplication on the two coefficient
            # basis vectors rather than the production transvection update.
            first, second = field_multiply(coefficient, 1), field_multiply(coefficient, 2)
            if first & 1:
                low ^= old_low
            if first & 2:
                high ^= old_low
            if second & 1:
                low ^= old_high
            if second & 2:
                high ^= old_high
        x, y, u, v = letters
        assert low.bit_count() % 2 == int(x == v and y == u) and high.bit_count() % 2 == 0
        assert max(low, high).bit_length() <= 19
        result.append((low, high))
    return tuple(result)


def profiles_by_frames(values, energy, claimed):
    half, zero = tuple(value//2 for value in energy), (0,)*5
    heights = (zero, zero, half, half, half, half, zero, zero)
    for side in range(2):
        matrix, scalar = literal_matrix(values, 5, 18, heights, side)
        certify_profile(matrix, scalar, claimed[side])


def controls():
    rng = random.Random(2026103212)
    cases = [((0, 0, 1, 2),), ((1, 4, 0, 1),),
             ((0, 0, 4, 2), (1, 3, 1, 3))]
    for length in (3, 5, 8, 12, 16):
        operations = []
        for step in range(length):
            target, source = rng.sample(range(5), 2)
            operations.append((rng.randrange(2), target, source, rng.randrange(1, 4)))
        cases.append(tuple(operations))
    base = five_letter_source()[0]
    coefficients = profile_count = frobenius_checks = 0
    for operations in cases:
        values = transform(base, 5, operations)
        assert values == literal_contraction(base, operations)
        coefficients += 625
        conjugated = tuple((pair, target, source, field_multiply(c, c)) for pair, target, source, c in operations)
        conjugate_values = literal_contraction(base, conjugated)
        assert conjugate_values == tuple((low ^ high, high) for low, high in values)
        frobenius_checks += 625
        profiles = complete_profiles(values, {'natural': GAUGES['natural']})['natural']
        assert complete_profiles(conjugate_values, {'natural': GAUGES['natural']})['natural'] == profiles
        profiles_by_frames(values, GAUGES['natural'], profiles)
        profile_count += 2
    return {'scope': 'Independent paired inverse identities, literal four-map contractions, field Frobenius, and polynomial Smith frames.',
            'status': 'passed', 'operation_chains': len(cases),
            'independent_literal_coefficient_identities': coefficients,
            'literal_Frobenius_coefficient_identities': frobenius_checks,
            'complete_profiles_certified_by_polynomial_frames': profile_count,
            'paired_inverse_matrix_entries': 50*len(cases),
            'operations_per_chain': [len(ops) for ops in cases]}


def replay(path):
    record = json.loads(path.read_text())
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    base = tuple(map(tuple, record['reference_values']))
    assert digest(base) == record['reference_digest'] and base == five_letter_source()[0]
    queries = record['query_records']
    counters = record['counters']
    assert counters['single_queries'] == sum(row['stage'] == 'single' for row in queries)
    assert counters['pair_queries'] == sum(row['stage'] == 'pair' for row in queries)
    assert counters['literal_endpoint_entries'] == 625*(len(queries)+1)
    expected = []
    moves = [(i, j) for i in range(5) for j in range(5) if i != j]
    for pair in (0, 1):
        for i, j in moves:
            for c in (1, 2):
                expected.append(([[pair, i, j, c]], 'single'))
    for i, j in moves:
        for k, ell in moves:
            for a, b in ((1, 1), (1, 2), (2, 1), (2, 2), (2, 3)):
                expected.append(([[0, i, j, a], [1, k, ell, b]], 'pair'))
    assert all(row['operations'] == operations and row['stage'] == stage
               for row, (operations, stage) in zip(queries, expected))
    if record['completed_requested_schedule']:
        assert len(queries) == len(expected) == 2080
    lookup = {tuple(map(tuple, row['operations'])): row for row in queries}
    literal_cache = {}
    endpoints = profiles = retained = 0
    gauges = {}
    for name, gauge in record['gauges'].items():
        energy = gauge['energies']
        profiles_by_frames(base, energy, gauge['reference_profiles'])
        same_cost(root_for_profiles(gauge['reference_profiles'], energy, energy), gauge['reference_cost'])
        profiles += 2
        endpoints += 625
        costs = [gauge['reference_cost']['rho']]
        for witness in gauge['retained']:
            key = tuple(map(tuple, witness['operations']))
            if key not in literal_cache:
                literal_cache[key] = literal_contraction(base, witness['operations'])
            values = literal_cache[key]
            assert values == tuple(map(tuple, witness['values'])) and digest(values) == witness['digest']
            row = lookup[key]
            assert row['digest'] == witness['digest'] and row['stage'] == witness['stage']
            assert row['gauges'][name]['profiles'] == witness['profiles']
            profiles_by_frames(values, energy, witness['profiles'])
            cost = root_for_profiles(witness['profiles'], energy, energy)
            same_cost(cost, witness)
            costs.append(cost['rho'])
            profiles += 2
            endpoints += 625
            retained += 1
        recorded_costs = [row['gauges'][name]['rho'] for row in queries]+[gauge['reference_cost']['rho']]
        assert abs(min(costs)-gauge['best_recorded_rho']) < 2e-8
        assert abs(min(recorded_costs)-gauge['best_recorded_rho']) < 2e-8
        gauges[name] = {'best_rho': gauge['best_recorded_rho'], 'retained': len(gauge['retained']),
                        'distinct_profile_classes': gauge['distinct_profile_classes']}
    assert record['strict_pressure_candidate_found'] == any(row['critical_pressure'] < -1e-9
        for gauge in record['gauges'].values() for row in gauge['retained'])
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(), 'status': 'passed',
            'source_provenance': provenance, 'completed_requested_schedule': record['completed_requested_schedule'],
            'recorded_query_counters': counters, 'retained_literal_witnesses': retained,
            'independent_four_map_contractions': len(literal_cache), 'literal_endpoint_entries': endpoints,
            'complete_profiles_certified_by_polynomial_frames': profiles, 'gauges': gauges,
            'strict_pressure_candidate_found': record['strict_pressure_candidate_found']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    result = replay(args.search) if args.search else controls()
    result.update(checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
