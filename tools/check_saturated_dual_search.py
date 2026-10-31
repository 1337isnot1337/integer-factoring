#!/usr/bin/env python3
"""Replay literal sparse saturation witnesses and their complete Smith profiles."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_block_smith_crossings import binary_rank
from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from replay_source_provenance import resolve_source_hashes
from saturated_dual_kernel_directions import literal_directions
from search_f4_recentered_crossings import digest, root_for_profiles


ROOT = Path(__file__).resolve().parent.parent


def replay(path):
    started = time.monotonic()
    record = json.loads(path.read_text())
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    direction_path = ROOT/record['arguments']['directions']
    assert sha256(direction_path.read_bytes()).hexdigest() == record['direction_record_sha256']
    direction_record = json.loads(direction_path.read_text())
    directions = literal_directions(direction_record)
    indices = record['independent_direction_indices']
    encoded = [sum(low << (19*entry) for entry, (low, high) in enumerate(directions[index])) for index in indices]
    assert binary_rank(encoded) == len(indices) == record['parameter_dimension_over_F4']
    base = record['reference_values']
    assert digest(base) == record['reference_digest']
    endpoints = profiles_checked = witnesses = 0
    counts = {}
    for name, gauge in record['gauges'].items():
        energy = gauge['energies']
        half, zero = tuple(value//2 for value in energy), (0,)*5
        heights = (zero, zero, half, half, half, half, zero, zero)
        for side in range(2):
            matrix, scalar = literal_matrix(base, 5, 18, heights, side)
            certify_profile(matrix, scalar, gauge['reference_profiles'][side])
            profiles_checked += 1
        same_cost(root_for_profiles(gauge['reference_profiles'], energy, energy), gauge['reference_cost'])
        costs = [gauge['reference_cost']['rho']]
        endpoints += 625
        for witness in gauge['retained']:
            values = [list(pair) for pair in base]
            assert len({index for index, coefficient in witness['assignment']}) == len(witness['assignment'])
            for index, coefficient in witness['assignment']:
                assert index in indices and coefficient in (1, 2, 3)
                for entry, (low, high) in enumerate(directions[index]):
                    assert high == 0
                    if coefficient & 1:
                        values[entry][0] ^= low
                    if coefficient >> 1:
                        values[entry][1] ^= low
            assert values == witness['values'] and digest(values) == witness['digest']
            for side in range(2):
                matrix, scalar = literal_matrix(values, 5, 18, heights, side)
                certify_profile(matrix, scalar, witness['profiles'][side])
                profiles_checked += 1
            cost = root_for_profiles(witness['profiles'], energy, energy)
            same_cost(cost, witness)
            costs.append(cost['rho'])
            endpoints += 625
            witnesses += 1
        assert abs(min(costs)-gauge['best_recorded_rho']) < 2e-8
        counts[name] = {'retained': len(gauge['retained']), 'best_rho': gauge['best_recorded_rho'],
                        'distinct_profile_classes': gauge['distinct_profile_classes']}
    counters = record['counters']
    assert counters['single_queries'] <= 2*len(indices)
    assert counters['pair_queries'] <= 5*len(indices)*(len(indices)-1)//2
    assert counters['random_queries'] <= record['arguments']['random_queries']
    if record['completed_requested_schedule']:
        assert counters['single_queries'] == 2*len(indices)
        assert counters['pair_queries'] == record['planned_pair_queries']
        assert counters['random_queries'] == record['arguments']['random_queries']
    assert record['strict_pressure_candidate_found'] == any(w['critical_pressure'] < -1e-9
        for gauge in record['gauges'].values() for w in gauge['retained'])
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(), 'status': 'passed',
            'source_provenance': provenance, 'parameter_dimension_over_F4': len(indices),
            'retained_literal_witnesses': witnesses, 'literal_endpoint_entries': endpoints,
            'complete_profiles_certified_by_polynomial_frames': profiles_checked, 'gauges': counts,
            'recorded_query_counters': counters, 'completed_requested_schedule': record['completed_requested_schedule'],
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
        args.output.write_text(json.dumps({'scope': 'Independent literal sparse constructions and complete polynomial Smith certificates; numerical roots remain exploratory.',
                                         'records': results, 'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2)+'\n')


if __name__ == '__main__':
    main()
