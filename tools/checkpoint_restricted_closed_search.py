#!/usr/bin/env python3
"""Freeze a completed closed-place search and replay every retained winner.

Complete weighted orders are recomputed by the independent exact minor
oracle at rationalized types, heights, and consistent Hermite allocations.
The final logarithmic costs remain floating diagnostic values.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import gzip
import json
import math
from pathlib import Path

from scipy.special import logsumexp

from f4_closed_place_oracle import F4ClosedPlaceOracle
from paired_mixed_place_pressure import PairedMixedPressure
from reflection_laurent_physical_frames import digest


def rational_vector(values, maximum_denominator):
    result = [Fraction(float(value)).limit_denominator(maximum_denominator) for value in values]
    total = sum(result)
    if min(result) <= 0 or total <= 0:
        raise ValueError('Use strictly positive normalized weights.')
    return [value / total for value in result]


def encoded(value):
    return [value.numerator, value.denominator]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    data = json.loads(raw)
    if not data['completed'] or data['stopping_reason'] != 'complete trial list':
        raise ValueError('Only a complete planned search can be archived by this checker.')
    for name, expected in data['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name
    for name, expected in data['input_sha256'].items():
        assert sha256(Path(data['arguments'][name]).read_bytes()).hexdigest() == expected, name
    selected = data['selected_cases']
    expected_keys = {(case, partition, restart) for case, row in enumerate(selected)
                     for partition in range(len(row['partitions']))
                     for restart in range(data['arguments']['starts'])}
    actual_keys = {(row['case_index'], row['partition_index'], row['restart']) for row in data['records']}
    assert expected_keys == actual_keys and len(actual_keys) == len(data['records'])
    assert data['attempted_trials'] == len(actual_keys)
    assert data['completed_trials'] == sum(row['completed'] for row in data['records'])
    assert data['unresolved_objective_calls'] == sum(row['unresolved_objective_calls'] for row in data['records'])
    assert all(row['completed'] for row in data['records'])
    archive = gzip.compress(raw, compresslevel=9, mtime=0)
    assert gzip.decompress(archive) == raw
    args.archive.write_bytes(archive)
    winners, replays = [], []
    totals = {'complete_profiles': 0, 'complete_channels': 0, 'endpoint_entries': 0,
              'allocation_identities': 0, 'maximum_saved_profile_error': 0.,
              'maximum_rational_profile_error': 0., 'maximum_saved_rho_error': 0.,
              'maximum_rational_rho_change': 0.}
    for case_index, case in enumerate(selected):
        records = [row for row in data['records'] if row['case_index'] == case_index and row['best'] is not None]
        minimal = min(records, key=lambda row: row['best']['rho'])
        winner = data['winners'][str(case_index)]
        assert {k: v for k, v in minimal.items() if k != 'best'} == {k: v for k, v in winner.items() if k != 'best'}
        assert minimal['best'] == {k: v for k, v in winner['best'].items() if k != 'place_records'}
        seed, best = case['seed'], winner['best']
        assert digest(seed['values']) == case['tensor_sha256'] == winner['tensor_sha256']
        oracle = F4ClosedPlaceOracle(seed)
        validated = PairedMixedPressure(oracle, winner['partition'])
        totals['endpoint_entries'] += len(seed['values'])
        px, py = [rational_vector(p, 10**12) for p in best['paired_types']]
        outer = Fraction(float(best['outer_log_scale'])).limit_denominator(10**8)
        assert outer > 0
        saved_costs, rational_costs, local_records = [], [], []
        for cut, side in enumerate(best['place_records']):
            assert len(side) == len(validated.groups[cut])
            allocation = rational_vector([row['group_allocation'] for row in side], 10**12)
            assert abs(sum(row['group_allocation'] for row in side) - 1) < 2e-10
            saved_side, rational_side, rational_records = [], [], []
            for group, alpha, expected_group in zip(side, allocation, validated.groups[cut]):
                point = None if group['point'] is None else tuple(group['point'])
                assert point == expected_group['point']
                assert group['conjugate_labels'] == expected_group['labels']
                multiplicity = len(group['conjugate_labels'])
                assert multiplicity > 0
                assert abs(group['allocation_per_root'] * multiplicity - group['group_allocation']) < 2e-10
                assert abs(group['allocation_per_root'] * group['local_log_scale'] - best['outer_log_scale']) < 2e-10
                totals['allocation_identities'] += 2
                saved_log, _, _, saved_orders = oracle.log_profile_with_gradient(
                    cut, point, group['heights'], group['local_log_scale'])
                assert len(saved_orders) == len(group['restored_orders']) == oracle.n
                saved_error = max(abs(float(a) - float(b)) for a, b in zip(saved_orders, group['restored_orders']))
                assert saved_error < 2e-8
                saved_mean = oracle.type_shift(group['heights'], *best['paired_types'])
                saved_side.append(saved_log + saved_mean)
                heights = [Fraction(float(value)).limit_denominator(10**8) for value in group['heights']]
                local_scale = outer * multiplicity / alpha
                orders = oracle.profile(cut, point, heights, local_scale)
                native_log, _, _, native_orders = oracle.log_profile_with_gradient(
                    cut, point, list(map(float, heights)), float(local_scale))
                assert len(orders) == len(native_orders) == oracle.n
                rational_error = max(abs(float(a) - float(b)) for a, b in zip(orders, native_orders))
                assert rational_error < 2e-7
                mean = oracle.type_shift(heights, px, py)
                logarithm = float(logsumexp([-float(value) for value in orders]))
                assert abs(logarithm - native_log) < 2e-7
                rational_side.append(logarithm + float(mean))
                rational_records.append({'point': point, 'conjugate_labels': group['conjugate_labels'],
                                         'group_allocation': encoded(alpha),
                                         'allocation_per_root': encoded(alpha / multiplicity),
                                         'local_log_scale': encoded(local_scale),
                                         'heights': list(map(encoded, heights)),
                                         'exact_orders': list(map(encoded, orders)),
                                         'exact_type_shift': encoded(mean),
                                         'floating_complete_cost': rational_side[-1]})
                totals['complete_profiles'] += 1
                totals['complete_channels'] += len(orders)
                totals['maximum_saved_profile_error'] = max(totals['maximum_saved_profile_error'], saved_error)
                totals['maximum_rational_profile_error'] = max(totals['maximum_rational_profile_error'], rational_error)
            saved_costs.append(max(saved_side))
            rational_costs.append(max(rational_side))
            local_records.append(rational_records)
        entropy = -sum(float(p) * math.log(float(p)) for p in px + py)
        saved_entropy = -sum(p * math.log(p) for side in best['paired_types'] for p in side)
        saved_rho = (seed['degree'] * best['outer_log_scale'] + sum(saved_costs)) / (2 * saved_entropy)
        rational_rho = (seed['degree'] * float(outer) + sum(rational_costs)) / (2 * entropy)
        assert abs(saved_rho - best['rho']) < 2e-8
        totals['maximum_saved_rho_error'] = max(totals['maximum_saved_rho_error'], abs(saved_rho - best['rho']))
        totals['maximum_rational_rho_change'] = max(totals['maximum_rational_rho_change'], abs(rational_rho - best['rho']))
        compact = {'case_index': case_index, 'source_index': case['source_index'],
                   'subset_x': case['subset_x'], 'subset_y': case['subset_y'],
                   'degree': seed['degree'], 'tensor_sha256': case['tensor_sha256'],
                   'partition_index': winner['partition_index'], 'restart': winner['restart'],
                   'rho': best['rho'], 'rational_replay_rho': rational_rho,
                   'old_two_end_weighted_rho': case['old_weighted_rho'],
                   'unweighted_adaptive_rho': case['adaptive_rho']}
        winners.append(compact)
        replays.append(compact | {'seed': seed, 'type_x': list(map(encoded, px)),
                                  'type_y': list(map(encoded, py)), 'outer_log_scale': encoded(outer),
                                  'local_profiles': local_records})
    common = {'scope': 'Completed bounded weighted closed-place search and independent exact local-profile replay; no global optimum or improved factoring exponent.',
              'raw_sha256': sha256(raw).hexdigest(), 'archive_sha256': sha256(archive).hexdigest(),
              'archive_bytes': len(archive), 'uncompressed_bytes': len(raw),
              'completed_trials': data['completed_trials'], 'attempted_trials': data['attempted_trials'],
              'selected_cases': len(selected), 'elapsed_seconds': data['elapsed_seconds'],
              'unresolved_objective_calls': data['unresolved_objective_calls'],
              'source_sha256': data['source_sha256'], 'input_sha256': data['input_sha256'],
              'replay_totals': totals}
    summary = common | {'winners': winners, 'best_selected_new_place_rho': min(row['rho'] for row in winners)}
    args.summary.write_text(json.dumps(summary, indent=2) + '\n')
    sources = (Path(__file__).name, 'f4_closed_place_oracle.py', 'paired_mixed_place_pressure.py',
               'f16_weighted_local_oracle.py', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    replay = common | {'status': 'pass', 'rational_winner_replays': replays,
                       'checker_source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                                for name in sources}}
    args.output.write_text(json.dumps(replay, indent=2) + '\n')
    print(json.dumps({'status': 'pass', 'best_selected_new_place_rho': summary['best_selected_new_place_rho'],
                      'completed_trials': data['completed_trials'], 'totals': totals}))


if __name__ == '__main__':
    main()
