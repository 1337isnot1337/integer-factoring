#!/usr/bin/env python3
"""Archive interrupted/resumed fused searches and replay their exact orders."""
from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from math import log
from pathlib import Path

from scipy.special import logsumexp

from fused_physical_basis_oracle import FusedPhysicalBasisOracle


def rational(value):
    return Fraction(str(value))


def archive_bytes(path):
    raw = path.read_bytes()
    return gzip.decompress(raw) if path.suffix == '.gz' else raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, default=Path('evidence/fused-physical-basis-prefix.json.gz'))
    parser.add_argument('--continuation', type=Path, default=Path('evidence/fused-physical-basis-continuation.json.gz'))
    parser.add_argument('--baseline', type=Path, default=Path('evidence/fused-physical-basis-baseline.json.gz'))
    parser.add_argument('--summary', type=Path, default=Path('evidence/fused-physical-basis-search-summary.json'))
    parser.add_argument('--replay', type=Path, default=Path('evidence/fused-physical-basis-search-replay.json'))
    args = parser.parse_args()
    raw = [archive_bytes(path) for path in (args.prefix, args.continuation)]
    reports = list(map(json.loads, raw))
    old, new = reports
    continuation = new['continuation']
    assert continuation['prefix_sha256'] == sha256(raw[0]).hexdigest()
    assert continuation['previous_checkpoint_seconds'] == old['elapsed_seconds']
    assert continuation['previous_saved_attempts'] == len(old['records'])
    assert new['budget_completed']
    assert old['source_sha256'] == new['source_sha256']
    assert old['baseline_sha256'] == new['baseline_sha256']
    records = old['records'] + new['records']
    assert [row['trial'] for row in records] == list(range(len(records)))
    assert continuation['first_continuation_trial'] == len(old['records'])
    assert abs(new['combined_checkpoint_seconds'] - sum(row['elapsed_seconds'] for row in reports)) < 1e-7
    assert new['combined_checkpoint_seconds'] >= old['arguments']['seconds']
    directory = Path(__file__).parent
    for name, expected in old['source_sha256'].items():
        assert sha256((directory/name).read_bytes()).hexdigest() == expected
    baseline = archive_bytes(args.baseline)
    assert sha256(baseline).hexdigest() == old['baseline_sha256']
    archives = []
    archive_names = ('fused-physical-basis-prefix.json.gz', 'fused-physical-basis-continuation.json.gz',
                     'fused-physical-basis-baseline.json.gz')
    for name, content in zip(archive_names, raw+[baseline]):
        path = args.summary.with_name(name)
        compressed = gzip.compress(content, mtime=0)
        path.write_bytes(compressed)
        archives.append({'path': str(path), 'sha256': sha256(compressed).hexdigest(),
                         'uncompressed_sha256': sha256(content).hexdigest(),
                         'uncompressed_bytes': len(content), 'compressed_bytes': len(compressed)})
    winners = {}
    for report in reports:
        for q, candidate in report['winners'].items():
            if q not in winners or candidate['best']['rho'] < winners[q]['best']['rho']:
                winners[q] = candidate
    trials = [[row['trial'], row['alphabet'], row['basis_name'], row['baseline_partition_index'],
               row['completed_trial'], row['queries'], row['elapsed_seconds'],
               row['best']['rho'] if row['best'] else None, row['interruption']] for row in records]
    summary = {'scope': 'Complete recorded two-hour fused physical-basis search, with interrupted prefix and separately timed continuation. No global optimum or factoring gain certified.',
               'requested_total_seconds': old['arguments']['seconds'],
               'prefix_seconds': old['elapsed_seconds'], 'continuation_seconds': new['elapsed_seconds'],
               'combined_checkpoint_seconds': new['combined_checkpoint_seconds'],
               'attempted_trials': len(records), 'completed_trials': sum(row['completed_trial'] for row in records),
               'profile_queries': sum(row['queries'] for row in records),
               'continuation': continuation, 'source_sha256': old['source_sha256'],
               'baseline_partitions': old['baseline_partitions'], 'archives': archives,
               'trial_fields': ['trial', 'alphabet', 'basis_name', 'partition_index', 'completed',
                                'profile_queries', 'seconds', 'best_rho', 'interruption'],
               'trials': trials, 'winners': winners,
               'checkpoint_source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    args.summary.write_text(json.dumps(summary, indent=2)+'\n')
    controls = []
    for q, winner in sorted(winners.items()):
        oracle = FusedPhysicalBasisOracle(winner['m'], winner['theta'], winner['first_basis'], winner['second_basis'])
        best = winner['best']
        types = []
        for row in best['paired_types']:
            values = list(map(rational, row))
            total = sum(values)
            types.append([value/total for value in values])
        entropy = sum(-sum(float(value)*log(float(value)) for value in row) for row in types)
        costs, profiles, point_sets = [], [], []
        outer_scale = rational(best['outer_log_scale'])
        for cut, places in enumerate(best['place_records']):
            current, points = [], set()
            allocations = [rational(place['group_allocation']) for place in places]
            allocation_sum = sum(allocations)
            assert all(value > 0 for value in allocations)
            assert abs(float(allocation_sum)-1) < 2e-12
            allocations = [value/allocation_sum for value in allocations]
            for place, allocation in zip(places, allocations):
                point = None if place['point'] is None else tuple(place['point'])
                assert point != (1, 1)
                points.add(point)
                heights = tuple(map(rational, place['heights']))
                multiplicity = len(place['conjugate_labels'])
                assert multiplicity > 0
                saved_per_root = rational(place['allocation_per_root'])
                assert abs(float(saved_per_root*multiplicity-rational(place['group_allocation']))) < 2e-12
                assert abs(float(saved_per_root*rational(place['local_log_scale'])-outer_scale)) < 2e-12
                per_root = allocation/multiplicity
                # Exact rational normalization restores the shared-scale
                # identities, instead of trusting independently rounded
                # floating allocations and local scales.
                scale = outer_scale/per_root
                orders = oracle.profile(cut, point, heights, scale)
                error = max(abs(float(actual)-saved) for actual, saved in zip(orders, place['restored_orders']))
                assert error < 2e-8
                shift = oracle.type_shift(heights, types[0], types[1])
                cost = float(logsumexp([-float(value) for value in orders])) + float(shift)
                assert abs(cost-place['complete_cost']) < 2e-8
                current.append(cost)
                profiles.append({'cut': cut, 'point': place['point'], 'channels': len(orders),
                                 'exact_group_allocation': str(allocation),
                                 'exact_allocation_per_root': str(per_root),
                                 'exact_local_scale': str(scale),
                                 'exact_orders': [str(value) for value in orders],
                                 'exact_type_mean': str(shift), 'complete_cost': cost,
                                 'saved_order_error': error})
            costs.append(max(current))
            point_sets.append(points)
        assert point_sets[0].isdisjoint(point_sets[1])
        rho = (oracle.degree*float(outer_scale)+sum(costs))/(2*entropy)
        assert abs(rho-best['rho']) < 2e-8
        controls.append({'alphabet': int(q), 'degree': oracle.degree, 'trial': winner['trial'],
                         'first_basis': winner['first_basis'], 'second_basis': winner['second_basis'],
                         'endpoint_entries': int(q)**4, 'paired_types': [[str(v) for v in row] for row in types],
                         'recorded_rho': best['rho'], 'recomputed_rho': rho, 'profiles': profiles})
    result = {'scope': 'Exact rational complete weighted-order replay using bordered polynomial minors; logarithmic pressure recomputation is floating. No threshold certificate.',
              'passed': True, 'summary_sha256': sha256(args.summary.read_bytes()).hexdigest(),
              'checkpoint_source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'source_sha256': old['source_sha256'], 'controls': controls,
              'complete_profiles': sum(len(row['profiles']) for row in controls),
              'complete_channels': sum(sum(p['channels'] for p in row['profiles']) for row in controls),
              'endpoint_entries': sum(row['endpoint_entries'] for row in controls)}
    args.replay.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'attempted': summary['attempted_trials'], 'completed': summary['completed_trials'],
                      'seconds': summary['combined_checkpoint_seconds'], 'profiles': result['complete_profiles'],
                      'channels': result['complete_channels'], 'rhos': [row['recomputed_rho'] for row in controls]}))


if __name__ == '__main__':
    main()
