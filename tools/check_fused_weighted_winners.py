#!/usr/bin/env python3
"""Exact rational local-order replay of the retained full-cost search winners."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import log
from pathlib import Path

from scipy.special import logsumexp

from fused_weighted_place_profiles import FusedPlaceOracle


def rational(value):
    return Fraction(str(value))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('summary', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.summary.read_bytes()
    summary = json.loads(raw)
    records = []
    for key, winner in sorted(summary['winners'].items()):
        oracle = FusedPlaceOracle(winner['m'], winner['theta'])
        best = winner['best']
        types = []
        for row in best['paired_types']:
            values = list(map(rational, row))
            total = sum(values)
            types.append([value / total for value in values])
        entropy = sum(-sum(float(value) * log(float(value)) for value in row) for row in types)
        costs, local_records = [], []
        point_sets = []
        for cut, places in enumerate(best['place_records']):
            current, points = [], set()
            for place in places:
                point = None if place['point'] is None else tuple(place['point'])
                assert point != (1, 1)
                points.add(point)
                heights = tuple(map(rational, place['heights']))
                scale = rational(place['local_log_scale'])
                orders = oracle.profile(cut, point, heights, scale)
                error = max(abs(float(actual) - saved) for actual, saved in zip(orders, place['restored_orders']))
                assert error < 2e-8
                shift = oracle.type_shift(heights, types[0], types[1])
                log_sum = float(logsumexp([-float(value) for value in orders]))
                cost = log_sum + float(shift)
                assert abs(cost - place['complete_cost']) < 2e-8
                current.append(cost)
                local_records.append({'cut': cut, 'point': place['point'],
                                      'complete_channel_count': len(orders),
                                      'exact_orders': [str(value) for value in orders],
                                      'exact_mean_shift': str(shift),
                                      'largest_saved_order_error': error,
                                      'recomputed_log_profile_sum': log_sum,
                                      'recomputed_complete_cost': cost})
            costs.append(max(current))
            point_sets.append(points)
        assert point_sets[0].isdisjoint(point_sets[1])
        rho = (oracle.degree * float(rational(best['outer_log_scale'])) + sum(costs)) / (2 * entropy)
        assert abs(rho - best['rho']) < 2e-8
        records.append({'m': winner['m'], 'theta': winner['theta'], 'alphabet': int(key),
                        'trial': winner['trial'], 'completed_trial': winner['completed_trial'],
                        'exact_normalized_paired_types': [[str(value) for value in row] for row in types],
                        'recorded_rho': best['rho'], 'recomputed_rho': rho,
                        'local_profiles': local_records})
    directory = Path(__file__).parent
    names = (Path(__file__).name, 'fused_weighted_place_profiles.py', 'check_fused_closed_profiles.py',
             'check_fused_multipoint_profiles.py', 'fused_six_vertex_crossings.py',
             'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    result = {'scope': 'Exact rational weighted-order replay with floating logarithmic cost recomputation; no global optimum or improved factoring exponent certified.',
              'summary_sha256': sha256(raw).hexdigest(),
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'cases': records,
              'total_complete_profiles': sum(len(row['local_profiles']) for row in records),
              'total_channels': sum(sum(p['complete_channel_count'] for p in row['local_profiles']) for row in records)}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'complete_profiles': result['total_complete_profiles'], 'channels': result['total_channels'],
                      'rhos': [row['recomputed_rho'] for row in records]}))


if __name__ == '__main__':
    main()
