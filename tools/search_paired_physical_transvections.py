#!/usr/bin/env python3
"""Finite q5 paired GL5(F4) basis-change screen at two fixed gauges.

The schedule contains every one-transvection coefficient Frobenius orbit,
then every pair consisting of one X/V and one Y/U transvection, again
modulo coefficient Frobenius. Profiles and tensors are exact; scale and
root optimization are floating-point search evidence only.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from paired_physical_transvections import transform
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles
from search_saturated_dual_directions import GAUGES, PAIR_ORBITS, check_endpoint, complete_profiles


def schedule():
    moves = [(target, source) for target in range(5) for source in range(5) if target != source]
    for pair in (0, 1):
        for target, source in moves:
            for coefficient in (1, 2):
                yield ((pair, target, source, coefficient),), 'single'
    for tx, sx in moves:
        for ty, sy in moves:
            for a, b in PAIR_ORBITS:
                yield ((0, tx, sx, a), (1, ty, sy, b)), 'pair'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=float, default=1800)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.seconds <= 0:
        parser.error('The time budget must be positive.')
    started = time.monotonic()
    base = five_letter_source()[0]
    initial_profiles = complete_profiles(base, GAUGES)
    cache, gauges = {}, {}
    for name, energy in GAUGES.items():
        cost = root_for_profiles(initial_profiles[name], energy, energy)
        cache[name] = {initial_profiles[name]: cost}
        gauges[name] = {'energies': energy, 'reference_profiles': initial_profiles[name],
                        'reference_cost': cost, 'best_recorded_rho': cost['rho'], 'retained': []}
    names = (Path(__file__).name, 'paired_physical_transvections.py', 'search_saturated_dual_directions.py',
             'search_f4_recentered_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'weighted_smith_oracle_f4.py')
    report = {'scope': 'All paired elementary physical-basis changes and all one-X/one-Y pairs, over F4 up to coefficient Frobenius; fixed gauges and numerical scale optimization.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'arguments': vars(args) | {'output': str(args.output)}, 'reference_values': base,
              'reference_digest': digest(base), 'alphabet': 5, 'degree': 18,
              'planned_single_queries': 80, 'planned_pair_queries': 2000,
              'single_coefficient_orbits': (1, 2), 'pair_coefficient_orbits': PAIR_ORBITS}
    counters = {'single_queries': 0, 'pair_queries': 0, 'literal_endpoint_entries': 625}
    results = []

    def save(completed=False):
        for name in gauges:
            gauges[name]['distinct_profile_classes'] = len(cache[name])
        report.update(counters=counters, gauges=gauges, query_records=results,
                      completed_requested_schedule=completed, elapsed_seconds=time.monotonic()-started,
                      strict_pressure_candidate_found=any(row['critical_pressure'] < -1e-9
                          for gauge in gauges.values() for row in gauge['retained']))
        args.output.write_text(json.dumps(report, indent=2)+'\n')

    for operations, stage in schedule():
        if time.monotonic()-started >= args.seconds:
            break
        values = transform(base, 5, operations)
        check_endpoint(values, 18)
        counters['literal_endpoint_entries'] += 625
        profiles = complete_profiles(values, GAUGES)
        recorded = {'operations': operations, 'stage': stage, 'digest': digest(values), 'gauges': {}}
        for name, energy in GAUGES.items():
            profile = profiles[name]
            if profile not in cache[name]:
                cache[name][profile] = root_for_profiles(profile, energy, energy)
            cost, gauge = cache[name][profile], gauges[name]
            recorded['gauges'][name] = {'profiles': profile, **cost}
            if cost['rho'] < gauge['best_recorded_rho']-1e-11:
                gauge['best_recorded_rho'] = cost['rho']
                print(json.dumps({'event': 'new_best', 'gauge': name, 'operations': operations, **cost}), flush=True)
            existing = {tuple(map(tuple, row['profiles'])) for row in gauge['retained']}
            if profile not in existing and (len(gauge['retained']) < 24 or cost['rho'] < max(row['rho'] for row in gauge['retained'])):
                gauge['retained'].append({'operations': operations, 'stage': stage, 'values': values,
                                          'digest': digest(values), 'profiles': profile, **cost})
                gauge['retained'][:] = sorted(gauge['retained'], key=lambda row: row['rho'])[:24]
        results.append(recorded)
        counters[stage+'_queries'] += 1
        if len(results) % 100 == 0:
            save()
            print(json.dumps({'event': 'checkpoint', **counters,
                              'best': {name: row['best_recorded_rho'] for name, row in gauges.items()},
                              'profile_classes': {name: len(row) for name, row in cache.items()},
                              'seconds': time.monotonic()-started}), flush=True)
    completed = counters['single_queries'] == 80 and counters['pair_queries'] == 2000
    save(completed)
    print(json.dumps({'event': 'finished', **counters, 'completed': completed,
                      'best': {name: row['best_recorded_rho'] for name, row in gauges.items()},
                      'seconds': time.monotonic()-started}), flush=True)


if __name__ == '__main__':
    main()
