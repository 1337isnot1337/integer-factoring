#!/usr/bin/env python3
"""Archive and independently replay the completed stabilized-code searches."""
from __future__ import annotations
import argparse
from hashlib import sha256
from itertools import combinations
import gzip
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
from affine_place_profile_oracle import primitive_seed
from check_implicit_reflection_codes import matrix
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import physical_bit_energies
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_physical_code_gauges import start_params


def profile(seed, cut, infinity, heights, scale):
    q = seed['q_x']
    row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
    rows = [heights[row_legs[0]][a] + heights[row_legs[1]][b] for a in range(q) for b in range(q)]
    columns = [heights[column_legs[0]][a] + heights[column_legs[1]][b] for a in range(q) for b in range(q)]
    return MatrixOracle(matrix(seed, cut), q * q).profile(rows, columns, scale, infinity)['orders']


def natural_pressure(profiles, energies, scale):
    numerator = sum(float(logsumexp(-scale * np.asarray(side))) for side in profiles)
    equation = lambda root: numerator - 2 * root * float(logsumexp(-scale * np.asarray(energies) / root))
    return brentq(equation, 1e-10, 4, xtol=1e-13) / 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--natural', type=Path, required=True)
    parser.add_argument('--gauge', type=Path, required=True)
    parser.add_argument('--gauge-input', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    args = parser.parse_args()
    inputs = {name: getattr(args, name.replace('-', '_')) for name in ('natural', 'gauge', 'gauge-input')}
    raw = {name: path.read_bytes() for name, path in inputs.items()}
    data = {name: json.loads(value) for name, value in raw.items()}
    natural, gauge = data['natural'], data['gauge']
    assert natural['completed'] and natural['case_list_exhausted'] and gauge['completed']
    assert natural['completed_cases'] == len(natural['records'])
    assert gauge['input_sha256'] == sha256(raw['gauge-input']).hexdigest()
    for report in (natural, gauge):
        for name, expected in report['source_sha256'].items():
            assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected, name
    weights = physical_bit_energies(natural['arguments']['interval'])
    energy = lambda symbol: sum(value for bit, value in enumerate(weights) if symbol >> bit & 1)
    all_labels = sorted(range(1 << len(weights)), key=lambda symbol: (energy(symbol), symbol))
    expected = []
    for size in natural['arguments']['sizes']:
        level = energy(all_labels[size - 1])
        lower = [symbol for symbol in all_labels if energy(symbol) < level]
        boundary = [symbol for symbol in all_labels if energy(symbol) == level]
        expected.extend(lower + list(chosen) for chosen in combinations(boundary, size - len(lower)))
    assert [row['letters'] for row in natural['records']] == expected
    assert len(expected) == natural['expected_cases']
    totals = {'reconstructed_coefficients': 0, 'complete_integer_profiles': 0,
              'complete_integer_channels': 0, 'natural_pressure_queries': 0,
              'full_gauge_profiles': 0, 'full_gauge_channels': 0,
              'maximum_natural_pressure_error': 0., 'maximum_full_gauge_order_error': 0.,
              'maximum_full_gauge_support_error': 0.}
    natural_records, seeds = [], {}
    for record in natural['records']:
        raw_seed = restricted_tensor(record['interval'], record['letters'])
        seed, common = primitive_seed(raw_seed)
        assert digest(seed['values']) == record['tensor_sha256']
        assert {str(k): v for k, v in common.items()} == record['scalar_content']
        assert seed['degree'] == record['primitive_degree']
        assert record['energies'] == [energy(symbol) for symbol in record['letters']]
        q, energies = record['size'], record['energies']
        profiles = []
        for cut, infinity in ((1, True), (0, False)):
            heights = [[0] * q for _ in range(4)]
            for leg in ((2, 3) if infinity else (0, 1)):
                heights[leg] = [value // 2 for value in energies]
            orders = profile(seed, cut, infinity, heights, 1)
            assert all(value == round(value) for value in orders)
            profiles.append([int(value) for value in sorted(orders)])
        assert record['best']['rho'] == min(row['rho'] for row in record['pressure_queries'])
        for query in record['pressure_queries']:
            observed = natural_pressure(profiles, energies, query['scale'])
            error = abs(observed - query['rho'])
            assert error < 3e-10
            totals['maximum_natural_pressure_error'] = max(totals['maximum_natural_pressure_error'], error)
            totals['natural_pressure_queries'] += 1
        assert record['best']['params'] == start_params(np.asarray(energies) * record['best']['scale'], record['best']['scale']).tolist()
        totals['reconstructed_coefficients'] += len(seed['values'])
        totals['complete_integer_profiles'] += 2
        totals['complete_integer_channels'] += 2 * q * q
        seeds[record['tensor_sha256']] = seed
        natural_records.append({key: record[key] for key in ('size', 'letters', 'energies', 'primitive_degree', 'tensor_sha256', 'completed', 'unresolved', 'elapsed_seconds')}
                               | {'rho': record['best']['rho'], 'scale': record['best']['scale'],
                                  'integer_profiles_infinity_then_zero': profiles})
    cases = {case['label']: case for case in data['gauge-input']['cases']}
    expected_keys = {(label, angle) for label in cases for angle in (.25, math.pi / 4, 1.32)}
    actual_keys = {(record['label'], record['initial_angle']) for record in gauge['records']}
    assert expected_keys == actual_keys and len(actual_keys) == len(gauge['records'])
    gauge_records = []
    for label, case in cases.items():
        assert digest(case['seed']['values']) == case['tensor_sha256']
        seed = seeds[case['tensor_sha256']]
        assert case['seed']['values'] == [list(value) for value in seed['values']]
        winner = min((row for row in gauge['records'] if row['label'] == label),
                     key=lambda row: row['best']['support_value'])
        assert winner == gauge['winners'][label]
        best, q = winner['best'], seed['q_x']
        parameters = best['parameters'][:-1]
        angle, scale = best['parameters'][-1], parameters[-1]
        branch = (math.cos(angle), math.sin(angle))
        types = [best['details'][name] for name in ('type_x', 'type_y')]
        assert all(min(p) > 0 and abs(sum(p) - 1) < 1e-12 for p in types)
        entropy = -sum(value * math.log(value) for p in types for value in p)
        cost, saved_profiles = 0., []
        for side, (cut, infinity) in enumerate(((1, True), (0, False))):
            heights = [[0.] + parameters[(4 * side + leg) * (q - 1):(4 * side + leg + 1) * (q - 1)]
                       for leg in range(4)]
            orders = profile(seed, cut, infinity, heights, scale)
            expected_orders = best['details']['scaled_orders'][side]
            assert len(orders) == len(expected_orders) == q * q
            error = max(abs(a - b) for a, b in zip(sorted(orders), sorted(expected_orders)))
            assert error < 3e-8
            totals['maximum_full_gauge_order_error'] = max(totals['maximum_full_gauge_order_error'], error)
            cost += branch[side] * float(logsumexp(-np.asarray(orders) / branch[side]))
            cost += sum(sum(p * h for p, h in zip(types[0 if leg in (0, 3) else 1], heights[leg]))
                        for leg in range(4))
            saved_profiles.append(orders)
            totals['full_gauge_profiles'] += 1
            totals['full_gauge_channels'] += len(orders)
        observed = cost / entropy
        error = abs(observed - best['support_value'])
        assert error < 3e-8
        totals['maximum_full_gauge_support_error'] = max(totals['maximum_full_gauge_support_error'], error)
        gauge_records.append({'label': label, 'tensor_sha256': case['tensor_sha256'],
                              'support_value': best['support_value'], 'recomputed_support_value': observed,
                              'paired_types': types, 'parameters': best['parameters'],
                              'complete_matrix_profiles': saved_profiles})
    archives = {}
    for name, value in raw.items():
        archive = gzip.compress(value, compresslevel=9, mtime=0)
        assert gzip.decompress(archive) == value
        path = Path(str(args.prefix) + '-' + name + '.json.gz')
        path.write_bytes(archive)
        archives[name] = {'path': str(path), 'raw_sha256': sha256(value).hexdigest(),
                          'archive_sha256': sha256(archive).hexdigest(),
                          'raw_bytes': len(value), 'archive_bytes': len(archive)}
    names = (Path(__file__).name, 'restrict_implicit_reflection.py', 'check_implicit_reflection_codes.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'affine_place_profile_oracle.py')
    common = {'scope': 'Completed bounded exact-code searches; natural integer profiles are replayed completely, full-gauge weights and final costs remain numerical.',
              'archives': archives,
              'checker_source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                        for name in names},
              'natural_timing': {key: natural[key] for key in ('started_utc', 'finished_utc', 'elapsed_seconds')},
              'gauge_timing': {key: gauge[key] for key in ('started_utc', 'finished_utc', 'elapsed_seconds')},
              'natural_cases': len(natural['records']),
              'natural_unresolved_queries': sum(row['unresolved'] for row in natural['records']),
              'gauge_attempted_trials': len(gauge['records']),
              'gauge_completed_optimizations': sum(row['completed'] for row in gauge['records']),
              'gauge_unresolved_queries': sum(row['unresolved'] for row in gauge['records']),
              'counts': totals, 'passed': True}
    summary = dict(common, natural_results=[{key: value for key, value in row.items()
                                           if key != 'integer_profiles_infinity_then_zero'} for row in natural_records],
                   gauge_results=[{key: value for key, value in row.items()
                                   if key not in ('paired_types', 'parameters', 'complete_matrix_profiles')}
                                  for row in gauge_records])
    replay = dict(common, natural_profiles=natural_records, gauge_profiles=gauge_records)
    Path(str(args.prefix) + '-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    Path(str(args.prefix) + '-replay.json').write_text(json.dumps(replay, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': totals, 'archives': archives}))


if __name__ == '__main__':
    main()
