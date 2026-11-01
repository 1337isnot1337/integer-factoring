#!/usr/bin/env python3
"""Literal Gaussian-minor and independent complete-matrix gauge replay.

Tensor coefficients are checked exactly. Gauge comparisons and Gibbs
roots remain floating point; this is not an exact optimum or pressure
sign certificate.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_aligned_gaussian_complements import expand_square, minor, multiply, strict_matrix
from native_matrix_smith_f4 import MatrixOracle
from replay_source_provenance import resolve_source_hashes
from search_f4_recentered_crossings import digest


ROOT = Path(__file__).resolve().parent.parent


def literal_seed(seed):
    xcode, ycode = seed['code_x'], seed['code_y']
    assert len(xcode) == len(ycode) == len(set(xcode)) == len(set(ycode)) == 5
    assert all(0 <= mask < 256 and mask.bit_count() == 4 for mask in xcode+ycode)
    common = seed['entry_gcd_in_z']
    assert common > 0 and common.bit_count() % 2 == 1
    matrix = strict_matrix()
    for index, (polynomial, high) in enumerate(seed['values']):
        x, y, u, v = ((index//5**leg) % 5 for leg in range(4))
        assert high == 0 and polynomial.bit_count() % 2 == int(x == v and y == u)
        actual = minor(matrix, xcode[v] | (ycode[y] << 8), ycode[u] | (xcode[x] << 8))
        assert actual == expand_square(multiply(common, polynomial))
    return 625


def matrix_oracles(seed):
    first = [[(0, 0)]*25 for _ in range(25)]
    second = [[(0, 0)]*25 for _ in range(25)]
    for index, value in enumerate(seed['values']):
        x, y, u, v = ((index//5**leg) % 5 for leg in range(4))
        first[5*x+u][5*y+v] = tuple(value)
        second[5*x+y][5*u+v] = tuple(value)
    return MatrixOracle(first, 25), MatrixOracle(second, 25)


def gauge_replay(oracles, parameters):
    import numpy as np
    from scipy.optimize import brentq
    from scipy.special import logsumexp
    assert len(parameters) == 33 and parameters[-1] > 0
    rows = [(0.,)+tuple(parameters[4*i:4*i+4]) for i in range(8)]
    first = oracles[0].profile(row_weights=[rows[0][x]+rows[2][u] for x in range(5) for u in range(5)],
        column_weights=[rows[1][y]+rows[3][v] for y in range(5) for v in range(5)],
        scale=parameters[-1], infinity=True)['orders']
    second = oracles[1].profile(row_weights=[rows[4][x]+rows[5][y] for x in range(5) for y in range(5)],
        column_weights=[rows[6][u]+rows[7][v] for u in range(5) for v in range(5)],
        scale=parameters[-1])['orders']
    gx = np.asarray([sum(rows[leg][i] for leg in (0, 3, 4, 7)) for i in range(5)])
    gy = np.asarray([sum(rows[leg][i] for leg in (1, 2, 5, 6)) for i in range(5)])
    constant = logsumexp(-np.asarray(first))+logsumexp(-np.asarray(second))

    def pressure(rate):
        return constant-rate*(logsumexp(-gx/rate)+logsumexp(-gy/rate))

    hi = 8.
    while pressure(hi) > 0:
        hi *= 2
    rate = brentq(pressure, 1e-9, hi, xtol=2e-13, rtol=2e-13)
    px = np.exp(-gx/rate-logsumexp(-gx/rate)).tolist()
    py = np.exp(-gy/rate-logsumexp(-gy/rate)).tolist()
    return {'rho': rate/2, 'type_x': px, 'type_y': py, 'orders': first+second,
            'energy_x_times_log_s': gx.tolist(), 'energy_y_times_log_s': gy.tolist()}


def close(actual, expected, tolerance=2e-7):
    if isinstance(actual, (list, tuple)):
        assert len(actual) == len(expected)
        for first, second in zip(actual, expected):
            close(first, second, tolerance)
    else:
        assert abs(actual-expected) <= tolerance*max(1., abs(actual), abs(expected)), (actual, expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    parser.add_argument('--refinement', type=Path, action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    record = json.loads(args.record.read_text())
    assert record['completed_requested_schedule']
    provenance, unresolved = resolve_source_hashes(record['source_sha256'])
    assert not unresolved, unresolved
    for path, expected in record['input_sha256'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == expected
    minors = profiles = states = 0
    rows = []
    for candidate in record['candidates']:
        seed = candidate['literal_seed']
        assert digest(seed['values']) == candidate['source_row']['digest']
        minors += literal_seed(seed)
        oracles = matrix_oracles(seed)
        costs = []
        for run in candidate['runs']:
            if run['status'] != 'successful_profile_retained':
                continue
            actual = gauge_replay(oracles, run['params'])
            close(actual['rho'], run['rho'])
            close(actual['orders'], run['orders_times_log_s'])
            for key in ('type_x', 'type_y', 'energy_x_times_log_s', 'energy_y_times_log_s'):
                close(actual[key], run[key])
            costs.append(actual['rho'])
            profiles += 2
            states += 1
        close(min(costs), candidate['best']['rho'])
        rows.append({'family': candidate['family'], 'code_x': seed['code_x'], 'code_y': seed['code_y'],
                     'rho': candidate['best']['rho'], 'type_x': candidate['best']['type_x'],
                     'type_y': candidate['best']['type_y'],
                     'optimizer_successes': sum(run['optimizer_success'] for run in candidate['runs']),
                     'recorded_failed_oracle_evaluations': sum(sum(run.get('errors', {}).values()) for run in candidate['runs'])})
    refinements = []
    for path in args.refinement:
        item = json.loads(path.read_text())
        assert item['status'] == 'completed' and item['returncode'] == 0
        resolved, unresolved = resolve_source_hashes(item['source_sha256'])
        assert not unresolved, unresolved
        assert sha256((ROOT/item['source_record']).read_bytes()).hexdigest() == item['source_record_sha256']
        for key in ('seed', 'start', 'output'):
            assert sha256((ROOT/item[key]).read_bytes()).hexdigest() == item[key+'_sha256']
        seed = json.loads((ROOT/item['seed']).read_text())
        assert any(seed == row['literal_seed'] for row in record['candidates'])
        data = json.loads((ROOT/item['output']).read_text())
        best = data['records'][0]
        actual = gauge_replay(matrix_oracles(seed), best['parameters'])
        close(actual['rho'], best['rho_floating'])
        close(actual['type_x'], best['type_x_v'])
        close(actual['type_y'], best['type_y_u'])
        profiles += 2
        states += 1
        refinements.append({'provenance': str(path), 'provenance_sha256': sha256(path.read_bytes()).hexdigest(),
                            'family': item['family'], 'rho': actual['rho'],
                            'type_x': actual['type_x'], 'type_y': actual['type_y'],
                            'profile_queries': data['profile_queries_attempted'],
                            'failed_profile_queries': data['failed_profiles'],
                            'completed_restarts': len(data['restarts']), 'source_provenance': resolved})
    report = {'scope': 'Exact complete Gaussian-minor reconstruction and independent complete-matrix replay of saved floating gauges; no exact numerical optimum or pressure-sign certificate.',
              'record': str(args.record), 'record_sha256': sha256(args.record.read_bytes()).hexdigest(),
              'status': 'passed', 'source_provenance': provenance, 'literal_gaussian_minors': minors,
              'saved_gauge_states_replayed': states, 'complete_matrix_gauge_profiles_replayed': profiles,
              'candidates': rows, 'refinements': refinements,
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'elapsed_seconds': time.monotonic()-started}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
