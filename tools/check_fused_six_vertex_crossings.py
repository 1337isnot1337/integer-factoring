#!/usr/bin/env python3
"""Literal endpoint, binary Smith, and recorded-gauge replays for finite fusion."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from fused_six_vertex_crossings import construct, transform, expanded, oracles
from f16_identity_operator_profiles import MUL, evaluate_one
from check_physical_split_crossings import scalar_restriction
from check_block_smith_crossings import reverse_polynomial, smith_profile, toeplitz_rank
from search_fused_six_vertex_gauges import Pressure


def check_instance(m, theta):
    source = construct(m, theta)
    assert source['invariant'] and source['endpoint_identity']
    q = m+1; degree = source['degree']
    endpoints = profiles = ranks = 0
    rows = []
    for a, b in [(0, None), (2, 3), (9, 15), (14, 10), (8, 10), (None, 0)]:
        matrix = transform(source['matrix'], degree, a, b)
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                assert evaluate_one(value) == int(i == j)
                endpoints += 1
        native = oracles(matrix, q)
        both = []
        for side in range(2):
            cut = [[{} for _ in range(q*q)] for _ in range(q*q)]
            for x, y, u, v in product(range(q), repeat=4):
                row, col = (q*x+y, q*u+v) if side == 0 else (q*x+u, q*y+v)
                cut[row][col] = matrix[q*y+x][q*u+v]
            binary = scalar_restriction(expanded(cut))
            if side:
                binary = [[reverse_polynomial(value, degree) for value in row] for row in binary]
            orders, precision = smith_profile(binary, degree)
            assert all(orders[i:i+4] == (orders[i],)*4 for i in range(0, len(orders), 4))
            restored = tuple(orders[i]-(degree if side else 0) for i in range(0, len(orders), 4))
            raw = native[side].profile(infinity=bool(side))['orders']
            assert all(abs(raw[i]-raw[i+1]) < 1e-9 for i in range(0, len(raw), 2))
            assert restored == tuple(round(value) for value in raw[::2])
            for length in sorted({1, orders[len(orders)//2]+1, orders[-1]+1}):
                assert toeplitz_rank(binary, length) == sum(max(0, length-order) for order in orders)
                ranks += 1
            both.append(restored); profiles += 1
        rows.append({'projective_pair': [a, b], 'complete_profiles_zero_then_infinity': both})
    return {'m': m, 'theta': theta, 'degree': degree, 'cancelled_common_degree': source['common_degree'],
        'literal_endpoint_entries': endpoints, 'complete_independent_binary_profiles': profiles,
        'independent_toeplitz_ranks': ranks, 'projective_cases': rows}


def replay(path):
    data = json.loads(path.read_text())
    assert data['completed_requested_starts']
    for name, digest in data['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest, name
    source = construct(data['arguments']['spin'], data['arguments']['theta'])
    assert source['invariant'] and source['endpoint_identity']
    cached = {}
    for row in data['records']:
        key = row['zero_image'], row['infinity_image']
        if key not in cached:
            cached[key] = Pressure(transform(source['matrix'], source['degree'], *key), data['arguments']['spin']+1)
        pressure = cached[key]
        pressure.evaluate(np.asarray(row['parameters'])[pressure.active])
        assert abs(pressure.last['rho']-row['rho']) < 1e-9
        assert abs(pressure.last['critical_pressure']-row['critical_pressure']) < 1e-8
        for actual, recorded in zip(pressure.last['profiles'], row['profiles']):
            assert np.max(np.abs(np.array(actual['orders'])-recorded['orders'])) < 1e-7
    return {'record': path.name, 'completed_starts': len(data['records']),
        'optimizer_successes': sum(row['optimizer_success'] for row in data['records']),
        'best_replayed_rho': data['best']['rho'],
        'every_recorded_full_profile_and_cost_replayed': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--search-records', type=Path, nargs='*', default=[])
    args = parser.parse_args()
    rejected = []
    for m, theta in [(3, 2), (4, 2)]:
        result = construct(m, theta)
        assert not (result['invariant'] and result['endpoint_identity'])
        rejected.append({'m': m, 'theta': theta, 'invariant': result['invariant'],
                         'endpoint_identity': result['endpoint_identity']})
    result = {'status': 'Exact finite fused crossings and independent full profiles; no asymptotic improvement.',
        'cases': [check_instance(2, 2), check_instance(3, 4)], 'rejected_fusion_cases': rejected,
        'search_replays': [replay(path) for path in args.search_records]}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'cases': [{key: value for key, value in row.items() if key != 'projective_cases'}
                                for row in result['cases']],
                      'rejected': rejected, 'search_replays': result['search_replays']}))


if __name__ == '__main__':
    main()
