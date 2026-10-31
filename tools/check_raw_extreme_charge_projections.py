#!/usr/bin/env python3
"""Typed endpoints and original weighted cuts for raw extreme charge masks."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from native_matrix_smith_f4 import MatrixOracle
from relabelled_gaussian_sums import gaussian_tensor, charge_projection, flatten, natural_heights
from strict_core_charge_profiles import affine_charge_representatives
from strict_core_charge_gauges import CoreChargePressure


def endpoint_controls():
    interval, q = 2, 4
    base = gaussian_tensor(tuple(range(interval)), interval, False)
    words = list(product(range(q), repeat=2))
    totals = {word: tuple(sum((state >> i) & 1 for state in word) for i in range(interval))
              for word in words}
    pairs = [(x, v) for x in words for v in words if totals[x] == totals[v]]
    checked = 0
    for charges in affine_charge_representatives(2):
        values = charge_projection(base, charges)
        endpoint = [value[0].bit_count() & 1 for value in values]
        for x, v in pairs:
            for y, u in pairs:
                actual = 1
                for copy in range(2):
                    index = x[copy] | (y[copy] << interval) | (u[copy] << (2 * interval)) | (v[copy] << (3 * interval))
                    actual *= endpoint[index]
                assert actual == int(x == v and y == u)
                checked += 1
    return checked


def direct(interval, charges, offsets, bulk):
    q = 1 << interval
    indices = (0, interval - 1)
    core_label = lambda state: (state & 1) | (((state >> (interval - 1)) & 1) << 1)
    values = charge_projection(gaussian_tensor(tuple(range(interval)), interval, False),
                               [charges[core_label(state)] for state in range(q)])
    heights = natural_heights(tuple(range(interval)), interval).copy()
    for place in range(2):
        for leg in range(4):
            heights[place, leg] += [offsets[place, leg, core_label(state)] for state in range(q)]
    for state in range(q):
        count = state.bit_count()
        heights[0, 0, state] += bulk[0] * count
        heights[0, 1, state] -= bulk[0] * count
        heights[1, 1, state] += bulk[1] * count
        heights[1, 3, state] += bulk[1] * count
    orders = []
    for place in range(2):
        first, second = ((0, 1), (2, 3)) if place == 0 else ((0, 2), (1, 3))
        rows = [heights[place, first[0], a] + heights[place, first[1], b]
                for a, b in product(range(q), repeat=2)]
        columns = [heights[place, second[0], a] + heights[place, second[1], b]
                   for a, b in product(range(q), repeat=2)]
        orders.append(MatrixOracle(flatten(values, interval, bool(place)), q * q).profile(
            rows, columns, infinity=bool(place))['orders'])
    log_b = sum(float(logsumexp(-np.asarray(values))) for values in orders)
    gx = heights[:, (0, 3), :].sum(axis=(0, 1))
    gy = heights[:, (1, 2), :].sum(axis=(0, 1))
    equation = lambda rate: log_b - rate * float(logsumexp(-gx / rate) + logsumexp(-gy / rate))
    upper = 4.0
    while equation(upper) > 0:
        upper *= 2
    return orders, brentq(equation, 1e-9, upper, xtol=1e-13) / 2


def verify():
    endpoint_entries = endpoint_controls()
    rng = np.random.default_rng(753680)
    records = []
    for interval in (2, 4):
        for charges in ((0, 0, 0, 1), (0, 1, 1, 1), (0, 1, 2, 0)):
            model = CoreChargePressure(interval, charges, family='raw_extremes')
            assert model.core_indices == [0, interval - 1]
            for trial in range(2):
                offsets = rng.integers(-4, 5, (2, 4, 4)) / 4
                offsets[:, :, 0] = 0
                bulk = rng.integers(-4, 5, 2) / 4
                expected, expected_rho = direct(interval, charges, offsets, bulk)
                actual_rho, detail = model.evaluate(1.0, offsets, bulk)
                assert abs(actual_rho - expected_rho) < 1e-11
                for place in range(2):
                    actual = np.asarray(detail['active_orders'][place])
                    data = model.records[place][0]
                    for order, charge in zip(data['spectator_orders'], data['spectator_charges']):
                        actual = np.r_[actual, actual + order + bulk[place] * charge]
                    assert np.array_equal(np.sort(actual), np.sort(expected[place]))
                    records.append({'interval': interval, 'charges': list(charges),
                                    'trial': trial, 'place': place, 'complete_rank': len(actual)})
    root = Path(__file__).resolve().parent.parent
    files = ['tools/strict_core_charge_gauges.py', 'tools/check_raw_extreme_charge_projections.py',
             'tools/search_raw_extreme_charge_projections.py', 'research/raw-extreme-charge-projections.md']
    return {'status': 'pass', 'typed_endpoint_entries': endpoint_entries,
            'original_weighted_cut_profiles': records,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in files},
            'scope': 'Raw Gaussian typed endpoints and exact weighted cut controls; no improved exponent claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'typed_endpoint_entries': result['typed_endpoint_entries'],
                      'original_weighted_cut_profiles': len(result['original_weighted_cut_profiles'])}))


if __name__ == '__main__':
    main()
