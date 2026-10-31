#!/usr/bin/env python3
"""Original weighted-cut controls for regular physical core coordinates."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from relabelled_gaussian_sums import RelabelingPressure, gaussian_tensor, charge_projection
from strict_core_charge_gauges import CoreChargePressure


def verify():
    rng = np.random.default_rng(659108)
    records = []
    control_model = None
    for interval, core, charges in ((2, 2, (0, 1, 2, 0)),
                                    (4, 2, (0, 1, 1, 1)), (4, 1, (0, 1))):
        model = CoreChargePressure(interval, charges, core)
        q, mask = 1 << interval, (1 << core) - 1
        original = RelabelingPressure(charge_projection(
            gaussian_tensor(tuple(range(interval)), interval),
            [charges[state & mask] for state in range(q)]), tuple(range(interval)), interval)
        for trial in range(3):
            offsets = rng.integers(-8, 9, (2, 4, 1 << core)) / 4
            offsets[:, :, 0] = 0
            bulk = rng.integers(-4, 5, 2) / 4
            full = np.asarray([[[offsets[p, leg, state & mask] for state in range(q)]
                                for leg in range(4)] for p in range(2)])
            for state in range(q):
                count = state.bit_count()
                full[0, 0, state] += bulk[0] * count
                full[0, 1, state] -= bulk[0] * count
                full[1, 1, state] += bulk[1] * count
                full[1, 3, state] += bulk[1] * count
            actual, detail = model.evaluate(1.0, offsets, bulk)
            expected, original_detail = original.evaluate(1.0, full)
            assert abs(actual - expected) < 1e-11
            for place in range(2):
                complete = np.asarray(detail['active_orders'][place])
                data = model.records[place][0]
                for order, charge in zip(data['spectator_orders'], data['spectator_charges']):
                    complete = np.r_[complete, complete + order + bulk[place] * charge]
                expected_orders = np.asarray(original_detail['orders'][place])
                assert np.array_equal(np.sort(complete), np.sort(expected_orders))
                records.append({'interval': interval, 'core': core, 'charges': list(charges),
                                'trial': trial, 'place': place, 'complete_rank': len(complete),
                                'active_rank': model.records[place][1]['matrix_order']})
        if interval == 4 and core == 2:
            control_model = model
    model = control_model
    offsets = rng.normal(0, 0.3, (2, 4, 4))
    offsets[:, :, 0] = 0
    bulk, ell = np.array([0.17, -0.29]), 0.43
    _, (do, db, ds), _ = model.evaluate(ell, offsets, bulk, True)
    epsilon = 1e-5
    errors = []
    for place in range(2):
        for leg in range(4):
            for state in range(1, 4):
                shifted = offsets.copy()
                shifted[place, leg, state] += epsilon
                plus = model.evaluate(ell, shifted, bulk)[0]
                shifted[place, leg, state] -= 2 * epsilon
                minus = model.evaluate(ell, shifted, bulk)[0]
                errors.append(abs((plus - minus) / (2 * epsilon) - do[place, leg, state]))
    for index in range(2):
        shifted = bulk.copy()
        shifted[index] += epsilon
        plus = model.evaluate(ell, offsets, shifted)[0]
        shifted[index] -= 2 * epsilon
        minus = model.evaluate(ell, offsets, shifted)[0]
        errors.append(abs((plus - minus) / (2 * epsilon) - db[index]))
    numerical = (model.evaluate(ell + epsilon, offsets, bulk)[0]
                 - model.evaluate(ell - epsilon, offsets, bulk)[0]) / (2 * epsilon)
    errors.append(abs(numerical - ds))
    assert max(errors) < 2e-7
    root = Path(__file__).resolve().parent.parent
    files = ['tools/strict_core_charge_profiles.py', 'tools/strict_core_charge_gauges.py',
             'tools/check_strict_core_charge_gauges.py', 'research/strict-core-charge-projection.md']
    return {'status': 'pass', 'original_weighted_cut_profiles': records,
            'gradient_components_checked': len(errors), 'gradient_errors': errors,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in files},
            'scope': 'Exact quarter-integer weighted profiles and floating derivative controls; no improved rate claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'original_weighted_cut_profiles': len(result['original_weighted_cut_profiles']),
                      'gradient_components_checked': result['gradient_components_checked'],
                      'maximum_gradient_error': max(result['gradient_errors'])}))


if __name__ == '__main__':
    main()
