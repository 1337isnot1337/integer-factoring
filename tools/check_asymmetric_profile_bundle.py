#!/usr/bin/env python3
"""Held-out plane and derivative controls for asymmetric fixed-type bundles."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from asymmetric_profile_bundle import AsymmetricBundleOracle, raw_from_scaled, scaled_from_raw
from mixed_graph_retractions import natural, tensor
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(60520261102)
    records, derivative_errors, planes, minimum = [], [], 0, 0.
    coefficients = tuple(int(x) for x in rng.integers(0, 4, 20))
    seed = tensor(coefficients)
    for logtan in (-1.1, 0., 1.4):
        oracle = AsymmetricBundleOracle(seed, logtan)
        points = []
        for _ in range(8):
            scaled = np.r_[natural()[:-1]+rng.normal(0, .3, 32), math.log(.21)+rng.normal(0, .1), logtan]
            raw = raw_from_scaled(scaled)
            assert max(abs(scaled_from_raw(raw, logtan)-scaled)) < 1e-12
            result = oracle.query(raw)
            rate, value, gradient, details, physical = result
            fixed_type_rate = (value+physical@raw)/details['entropy_sum']
            assert abs(fixed_type_rate-rate) < 1e-10
            points.append((raw, value, gradient))
            if len(points) == 1:
                for coordinate in range(33):
                    delta = np.zeros(33)
                    delta[coordinate] = 1e-6
                    finite = (oracle.query(raw+delta)[1]-oracle.query(raw-delta)[1])/2e-6
                    error = abs(float(finite-gradient[coordinate]))
                    assert error < 2e-5, (logtan, coordinate, error)
                    derivative_errors.append(error)
            records.append({'log_tangent': logtan, 'raw_parameters': raw.tolist(),
                            'log_profile_value': value, 'subgradient': gradient.tolist(),
                            'support_value': rate})
        for x, f, g in points:
            for y, fy, _ in points:
                if x is y:
                    continue
                slack = float(fy-f-g@(y-x))
                minimum = min(minimum, slack)
                assert slack >= -2e-7, slack
                planes += 1
        # Construct the scaled native model last: the native backend has one active seed.
        scaled_model = ScaledAsymmetricPressure(seed)
        for raw, f, _ in points:
            rate, _ = scaled_model.evaluate(scaled_from_raw(raw, logtan), False)
            # Recover the root from the corresponding saved query.
            record = next(v for v in records if v['log_tangent'] == logtan and v['raw_parameters'] == raw.tolist())
            assert abs(rate-record['support_value']) < 1e-10
    names = ['asymmetric_profile_bundle.py', 'check_asymmetric_profile_bundle.py',
             'asymmetric_crossing_pressure.py', 'mixed_graph_retractions.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp']
    report = {'status': 'pass', 'nonbinary_coefficients': coefficients,
              'finite_log_profile_derivatives': len(derivative_errors),
              'maximum_derivative_error': max(derivative_errors),
              'directed_heldout_plane_inequalities': planes, 'minimum_slack': minimum,
              'records': records,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Numerical checks of complete-profile bundle ingredients; not an exact inequality '
                       'certificate, a global minimum or an improving recursion.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'finite_log_profile_derivatives',
                     'directed_heldout_plane_inequalities', 'maximum_derivative_error')}))


if __name__ == '__main__':
    main()
