#!/usr/bin/env python3
"""Independent normalization and derivative controls for mixed support costs."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from affine_place_profile_oracle import AffinePlaceOracle
from asymmetric_crossing_pressure import AsymmetricPressure
from asymmetric_mixed_place_pressure import MixedSupportingPressure, ScaledMixedSupportingPressure
from reflection_dual_rail_gaussian import tensor
from search_affine_coupled_reflection_frames import affine_seed, catalogue
from search_fused_mixed_weighted import MixedPressure
from search_physical_code_gauges import start_params


def baseline(row, column):
    return {'row_places': row, 'column_places': column,
            'row_allocations': [1 / len(row)] * len(row),
            'column_allocations': [1 / len(column)] * len(column),
            'outer_log_scale': .13}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(312431)
    base = tensor(4)
    frames = catalogue(base)
    cases = [('base', base), ('finite_F4', affine_seed(base, frames, [(56, 2), (58, 2)])),
             ('quadratic', affine_seed(base, frames, [(19, 3), (89, 2)]))]
    partitions = [baseline([0], [None]), baseline([0, 2], [None, 3]),
                  baseline([0, 'root:(6,1):a', 'root:(6,1):b'],
                           [None, 'root:(4,1):a', 'root:(4,1):b'])]
    controls, derivatives, native_controls = [], [], []
    for label, seed in cases:
        oracle = AffinePlaceOracle(seed)
        for partition_index, partition in enumerate(partitions):
            old = MixedPressure(oracle, partition, 'full', False, True, math.inf)
            new = MixedSupportingPressure(oracle, partition, 'full', False, True, math.inf)
            values = old.initial(partition, rng, 0)
            for slices in old.height_slices:
                for current in slices:
                    values[current] = rng.normal(0, .2, old.height_size)
            values[old.type_slice] = rng.normal(0, .2, 6)
            for current in old.allocation_slices:
                values[current] = rng.normal(0, .15, current.stop - current.start)
            converted = np.r_[values, math.pi / 4]
            converted[new.scale_index] /= math.sqrt(2)
            for slices in new.height_slices:
                for current in slices:
                    converted[current] /= math.sqrt(2)
            expected = math.sqrt(2) * old.evaluate(values, smoothing=.007)[0]
            actual = new.evaluate(converted, smoothing=.007 / math.sqrt(2))[0]
            assert abs(actual - expected) < 2e-10, (label, partition_index, actual, expected)
            controls.append({'seed': label, 'partition': partition_index,
                             'old_balanced_rho': expected / math.sqrt(2),
                             'balanced_support': actual, 'absolute_error': abs(actual - expected)})
            # Each finite-difference point is away from integer-weight ties.
            if partition_index == 1:
                converted[new.angle_index] = .59
                _, gradient = new.evaluate(converted, smoothing=.013, record=False)
                epsilon = 2e-6
                observed = []
                for index in range(new.dimension):
                    plus, minus = converted.copy(), converted.copy()
                    plus[index] += epsilon
                    minus[index] -= epsilon
                    observed.append((new.evaluate(plus, smoothing=.013, record=False)[0]
                                     - new.evaluate(minus, smoothing=.013, record=False)[0]) / (2 * epsilon))
                error = float(np.max(np.abs(gradient - observed)))
                assert error < 2e-6, (label, 'unscaled', error)
                derivatives.append({'seed': label, 'parameterization': 'branch-weighted heights',
                                    'coordinates': new.dimension, 'maximum_error': error})
                scaled = ScaledMixedSupportingPressure(oracle, partition, 'full', False, True, math.inf)
                scaled_values = converted.copy()
                for slices in new.height_slices:
                    for current in slices:
                        scaled_values[current] /= converted[new.scale_index]
                scaled_value, gradient = scaled.evaluate(scaled_values, smoothing=.013, record=False)
                assert abs(scaled_value - new.evaluate(converted, smoothing=.013, record=False)[0]) < 1e-10
                observed = []
                for index in range(scaled.dimension):
                    plus, minus = scaled_values.copy(), scaled_values.copy()
                    plus[index] += epsilon
                    minus[index] -= epsilon
                    observed.append((scaled.evaluate(plus, smoothing=.013, record=False)[0]
                                     - scaled.evaluate(minus, smoothing=.013, record=False)[0]) / (2 * epsilon))
                error = float(np.max(np.abs(gradient - observed)))
                assert error < 2e-6, (label, 'scaled', error)
                derivatives.append({'seed': label, 'parameterization': 'b=t*B',
                                    'coordinates': scaled.dimension, 'maximum_error': error})
        # The separate two-place native implementation optimizes types via
        # its Gibbs root. Insert those types into the explicit mixed formula.
        native = AsymmetricPressure(seed, cut=1)
        partition = partitions[0]
        mixed = MixedSupportingPressure(oracle, partition, 'full', False, True, math.inf)
        for angle in (.24, math.pi / 4, 1.21):
            params = start_params(np.asarray([0, 12, 4, 16]) * .12, .12) / math.sqrt(2)
            params[:-1] += rng.normal(0, .017, len(params) - 1)
            native.oracle.activate()
            expected, data = native.evaluate(np.r_[params, angle], gradient=False)
            values = mixed.initial(partition, rng, 0)
            values[mixed.height_slices[0][0]] = params[12:24]
            values[mixed.height_slices[1][0]] = params[:12]
            values[mixed.scale_index] = params[-1]
            values[mixed.angle_index] = math.pi / 2 - angle
            values[mixed.type_slice] = np.concatenate([np.log(np.asarray(p)[1:] / p[0])
                                                       for p in (data['type_x'], data['type_y'])])
            actual = mixed.evaluate(values)[0]
            assert abs(actual - expected) < 2e-9, (label, angle, actual, expected)
            native_controls.append({'seed': label, 'native_angle': angle,
                                    'native_support': expected, 'mixed_support': actual,
                                    'absolute_error': abs(actual - expected)})
    names = (Path(__file__).name, 'asymmetric_mixed_place_pressure.py',
             'asymmetric_crossing_pressure.py', 'search_fused_mixed_weighted.py',
             'affine_place_profile_oracle.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Floating objective equivalence and derivative controls. No numerical optimum or factoring improvement certified.',
              'balanced_controls': controls, 'derivative_controls': derivatives,
              'native_two_place_controls': native_controls,
              'total_derivative_coordinates': sum(row['coordinates'] for row in derivatives),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: len(value) if isinstance(value, list) else value
                      for key, value in report.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
