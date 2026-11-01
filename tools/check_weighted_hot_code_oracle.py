#!/usr/bin/env python3
"""Bind the block-gradient interface to complete literal-profile controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from affine_place_profile_oracle import primitive_seed
from asymmetric_crossing_pressure import AsymmetricPressure
from check_weighted_hot_compounds import logical
from hot_compound_pressure import evaluate as mirrored_pressure
from restrict_implicit_reflection import restricted_tensor
from search_affine_coupled_reflection_frames import balanced
from weighted_hot_code_oracle import HotWeightedOracle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    evidence = Path(__file__).resolve().parents[1] / 'evidence/weighted-hot-compound-profiles.json'
    raw = evidence.read_bytes()
    controls = json.loads(raw)
    assert controls['passed']
    cache, channels, inequalities = {}, 0, 0
    maximum_error = maximum_gradient_error = 0.0
    for record in controls['records']:
        key = record['interval'], tuple(record['x_letters']), tuple(record['y_letters'])
        if key not in cache:
            interval, first, second = key
            literal = restricted_tensor(interval, [logical(x, interval // 2) for x in first],
                                         [logical(y, interval // 2) for y in second])
            _, common = primitive_seed(literal)
            assert len(common) == 1
            shift = interval * interval // 2 - next(iter(common))
            cache[key] = HotWeightedOracle(first, second), shift
        oracle, shift = cache[key]
        for side, field in ((0, 'additional_infinity_weights'), (1, 'additional_zero_weights')):
            legs = [[float(Fraction(value)) for value in leg] for leg in record[field]]
            if side == 0:
                legs = [legs[0], legs[2], legs[1], legs[3]]
            heights = np.concatenate(legs)
            cone, bounds = oracle.cone(1)
            assert np.all(cone @ heights >= bounds)
            inequalities += len(bounds)
            actual = oracle.log_profile_with_gradient(heights, 1)[3] + (shift if side else -shift)
            expected = np.asarray(record['complete_orders_times_four'][side]) / 4
            error = float(np.max(np.abs(actual - expected)))
            assert error == 0
            maximum_error = max(maximum_error, error)
            channels += len(actual)
    rng = np.random.default_rng(387612)
    gradient_coordinates = 0
    pressure_controls = []
    for key, (oracle, _) in cache.items():
        heights = rng.uniform(-.1, .1, oracle.height_count)
        scale = .173846
        _, gradient, scale_gradient, _ = oracle.log_profile_with_gradient(heights, scale)
        step = 2e-7
        for index in range(oracle.height_count):
            plus, minus = heights.copy(), heights.copy()
            plus[index] += step; minus[index] -= step
            numeric = (oracle.log_profile_with_gradient(plus, scale)[0]
                       - oracle.log_profile_with_gradient(minus, scale)[0]) / (2 * step)
            error = abs(numeric - gradient[index])
            assert error < 3e-7
            maximum_gradient_error = max(maximum_gradient_error, float(error))
            gradient_coordinates += 1
        numeric = (oracle.log_profile_with_gradient(heights, scale + step)[0]
                   - oracle.log_profile_with_gradient(heights, scale - step)[0]) / (2 * step)
        error = abs(numeric - scale_gradient)
        assert error < 3e-7
        maximum_gradient_error = max(maximum_gradient_error, float(error))
        gradient_coordinates += 1
        interval, first, second = key
        legs = []
        for code in (first, second, second, first):
            potentials = rng.uniform(.01, .2, len(code))
            legs.append([sum(potentials[j] for j, subset in enumerate(code)
                             if subset and subset & value == subset) for value in code])
        heights = np.concatenate(legs)
        rate, derivative, scale_derivative, details = mirrored_pressure(oracle, heights, scale)
        literal = restricted_tensor(interval, [logical(x, interval // 2) for x in first],
                                     [logical(y, interval // 2) for y in second])
        original = AsymmetricPressure(literal)
        params = []
        for side, extra in enumerate(([legs[0], legs[2], legs[1], legs[3]], legs)):
            for leg, values in enumerate(extra):
                energies = oracle.energies_x if leg in (0, 3) else oracle.energies_y
                for index in range(1, len(values)):
                    natural = scale * energies[index] / 2 if leg in ((2, 3) if side == 0 else (0, 1)) else 0
                    params.append(values[index] + natural)
        actual = balanced(original, np.r_[params, scale])[0]
        assert abs(actual - rate) < 2e-11
        pressure_controls.append({'interval': interval, 'x_letters': first, 'y_letters': second,
                                  'rho': rate, 'literal_rho': actual})
        for index in rng.choice(oracle.height_count, min(8, oracle.height_count), replace=False):
            plus, minus = heights.copy(), heights.copy()
            plus[index] += step; minus[index] -= step
            numeric = (mirrored_pressure(oracle, plus, scale)[0]
                       - mirrored_pressure(oracle, minus, scale)[0]) / (2 * step)
            error = abs(numeric - derivative[index])
            assert error < 3e-7
            maximum_gradient_error = max(maximum_gradient_error, float(error))
            gradient_coordinates += 1
        numeric = (mirrored_pressure(oracle, heights, scale + step)[0]
                   - mirrored_pressure(oracle, heights, scale - step)[0]) / (2 * step)
        error = abs(numeric - scale_derivative)
        assert error < 3e-7
        maximum_gradient_error = max(maximum_gradient_error, float(error))
        gradient_coordinates += 1
    names = (Path(__file__).name, 'weighted_hot_code_oracle.py', 'hot_code_charge_profiles.py',
             'hot_code_charge_profiles.cpp', 'selected_gaussian_coefficients.cpp',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp', 'restrict_implicit_reflection.py',
             'hot_compound_pressure.py', 'asymmetric_crossing_pressure.py',
             'search_affine_coupled_reflection_frames.py')
    report = {'scope': 'Complete block profiles bound to the independently checked literal reflection records; derivatives checked inside generic Smith chambers.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'literal_control_sha256': sha256(raw).hexdigest(),
              'counts': {'literal_profile_channels': channels, 'cone_inequalities': inequalities,
                         'gradient_coordinates': gradient_coordinates,
                         'literal_common_type_pressure_controls': len(pressure_controls)},
              'pressure_controls': pressure_controls,
              'maximum_order_error': maximum_error, 'maximum_gradient_error': maximum_gradient_error,
              'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': report['counts'],
                      'maximum_gradient_error': maximum_gradient_error}))


if __name__ == '__main__':
    main()
