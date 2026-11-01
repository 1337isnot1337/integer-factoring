#!/usr/bin/env python3
"""Independent selected-minor and stabilized-code controls for large reflections."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from affine_place_profile_oracle import primitive_seed
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from reflection_dual_rail_gaussian import physical_bit_energies, reflection_letters, tensor
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor
from search_physical_code_gauges import start_params
from selected_gaussian_coefficients import selected_tensor


def direct_restriction(full, xletters, yletters):
    q = full['q_x']
    return [full['values'][x + q * (y + q * (u + q * v))]
            for v in xletters for u in yletters for y in yletters for x in xletters]


def matrix(seed, cut):
    q = seed['q_x']
    result = [[(0, 0)] * (q * q) for _ in range(q * q)]
    for index, value in enumerate(seed['values']):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        row, column = (q * x + y, q * u + v) if cut == 0 else (q * x + u, q * y + v)
        result[row][column] = value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    existing_coefficients = original_minor_coefficients = endpoint_entries = 0
    for interval in (4, 8):
        full = tensor(interval)
        q = full['q_x']
        for first, second in ((tuple(range(q)), tuple(range(q))),
                              ((0, 1), (0, 2, 3)), ((0, 2, 3), (0, 1))):
            actual = restricted_tensor(interval, first, second)
            assert actual['values'] == direct_restriction(full, first, second)
            existing_coefficients += len(actual['values'])
            endpoint_entries += len(actual['values'])
    weights = physical_bit_energies(16)
    energy = lambda symbol: sum(value for bit, value in enumerate(weights) if symbol >> bit & 1)
    letters = sorted(range(256), key=lambda symbol: (energy(symbol), symbol))[:16]
    raw = restricted_tensor(16, letters)
    larger = restricted_tensor(32, [symbol << 8 for symbol in letters])
    scalar_difference = (32 * 32 - 16 * 16) // 2
    assert larger['values'] == [(first << scalar_difference, second << scalar_difference)
                                for first, second in raw['values']]
    assert raw['construction_counts']['cached_pfaffians'] == len(letters) ** 4
    assert larger['construction_counts']['cached_pfaffians'] == len(letters) ** 4
    endpoint_entries += len(raw['values']) + len(larger['values'])
    ambient = reflection_letters(16)
    minor_controls = []
    for indices in ((0, 1, 15), (2, 7, 13), (3, 9, 14), (4, 10, 12)):
        logical = [letters[index] for index in indices]
        expected = selected_tensor(tuple(range(16)), 16, [ambient[symbol] for symbol in logical])
        actual = restricted_tensor(16, logical)
        assert actual['values'] == expected['values']
        original_minor_coefficients += len(actual['values'])
        endpoint_entries += len(actual['values'])
        minor_controls.append({'logical_labels': logical, 'original_selected_labels': [ambient[s] for s in logical],
                               'coefficient_count': len(actual['values']), 'tensor_sha256': digest(actual['values'])})
    seed, common = primitive_seed(raw)
    q, energies = len(letters), np.asarray([energy(symbol) for symbol in letters], dtype=float)
    rectangular = RectangularOracle(seed['values'], q, q)
    combined = rectangular.prefixes(start_params(energies, 1), 1)['orders']
    profiles = []
    for cut, infinity in ((1, True), (0, False)):
        if infinity:
            row_weights = [energies[second] / 2 for first in range(q) for second in range(q)]
            column_weights = [energies[second] / 2 for first in range(q) for second in range(q)]
        else:
            row_weights = [(energies[first] + energies[second]) / 2
                           for first in range(q) for second in range(q)]
            column_weights = [0] * (q * q)
        independent = MatrixOracle(matrix(seed, cut), q * q).profile(row_weights, column_weights, 1, infinity)['orders']
        assert all(value == round(value) for value in independent)
        expected = combined[:q * q] if infinity else combined[q * q:]
        assert sorted(independent) == sorted(expected)
        profiles.append([int(value) for value in sorted(independent)])

    def rho(scale):
        numerator = sum(float(logsumexp(-scale * np.asarray(profile))) for profile in profiles)
        equation = lambda root: numerator - 2 * root * float(logsumexp(-scale * energies / root))
        root = brentq(equation, 1e-10, 4, xtol=1e-13)
        return root / 2

    fit = minimize_scalar(rho, bounds=(.0001, 1), method='bounded', options={'xatol': 1e-11})
    assert abs(float(fit.fun) - .7147366313237278) < 2e-11
    names = (Path(__file__).name, 'restrict_implicit_reflection.py', 'reflection_dual_rail_gaussian.py',
             'selected_gaussian_coefficients.py', 'selected_gaussian_coefficients.cpp',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Exact restricted coefficients, interval stabilization, and complete integer natural profiles; final pressure optimization is numerical.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'counts': {'existing_full_tensor_coefficients': existing_coefficients,
                         'independent_original_minor_coefficients': original_minor_coefficients,
                         'stabilized_complete_coefficients': len(raw['values']),
                         'strict_endpoint_entries': endpoint_entries,
                         'complete_integer_profile_channels': 2 * q * q},
              'large_code': {'interval': 16, 'labels': letters, 'energies': energies.astype(int).tolist(),
                             'scalar_content': common, 'degree': seed['degree'],
                             'tensor_sha256': digest(seed['values']),
                             'raw_construction_counts': raw['construction_counts'],
                             'stabilized_interval': 32, 'common_monomial_difference': scalar_difference,
                             'natural_profiles_infinity_then_zero': profiles,
                             'numerical_rho': float(fit.fun), 'numerical_scale': float(fit.x)},
              'original_minor_controls': minor_controls, 'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': report['counts'], 'numerical_rho': float(fit.fun)}))


if __name__ == '__main__':
    main()
