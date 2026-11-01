#!/usr/bin/env python3
"""Independent normalization and derivative controls for numerical pressure.

These tests validate the search interface; they do not certify a pressure
inequality or a new factoring exponent.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from physical_split_crossings import identity_split
from reflection_dual_rail_gaussian import tensor
from search_physical_code_gauges import start_params
from search_reflection_physical_bases import change_basis, row_shear


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    base = tensor(4)
    eye = identity_split(4)
    left = row_shear(row_shear(eye, 0, 3, 2), 1, 0, 3)
    right = row_shear(eye, 3, 0, 1)
    variants = [('natural', base), ('nontrivial-basis', change_basis(base, left, right))]
    rng = np.random.default_rng(312108)
    records = []
    for name, seed in variants:
        for cut in (0, 1):
            oracle = AsymmetricPressure(seed, cut)
            params = start_params(np.array([0, 12, 4, 16]) * .17, .17)
            params[:-1] += rng.normal(0, .12, len(params) - 1)
            ordinary = oracle.oracle.types(params, cut)[0]
            balanced = oracle.evaluate(np.r_[params / math.sqrt(2), math.pi / 4], gradient=False)[0]
            normalized_error = abs(balanced - math.sqrt(2) * ordinary)
            assert normalized_error < 2e-10
            point = np.r_[params, .41 if cut else 1.07]
            value, derivative, details = oracle.evaluate(point)
            errors = []
            epsilon = 1e-6
            for j in range(len(point)):
                direction = np.zeros(len(point))
                direction[j] = epsilon
                positive = oracle.evaluate(point + direction, gradient=False)[0]
                negative = oracle.evaluate(point - direction, gradient=False)[0]
                errors.append(abs((positive - negative) / (2 * epsilon) - derivative[j]))
            assert max(errors) < 2e-6
            shift_controls = []
            for shift in (1, 7):
                shifted = dict(seed, values=[(lo << shift, hi << shift) for lo, hi in seed['values']])
                shifted_oracle = AsymmetricPressure(shifted, cut)
                changed = shifted_oracle.evaluate(point, gradient=False)[0]
                error = abs(changed - value)
                assert error < 2e-10
                shift_controls.append({'shift': shift, 'absolute_error': error})
            records.append({'basis': name, 'cut_at_infinity': cut,
                            'balanced_normalization_error': normalized_error,
                            'finite_difference_coordinates': len(point),
                            'maximum_derivative_error': max(errors),
                            'support_value_at_test': value, 'monomial_shift_controls': shift_controls})
    names = ('check_asymmetric_pressure_oracle.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'search_reflection_physical_bases.py')
    result = {'status': 'pass', 'scope': 'numerical search-interface controls; no rate certificate',
              'records': records,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'complete_normalization_controls': len(records),
                      'derivative_coordinates': sum(r['finite_difference_coordinates'] for r in records),
                      'monomial_controls': sum(len(r['monomial_shift_controls']) for r in records),
                      'maximum_derivative_error': max(r['maximum_derivative_error'] for r in records)}))


if __name__ == '__main__':
    main()
