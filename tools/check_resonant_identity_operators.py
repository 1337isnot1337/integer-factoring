#!/usr/bin/env python3
"""Full Weyl, pressure, and original physical controls for resonant blocks."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import random

import numpy as np

from resonant_identity_operators import determinant, kernel, pressure, profiles, signature
from f4_repaired_identity_filters import add, evaluate_one
from f4_identity_operator_filters import zero_profiles, infinity_profiles, polynomial_seed, profile_data
from check_f4_repaired_identity_filters import expanded, physical_profiles
from search_coupled_identity_operator_filters import Pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(4201337)
    controls, pressure_checks = [], 0
    for n, P, selected in ((4, 1, (0, 1, 2, 3)), (8, 1, (0, 1, 6, 7)),
                           (8, 3, (4, 5, 6, 7)), (32, 1, (0, 7, 24, 31)),
                           (32, 1, (28, 29, 30, 31))):
        geometry = 'selected:'+','.join(map(str, selected))
        for case in range(8):
            while True:
                block = [[{degree: rng.randrange(4) for degree in (-1, 0, 1)} for _ in range(2)] for _ in range(2)]
                block = [[{degree: value for degree, value in entry.items() if value} for entry in row] for row in block]
                if evaluate_one(determinant(block)):
                    break
            values = signature(n, P, block)
            compressed = profiles(n, P, selected, values)
            operator = kernel(block)
            expected = (zero_profiles(n, P, 4, operator, geometry), infinity_profiles(n, P, 4, operator, geometry))
            assert compressed == expected
            scales = math.log(2*n)/n*np.array((.5, 1., 2.))
            charges = np.array((-2., 0., float(n)))
            direct = Pressure(n, P, 4, geometry, *expected).grid(scales, charges)
            actual = pressure(n, P, selected, values, scales, charges)
            assert np.max(np.abs(direct-actual)) < 1e-10
            pressure_checks += len(scales)
            controls.append({'n': n, 'ramification': P, 'selected': selected,
                             'matrix': [[sorted(value.items()) for value in row] for row in block],
                             'signature': values})
    late = [[{0: 1}, {0: 2, 600: 2}], [{}, {0: 1}]]
    values = signature(4, 1, late)
    assert profiles(4, 1, (0, 1, 2, 3), values) == (
        zero_profiles(4, 1, 4, kernel(late), 'selected:0,1,2,3'),
        infinity_profiles(4, 1, 4, kernel(late), 'selected:0,1,2,3'))
    endpoint = ((1, 2), (1, 1))
    correction = ((1, 3), (2, 0))
    block = [[add({0: endpoint[i][j]}, {0: correction[i][j], -1: correction[i][j]})
              for j in range(2)] for i in range(2)]
    operator = kernel(block)
    tensor, clearing, gaussian = polynomial_seed(4, 1, 4, operator, 'selected:0,1,2,3')
    values = signature(4, 1, block)
    zero, infinity = profiles(4, 1, (0, 1, 2, 3), values)
    data = profile_data(4, 1, 0, 4, operator, 'selected:0,1,2,3', zero)
    expected = (
        tuple(sorted(a+b+2*clearing for a in data['zero_active'] for b in expanded(data['zero_spectators']))),
        tuple(sorted(a+b-2*clearing for a in data['infinity_active'] for b in expanded(data['infinity_spectators']))),
    )
    actual, precision = physical_profiles(tensor, 4, 1, 0)
    assert actual == expected
    result = {'status': 'Exact closed resonance profiles, full Weyl comparisons, and original physical controls; no pressure claim.',
              'complete_two_cut_controls': controls, 'independent_pressure_values': pressure_checks,
              'late_extension_coefficient_degree': 600,
              'original_physical_control': {'n': 4, 'ramification': 1, 'selected': (0, 1, 2, 3),
                                            'endpoint': endpoint, 'correction': correction,
                                            'matrix_order': 256, 'cut_profiles': 2,
                                            'endpoint_coefficients_checked': len(set(tensor)|set(gaussian)),
                                            'maximum_native_precision': precision}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'complete_two_cut_comparisons': len(controls)+1,
                      'independent_pressure_values': pressure_checks,
                      'original_physical_cut_profiles': 2,
                      'original_endpoint_coefficients_checked': len(set(tensor)|set(gaussian))}))


if __name__ == '__main__':
    main()
