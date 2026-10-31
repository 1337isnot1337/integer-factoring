#!/usr/bin/env python3
"""Independent support, original-cut, and pressure-gradient controls."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_rectangular_defect_support import ONE, ZERO, Z, UNIT, T, Rat, mm, pairing, car_checks
from native_matrix_smith_f4 import MatrixOracle
from relabelled_gaussian_sums import gaussian_tensor, charge_projection, flatten, natural_heights
from strict_core_charge_profiles import (
    core_data, strict_images, affine_charge_representatives, profile, pressure,
)


def complete_profile(record):
    counts = Counter(record['active_orders'])
    for order in record['spectator_orders']:
        shifted = Counter({value + order: count for value, count in counts.items()})
        counts.update(shifted)
    return counts


def direct_profile(interval, core, charges, place):
    q = 1 << interval
    mask = (1 << core) - 1
    values = charge_projection(gaussian_tensor(tuple(range(interval)), interval),
                               [charges[state & mask] for state in range(q)])
    index = int(place == 'infinity')
    heights = natural_heights(tuple(range(interval)), interval)[index]
    first, second = ((0, 1), (2, 3)) if index == 0 else ((0, 2), (1, 3))
    rows = [heights[first[0]][a] + heights[first[1]][b] for a, b in product(range(q), repeat=2)]
    columns = [heights[second[0]][a] + heights[second[1]][b] for a, b in product(range(q), repeat=2)]
    result = MatrixOracle(flatten(values, interval, bool(index)), q * q).profile(
        rows, columns, infinity=bool(index))
    return Counter(int(value) for value in result['orders'])


def strict_infinity_identity(interval):
    m = interval
    a = [[UNIT * Z ** (1 - 2 * m + 2 * (i | j)) for j in range(m)] for i in range(m)]
    diagonal = [Z ** (1 - 2 * m + 2 * i) for i in range(m)]
    middle = [[a[i][j] / diagonal[j] for j in range(m)] for i in range(m)]
    square = mm(middle, a)
    for i, row in enumerate(square):
        for j, value in enumerate(row):
            assert value == ((ONE + Z ** (2 * m)) * diagonal[i] if i == j else ZERO)
    return {'interval': m, 'exact_bilinear_entries': m * m}


def support_control(interval, core, place):
    initial, _, names = strict_images(place, interval, core)
    data = core_data(place, interval, core)
    final = {name: ([Rat(*entry) for entry in image['creation']],
                    [Rat(*entry) for entry in image['annihilation']])
             for name, image in data['active_images'].items()}
    checks = 0
    for group in names.values():
        checks += car_checks(initial, group)
        checks += car_checks(final, group)
    for first in initial:
        for second in initial:
            assert pairing(initial[first], initial[second]) == pairing(final[first], final[second])
    assert len(data['active_indices']) <= 4 * core
    return {'interval': interval, 'core': core, 'place': place,
            'active_modes': len(data['active_indices']),
            'spectator_modes': len(data['spectator_orders']),
            'regular_shears': len(data['regular_shears']),
            'canonical_pairings': checks, 'all_pairings': len(initial) ** 2,
            'both_conjugate_supports_zero': True,
            'left_and_untransported_right_generators_integral': True}


def verify():
    identities = [strict_infinity_identity(m) for m in (1, 2, 4, 8)]
    support = [support_control(m, k, place)
               for m in (2, 4, 8, 16, 32) for k in (1, 2, 4) if k <= m
               for place in ('zero', 'infinity')]
    cases = [(2, 2, charge) for charge in affine_charge_representatives(2)]
    cases += [(4, 1, (0, 1)), (4, 2, (0, 0, 0, 1)), (4, 2, (0, 1, 2, 0))]
    profiles = []
    gradient_profiles = None
    for m, core, charges in cases:
        current = []
        for place in ('zero', 'infinity'):
            record = profile(place, m, charges, core)
            assert complete_profile(record) == direct_profile(m, core, charges, place)
            profiles.append({'interval': m, 'core': core, 'charges': list(charges), 'place': place,
                             'active_rank': len(record['active_orders']),
                             'complete_rank': sum(complete_profile(record).values())})
            current.append(record)
        if m == 4 and core == 2:
            gradient_profiles = current
    comparisons = []
    for place in ('zero', 'infinity'):
        first = profile(place, 8, (0, 1), 1)
        second = profile(place, 8, (0, 1, 0, 1), 2)
        assert complete_profile(first) == complete_profile(second)
        comparisons.append({'interval': 8, 'place': place,
                            'first_active_modes': first['active_modes'],
                            'second_active_modes': second['active_modes'],
                            'full_rank': sum(complete_profile(first).values())})
    point = np.array([0.37, 0.23, -0.41])
    value, derivative = pressure(4, gradient_profiles, *point, gradient=True)
    errors = []
    for coordinate in range(3):
        shift = np.zeros(3)
        shift[coordinate] = 1e-5
        numerical = (pressure(4, gradient_profiles, *(point + shift))
                     - pressure(4, gradient_profiles, *(point - shift))) / (2e-5)
        errors.append(abs(numerical - derivative[coordinate]))
    assert max(errors) < 1e-7
    root = Path(__file__).resolve().parent.parent
    files = ['tools/strict_core_charge_profiles.py', 'tools/check_strict_core_charge_profiles.py']
    return {'status': 'pass', 'strict_infinity_identities': identities,
            'support_cases': support, 'original_full_cut_profiles': profiles,
            'different_core_spectator_controls': comparisons,
            'pressure_gradient_errors': errors,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in files},
            'scope': 'Exact bounded profiles and independent finite controls; no improved pressure claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'support_cases': len(result['support_cases']),
                      'original_full_cut_profiles': len(result['original_full_cut_profiles']),
                      'different_core_spectator_controls': len(result['different_core_spectator_controls']),
                      'maximum_gradient_error': max(result['pressure_gradient_errors'])}))


if __name__ == '__main__':
    main()
