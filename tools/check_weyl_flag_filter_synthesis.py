#!/usr/bin/env python3
"""Exact controls for constant F4 flags and complete F16 minor profiles."""
from itertools import product
import json
from pathlib import Path
import random

import f16_diagonal_filter_profiles as profile
from f4_repaired_identity_filters import reduced_zero_profiles
from weyl_flag_filter_synthesis import (
    coefficient_tables, flag_basis, full_flags, identity_flag, synthesize,
)
from weyl_lattice_filter_synthesis import synthesize as old_synthesize


def main():
    rng = random.Random(0xF416F1A6)
    matrix_controls = 0
    for size in (1, 2, 3):
        for trial in range(12):
            orders = sorted(rng.randrange(16) for _ in range(size))
            matrix = [[profile.encode({orders[i]: 1}) if i == j else profile.ZERO
                       for j in range(size)] for i in range(size)]
            if size > 1:
                for step in range(6):
                    i, j = rng.sample(range(size), 2)
                    degree = 600 if trial == 0 and step == 0 else rng.randrange(8)
                    scalar = profile.encode({degree: rng.randrange(1, 16)})
                    if step & 1:
                        for row in matrix:
                            row[i] = profile.add(row[i], profile.multiply(scalar, row[j]))
                    else:
                        matrix[i] = [profile.add(a, profile.multiply(scalar, b))
                                     for a, b in zip(matrix[i], matrix[j])]
            assert profile.smith_orders(matrix) == tuple(orders)
            matrix_controls += 1
    native_controls = 0
    for n, P, geometry, modes in ((2, 1, 'initial', 2), (4, 1, 'mixed', 2),
                                   (4, 3, 'terminal', 2), (4, 1, 'initial', 3),
                                   (4, 1, 'split:1:2', 3)):
        states = tuple(product(range(2), repeat=modes))
        for trial in range(6):
            alpha = {state: {rng.randrange(6): rng.randrange(1, 4)} for state in states}
            actual = profile.zero_profiles(n, P, modes, alpha, geometry)
            control = reduced_zero_profiles(n, P, modes, alpha, geometry)
            assert tuple((b['charge'], tuple(b['orders'])) for b in actual) == tuple(
                (b['charge'], tuple(b['orders'])) for b in control)
            native_controls += 1
    comparison_controls = 0
    for n, geometry in ((2, 'initial'), (4, 'mixed')):
        identity = (0, identity_flag(2), 0)
        base = flag_basis(n, 1, geometry, 2, identity)
        for trial in range(80):
            degrees = tuple(rng.randrange(3) for _ in range(4))
            targets = tuple(value + rng.randrange(-2, 3) for value in base['base']['gaussian_orders'])
            old = old_synthesize(n, 1, geometry, degrees, targets, construct=False, modes=2)
            new = synthesize(n, 1, geometry, degrees, targets, identity, modes=2,
                             construct=False, field_size=4)
            assert bool(old) == bool(new)
            comparison_controls += 1
    witness_controls = []
    for n, P, geometry, modes in ((3, 1, 'initial', 2), (4, 3, 'mixed', 2),
                                   (4, 1, 'initial', 3), (4, 1, 'split:1:2', 3)):
        for trial in range(4):
            if modes == 2:
                flags = (0, rng.randrange(5), 0)
            else:
                flags = (0, rng.randrange(105), rng.randrange(105), 0)
            data = flag_basis(n, P, geometry, modes, flags)
            degrees = tuple(rng.choice((0, 1, n, n + 1)) for _ in data['base']['states'])
            _, _, equations = coefficient_tables(n, P, geometry, modes, flags, degrees)
            target = list(data['minimum_rows'])
            priority = list(range(len(target)))
            rng.shuffle(priority)
            for index in priority:
                low = target[index]
                high = max(power for power, equation in equations[index]) + 1
                while high - low > 1:
                    middle = (low + high) // 2
                    target[index] = middle
                    if synthesize(n, P, geometry, degrees, target, flags, modes, False):
                        low = middle
                    else:
                        high = middle
                target[index] = low
            witness = synthesize(n, P, geometry, degrees, target, flags, modes)
            verification = profile.verify_flag_witness(witness)
            witness_controls.append({'n': n, 'ramification': P, 'geometry': geometry,
                                     'modes': modes, 'flags': flags, 'degree_bounds': degrees,
                                     'row_targets': target, 'endpoint_values': witness['endpoint_values'],
                                     'actual_degrees': witness['actual_degrees'],
                                     'verification': verification})
    result = {'status': 'All exact flag/profile controls passed.',
              'full_flag_counts': {'dimension2': len(full_flags(2)), 'dimension3': len(full_flags(3))},
              'known_Smith_matrices_over_F16': matrix_controls,
              'late_coefficient_degree': 600, 'independent_native_F4_profiles': native_controls,
              'old_and_new_identity_lattice_comparisons': comparison_controls,
              'polynomial_flag_witnesses': witness_controls}
    path = Path(__file__).resolve().parents[1] / 'evidence' / 'weyl-flag-filter-synthesis.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'polynomial_flag_witnesses'}))
    print('Verified polynomial witnesses:', len(witness_controls))


if __name__ == '__main__':
    main()
