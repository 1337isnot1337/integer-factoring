#!/usr/bin/env python3
"""Audit extra-quartic CAR programs against independent original tensors."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import point_control, weighted_controls
from check_rectangular_extra_quartics import PATTERNS, extra_quartic_seed
from rectangular_coordinate_transforms import linear_combination_coordinate_matrices
from rectangular_defect_coordinates import coordinate_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    for label, pattern in enumerate(PATTERNS):
        coefficient, ell = ((2, 0), (3, -2), (1, 3), (2, -1))[label]
        extra = (*(mode*3 for mode in pattern), [(coefficient, ell)], 0)
        cases.append((f'pattern_{label}', [extra], label == 2))
    extras = [(*(mode*3 for mode in pattern), [(1+label % 3, label-2)], 0)
              for label, pattern in enumerate(PATTERNS)]
    # The numerator's zero at one cancels this marginal c denominator.
    extras.append((0, 3, 0, 0, [(2, 0), (2, 1)], 1))
    cases.append(('mixed_with_c_cancellation', extras, True))
    seeds = {name: extra_quartic_seed(extras, tuple(range(4)), 4, check_types=False)
             for name, extras, _ in cases}
    records = []
    for place in ('infinity', 'zero'):
        for name, extras, full in cases:
            data = data_for(place, 4, full_core=full)
            compiled = coordinate_matrix(data, extra_quartics=extras)
            assert compiled[2]['exact_polynomial_matrix']
            seed = seeds[name]
            power = seed['global_w_power']*(-1 if place == 'infinity' else 1)
            if power:
                compiled = linear_combination_coordinate_matrices([(compiled, [(1, power)])])
            points = [point_control(data, compiled, z, seed) for z in (2, 3)]
            profiles = weighted_controls(data, compiled, 2, seed)
            records.append({'case': name, 'place': place, 'M': 4, 'full_core': full,
                            'extra_quartics': extras, 'global_w_power': seed['global_w_power'],
                            'entrywise_checks': points, 'weighted_profiles': profiles})
            print(json.dumps({'case': name, 'place': place,
                              'entry_checks': sum(point['entries_checked'] for point in points)}), flush=True)
    result = {'status': 'independent original-tensor extra-quartic CAR and gauge controls',
              'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
