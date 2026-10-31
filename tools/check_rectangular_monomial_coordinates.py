#!/usr/bin/env python3
"""Independent dense controls for mixed/full-alphabet odd core monomials."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_rectangular_core_monomials import core_monomial_seed
from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import point_control, weighted_controls
from check_full_alphabet_quartic_profiles import profile_controls
from rectangular_coordinate_transforms import linear_combination_coordinate_matrices
from rectangular_defect_coordinates import coordinate_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    minus = (0, 3, 0, None, [(2, -1)], 0)
    pair = [(0, 3, 0, None, [(1, 0)], 0), (0, 0, None, 3, [(2, -2)], 0)]
    rational = (None, 3, 0, 3, [(3, 1)], 1)
    full_pair = [(0, None, 0, 3, [(1, 0)], 0), (3, 3, None, 3, [(2, 2)], 0)]
    cases = [('mixed_single', None, False, [minus]),
             ('mixed_pair_and_pole', None, False, pair+[rational]),
             ('full_gaussian_pair', 0, False, pair),
             ('full_gamma2_single', 2, False, [minus]),
             ('full_gamma2_pair_and_pole', 2, False, pair+[rational]),
             ('full_gamma2_xT_pair', 2, True, full_pair)]
    seeds = {name: core_monomial_seed(extras, tuple(range(4)), 4,
                                     full_alphabet_gamma=gamma, check_types=False)
             for name, gamma, _, extras in cases}
    records = []
    for place in ('infinity', 'zero'):
        for name, gamma, full, extras in cases:
            data = data_for(place, 4, full_core=full)
            compiled = coordinate_matrix(data, full_alphabet_gamma=gamma, extra_monomials=extras)
            assert compiled[2]['exact_polynomial_matrix']
            seed = seeds[name]
            power = seed['global_w_power']*(-1 if place == 'infinity' else 1)
            if power:
                compiled = linear_combination_coordinate_matrices([(compiled, [(1, power)])])
            points = [point_control(data, compiled, z, seed) for z in (2, 3)]
            profiles = (weighted_controls(data, compiled, 2, seed)['cases'] if gamma is None
                        else profile_controls(data, compiled, seed, cases=2))
            records.append({'case': name, 'place': place, 'M': 4, 'full_core': full,
                            'full_alphabet_gamma': gamma, 'extra_monomials': extras,
                            'global_w_power': seed['global_w_power'],
                            'entrywise_checks': points, 'weighted_profiles': profiles})
            print(json.dumps({'case': name, 'place': place,
                              'entries': sum(point['entries_checked'] for point in points)}), flush=True)
    result = {'status': 'independent dense CAR and gauge controls for odd core perturbations',
              'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
