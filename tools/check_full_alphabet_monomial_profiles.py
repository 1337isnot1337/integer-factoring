#!/usr/bin/env python3
"""Independent dense controls of full-alphabet rank-restoring cubic repairs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_full_alphabet_monomial_residual import (
    full_residual_monomial_seed, cubic_pair, prepare_residual_repair,
)
from check_full_alphabet_quartic_repair import partial_repair
from check_full_alphabet_quartic_profiles import profile_controls
from check_rectangular_defect_coordinates import point_control
from check_rectangular_defect_core_basis import data_for
from rectangular_coordinate_transforms import linear_combination_coordinate_matrices
from rectangular_defect_coordinates import coordinate_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    for gamma in (0, 2, 3):
        assert prepare_residual_repair(gamma, [])['repair_matrix'] == partial_repair(gamma)[0]
    cases = [('rank_restore_02_binary', 1, 0, 2, 1, 'partial', False, 0, 0),
             ('rank_restore_13_nonbinary', 1, 1, 3, 2, 'partial', False, 0, 0),
             ('rank_restore_07_full_core', 1, 0, 7, 3, 'partial', True, 0, 0),
             ('rank_restore_16_inverse', 1, 1, 6, 2, 'inverse', False, -2, 1),
             ('critical_gaussian_base', 0, 0, 2, 2, 'partial', False, 1, -1),
             ('critical_nonbinary_base', 2, 0, 2, 1, 'partial', False, -1, 2)]
    prepared = []
    for name, gamma, minus, plus, coefficient, repair, full, a, b in cases:
        extras = cubic_pair(3, minus, plus, coefficient, a, b)
        if name == 'critical_nonbinary_base':
            # The extra numerator cancels one c pole exactly; the residual is unchanged.
            modes, terms, pole = extras[0][:4], extras[0][4], extras[0][5]
            extras[0] = (*modes, [(1, a), (1, a+1)], pole+1)
            extras += [(None, 3, 0, 3, [(2, 1)], 0)]
        seed = full_residual_monomial_seed(extras, tuple(range(4)), 4, gamma, repair, check_types=False)
        # Represent every canonical strength, including one, on the gamma-zero native base.
        native_spec = ([(0, 3, 0, 3, [(gamma, 0)], 2)] if gamma else [])+extras
        prepared.append((name, repair, full, native_spec, seed))
    records = []
    for place in ('infinity', 'zero'):
        for name, repair, full, specifications, seed in prepared:
            data = data_for(place, 4, full_core=full)
            compiled = coordinate_matrix(data, full_alphabet_gamma=0, repair_kind=repair,
                                         residual_monomials=specifications)
            assert compiled[2]['exact_polynomial_matrix']
            power = seed['global_w_power']*(-1 if place == 'infinity' else 1)
            if power:
                compiled = linear_combination_coordinate_matrices([(compiled, [(1, power)])])
            points = [point_control(data, compiled, z, seed) for z in (2, 3)]
            profiles = profile_controls(data, compiled, seed, cases=2)
            records.append({'case': name, 'place': place, 'M': 4, 'full_core': full,
                            'repair_kind': repair, 'residual_monomials': specifications,
                            'global_w_power': seed['global_w_power'],
                            'entrywise_checks': points, 'weighted_profiles': profiles})
            print(json.dumps({'case': name, 'place': place,
                              'entries': sum(point['entries_checked'] for point in points)}), flush=True)
    result = {'status': 'independent full tensor/CAR controls for critical odd residual corrections',
              'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
