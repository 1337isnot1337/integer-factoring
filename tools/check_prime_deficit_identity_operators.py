#!/usr/bin/env python3
"""Independent original, rectangular, and full-cost checks of prime deficit operators."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from coupled_identity_operator_filters import identity_kernel, states
from coupled_identity_weyl_profiles import zero_profiles as existing_profiles
from check_prime_state_filter_profiles import physical_cut_profiles, expand_spectators
from prime_deficit_identity_operators import (
    geometry_data, zero_profiles, full_zero_profiles, infinity_profiles,
    profile_data, polynomial_seed, Pressure,
)
from search_qutrit_deficit_operators import projectors, rank_one_catalogue


def example(modes, prime=3, filtered=False):
    result = identity_kernel(prime, modes)
    if modes == 1:
        result[(prime-1,), (prime-1,)] = {-1: 1}
        if filtered:
            result[(0,), (0,)] = {0: 2}
    else:
        first, second = (0, 2), (1, 1)
        result[first, second] = {0: 1, -1: prime-1}
        result[second, first] = {0: 1, -1: prime-1}
        if filtered:
            result[first, second] = {-1: 1}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    full_cases = [(3, 3, 1, 0, 1, 1, False), (3, 3, 1, 0, 1, 3, True)]
    if args.extended:
        full_cases += [(3, 9, 1, 0, 2, 1, False), (3, 9, 1, 1, 1, 1, False)]
    full = []
    for prime, parent, deficit, initial, terminal, P, filtered in full_cases:
        F = example(initial+terminal, prime, filtered)
        endpoint = 'filtered' if filtered else 'identity'
        options = (prime, parent, deficit, initial, terminal, P)
        actual = zero_profiles(*options, F, endpoint)
        expected = full_zero_profiles(*options, F, endpoint)
        assert actual == expected
        full.append({'prime': prime, 'parent': parent, 'n': parent-deficit,
                      'deficit': deficit, 'initial_modes': initial, 'terminal_modes': terminal,
                      'ramification': P, 'endpoint': endpoint,
                      'full_active_matrix_order': prime**(2*(initial+terminal)+deficit),
                      'reduced_matrix_order': prime**(initial+terminal+deficit)})
    physical = []
    for P, J, filtered in ((1, 0, False), (3, -1, False), (1, 2, True)):
        F = example(1, 3, filtered)
        endpoint = 'filtered' if filtered else 'identity'
        tensor, clearing, gaussian = polynomial_seed(3, 3, 1, 0, 1, P, F, endpoint)
        dense = {powers: tuple(value.get(i, 0) for i in range(max(value)+1)) for powers, value in tensor.items()}
        actual = physical_cut_profiles(dense, 3, 2, P, J)
        data = profile_data(3, 3, 1, 0, 1, P, J, F, endpoint=endpoint)
        expected = (
            tuple(sorted(a+b+2*clearing for a in data['zero_active'] for b in expand_spectators(data['zero_spectators'], 3))),
            tuple(sorted(a+b-2*clearing for a in data['infinity_active'] for b in expand_spectators(data['infinity_spectators'], 3))),
        )
        assert actual == expected
        degree = max(power for value in tensor.values() for power in value)
        physical.append({'n': 2, 'ramification': P, 'charge': J, 'endpoint': endpoint,
                          'original_matrix_order': 81, 'clearing': clearing, 'polynomial_degree': degree,
                          'endpoint_coefficients_checked': len(set(tensor)|set(gaussian))})
    old_controls = 0
    for initial, terminal in ((0, 2), (1, 1)):
        F = example(2)
        selected = geometry_data(3, 3, 0, initial, terminal)['selected']
        geometry = 'selected:'+','.join(map(str, selected))
        assert zero_profiles(3, 3, 0, initial, terminal, 1, F) == existing_profiles(3, 3, 1, 2, F, geometry)
        old_controls += 1
    gaussian_controls, pressure_controls, pressure_error = 0, 0, 0.0
    for initial in (0, 1):
        options = (3, 27, 1, initial, 2-initial, 1)
        geometry = geometry_data(*options[:-1])
        F = example(2)
        blocks = zero_profiles(*options, F)
        infinity = infinity_profiles(*options, F)
        gaussian = zero_profiles(*options, identity_kernel(3, 2))
        expected = {}
        for first, second in product(states(3, 2), states(3, 3)):
            charge = sum(first)-sum(second)
            order = 2*(sum(j*digit for j, digit in zip(geometry['selected'], first))
                       +sum((i+1)*digit for i, digit in zip(geometry['retained'], second)))
            expected.setdefault(charge, []).append(order)
        assert gaussian == tuple({'charge': charge, 'orders': tuple(sorted(values))}
                                 for charge, values in sorted(expected.items()))
        gaussian_controls += 1
        evaluator = Pressure(*options, blocks, infinity)
        for J, scale in product((-3, 0, 2), (.03, .2, .7)):
            data = profile_data(*options, J, F, blocks)
            value = float(logsumexp(-scale*np.array(data['zero_active']))
                          +logsumexp(-scale*np.array(data['infinity_active'])))
            for order in data['zero_spectators']+data['infinity_spectators']:
                value += float(logsumexp(-scale*order*np.arange(3)))
            for order in data['paired_physical_first']+data['paired_physical_second']:
                value -= math.sqrt(2)*float(logsumexp(-scale*order*np.arange(3)/math.sqrt(2)))
            actual = float(evaluator.grid([scale], [J])[0])
            pressure_error = max(pressure_error, abs(actual-value))
            assert abs(actual-value) < 2e-10
            pressure_controls += 1
    F = identity_kernel(3, 1)
    F[(2,), (2,)] = {0: 1, 600: 2, 601: 1}
    late = infinity_profiles(3, 3, 1, 0, 1, 1, F)
    assert late == (-1204, -1, 0)
    assert len(rank_one_catalogue(26, 1)) == 1129
    counts = {size: len(projectors(size)) for size in (2, 3)}
    result = {'status': 'Exact prime-field deficit operator profiles with independent original tensors and full physical costs; no improved bound claim.',
              'full_rectangular_controls': full, 'original_physical_controls': physical,
              'independent_prime_power_specializations': old_controls,
              'n26_Gaussian_recovery_controls': gaussian_controls,
              'complete_profile_pressure_controls': pressure_controls,
              'maximum_pressure_replay_error': pressure_error,
              'degree_601_complete_source_infinity_profile': late,
              'constant_idempotent_counts': counts,
              'rank_one_projectors_by_occupation': {1: 12, 2: 117, 3: 12},
              'bounded_catalogue_kernel_count': 1129, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'full_rectangular_controls': len(full), 'original_cut_profiles': 2*len(physical),
                      'original_endpoint_coefficients': sum(row['endpoint_coefficients_checked'] for row in physical),
                      'n26_Gaussian_controls': gaussian_controls, 'pressure_controls': pressure_controls,
                      'maximum_pressure_error': pressure_error, 'catalogue_kernels': 1129}))


if __name__ == '__main__':
    main()
