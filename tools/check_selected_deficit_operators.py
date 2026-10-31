#!/usr/bin/env python3
"""Original physical and rectangular controls for arbitrary one-deficit selections."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from coupled_identity_operator_filters import identity_kernel, states
from check_prime_state_filter_profiles import physical_cut_profiles, expand_spectators
from prime_deficit_identity_operators import (
    operator_geometry, zero_profiles, full_zero_profiles, infinity_profiles,
    profile_data, polynomial_seed, Pressure,
)


def example(prime, modes, filtered=False):
    kernel = identity_kernel(prime, modes)
    if modes == 1:
        kernel[(prime-1,), (prime-1,)] = {-1: 1}
        if filtered:
            kernel[(0,), (0,)] = {0: prime-1}
    else:
        a, b = (0, 1), (1, 0)
        kernel[a, b] = {-1: 1} if filtered else {0: 1, -1: prime-1}
        kernel[b, a] = {0: 1, -1: prime-1}
    return kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    physical = []
    for prime, parent, selected, P, J, filtered in (
        (2, 4, (0, 1), 1, 0, False),
        (2, 4, (0, 2), 3, -1, False),
        (2, 4, (1, 2), 1, 2, True),
        (3, 3, (0,), 1, -2, True),
        (3, 3, (1,), 3, 1, False),
    ):
        options = (prime, parent, 1, 0, len(selected), P)
        endpoint = 'filtered' if filtered else 'identity'
        kernel = example(prime, len(selected), filtered)
        tensor, clearing, gaussian = polynomial_seed(*options, kernel, endpoint, selected_indices=selected)
        dense = {powers: tuple(value.get(i, 0) for i in range(max(value)+1))
                 for powers, value in tensor.items()}
        actual = physical_cut_profiles(dense, prime, parent-1, P, J)
        data = profile_data(*options, J, kernel, endpoint=endpoint, selected_indices=selected)
        expected = (
            tuple(sorted(a+b+2*clearing for a in data['zero_active']
                         for b in expand_spectators(data['zero_spectators'], prime))),
            tuple(sorted(a+b-2*clearing for a in data['infinity_active']
                         for b in expand_spectators(data['infinity_spectators'], prime))),
        )
        assert actual == expected
        physical.append({'prime': prime, 'parent': parent, 'selected': selected,
                         'ramification': P, 'charge': J, 'endpoint': endpoint,
                         'original_matrix_order': prime**(2*(parent-1)),
                         'endpoint_coefficients_checked': len(set(tensor)|set(gaussian)),
                         'polynomial_degree': max(degree for value in tensor.values() for degree in value)})
    rectangular = []
    cases = [(2, 8, (1, 4), 1), (2, 8, (0, 1), 3)]
    if args.extended:
        cases.append((3, 9, (2, 5), 1))
    for prime, parent, selected, P in cases:
        options = (prime, parent, 1, 0, len(selected), P)
        kernel = example(prime, len(selected))
        assert zero_profiles(*options, kernel, selected_indices=selected) == full_zero_profiles(
            *options, kernel, selected_indices=selected)
        rectangular.append({'prime': prime, 'parent': parent, 'selected': selected,
                            'ramification': P, 'full_active_matrix_order': prime**(2*len(selected)+1),
                            'reduced_matrix_order': prime**(len(selected)+1)})
    compatibility = []
    for prime, parent, initial, terminal, P in (
        (2, 4, 2, 0, 1), (2, 8, 0, 2, 3), (2, 8, 1, 1, 1),
        (3, 9, 0, 2, 1), (3, 9, 1, 1, 3),
    ):
        options = (prime, parent, 1, initial, terminal, P)
        selected = operator_geometry(*options[:-1])['selected']
        kernel = example(prime, initial+terminal)
        assert zero_profiles(*options, kernel) == zero_profiles(*options, kernel, selected_indices=selected)
        assert infinity_profiles(*options, kernel) == infinity_profiles(*options, kernel, selected_indices=selected)
        compatibility.append({'prime': prime, 'parent': parent, 'selected': selected, 'ramification': P})
    pressure_controls, pressure_error, gaussian_controls = 0, 0.0, 0
    for selected in ((0, 1), (0, 13), (7, 18), (24, 25)):
        options = (3, 27, 1, 0, 2, 1)
        geometry = operator_geometry(*options[:-1], selected_indices=selected)
        gaussian = zero_profiles(*options, identity_kernel(3, 2), selected_indices=selected)
        expected = {}
        for first, second in product(states(3, 2), states(3, 3)):
            charge = sum(first)-sum(second)
            order = 2*(sum(j*digit for j, digit in zip(selected, first))
                       +sum((i+1)*digit for i, digit in zip(geometry['retained'], second)))
            expected.setdefault(charge, []).append(order)
        assert gaussian == tuple({'charge': charge, 'orders': tuple(sorted(values))}
                                 for charge, values in sorted(expected.items()))
        gaussian_controls += 1
        kernel = example(3, 2)
        blocks = zero_profiles(*options, kernel, selected_indices=selected)
        infinity = infinity_profiles(*options, kernel, selected_indices=selected)
        pressure = Pressure(*options, blocks, infinity, selected_indices=selected)
        for J, scale in product((-3, 0, 2), (.03, .2, .7)):
            data = profile_data(*options, J, kernel, blocks, selected_indices=selected)
            value = float(logsumexp(-scale*np.array(data['zero_active']))
                          +logsumexp(-scale*np.array(data['infinity_active'])))
            for order in data['zero_spectators']+data['infinity_spectators']:
                value += float(logsumexp(-scale*order*np.arange(3)))
            for order in data['paired_physical_first']+data['paired_physical_second']:
                value -= math.sqrt(2)*float(logsumexp(-scale*order*np.arange(3)/math.sqrt(2)))
            actual = float(pressure.grid([scale], [J])[0])
            pressure_error = max(pressure_error, abs(actual-value))
            assert abs(actual-value) < 2e-10
            pressure_controls += 1
    rejections = 0
    for deficit, selected in ((2, (0, 1)), (1, (1, 0)), (1, (1, 1)), (1, (0,)), (1, (0, 7))):
        try:
            operator_geometry(2, 8, deficit, 0, 2, selected)
        except ValueError:
            rejections += 1
        else:
            raise AssertionError('An invalid selected geometry was accepted.')
    result = {'status': 'Arbitrary one-deficit selections with independent full physical and rectangular controls; no improved bound claim.',
              'original_physical_controls': physical, 'full_rectangular_controls': rectangular,
              'backward_compatible_end_selections': compatibility, 'n26_Gaussian_recovery_controls': gaussian_controls,
              'complete_profile_pressure_controls': pressure_controls, 'maximum_pressure_replay_error': pressure_error,
              'invalid_selection_rejections': rejections, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'original_cut_profiles': 2*len(physical),
                      'original_endpoint_coefficients': sum(row['endpoint_coefficients_checked'] for row in physical),
                      'rectangular_controls': len(rectangular), 'compatibility_controls': len(compatibility),
                      'n26_Gaussian_controls': gaussian_controls, 'pressure_controls': pressure_controls,
                      'maximum_pressure_error': pressure_error, 'invalid_selection_rejections': rejections}))


if __name__ == '__main__':
    main()
