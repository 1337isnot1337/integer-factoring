#!/usr/bin/env python3
"""Independent exact controls for the n=2^k+1 terminal operator core."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from above_power_terminal_operators import (
    geometry_data, weyl_data, zero_profiles, infinity_profiles,
    compile_family, profile_data, polynomial_seed, AbovePowerPressure,
)
from check_f4_identity_operator_filters import colored_kernel
from check_f4_repaired_identity_filters import expanded, physical_profiles
from compiled_f4_identity_operators import combine_kernel
from coupled_identity_mode_filters import add as add2, multiply as multiply2, shift, algebra_multiply
from coupled_identity_operator_filters import _components, identity_kernel, states
from coupled_identity_weyl_profiles import determinant
from f4_repaired_identity_filters import add, multiply, laurent_orders


def literal_active_blocks(n, P, modes, kernel):
    """Expand the complete square exterior tensor before the Weyl transform."""
    data = geometry_data(n, modes)
    size = modes+1
    zero = (0,)*(4*size)
    gaussian = {zero: {0: 1}}
    for i, j in enumerate(data['augmented']):
        powers = list(zero)
        powers[size+i] = powers[2*size+i] = 1
        gaussian = algebra_multiply(gaussian,
            {zero: {0: 1}, tuple(powers): {2*P*j: 1}}, 2)
    for i, j in enumerate(data['retained']):
        powers = list(zero)
        powers[i] = powers[3*size+i] = 1
        gaussian = algebra_multiply(gaussian,
            {zero: {0: 1}, tuple(powers): {2*P*(j+1): 1}}, 2)
    left, right = [], []
    for a, j in enumerate(data['selected']):
        powers = list(zero)
        powers[size+a] = 1
        first = {tuple(powers): {2*P*j: 1}}
        powers = list(zero)
        powers[2*size+a] = 1
        second = {tuple(powers): {0: 1}}
        for i in range(size):
            powers = list(zero)
            powers[i] = 1
            value = multiply2({0: 1, 1: 1}, data['augmented_v'][i][a], 2)
            if value:
                first[tuple(powers)] = {2*P*degree: coefficient for degree, coefficient in value.items()}
            powers = list(zero)
            powers[3*size+i] = 1
            value = data['augmented_x'][i][a]
            if value:
                second[tuple(powers)] = {2*P*(degree-j): coefficient for degree, coefficient in value.items()}
        left.append(first)
        right.append(second)
    components = _components(gaussian, left, right, 2)
    tensor = {}
    for index, scalar in kernel.items():
        for powers, polynomial in components[index].items():
            tensor[powers] = add(tensor.get(powers, {}), multiply(polynomial, {2*d: v for d, v in scalar.items()}))
    tensor = {powers: polynomial for powers, polynomial in tensor.items() if polynomial}
    result = []
    for charge in range(-size, size+1):
        labels = [(y, x) for y, x in product(states(2, size), repeat=2) if sum(y)-sum(x) == charge]
        matrix = [[tensor.get(v+y+u+x, {}) for u, v in labels] for y, x in labels]
        result.append({'charge': charge, 'orders': laurent_orders(matrix)[0]})
    return tuple(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    reflected, gram, regularity, gaussian_cases = 0, 0, 0, []
    for n, modes in ((3, 1), (3, 2), (5, 1), (5, 3), (9, 2),
                      (17, 3), (33, 3), (65, 4)):
        data = geometry_data(n, modes)
        N = n-1
        old_unit = {degree: 1 for degree in range(N)}
        for i in range(N):
            for j in range(N):
                # B H = c^(N−1) R, checked as complete polynomials.
                left = {}
                right = {}
                for k in range(N):
                    left = add2(left, { (i & (N-1-k))+(k & j): 1}, 2)
                    right = add2(right, { (k & i)+(k & (N-1-j))+j-k: 1}, 2)
                expected = old_unit if i+j == N-1 else {}
                assert left == expected == right
                reflected += 2
                # Both B and D^-1 B D are regular; its inverse is B / unit.
                assert i & (N-1-j) >= 0
                assert (i & (N-1-j))+j-i >= 0
                regularity += 2
        for a, j in enumerate(data['augmented']):
            for b, k in enumerate(data['augmented']):
                original, reduced = {}, {}
                for i in range(n):
                    original = add2(original, {(i & j)+(i & k)-i-1: 1}, 2)
                for h, i in enumerate(data['retained']):
                    reduced = add2(reduced, shift(multiply2(
                        data['augmented_v'][h][a], data['augmented_x'][h][b], 2), -i-1, 2), 2)
                assert original == reduced
                gram += 1
        parameter = weyl_data(n, modes)
        for i, row in enumerate(parameter['inverse_parameter_numerator']):
            for j, value in enumerate(row):
                if value:
                    assert min(value) >= max(0, data['augmented'][i]-data['augmented'][j])
                regularity += 1
        gaussian = zero_profiles(n, 1, modes, identity_kernel(2, modes))
        expected = {}
        for y, x in product(states(2, modes+1), repeat=2):
            charge = sum(y)-sum(x)
            order = 2*(sum(j*a for j, a in zip(data['augmented'], y))
                       +sum((j+1)*a for j, a in zip(data['retained'], x)))
            expected.setdefault(charge, []).append(order)
        assert gaussian == tuple({'charge': q, 'orders': tuple(sorted(values))}
                                  for q, values in sorted(expected.items()))
        gaussian_cases.append({'n': n, 'modes': modes, 'selected': data['selected'],
                               'retained': data['retained'], 'active_zero_modes': 2*(modes+1)})
    active_cases = [(3, 1, 1), (3, 2, 1), (5, 2, 3), (9, 3, 1)]
    if args.extended:
        active_cases += [(17, 3, 3), (33, 3, 1)]
    active, compiled_count, pressure_count = [], 0, 0
    for n, modes, P in active_cases:
        selected = geometry_data(n, modes)['selected']
        kernel = colored_kernel(modes, selected)
        compact = zero_profiles(n, P, modes, kernel)
        assert compact == literal_active_blocks(n, P, modes, kernel)
        constant = identity_kernel(2, modes)
        difference = dict(kernel)
        for index, value in constant.items():
            difference[index] = add(difference.get(index, {}), value)
        perturbations = [difference]
        if modes > 1:
            sector = [state for state in states(2, modes) if sum(state) == 1]
            perturbations.append({(sector[0], sector[-1]): {0: 1, -2: 1}})
        compiled = compile_family(n, P, modes, constant, perturbations)
        for colors in product(range(4), repeat=len(perturbations)):
            F = combine_kernel(constant, perturbations, colors)
            try:
                actual = compiled(colors)
            except ValueError:
                continue
            assert actual == (zero_profiles(n, P, modes, F), infinity_profiles(n, P, modes, F))
            compiled_count += 1
        for J, scale in product((-3, 0, 2), (.03, .2, .7)):
            data = profile_data(n, P, J, modes, kernel, compact)
            direct = float(logsumexp(-scale*np.array(data['zero_active'])/P)
                           +logsumexp(-scale*np.array(data['infinity_active'])/P))
            direct += sum(np.logaddexp(0, -scale*value/P) for value in data['zero_spectators'])
            direct += sum(np.logaddexp(0, -scale*value/P) for value in data['infinity_spectators'])
            direct -= math.sqrt(2)*sum(np.logaddexp(0, -scale*value/(P*math.sqrt(2)))
                                      for value in data['paired_physical_first']+data['paired_physical_second'])
            actual = AbovePowerPressure(n, P, modes, compact, data['infinity_active'])(math.log(scale), J)
            assert abs(actual-direct) < 2e-10, (n, modes, P, J, actual, direct)
            pressure_count += 1
        active.append({'n': n, 'ramification': P, 'modes': modes,
                       'square_matrix_order': 2**(2*(modes+1)), 'Weyl_matrix_order': 2**(modes+1),
                       'zero_blocks': compact})
    physical = []
    for n, modes, P, J in ((3, 1, 1, 0), (3, 1, 3, -1),
                           (3, 2, 1, 0), (3, 2, 3, -2)):
        kernel = colored_kernel(modes, geometry_data(n, modes)['selected'])
        tensor, clearing, gaussian = polynomial_seed(n, P, modes, kernel)
        data = profile_data(n, P, J, modes, kernel)
        expected = (
            tuple(sorted(a+b+2*clearing for a in data['zero_active'] for b in expanded(data['zero_spectators']))),
            tuple(sorted(a+b-2*clearing for a in data['infinity_active'] for b in expanded(data['infinity_spectators']))))
        actual, precision = physical_profiles(tensor, n, P, J)
        assert actual == expected
        physical.append({'n': n, 'modes': modes, 'ramification': P, 'charge': J,
                         'matrix_order': 2**(2*n), 'endpoint_coefficients': len(set(tensor)|set(gaussian)),
                         'maximum_native_precision': precision})
    constant = identity_kernel(2, 2)
    source = {((1, 0), (0, 1)): {0: 1, 600: 1}}
    compiled = compile_family(9, 1, 2, constant, (source,))
    F = combine_kernel(constant, (source,), (2,))
    assert compiled((2,)) == (zero_profiles(9, 1, 2, F), infinity_profiles(9, 1, 2, F))
    rejected = 0
    for n, modes in ((4, 2), (3, 3), (3, 0)):
        try:
            geometry_data(n, modes)
        except ValueError:
            rejected += 1
    assert rejected == 3
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/above_power_terminal_operators.py', 'tools/check_above_power_terminal_operators.py']
    result = {'status': 'Exact above-power terminal core and complete Weyl profiles; no improved-bound claim.',
              'reflected_Pascal_identities': reflected, 'full_vs_reduced_Gram_identities': gram,
              'weighted_regularity_entries': regularity, 'gaussian_recovery_cases': gaussian_cases,
              'literal_active_vs_Weyl_cases': active, 'complete_compiler_cases': compiled_count+1,
              'full_pressure_cases': pressure_count, 'original_physical_cases': physical,
              'degree_600_second_plane_control': True, 'invalid_geometries_rejected': rejected,
              'extended': args.extended,
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in sources}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'active_vs_Weyl': len(active), 'complete_compiler': compiled_count+1,
                      'original_cut_profiles': 2*len(physical),
                      'endpoint_coefficients': sum(row['endpoint_coefficients'] for row in physical),
                      'pressure': pressure_count, 'Pascal': reflected, 'Gram': gram,
                      'weighted_regularity': regularity}))


if __name__ == '__main__':
    main()
