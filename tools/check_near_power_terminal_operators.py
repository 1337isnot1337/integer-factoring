#!/usr/bin/env python3
"""Independent rectangular, dummy-Weyl, compiler, and original-cut controls."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from math import comb
from pathlib import Path
import random

import numpy as np
from scipy.special import logsumexp

from check_f4_identity_operator_filters import colored_kernel
from check_f4_repaired_identity_filters import expanded, physical_profiles
from compiled_f4_identity_operators import combine_kernel
from coupled_identity_mode_filters import add as add2, multiply as multiply2, shift
from coupled_identity_operator_filters import identity_kernel, states
from f4_repaired_identity_filters import add, MUL
from near_power_terminal_operators import (
    geometry_data, zero_profiles as full_profiles, infinity_profiles,
    profile_data, polynomial_seed, TerminalPressure,
)
from near_power_terminal_weyl import zero_profiles, compile_family, remove_dummy, weyl_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    direct_forms, interpolation, geometries = 0, 0, []
    for n, modes, initial in ((3, 1, 0), (3, 2, 0), (7, 2, 0), (7, 3, 0),
                              (15, 3, 0), (31, 3, 0), (63, 4, 0),
                              (3, 2, 1), (7, 3, 1), (31, 3, 1), (63, 4, 1)):
        data = geometry_data(n, modes, initial)
        for i in range(n):
            for k, j in enumerate(data['augmented']):
                direct = {}
                for h in range(min(i, j)+1):
                    for d in range(h+1):
                        if comb(j, h)*comb(i, h)*comb(h, d) % 2:
                            direct[d] = direct.get(d, 0)^1
                direct = {d: value for d, value in direct.items() if value}
                assert direct == data['forms'][i][k]
                direct_forms += 1
        for a, i in enumerate(data['discarded']):
            for j in range(modes+1):
                reconstructed = {}
                for b, t in enumerate(data['retained']):
                    reconstructed = add2(reconstructed, multiply2(
                        data['interpolation_numerators'][a][b], data['augmented_x'][b][j], 2), 2)
                assert reconstructed == multiply2(data['unit'], data['forms'][i][j], 2)
                interpolation += 1
        geometries.append({'n': n, 'initial_modes': initial, 'selected': data['selected'], 'retained': data['retained'],
                           'weighted_order_margins_nonnegative': True})
    rank_one_entries, rank_one_geometries = 0, 0
    for n in (3, 7, 15, 31, 63):
        for modes in range(1, min(4, n-1)+1):
            for initial in (0, 1):
                if initial >= modes:
                    continue
                data = weyl_data(n, modes, initial)
                selected = data['selected']
                c, d = {0: 1, 1: 1}, {0: 1, n+1: 1}
                e = add2(d, multiply2(c, {j: 1 for j in selected}, 2), 2)
                h = min(e)
                assert h == initial
                e0 = shift(e, -h, 2)
                delta = multiply2(d, e0, 2)
                for i in range(modes+1):
                    for j in range(modes+1):
                        numerator = add2(e0 if i == j else {}, shift(c, n+1-h+selected[i], 2), 2)
                        assert multiply2(numerator, data['denominator_unit'], 2) == multiply2(
                            data['inverse_parameter_numerator'][i][j], delta, 2)
                        rank_one_entries += 1
                rank_one_geometries += 1
    cases = [(3, 1, 1, 0), (3, 2, 1, 0), (7, 2, 3, 0), (7, 3, 1, 0),
             (3, 2, 1, 1), (7, 3, 1, 1)]
    if args.extended:
        cases += [(15, 3, 3, 0), (31, 3, 1, 0), (31, 3, 1, 1)]
    records, compiler_count, pressure_count = [], 0, 0
    for n, modes, P, initial in cases:
        selected = geometry_data(n, modes, initial)['selected']
        kernel = colored_kernel(modes, selected)
        compact = zero_profiles(n, P, modes, kernel, initial)
        assert compact == full_profiles(n, P, modes, kernel, initial)
        gaussian = zero_profiles(n, P, modes, identity_kernel(2, modes), initial)
        expected = {}
        for first, second in product(states(2, modes), states(2, modes+1)):
            charge = sum(first)-sum(second)
            order = 2*P*(sum(j*a for j, a in zip(selected, first))
                         +sum((i+1)*b for i, b in zip(geometry_data(n, modes, initial)['retained'], second)))
            expected.setdefault(charge, []).append(order)
        assert gaussian == tuple({'charge': charge, 'orders': tuple(sorted(values))}
                                 for charge, values in sorted(expected.items()))
        constant = identity_kernel(2, modes)
        difference = dict(kernel)
        for index, value in constant.items():
            difference[index] = add(difference.get(index, {}), value)
        perturbations = [difference]
        alphabet = states(2, modes)
        if modes > 1:
            sector = [state for state in alphabet if sum(state) == 1]
            perturbations.append({(sector[0], sector[-1]): {0: 1, -2: 1}})
        compiled = compile_family(n, P, modes, constant, perturbations, initial)
        for colors in product(range(4), repeat=len(perturbations)):
            F = combine_kernel(constant, perturbations, colors)
            try:
                actual = compiled(colors)
            except ValueError:
                continue
            assert actual == (zero_profiles(n, P, modes, F, initial), infinity_profiles(n, P, modes, F, initial))
            compiler_count += 1
        for J, scale in product((-3, 0, 2), (.03, .2, .7)):
            data = profile_data(n, P, J, modes, kernel, compact, initial)
            zero = np.array(data['zero_active'], dtype=float)/P
            infinity = np.array(data['infinity_active'], dtype=float)/P
            direct = float(logsumexp(-scale*zero)+logsumexp(-scale*infinity))
            direct += sum(math.log1p(math.exp(-scale*value/P)) for value in data['zero_spectators'])
            direct += sum(np.logaddexp(0, -scale*value/P) for value in data['infinity_spectators'])
            ratio = math.sqrt(2)
            direct -= ratio*sum(np.logaddexp(0, -scale*value/(P*ratio))
                                for value in data['paired_physical_first']+data['paired_physical_second'])
            actual = TerminalPressure(n, P, modes, compact, data['infinity_active'], initial)(math.log(scale), J)
            assert abs(direct-actual) < 2e-10
            pressure_count += 1
        records.append({'n': n, 'ramification': P, 'initial_modes': initial, 'modes': modes,
                        'rectangular_matrix_order': 2**(2*modes+1),
                        'Weyl_matrix_order': 2**(modes+1), 'zero_blocks': compact})
    physical = []
    for n, modes, P, J, initial in ((3, 1, 3, -1, 0), (3, 2, 1, 0, 0),
                                   (3, 2, 3, 2, 0), (3, 2, 1, -2, 1), (3, 2, 3, 1, 1)):
        kernel = colored_kernel(modes, geometry_data(n, modes, initial)['selected'])
        tensor, clearing, gaussian = polynomial_seed(n, P, modes, kernel, initial)
        data = profile_data(n, P, J, modes, kernel, initial_modes=initial)
        expected = (
            tuple(sorted(a+b+2*clearing for a in data['zero_active'] for b in expanded(data['zero_spectators']))),
            tuple(sorted(a+b-2*clearing for a in data['infinity_active'] for b in expanded(data['infinity_spectators']))),
        )
        actual, precision = physical_profiles(tensor, n, P, J)
        assert actual == expected
        physical.append({'n': n, 'ramification': P, 'initial_modes': initial, 'modes': modes, 'charge': J,
                          'original_matrix_order': 2**(2*n),
                          'endpoint_coefficients_checked': len(set(tensor)|set(gaussian)),
                          'maximum_native_precision': precision})
    # A high-degree second-plane coefficient must survive complete compilation.
    constant = identity_kernel(2, 2)
    source = {((1, 0), (0, 1)): {0: 1, 600: 1}}
    compiled = compile_family(7, 1, 2, constant, (source,))
    F = combine_kernel(constant, (source,), (2,))
    assert compiled((2,)) == (zero_profiles(7, 1, 2, F), infinity_profiles(7, 1, 2, F))
    rejected = 0
    try:
        remove_dummy({-2: [0], -1: [0], 0: [0], 1: [0]}, 6, 1)
    except ArithmeticError:
        rejected += 1
    for n, modes in ((4, 2), (7, 7), (3, 0)):
        try:
            geometry_data(n, modes)
        except ValueError:
            rejected += 1
    assert rejected == 4
    result = {'status': 'Exact near-power terminal reduction, dummy Weyl division, complete compiler, and original physical cuts; no improved bound claim.',
              'direct_Pascal_reversal_identities': direct_forms,
              'exact_interpolation_identities': interpolation, 'weighted_geometries': geometries,
              'closed_rank_one_parameter_geometries': rank_one_geometries,
              'closed_rank_one_parameter_entry_identities': rank_one_entries,
              'full_rectangular_vs_reduced': records,
              'compiled_two_cut_controls': compiler_count+1,
              'complete_profile_physical_partition_pressure_controls': pressure_count,
              'original_two_cut_cases': physical,
              'degree_600_second_plane_control': True,
              'invalid_geometry_or_dummy_rejections': rejected, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'rectangular_vs_Weyl_profiles': len(records),
                      'compiled_two_cut_controls': compiler_count+1,
                      'original_cut_profiles': 2*len(physical),
                      'endpoint_coefficients_checked': sum(row['endpoint_coefficients_checked'] for row in physical),
                      'pressure_controls': pressure_count,
                      'reversal_identities': direct_forms, 'interpolation_identities': interpolation,
                      'rank_one_parameter_geometries': rank_one_geometries,
                      'rank_one_parameter_entries': rank_one_entries}))


if __name__ == '__main__':
    main()
