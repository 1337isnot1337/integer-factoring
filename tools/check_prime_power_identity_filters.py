#!/usr/bin/env python3
"""Independent reflected-Pascal, joint-profile and physical-cut controls."""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import comb
from pathlib import Path

from coupled_identity_mode_filters import add, multiply, shift, polynomial_seed
from check_coupled_identity_extreme_filters import raw_form
from check_coupled_identity_mode_filters import direct_profiles, expand_modes
from prime_power_identity_filters import beta_polynomials, profile_data


def gram_controls():
    rows = []
    identities = weighted_entries = 0
    for prime, n in ((2, 1), (2, 2), (2, 4), (2, 8), (2, 16),
                     (3, 3), (3, 9), (5, 5), (7, 7)):
        h = tuple(tuple(raw_form(i, j, prime) for j in range(n)) for i in range(n))
        unit = {power: comb(n - 1, power) * (-1) ** power % prime
                for power in range(n) if comb(n - 1, power) % prime}
        for j, k in product(range(n), repeat=2):
            value = {}
            for i in range(n):
                value = add(value, shift(multiply(h[i][j], h[i][k], prime), -i, prime), prime)
            expected = shift(unit, j + 1 - n, prime) if j == k else {}
            assert value == expected
            identities += 1
        for i, k in product(range(n), repeat=2):
            value = h[i][n - 1 - k]
            assert value and min(value) >= max(0, i - k)
            # Both B and its weighted companion have invertible constant
            # upper-triangular residues; the diagonal entries are nonzero.
            if i > k:
                assert value.get(0, 0) == 0
            if i == k:
                assert value.get(0, 0) != 0
            weighted_entries += 1
        rows.append({'prime': prime, 'n': n, 'Gram_dimension': n,
                     'minimum_weighted_entry_order': min(min(h[i][n - 1 - k]) + k - i
                                                        for i, k in product(range(n), repeat=2))})
    return {'full_polynomial_Gram_entries': identities,
            'weighted_regular_change_entries': weighted_entries, 'cases': rows}


def old_core_controls():
    from f4_repaired_identity_filters import reduced_zero_profiles as f4_profiles
    from coupled_identity_extreme_filters import zero_profiles as prime_profiles
    records = []
    for prime, n, P, geometry, selected, field in (
            (2, 4, 1, 'initial', (0, 1), 4),
            (2, 4, 3, 'mixed', (0, 3), 4),
            (2, 4, 1, 'terminal', (1, 2, 3), 4),
            (2, 8, 1, 'split:1:2', (0, 6, 7), 4),
            (3, 3, 1, 'terminal', (1, 2), 3)):
        modes = len(selected)
        states = tuple(product(range(prime), repeat=modes))
        for trial in range(3):
            alpha = {state: {(index * index + trial * index + 1) % 5 - 2:
                             (index + trial) % (field - 1) + 1}
                     for index, state in enumerate(states)}
            formula = profile_data(prime, n, P, 1, selected, alpha, field)
            if prime == 2:
                actual = f4_profiles(n, P, modes, alpha, geometry)
            else:
                actual = prime_profiles(prime, n, P, modes, alpha, geometry, strict_endpoint=False)
            assert tuple((b['charge'], tuple(b['orders'])) for b in actual) == tuple(
                (b['charge'], tuple(b['orders'])) for b in formula['zero_charge_blocks_at_J1'])
            records.append({'prime': prime, 'n': n, 'ramification': P,
                            'geometry': geometry, 'selected': selected, 'coefficient_field': field,
                            'trial': trial})
    return records


def original_controls(include_qutrit=False):
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel
    records = []
    endpoint_equalities = 0
    cases = [(2, 2, 1, (0, 1)), (2, 4, 1, (1, 2)),
             (2, 4, 1, (0, 2)), (2, 4, 1, (0, 1, 3)),
             (2, 4, 3, (0, 2, 3)), (2, 4, 1, (0, 1, 2, 3))]
    if include_qutrit:
        cases.append((3, 3, 1, (0, 2)))
    for prime, n, P, selected in cases:
        modes = len(selected)
        states = tuple(product(range(prime), repeat=modes))
        alpha = {state: {(index * index + index) % 3 - 1: 1}
                 for index, state in enumerate(states)}
        # A truly joint polynomial weight, beyond a product of mode weights.
        if prime == 2:
            alpha[states[len(states) // 2 - 1]] = {-1: 1, 0: 1, 1: 1}
        else:
            alpha[states[len(states) // 2 - 1]] = {0: 1, 1: 1}
            alpha[states[-1]] = {0: 2}
        tensor, clearing = polynomial_seed(prime, n, P, modes, alpha, selected,
                                           strict_endpoint=False, repair_endpoint=True)
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, prime), prime)
        evaluate = lambda value: sum(coefficient * (-1) ** degree
                                     for degree, coefficient in enumerate(value)) % prime
        for powers in product(range(prime), repeat=4 * n):
            assert evaluate(tensor.get(powers, ())) == evaluate(gaussian.get(powers, ()))
            endpoint_equalities += 1
        for charge in (0, 2):
            actual, precision = direct_profiles(prime, tensor, n, P, charge)
            formula = profile_data(prime, n, P, charge, selected, alpha)
            expected = []
            for cut, scalar in (('zero', 2 * clearing), ('infinity', -2 * clearing)):
                spectators = expand_modes(prime, formula[cut + '_spectators'])
                expected.append(tuple(sorted(a + b + scalar for a in formula[cut + '_active']
                                             for b in spectators)))
            assert actual == expected, (prime, n, P, selected, charge)
            records.append({'prime': prime, 'n': n, 'ramification': P, 'selected': selected,
                            'charge': charge, 'full_matrix_order': prime ** (2 * n),
                            'common_Laurent_clearing': clearing, 'maximum_jet_precision': precision})
        print('PHYSICAL_PASS', prime, n, P, selected, flush=True)
    return {'physical_endpoint_coefficient_equalities': endpoint_equalities,
            'complete_original_both_cut_comparisons': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-qutrit', action='store_true')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1]
                        / 'evidence' / 'prime-power-identity-filters.json')
    args = parser.parse_args()
    result = {'status': 'Exact all-index prime-power profile controls passed; no exponent improvement claimed.',
              'reflected_Pascal': gram_controls(), 'old_active_core_comparisons': old_core_controls(),
              'original_physical_tensors': original_controls(args.include_qutrit)}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS', args.output)


if __name__ == '__main__':
    main()
