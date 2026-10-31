#!/usr/bin/env python3
"""Exact symbolic identities and independent profiles for the Weyl reduction."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from math import comb, factorial
from pathlib import Path
import random
import time

from coupled_identity_mode_filters import add, multiply
from coupled_identity_operator_filters import identity_kernel, states, zero_profiles as full_profiles
from coupled_identity_weyl_profiles import zero_profiles, weyl_data
from check_coupled_identity_operator_filters import example_kernel


def matrix_multiply(first, second, prime):
    rows = {}
    for (i, j), value in second.items():
        rows.setdefault(i, []).append((j, value))
    result = {}
    for (i, k), value in first.items():
        for j, other in rows.get(k, ()):
            entry = i, j
            result[entry] = add(result.get(entry, {}), multiply(value, other, prime), prime)
    return {entry: value for entry, value in result.items() if value}


def one_mode_kernel(prime, a, b):
    alphabet = tuple(product(range(prime), repeat=2))
    result = {}
    for row, (i, j) in enumerate(alphabet):
        for column, (k, ell) in enumerate(alphabet):
            tail = ell - a + i
            if i <= a and k <= b and 0 <= tail < prime and j == b - k + tail:
                coefficient = (comb(a, i) * comb(b, k) * (-1)**(b-k)
                               * pow(factorial(b) * factorial(tail), -1, prime)) % prime
                if coefficient:
                    result[row, column] = {a-i: coefficient}
    return result


def one_mode_weyl_unit(prime, sign):
    alphabet = tuple(product(range(prime), repeat=2))
    lookup = {state: i for i, state in enumerate(alphabet)}
    result = {}
    for column, (i, j) in enumerate(alphabet):
        for digit in range(min(i, j) + 1):
            coefficient = ((factorial(i) // factorial(i-digit))
                           * (factorial(j) // factorial(j-digit))
                           * pow(factorial(digit), -1, prime) * sign**digit) % prime
            row = lookup[i-digit, j-digit]
            assert (i-digit) - (j-digit) == i-j
            result[row, column] = {digit: coefficient}
    return result


def one_mode_reduced_lift(prime, a, b):
    alphabet = tuple(product(range(prime), repeat=2))
    lookup = {state: i for i, state in enumerate(alphabet)}
    derivative, power = prime-1-a, prime-1-b
    result = {}
    for column, (i, j) in enumerate(alphabet):
        degree = i + power
        if degree < prime and degree >= derivative:
            coefficient = (-factorial(a) * (-1)**power
                           * (factorial(degree) // factorial(degree-derivative))) % prime
            row = lookup[degree-derivative, j]
            result[row, column] = {-power: coefficient}
    return result


def symbolic_controls(primes=(2, 3, 5, 7)):
    records = []
    for prime in primes:
        start = time.monotonic()
        unit = one_mode_weyl_unit(prime, 1)
        inverse = one_mode_weyl_unit(prime, -1)
        identity = {(i, i): {0: 1} for i in range(prime**2)}
        assert matrix_multiply(unit, inverse, prime) == identity
        reference = one_mode_kernel(prime, prime-1, prime-1)
        right = matrix_multiply(inverse, reference, prime)
        checked = diagonal = 0
        for a, b in product(range(prime), repeat=2):
            lifted = one_mode_reduced_lift(prime, a, b)
            actual = matrix_multiply(unit, matrix_multiply(lifted, right, prime), prime)
            expected = one_mode_kernel(prime, a, b)
            assert actual == expected, (prime, a, b, actual, expected)
            checked += 1
            if a == b:
                h = prime-1-a
                for i in range(prime):
                    k = prime-1-i
                    expected_value = { -h: (-1)**h * comb(k, h) % prime } if h <= k else {}
                    expected_value = {degree: value for degree, value in expected_value.items() if value}
                    assert lifted.get((i*prime, i*prime), {}) == expected_value
                    diagonal += 1
        records.append({'prime': prime, 'symbolic_matrix_identities': checked,
                        'diagonal_transform_coefficients': diagonal,
                        'coefficient_ring': 'F_p[c,c^-1]',
                        'matrix_order': prime**2, 'seconds': time.monotonic()-start})
    return records


def random_kernel(prime, modes, seed, selected):
    rng = random.Random(seed)
    kernel = identity_kernel(prime, modes)
    alphabet = states(prime, modes)
    for row in alphabet:
        for column in alphabet:
            if sum(row) != sum(column):
                continue
            if rng.randrange(3) == 0:
                degree = rng.choice((-2, -1, 0, 1, 2))
                scalar = rng.randrange(1, prime)
                value = add(kernel.get((row, column), {}),
                            {degree: scalar, degree+1: -scalar % prime}, prime)
                if value:
                    kernel[row, column] = value
                else:
                    kernel.pop((row, column), None)
    # The random corrections vanish at q=1, so this kernel needs no repair.
    return kernel


def profile_control(prime, n, ramification, modes, geometry, seed):
    start = time.monotonic()
    data = weyl_data(prime, n, modes, geometry)
    if seed < 0:
        kernel = example_kernel(prime, modes, True, 'filtered')
        endpoint = 'filtered'
    else:
        kernel = random_kernel(prime, modes, seed, data['selected'])
        endpoint = 'identity'
    reduced = zero_profiles(prime, n, ramification, modes, kernel, geometry, endpoint)
    full = full_profiles(prime, n, ramification, modes, kernel, geometry, endpoint)
    assert reduced == full, (prime, n, ramification, modes, geometry, seed, reduced, full)
    assert sum(len(block['orders']) for block in reduced) == prime**(2*modes)
    return {'prime': prime, 'n': n, 'ramification': ramification, 'modes': modes,
            'geometry': geometry, 'seed': seed, 'endpoint_policy': endpoint,
            'full_matrix_order': prime**(2*modes), 'reduced_matrix_order': prime**modes,
            'zero_unit_constant': data['denominator_unit'][0],
            'charge_profiles': reduced, 'seconds': time.monotonic()-start}


def replay_screen(path):
    from search_coupled_identity_operator_filters import decode_kernel, candidates
    record = json.loads(path.read_text())
    if 'runs' in record:
        count = 0
        for run in record['runs']:
            selected = weyl_data(2, run['n'], run['modes'], run['geometry'])['selected']
            catalogue = candidates(run['modes'], run['ramification'], selected, run['family'])
            assert len(catalogue) == run['candidate_count']
            for row in run['profiles']:
                kernel = catalogue[row['candidate_index']]
                actual = zero_profiles(2, run['n'], run['ramification'], run['modes'],
                                       kernel, run['geometry'], run['endpoint_policy'])
                expected = tuple({'charge': block['charge'], 'orders': tuple(block['orders'])}
                                 for block in row['zero_blocks'])
                assert actual == expected
                count += 1
        return {'source': str(path), 'complete_charge_profile_comparisons': count,
                'run_count': len(record['runs'])}
    count = 0
    for row in record['rows']:
        kernel = decode_kernel(row['kernel'])
        actual = zero_profiles(2, record['n'], record['ramification'], record['modes'],
                               kernel, record['geometry'], record['endpoint_policy'])
        expected = tuple({'charge': block['charge'], 'orders': tuple(block['orders'])}
                         for block in row['zero_blocks'])
        assert actual == expected
        count += 1
    return {'source': str(path), 'n': record['n'], 'modes': record['modes'],
            'geometry': record['geometry'], 'complete_charge_profile_comparisons': count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--screen', type=Path, action='append', default=[])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    symbolic = symbolic_controls()
    cases = [
        (2, 2, 1, 1, 'initial', 1), (3, 2, 3, 1, 'initial', 2),
        (5, 2, 1, 1, 'initial', 3), (2, 2, 1, 2, 'initial', 4),
        (2, 3, 3, 2, 'initial', 5), (2, 3, 1, 3, 'initial', 6),
        (3, 2, 1, 2, 'initial', 7), (2, 4, 3, 2, 'mixed', 8),
        (2, 4, 1, 2, 'terminal', 9), (2, 4, 1, 3, 'terminal', 10),
        (2, 3, 1, 3, 'initial', -1), (3, 2, 1, 2, 'initial', -1),
    ]
    if args.extended:
        cases += [(5, 5, 1, 1, 'terminal', 11),
                  (3, 3, 1, 2, 'mixed', 12),
                  (2, 16, 3, 3, 'initial', 13),
                  (2, 16, 1, 3, 'terminal', 14)]
    profiles = [profile_control(*case) for case in cases]
    screens = [replay_screen(path) for path in args.screen]
    result = {'status': 'Exact Weyl identities and independent complete-profile controls; no pressure claim.',
              'symbolic_controls': symbolic, 'profile_controls': profiles,
              'complete_screen_replays': screens, 'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'symbolic_matrix_identities': sum(row['symbolic_matrix_identities'] for row in symbolic),
                      'diagonal_transform_coefficients': sum(row['diagonal_transform_coefficients'] for row in symbolic),
                      'independent_full_profile_comparisons': len(profiles),
                      'full_zero_orders_compared': sum(row['full_matrix_order'] for row in profiles),
                      'screen_profile_comparisons': sum(row['complete_charge_profile_comparisons'] for row in screens)}))


if __name__ == '__main__':
    main()
