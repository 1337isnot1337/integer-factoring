#!/usr/bin/env python3
"""Original-tensor checks of the endpoint-identity operator kernel family."""
from __future__ import annotations
import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import time

from coupled_identity_operator_filters import (
    identity_kernel, normalized_kernel, polynomial_seed, profile_data,
    states, active_tensor, original_components, _laurent_orders,
)
from coupled_identity_mode_filters import active_tensor as diagonal_tensor, monomial_filters
from check_prime_state_filter_profiles import physical_cut_profiles, expand_spectators


def entry_list(kernel):
    return [{'row': row, 'column': column, 'polynomial': value}
            for (row, column), value in sorted(kernel.items())]


def example_kernel(prime, modes, variant, endpoint='identity'):
    kernel = identity_kernel(prime, modes)
    first = (0,) * (modes - 2) + ((0, 1) if prime == 2 else (0, 2))
    second = (0,) * (modes - 2) + ((1, 0) if prime == 2 else (1, 1))
    kernel[first, second] = {0: 1, -1: -1 % prime}
    kernel[second, first] = {0: 1, -1: -1 % prime}
    if variant:
        top = (prime - 1,) * modes
        kernel[top, top] = {-1: 1}
    if endpoint == 'filtered':
        kernel[first, second] = {-1: 1}
        if prime > 2:
            kernel[first, first] = {0: 2}
    return kernel


def physical_control(prime, n, ramification, charge, modes, geometry, variant, endpoint='identity'):
    start = time.monotonic()
    kernel = example_kernel(prime, modes, variant, endpoint)
    tensor, clearing = polynomial_seed(prime, n, ramification, modes, kernel, geometry, endpoint)
    _, gaussian = original_components(prime, n, modes, geometry)
    degree_bound = ramification * (
        max(degree for value in gaussian.values() for degree in value)
        + 2 * n * modes * (prime - 1)
    ) + max(degree + clearing for value in kernel.values() for degree in value)
    actual_degree = max(degree for value in tensor.values() for degree in value)
    assert actual_degree <= degree_bound
    dense = {variables: tuple(value.get(degree, 0) for degree in range(max(value) + 1))
             for variables, value in tensor.items()}
    actual = physical_cut_profiles(dense, prime, n, ramification, charge)
    compact = profile_data(prime, n, ramification, charge, modes, kernel, geometry, endpoint=endpoint)
    expected = (
        tuple(sorted(a + b + 2 * clearing for a in compact['zero_active']
                     for b in expand_spectators(compact['zero_spectators'], prime))),
        tuple(sorted(a + b - 2 * clearing for a in compact['infinity_active']
                     for b in expand_spectators(compact['infinity_spectators'], prime))),
    )
    assert actual == expected, (prime, n, ramification, charge, modes, geometry, variant, actual, expected)
    return {'prime': prime, 'n': n, 'ramification': ramification, 'charge': charge,
            'modes': modes, 'geometry': geometry, 'top_digit_pole': variant,
            'endpoint_policy': endpoint,
            'kernel': entry_list(kernel), 'clearing_power': clearing,
            'original_matrix_size': prime**(2 * n),
            'endpoint_coefficients_checked': len(set(tensor) | set(gaussian)),
            'polynomial_degree': actual_degree, 'proved_degree_bound': degree_bound,
            'zero_profile': dict(sorted(Counter(actual[0]).items())),
            'infinity_profile': dict(sorted(Counter(actual[1]).items())),
            'seconds': time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--extended', action='store_true')
    args = parser.parse_args()
    diagonal = []
    for prime, n, modes in ((2, 3, 2), (2, 3, 3), (3, 2, 2)):
        alpha = monomial_filters(prime, modes, tuple(sum(state)**2 - 2 * sum(state)
                                                    for state in states(prime, modes)))
        kernel = {(state, state): value for state, value in alpha.items()}
        actual = active_tensor(prime, n, 1, modes, kernel)
        expected = diagonal_tensor(prime, n, 1, modes, alpha)
        assert actual == expected
        diagonal.append({'prime': prime, 'n': n, 'modes': modes,
                         'exact_active_coefficients': len(actual)})
    cases = [(2, 2, 1, 0, 2, 'initial', False),
             (2, 2, 3, 1, 2, 'initial', True),
             (2, 3, 1, -1, 2, 'initial', True),
             (2, 3, 1, 2, 3, 'initial', True),
             (3, 2, 1, 0, 2, 'initial', True),
             (2, 2, 1, 0, 2, 'initial', True, 'filtered'),
             (2, 3, 1, 2, 3, 'initial', True, 'filtered'),
             (3, 2, 1, -1, 2, 'initial', True, 'filtered')]
    if args.extended:
        cases.append((2, 4, 1, 1, 2, 'mixed', True))
    physical = [physical_control(*case) for case in cases]
    invalid = identity_kernel(2, 2)
    invalid[(0, 1), (1, 0)] = {0: 1}
    try:
        normalized_kernel(2, 2, invalid)
    except ValueError:
        pass
    else:
        raise AssertionError('A nonidentity endpoint was accepted.')
    invalid = identity_kernel(2, 2)
    invalid[(1, 0), (0, 1)] = {0: 1}
    try:
        normalized_kernel(2, 2, invalid, endpoint='filtered')
    except ValueError:
        pass
    else:
        raise AssertionError('A descending physical-energy endpoint transfer was accepted.')
    invalid = identity_kernel(2, 2)
    del invalid[(1, 1), (1, 1)]
    try:
        normalized_kernel(2, 2, invalid, endpoint='filtered')
    except ValueError:
        pass
    else:
        raise AssertionError('A singular endpoint was accepted.')
    # The native path must retain a coefficient past its initial 512-term
    # window, and it must reject a source beyond the supported maximum.
    late = [[{0: 1}, {0: 1}], [{0: 1}, {0: 1, 600: 1}]]
    assert _laurent_orders(late, 2) == (0, 600)
    try:
        _laurent_orders([[{0: 1, 8192: 1}]], 2)
    except ArithmeticError:
        pass
    else:
        raise AssertionError('An oversized complete polynomial was silently truncated.')
    result = {'status': 'Exact operator-kernel, endpoint, and physical-profile controls; no pressure claim.',
              'diagonal_specializations': diagonal, 'original_physical_profiles': physical,
              'invalid_endpoint_rejected': True,
              'nonregular_endpoint_repair_rejected': True,
              'singular_endpoint_rejected': True,
              'late_coefficient_profile': [0, 600],
              'unsupported_source_precision_rejected': True,
              'extended': args.extended}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'diagonal_specializations': len(diagonal),
                      'original_cut_profiles': 2 * len(physical),
                      'endpoint_coefficients_checked': sum(row['endpoint_coefficients_checked'] for row in physical)}))


if __name__ == '__main__':
    main()
