#!/usr/bin/env python3
"""Exact controls for residue-compressed polynomial state-filter synthesis."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
from random import Random

from polynomial_state_filter_synthesis import (
    difference_matrix, polynomial_seed, profile_data, residue_patterns,
    synthesize, transform,
)


def direct_feasible(prime, m, degrees, orders):
    """Independent full coefficient system, without residue compression."""
    variables = [(digit, degree) for digit, bound in enumerate(degrees)
                 for degree in range(bound + 1)]
    matrix = difference_matrix(prime)
    equations = []
    for k, order in enumerate(orders):
        for degree in range(order):
            row = []
            for digit, power in variables:
                remainder = degree - power
                row.append(matrix[k][digit].get(remainder // m, 0)
                           if remainder >= 0 and remainder % m == 0 else 0)
            equations.append(row + [0])
    for digit in range(prime):
        equations.append([int(j == digit) for j, degree in variables] + [1])
    # Only forward elimination is needed to detect an inconsistent row.
    pivots = {}
    for row in equations:
        for column in range(len(variables)):
            coefficient = row[column] % prime
            if not coefficient:
                continue
            if column in pivots:
                old = pivots[column]
                row = [(a - coefficient * b) % prime for a, b in zip(row, old)]
            else:
                inverse = pow(coefficient, -1, prime)
                pivots[column] = [value * inverse % prime for value in row]
                break
        else:
            if row[-1] % prime:
                return False
    return True


def coefficient_controls():
    rng = Random(13371337)
    count = feasible = 0
    for prime in (2, 3, 5, 7):
        for _ in range(200):
            m = rng.randrange(1, 7)
            degrees = tuple(rng.randrange(0, 9) for _ in range(prime))
            orders = tuple(rng.randrange(0, 9 + (prime - 1) * m) for _ in range(prime))
            expected = direct_feasible(prime, m, degrees, orders)
            actual = synthesize(prime, m, degrees, orders)
            assert (actual is not None) == expected
            assert (synthesize(prime, m, degrees, orders, construct=False) is True) == expected
            assert len(residue_patterns(m, degrees, orders)) <= 2 * prime + 1
            feasible += expected
            count += 1
    return {'full_affine_system_comparisons': count, 'feasible_controls': feasible,
            'characteristics': [2, 3, 5, 7]}


def enumerate_endpoint_polynomials(prime, bound):
    for high in product(range(prime), repeat=bound):
        coefficients = ((1 - sum(high)) % prime,) + high
        yield {degree: value for degree, value in enumerate(coefficients) if value}


def brute_force_controls():
    prime = 3
    polynomials = tuple(enumerate_endpoint_polynomials(prime, 2))
    total = queries = 0
    rng = Random(883)
    for m in (1, 2, 3, 5):
        profiles = []
        for alpha in product(polynomials, repeat=prime):
            beta = transform(alpha, m, prime)
            profiles.append((tuple(max(value) for value in alpha),
                             tuple(min(value) for value in beta)))
            total += 1
        bounds = list(product(range(3), repeat=prime))
        targets = [tuple(rng.randrange(2 * m + 4) for _ in range(prime)) for _ in range(80)]
        targets += [(0, m, 2 * m), (0, 0, 0), (1, 1, 1)]
        for degree_bounds, order_bounds in product(bounds, targets):
            expected = any(all(a <= b for a, b in zip(degrees, degree_bounds))
                           and all(a >= b for a, b in zip(orders, order_bounds))
                           for degrees, orders in profiles)
            actual = synthesize(prime, m, degree_bounds, order_bounds, construct=False)
            assert (actual is True) == expected
            queries += 1
    return {'enumerated_endpoint_polynomial_vectors': total,
            'brute_force_feasibility_comparisons': queries,
            'maximum_individual_degree': 2, 'characteristic': prime}


def large_scale_controls():
    prime = 3
    scale = 2 ** 209 + 3 ** 117
    controls = []
    for m, degrees, orders in (
            (1, (0, 0, 0), (0, 1, 2)),
            (3, (2, 4, 5), (2, 5, 2)),
            (4, (4, 0, 3), (0, 0, 8)),
            (2, (3, 2, 4), (3, 4, 2))):
        small = synthesize(prime, m, degrees, orders)
        large = synthesize(prime, m * scale,
                           tuple(value * scale for value in degrees),
                           tuple(value * scale for value in orders))
        assert (small is None) == (large is None)
        controls.append({'m_bits': (m * scale).bit_length(),
                         'feasible': large is not None,
                         'residue_patterns': len(residue_patterns(
                             m * scale, tuple(value * scale for value in degrees),
                             tuple(value * scale for value in orders)))})
    assert synthesize(3, scale, (0, 0, 0), (scale ** 2,) * 3) is None
    return {'scaled_instances': controls,
            'very_large_impossible_order_bit_length': (scale ** 2).bit_length()}


def prime_five_brute_force_controls():
    """Independent degree-one enumeration supplied by the audit agent."""
    prime = 5
    result = []
    for m in (1, 2, 7):
        catalog = set()
        for values in product(range(prime), repeat=prime):
            alpha = tuple({degree: value for degree, value in
                           ((0, constant), (1, (1 - constant) % prime)) if value}
                          for constant in values)
            beta = transform(alpha, m, prime)
            catalog.add((tuple(max(value) for value in alpha),
                         tuple(min(value) for value in beta)))
        cases = set()
        for degrees, orders in catalog:
            cases.add((degrees, orders))
            for k in range(prime):
                tightened = list(orders)
                tightened[k] += 1
                cases.add((degrees, tuple(tightened)))
                if degrees[k]:
                    smaller = list(degrees)
                    smaller[k] = 0
                    cases.add((tuple(smaller), orders))
        rng = Random(730919 + m)
        for _ in range(1000):
            degrees = tuple(rng.randrange(2) for _ in range(prime))
            orders = tuple(rng.randrange(4 * m + 2) for _ in range(prime))
            cases.add((degrees, orders))
        feasible = 0
        for degrees, orders in sorted(cases):
            expected = any(all(a <= b for a, b in zip(old_degrees, degrees))
                           and all(a >= b for a, b in zip(old_orders, orders))
                           for old_degrees, old_orders in catalog)
            actual = synthesize(prime, m, degrees, orders)
            assert (actual is not None) == expected
            feasible += expected
        result.append({'m': m, 'endpoint_one_vectors_enumerated': prime ** prime,
                       'distinct_degree_order_pairs': len(catalog),
                       'feasibility_cases': len(cases), 'feasible_cases': feasible})
    return result


def full_tensor_controls():
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel
    from check_truncated_symmetric_identity_filters import direct_profiles, expand_qutrit_modes

    rng = Random(441)
    choices = []
    for n, ramification in ((1, 1), (1, 3), (2, 1)):
        for _ in range(2 if n == 1 else 1):
            alpha = []
            for digit in range(3):
                high = [rng.randrange(3) for _ in range(3)]
                coefficients = [(1 - sum(high)) % 3] + high
                alpha.append({degree: value for degree, value in enumerate(coefficients) if value})
            choices.append((n, ramification, tuple(alpha)))
    endpoints = 0
    controls = []
    evaluate = lambda polynomial: sum(value * (-1) ** degree
                                      for degree, value in enumerate(polynomial)) % 3
    for n, ramification, alpha in choices:
        tensor = polynomial_seed(3, n, ramification, alpha)
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, 3), 3)
        for powers in product(range(3), repeat=4 * n):
            assert evaluate(tensor.get(powers, ())) == evaluate(gaussian.get(powers, ()))
            endpoints += 1
        for charge in (0, 2):
            actual, precision = direct_profiles(tensor, n, ramification, charge)
            data = profile_data(3, n, ramification, charge, alpha)
            expected = []
            for cut in ('zero', 'infinity'):
                spectators = expand_qutrit_modes(data[cut + '_spectators'])
                expected.append(tuple(sorted(a + b for a in data[cut + '_active']
                                             for b in spectators)))
            assert actual == expected
            controls.append({'n': n, 'ramification': ramification, 'charge': charge,
                             'alpha': alpha, 'matrix_order': 3 ** (2 * n),
                             'maximum_series_precision': precision})
    return {'endpoint_coefficient_equalities': endpoints,
            'complete_two_cut_controls': controls}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {
        'status': 'Exact endpoint-constrained synthesis and profile controls; no pressure claim.',
        'uncompressed_coefficient_systems': coefficient_controls(),
        'brute_force_polynomials': brute_force_controls(),
        'large_integer_scales': large_scale_controls(),
        'prime_five_independent_enumeration': prime_five_brute_force_controls(),
        'literal_polynomial_tensors': full_tensor_controls(),
    }
    encoded = json.dumps(result, indent=2, default=str) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
