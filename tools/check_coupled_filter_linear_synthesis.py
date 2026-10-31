#!/usr/bin/env python3
"""Brute-force and exact rank controls for polynomial coupled-filter synthesis."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

import coupled_identity_extreme_filters as extreme
import coupled_identity_mode_filters as initial
from coupled_filter_linear_synthesis import json_ready, synthesize, verify_witness


def endpoint_polynomials(prime, degree):
    for high in product(range(prime), repeat=degree):
        coefficients = ((1 - sum(high)) % prime,) + high
        yield {power: value for power, value in enumerate(coefficients) if value}


def literal_margin(prime, n, modes, geometry, alpha, balance):
    if geometry == 'initial':
        selected, retained = tuple(range(modes)), tuple(range(n - modes, n))
        tensor = initial.active_tensor(prime, n, 1, modes, alpha)
    else:
        data = extreme.geometry_data(prime, n, modes, geometry)
        selected, retained = data['selected'], data['retained']
        tensor = extreme.active_tensor(prime, n, 1, modes, alpha, geometry)
    margins = []
    for variables, polynomial in tensor.items():
        v, y = variables[:modes], variables[modes:2 * modes]
        u, x = variables[2 * modes:3 * modes], variables[3 * modes:]
        row = sum(j * digit for j, digit in zip(selected, y)) + sum(
            (j + 1) * digit for j, digit in zip(retained, x))
        column = sum(j * digit for j, digit in zip(selected, u)) + sum(
            (j + 1) * digit for j, digit in zip(retained, v))
        margins.append(Fraction(min(polynomial), 2) - balance * row - (1 - balance) * column)
    return min(margins)


def brute_force_controls():
    states = tuple(product(range(2), repeat=2))
    values = tuple(endpoint_polynomials(2, 2))
    balances = (Fraction(0), Fraction(1, 2), Fraction(1))
    count = vectors = 0
    for geometry in ('initial', 'mixed'):
        catalog = []
        for polynomials in product(values, repeat=4):
            alpha = dict(zip(states, polynomials))
            degrees = tuple(max(value) for value in polynomials)
            margins = tuple(literal_margin(2, 2, 2, geometry, alpha, balance)
                            for balance in balances)
            catalog.append((degrees, margins))
            vectors += 1
        for degrees, h, index in product(product(range(3), repeat=4), range(3), range(3)):
            expected = any(all(a <= b for a, b in zip(old_degrees, degrees)) and margins[index] >= h
                           for old_degrees, margins in catalog)
            result = synthesize(2, 2, 1, 2, geometry, h, degrees, balances[index])
            assert result['feasible'] == expected
            count += 1
    return {'literal_endpoint_polynomial_vectors': vectors,
            'exact_brute_force_feasibility_comparisons': count,
            'maximum_individual_degree': 2}


def scalar_only_boxes():
    records = []
    for prime, n, geometry in ((2, 8, 'initial'), (2, 8, 'terminal'), (2, 8, 'mixed'),
                               (3, 3, 'initial'), (3, 3, 'terminal'), (3, 3, 'mixed'),
                               (3, 9, 'mixed')):
        for h, balance in product(sorted({0, 1, 2, n // 2, n}),
                                  (Fraction(0), Fraction(1, 2))):
            result = synthesize(prime, n, 1, 2, geometry, h,
                                (h,) * (prime ** 2), balance)
            assert result['feasible'] and result['homogeneous_nullity'] == 1
            assert all(value == {h: 1} for value in result['alpha'].values())
            result.pop('alpha')
            records.append(result)
    return {'scope': 'Only the displayed finite degree boxes and sufficient energy-split certificate.',
            'exact_scalar_only_boxes': records}


def rare_state_compensation_controls():
    records = []
    for prime, n, geometry in ((2, 8, 'initial'), (2, 8, 'mixed'), (3, 3, 'mixed')):
        for h, large_degree, balance in product((1, 2), (n, 2 * n),
                                                (Fraction(0), Fraction(1, 2))):
            degrees = [large_degree] * (prime ** 2)
            for j in range(prime):
                degrees[(prime - 1) * prime + j] = h if j == prime - 1 else h - 1
            result = synthesize(prime, n, 1, 2, geometry, h, degrees, balance)
            assert not result['feasible']
            records.append(result)
    return {'scope': 'The specified degree caps; arbitrary polynomial filters are not excluded.',
            'infeasible_occupied_state_improvement_boxes': records}


def nonidentity_positive_control():
    prime, n, modes, geometry = 2, 4, 2, 'mixed'
    states = tuple(product(range(prime), repeat=modes))
    alpha = {state: {0: 1} for state in states}
    # A high-order endpoint-vanishing correction is regular at zero, but
    # its full degree makes it expensive at infinity. This confirms that
    # the sufficient regularity test can admit nonidentity filters.
    alpha[(0, 0)] = {0: 1, 32: 1, 33: 1}
    result = verify_witness(prime, n, 1, modes, geometry, 0,
                            (33, 0, 0, 0), Fraction(1, 2), alpha)
    assert literal_margin(prime, n, modes, geometry, alpha, Fraction(1, 2)) >= 0
    return {'prime': prime, 'n': n, 'geometry': geometry,
            'full_extra_polynomial_degree': 33, 'verification': result,
            'pressure_improvement_claimed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {
        'status': 'Exact sufficient-certificate synthesis; no pressure witness or global exclusion.',
        'literal_polynomial_enumeration': brute_force_controls(),
        'finite_scalar_only_regions': scalar_only_boxes(),
        'finite_compensation_regions': rare_state_compensation_controls(),
        'nonidentity_regular_control': nonidentity_positive_control(),
    }
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
