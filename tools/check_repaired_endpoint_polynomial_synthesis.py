#!/usr/bin/env python3
"""Independent exact checks for nonzero-endpoint polynomial synthesis."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
from random import Random

from polynomial_state_filter_synthesis import difference_matrix, transform, synthesize as strict_synthesize
from repaired_endpoint_polynomial_synthesis import synthesize


def affine_coordinate_feasible(prime, m, degrees, orders, coordinate):
    variables = tuple((j, degree) for j, bound in enumerate(degrees) for degree in range(bound + 1))
    matrix = difference_matrix(prime)
    equations = []
    for k, order in enumerate(orders):
        for degree in range(order):
            row = []
            for j, power in variables:
                remainder = degree - power
                row.append(matrix[k][j].get(remainder // m, 0)
                           if remainder >= 0 and remainder % m == 0 else 0)
            equations.append(row + [0])
    equations.append([int(j == coordinate) for j, power in variables] + [1])
    pivots = {}
    for equation in equations:
        row = list(equation)
        for column in range(len(variables)):
            coefficient = row[column] % prime
            if not coefficient:
                continue
            if column in pivots:
                row = [(a - coefficient * b) % prime for a, b in zip(row, pivots[column])]
            else:
                inverse = pow(coefficient, -1, prime)
                pivots[column] = [value * inverse % prime for value in row]
                break
        else:
            if row[-1] % prime:
                return False
    return True


def full_system_controls():
    rng = Random(573911)
    count = feasible = 0
    for prime in (2, 3, 5, 7):
        for _ in range(100):
            m = rng.randrange(1, 6)
            degrees = tuple(rng.randrange(0, 8) for _ in range(prime))
            orders = tuple(rng.randrange(0, 8 + (prime - 1) * m) for _ in range(prime))
            possible = [affine_coordinate_feasible(prime, m, degrees, orders, coordinate)
                        for coordinate in range(prime)]
            actual = synthesize(prime, m, degrees, orders)
            assert (actual is not None) == all(possible)
            if actual:
                assert actual['extension_degree'] == 1
                assert all(value[0] for value in actual['endpoint_values_in_extension_basis'])
                feasible += 1
            count += prime
    return {'independent_full_affine_coordinate_systems': count,
            'joint_bound_queries': 400, 'feasible_queries': feasible}


def brute_force_controls():
    prime = 3
    polynomials = tuple({degree: value for degree, value in enumerate(coefficients) if value}
                        for coefficients in product(range(prime), repeat=3)
                        if sum(coefficients) % prime)
    rng = Random(1458)
    vectors = queries = 0
    for m in (1, 2, 4):
        catalog = set()
        for alpha in product(polynomials, repeat=prime):
            beta = transform(alpha, m, prime)
            catalog.add((tuple(max(value) for value in alpha),
                         tuple(min(value) for value in beta)))
            vectors += 1
        targets = [tuple(rng.randrange(2 * m + 4) for _ in range(prime)) for _ in range(60)]
        targets += [(0, m, 2 * m), (1, 0, 1), (0, 0, 0)]
        for degrees, orders in product(product(range(3), repeat=prime), targets):
            expected = any(all(a <= b for a, b in zip(old_degrees, degrees))
                           and all(a >= b for a, b in zip(old_orders, orders))
                           for old_degrees, old_orders in catalog)
            actual = synthesize(prime, m, degrees, orders, construct=False)
            assert (actual is True) == expected
            queries += 1
    return {'nonzero_endpoint_polynomial_vectors_enumerated': vectors,
            'brute_force_bound_comparisons': queries,
            'individual_polynomial_degree_bound': 2}


def strict_extension_control():
    for m in (2, 3, 16, 2 ** 201 + 1):
        degrees, orders = (0, 0, 1), (1, 0, 1)
        assert strict_synthesize(3, m, degrees, orders) is None
        result = synthesize(3, m, degrees, orders)
        assert result and result['actual_degrees'] == degrees and result['actual_orders'] == orders
        alpha = result['alpha_basis'][0]
        assert alpha[0] == {0: 1} and alpha[1] == {0: 2} and alpha[2] == {1: 1}
    return {'strictly_broader_feasibility_example': 'alpha=(1,-1,x)',
            'degrees': [0, 0, 1], 'zero_orders': [1, 0, 1],
            'mode_scales_tested': [2, 3, 16, str(2 ** 201 + 1)]}


def polynomial_tensor_controls():
    from coupled_identity_mode_filters import polynomial_seed
    from check_coupled_identity_mode_filters import direct_profiles, expand_modes
    from check_all_field_pascal_crossings import pascal_kernel
    from check_truncated_symmetric_crossings import direct_kernel

    records = []
    endpoint_entries = 0
    for n, ramification, degrees, orders in ((2, 1, (0, 0, 1), (1, 0, 1)),
                                             (1, 3, (7, 4, 0), (0, 4, 7))):
        result = synthesize(3, n * ramification, degrees, orders)
        alpha = {(j,): value for j, value in enumerate(result['alpha_basis'][0])}
        tensor, clearing = polynomial_seed(3, n, ramification, 1, alpha,
                                            strict_endpoint=False, repair_endpoint=True)
        assert clearing == 0
        gaussian = direct_kernel(pascal_kernel(tuple(range(n)), 2 * n, 3), 3)
        evaluate = lambda value: sum(coefficient * (-1) ** i
                                     for i, coefficient in enumerate(value)) % 3
        for powers in product(range(3), repeat=4 * n):
            assert evaluate(tensor.get(powers, ())) == evaluate(gaussian.get(powers, ()))
            endpoint_entries += 1
        actual, precision = direct_profiles(3, tensor, n, ramification, 0)
        a = ramification
        beta_orders = [2 * value for value in result['actual_orders']]
        zero = [beta_orders[k] + a * (digit - k) for digit, k in product(range(3), repeat=2)]
        infinity = [ramification * (1 - 2 * n) * digit - 2 * degree
                    for digit, degree in enumerate(result['actual_degrees'])]
        zero_spectators = [ramification * (2 * i + 1) for i in range(1, n)]
        zero_spectators += [ramification * (2 * i + 1) for i in range(n - 1)]
        infinity_spectators = [ramification * (2 * i + 1 - 2 * n) for i in range(1, 2 * n)]
        expected = [tuple(sorted(x + y for x in active for y in expand_modes(3, spectators)))
                    for active, spectators in ((zero, zero_spectators), (infinity, infinity_spectators))]
        assert actual == expected
        records.append({'n': n, 'ramification': ramification, 'degree_bounds': degrees,
                        'order_bounds': orders, 'alpha': list(alpha.values()),
                        'full_matrix_order': 3 ** (2 * n),
                        'maximum_exact_jet_precision': precision})
    return {'physical_endpoint_equalities': endpoint_entries,
            'complete_polynomial_witness_cut_replays': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {
        'status': 'Exact nonzero-endpoint synthesis with diagonal repair; no pressure gain claimed.',
        'uncompressed_linear_systems': full_system_controls(),
        'brute_force_polynomials': brute_force_controls(),
        'strict_endpoint_comparison': strict_extension_control(),
        'literal_repaired_polynomial_witnesses': polynomial_tensor_controls(),
    }
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
