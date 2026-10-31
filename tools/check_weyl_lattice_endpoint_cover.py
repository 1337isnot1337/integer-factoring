#!/usr/bin/env python3
"""Exact F4 endpoint-cover tests and three-direction polynomial controls."""
from __future__ import annotations

import argparse
from itertools import accumulate, product
import json
from pathlib import Path
from random import Random

from f4_repaired_identity_filters import (
    add, evaluate_one, multiply, polynomial_seed, profile_data, reduced_zero_profiles,
)
from weyl_lattice_filter_synthesis import (
    binary_image_basis, bound_blocks, covering_pair, endpoint_image, lattice_basis,
    nonzero_endpoint_lift, ones_endpoint_possible, polynomial_degrees, synthesize,
)


def evaluate_endpoint(functionals, coefficients):
    result = []
    for functional in functionals:
        value = 0
        for index, coefficient in enumerate(coefficients):
            if functional & (1 << index):
                value ^= coefficient
        result.append(value)
    return result


def endpoint_controls():
    rng = Random(490011)
    count = assignments = feasible = 0
    for size in (4, 8):
        for _ in range(120):
            variables = rng.randrange(1, 6)
            functionals = tuple(rng.randrange(1 << variables) for _ in range(size))
            expected = False
            for coefficients in product(range(4), repeat=variables):
                assignments += 1
                if all(evaluate_endpoint(functionals, coefficients)):
                    expected = True
            basis = binary_image_basis(functionals)
            pair = covering_pair(basis, size)
            assert (pair is not None) == expected
            binary_expected = any(all(evaluate_endpoint(functionals, coefficients))
                                  for coefficients in product(range(2), repeat=variables))
            assert ones_endpoint_possible(functionals) == binary_expected
            if size == 4:
                assert bool(all(functionals)) == expected
            if expected:
                image = endpoint_image(tuple(range(variables)), {}, functionals)
                endpoint, lift = nonzero_endpoint_lift(image, size)
                coefficients = tuple(lift.get(j, 0) for j in range(variables))
                assert evaluate_endpoint(functionals, coefficients) == endpoint and all(endpoint)
                feasible += 1
            count += 1
    # Seven projective F2 functionals on three variables, plus one repeat.
    # Every projection is nonzero, but any two binary images miss a coordinate.
    fano = tuple(range(1, 8)) + (1,)
    assert all(fano) and covering_pair(binary_image_basis(fano), 8) is None
    assert not any(all(evaluate_endpoint(fano, coefficients))
                   for coefficients in product(range(4), repeat=3))
    return {'independent_endpoint_maps': count, 'complete_F4_coefficient_assignments': assignments,
            'feasible_endpoint_maps': feasible,
            'eight_state_full_projection_counterexample': list(fano),
            'counterexample_tested_over_all_F4_assignments': 64}


def transformed_row_orders(data, alpha):
    entries = {}
    for state, basis in zip(data['states'], data['coefficients']):
        for entry, coefficient in basis.items():
            polynomial = {degree: 1 for degree in polynomial_degrees(coefficient)}
            entries[entry] = add(entries.get(entry, {}), multiply(polynomial, alpha[state]))
    return tuple(min(min(value) for (row, column), value in entries.items() if row == i and value)
                 for i in range(len(data['states'])))


def physical_controls():
    from check_f4_repaired_identity_filters import expanded, physical_profiles
    from check_weyl_lattice_filter_synthesis import binary_rank

    records, endpoints = [], 0
    for n, ramification, charge, geometry in ((3, 1, 0, 'initial'), (3, 1, 2, 'initial'),
                                              (3, 3, -1, 'initial'), (4, 1, 1, 'terminal'),
                                              (4, 1, 0, 'split:1:2')):
        data = lattice_basis(n, ramification, geometry, modes=3)
        for block in data['blocks']:
            size = len(block['indices'])
            for matrix in (block['left'], block['right']):
                assert binary_rank([sum((value & 1) << j for j, value in enumerate(row))
                                    for row in matrix]) == size
        alpha = {state: {0: 1 + (i + charge) % 3} for i, state in enumerate(data['states'])}
        alpha[(1, 1, 1)] = {1: 2, 2: 1}
        alpha[(1, 1, 0)] = {0: 1, 2: 1, 3: 1}
        degrees = tuple(max(alpha[state]) for state in data['states'])
        targets = transformed_row_orders(data, alpha)
        witness = synthesize(n, ramification, geometry, degrees, targets, modes=3)
        assert witness and any(max(value) for value in witness['alpha'].values())
        actual = reduced_zero_profiles(n, ramification, 3, witness['alpha'], geometry)
        bounds = bound_blocks(n, ramification, geometry, targets, modes=3)
        for block, bound in zip(actual, bounds):
            assert block['charge'] == bound['charge']
            assert all(a >= b for a, b in zip(accumulate(block['orders']), accumulate(bound['orders'])))
        tensor, clearing, gaussian = polynomial_seed(n, ramification, 3, witness['alpha'], geometry)
        assert clearing == 0
        for powers in product(range(2), repeat=4 * n):
            assert evaluate_one(tensor.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
            endpoints += 1
        physical, precision = physical_profiles(tensor, n, ramification, charge)
        compact = profile_data(n, ramification, charge, 3, witness['alpha'], geometry, actual)
        expected = (
            tuple(sorted(a + b for a in compact['zero_active'] for b in expanded(compact['zero_spectators']))),
            tuple(sorted(a + b for a in compact['infinity_active'] for b in expanded(compact['infinity_spectators']))),
        )
        assert physical == expected
        records.append({'n': n, 'ramification': ramification, 'charge': charge, 'geometry': geometry,
                        'row_targets': targets, 'actual_degrees': witness['actual_degrees'],
                        'endpoint_values': witness['endpoint_values'],
                        'endpoint_image_dimension': witness['endpoint_image_dimension'],
                        'full_matrix_order': 2 ** (2 * n), 'maximum_exact_jet_precision': precision})
    return {'complete_physical_witnesses': records,
            'physical_endpoint_coefficients_checked': endpoints}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'Exact F4 endpoint-cover criterion and repaired three-mode polynomial witnesses.',
              'finite_field_endpoint_cover': endpoint_controls(),
              'three_mode_physical_controls': physical_controls()}
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
