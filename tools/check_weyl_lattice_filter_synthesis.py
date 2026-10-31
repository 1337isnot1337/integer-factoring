#!/usr/bin/env python3
"""Independent finite controls for exact reduced-lattice filter synthesis."""
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
    bound_blocks, endpoint_image, homogeneous_system, json_ready, lattice_basis,
    matrix_multiply, polynomial_degrees, polynomial_order, synthesize,
)


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def basis_controls():
    records = []
    for n, ramification, geometry in ((2, 1, 'initial'), (3, 3, 'initial'),
                                       (8, 1, 'initial'), (8, 3, 'terminal'),
                                       (8, 1, 'mixed'), (32, 3, 'mixed')):
        data = lattice_basis(n, ramification, geometry)
        for block in data['blocks']:
            size = len(block['indices'])
            for matrix in (block['left'], block['right']):
                assert binary_rank([sum((value & 1) << j for j, value in enumerate(row))
                                    for row in matrix]) == size
            actual = matrix_multiply(matrix_multiply(block['left'], block['gaussian']), block['right'])
            assert actual == block['diagonal']
            assert all(not actual[i][j] for i in range(size) for j in range(size) if i != j)
            assert tuple(polynomial_order(actual[i][i]) for i in range(size)) == block['orders']
        witness = synthesize(n, ramification, geometry, (0,) * 4, data['gaussian_orders'])
        assert witness and witness['endpoint_values'] == [1] * 4
        records.append({'n': n, 'ramification': ramification, 'geometry': geometry,
                        'gaussian_row_orders': data['gaussian_orders'],
                        'exact_unit_basis_blocks': len(data['blocks'])})
    return records


def row_orders(data, alpha):
    entries = {}
    for state, basis in zip(data['states'], data['coefficients']):
        for entry, coefficient in basis.items():
            polynomial = {degree: 1 for degree in polynomial_degrees(coefficient)}
            entries[entry] = add(entries.get(entry, {}), multiply(polynomial, alpha[state]))
    return tuple(min(min(value) for (row, column), value in entries.items() if row == i and value)
                 for i in range(4))


def brute_force_controls():
    rng = Random(918174)
    polynomials = tuple((a, b) for a, b in product(range(4), repeat=2) if a ^ b)
    enumerated = comparisons = 0
    strict_examples = []
    for n, geometry in ((2, 'initial'), (4, 'mixed')):
        data = lattice_basis(n, 1, geometry)
        entry_keys = tuple(sorted(set().union(*(basis.keys() for basis in data['coefficients']))))
        contributions = []
        for basis in data['coefficients']:
            contributions.append({poly: tuple(tuple(
                (coefficient if a & (1 << plane) else 0)
                ^ (coefficient << 1 if b & (1 << plane) else 0)
                for plane in range(2)) for key in entry_keys
                for coefficient in (basis.get(key, 0),))
                for poly in polynomials for a, b in (poly,)})
        catalog = set()
        for chosen in product(polynomials, repeat=4):
            values = [contributions[j][polynomial] for j, polynomial in enumerate(chosen)]
            row_minima = [None] * 4
            for entry_index, (row, column) in enumerate(entry_keys):
                first = second = 0
                for value in values:
                    first ^= value[entry_index][0]
                    second ^= value[entry_index][1]
                if first | second:
                    order = polynomial_order(first | second)
                    row_minima[row] = order if row_minima[row] is None else min(row_minima[row], order)
            assert all(value is not None for value in row_minima)
            degrees = tuple(int(bool(b)) for a, b in chosen)
            catalog.add((degrees, tuple(row_minima)))
            enumerated += 1
        targets = [tuple(rng.randrange(2 * n + 4) for _ in range(4)) for _ in range(96)]
        targets += [orders for degrees, orders in sorted(catalog)]
        targets = tuple(dict.fromkeys(targets))
        for degrees, target in product(product(range(2), repeat=4), targets):
            expected = any(all(a <= b for a, b in zip(old_degrees, degrees))
                           and all(a >= b for a, b in zip(old_orders, target))
                           for old_degrees, old_orders in catalog)
            actual = synthesize(n, 1, geometry, degrees, target, construct=False)
            assert (actual is True) == expected
            if actual and len(strict_examples) < 3:
                variables, pivots, endpoints, count = homogeneous_system(n, 1, geometry, degrees, target)
                image = endpoint_image(variables, pivots, endpoints)
                old = {endpoint: lift for endpoint, lift in image}
                strict = 15
                for endpoint in sorted(old, key=int.bit_length, reverse=True):
                    if strict & (1 << (endpoint.bit_length() - 1)):
                        strict ^= endpoint
                if strict:
                    witness = synthesize(n, 1, geometry, degrees, target)
                    strict_examples.append(json_ready(witness))
            comparisons += 1
    return {'complete_F4_polynomial_vectors_enumerated': enumerated,
            'individual_degree_bound': 1,
            'exact_bound_queries_against_full_catalog': comparisons,
            'nonzero_endpoint_feasible_but_endpoint_one_impossible_examples': strict_examples}


def physical_witness_controls():
    from check_f4_repaired_identity_filters import expanded, physical_profiles

    rng = Random(549111)
    records = []
    endpoint_entries = 0
    n, ramification, geometry = 2, 1, 'initial'
    data = lattice_basis(n, ramification, geometry)
    for case in range(8):
        supplied = {}
        for state in data['states']:
            while True:
                value = {degree: rng.randrange(4) for degree in range(4)}
                value = {degree: coefficient for degree, coefficient in value.items() if coefficient}
                if evaluate_one(value):
                    break
            supplied[state] = value
        degrees = tuple(max(supplied[state]) for state in data['states'])
        targets = row_orders(data, supplied)
        witness = synthesize(n, ramification, geometry, degrees, targets)
        assert witness
        alpha = witness['alpha']
        actual_blocks = reduced_zero_profiles(n, ramification, 2, alpha, geometry)
        certified = bound_blocks(n, ramification, geometry, targets)
        for actual, bound in zip(actual_blocks, certified):
            assert actual['charge'] == bound['charge']
            assert all(a >= b for a, b in zip(accumulate(actual['orders']), accumulate(bound['orders'])))
        tensor, clearing, gaussian = polynomial_seed(n, ramification, 2, alpha, geometry)
        assert clearing == 0
        for powers in product(range(2), repeat=4 * n):
            assert evaluate_one(tensor.get(powers, {})) == evaluate_one(gaussian.get(powers, {}))
            endpoint_entries += 1
        charge = case % 4 - 1
        physical, precision = physical_profiles(tensor, n, ramification, charge)
        compact = profile_data(n, ramification, charge, 2, alpha, geometry, actual_blocks)
        expected = (tuple(sorted(a + b for a in compact['zero_active']
                                 for b in expanded(compact['zero_spectators']))),
                    tuple(sorted(a + b for a in compact['infinity_active']
                                 for b in expanded(compact['infinity_spectators']))))
        assert physical == expected
        records.append({'degree_bounds': degrees, 'row_targets': targets,
                        'endpoint_values': witness['endpoint_values'],
                        'endpoint_image_dimension': witness['endpoint_image_dimension'],
                        'actual_degrees': witness['actual_degrees'],
                        'maximum_exact_jet_precision': precision})
    return {'complete_physical_witnesses': records,
            'physical_endpoint_coefficients_checked': endpoint_entries,
            'full_matrix_order': 16}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'Exact reduced Smith-basis synthesis and complete finite controls; no pressure claim.',
              'regular_polynomial_bases': basis_controls(),
              'complete_small_coefficient_search': brute_force_controls(),
              'literal_repaired_witnesses': physical_witness_controls()}
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
