#!/usr/bin/env python3
"""Independent coefficient, exhaustive affine, and physical synthesis controls."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import random

import numpy as np

from coupled_identity_mode_filters import add, shift
from coupled_identity_operator_filters import _laurent_orders
from prime_deficit_identity_operators import reduced_matrix, polynomial_seed, profile_data
from prime_deficit_operator_synthesis import (
    row_source, coefficient_system, constraints, synthesize, matrix_product_below, flag_kernel,
)
from check_prime_state_filter_profiles import physical_cut_profiles, expand_spectators


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(31297)
    coefficient_equalities = prefix_controls = 0
    source_cases = (((7, 18), (-1, 0, 1)), ((0, 13), (-1, 0, 1)),
                    ((0, 13), (-27, -1, 0, 1, 27)))
    for selected, powers in source_cases:
        data = coefficient_system(3, 27, selected, 1, min(powers), max(powers), powers)['data']
        for _ in range(4):
            kernel = {}
            for index in data['sources']:
                value = {degree: rng.randrange(3) for degree in powers if degree != 0}
                value[0] = (int(index[0] == index[1])-sum(value.values())) % 3
                kernel[index] = {degree: coefficient for degree, coefficient in value.items() if coefficient}
            kernel = {index: value for index, value in kernel.items() if value}
            literal = reduced_matrix(3, 27, 1, 0, 2, 1, kernel, selected_indices=selected)
            accumulated = {}
            cutoff = data['known_below']+min(powers)
            for index, scalar in kernel.items():
                for entry, polynomial in data['sources'][index].items():
                    for degree, coefficient in scalar.items():
                        value = {power: c for power, c in shift(polynomial, degree, 3, coefficient).items()
                                 if power < cutoff}
                        accumulated[entry] = add(accumulated.get(entry, {}), value, 3)
            for frame in data['frames']:
                indices, offset = frame['indices'], frame['row_offset']
                block = [[literal[i][j] for j in indices] for i in indices]
                transformed = matrix_product_below(frame['left'], block, 3, cutoff)
                floors = []
                for i, row in enumerate(transformed):
                    floors.append(min(degree for value in row for degree in value))
                    for j, value in enumerate(row):
                        assert value == accumulated.get((offset+i, j), {})
                        coefficient_equalities += 1
                actual = _laurent_orders(block, 3)
                floors.sort()
                assert all(sum(actual[:k]) >= sum(floors[:k]) for k in range(1, len(actual)+1))
                prefix_controls += 1
    system = coefficient_system(3, 3, (0,), 1, -1, 1)
    data = system['data']
    count = len(system['variables'])
    assignments = np.array(list(product(range(3), repeat=count)), dtype=np.int64)
    queries = []
    for shift_value in (-1, 0, 1):
        queries.append((tuple(value+shift_value for value in data['zero_orders']),
                        tuple(value-2*shift_value for value in data['infinity_orders'])))
    for _ in range(24):
        queries.append((tuple(value+rng.choice((-1, 0, 1)) for value in data['zero_orders']),
                        tuple(value+rng.choice((-2, 0, 2)) for value in data['infinity_orders'])))
    feasible = 0
    for zero, infinity in queries:
        equations = np.array(constraints(system, zero, infinity))
        allowed = np.all((assignments@equations[:, :-1].T-equations[:, -1]) % 3 == 0, axis=1)
        actual = synthesize(3, 3, (0,), 1, -1, 1, zero, infinity)
        assert (actual is not None) == bool(np.any(allowed))
        feasible += actual is not None
    flags = 0
    flag_system = coefficient_system(3, 27, (0, 13), 1, -27, 27, (-27, -1, 0, 1, 27))
    flag_data = flag_system['data']
    alphabet = flag_data['states']
    for _ in range(6):
        basis = [[int(i == j) for j in range(9)] for i in range(9)]
        for i in range(9):
            for j in range(i+1, 9):
                if sum(alphabet[i]) == sum(alphabet[j]):
                    basis[i][j] = rng.randrange(3)
        gains = [0]*9
        gains[rng.randrange(9)] = 1
        kernel = flag_kernel(3, (0, 13), basis, 27, gains)
        accumulated = {}
        for index, scalar in kernel.items():
            for entry, polynomial in flag_data['sources'][index].items():
                for degree, coefficient in scalar.items():
                    value = {power: c for power, c in shift(polynomial, degree, 3, coefficient).items()
                             if power < flag_data['known_below']-27}
                    accumulated[entry] = add(accumulated.get(entry, {}), value, 3)
        zero = tuple(min(degree for (row, column), value in accumulated.items() if row == i
                         for degree in value) for i in range(27))
        infinity = tuple(value+2*gain for value, gain in zip(flag_data['infinity_orders'], gains))
        witness = synthesize(3, 27, (0, 13), 1, -27, 27, zero, infinity,
                             exponents=(-27, -1, 0, 1, 27), infinity_basis=basis)
        assert witness is not None
        raw = [[{flag_data['infinity_orders'][j]-2*degree: coefficient
                 for degree, coefficient in witness['kernel'].get((row, column), {}).items()}
                for j, column in enumerate(alphabet)] for row in alphabet]
        for i in range(9):
            for j in range(9):
                value = {}
                for k in range(9):
                    value = add(value, raw[k][j], 3, basis[i][k])
                assert not value or min(value) >= infinity[i]
        actual = _laurent_orders(raw, 3)
        floors = sorted(infinity)
        assert all(sum(actual[:k]) >= sum(floors[:k]) for k in range(1, 10))
        flags += 1
    physical = []
    for prime, parent, selected, P, J in ((3, 3, (0,), 1, -1), (2, 4, (0, 2), 3, 2)):
        local = coefficient_system(prime, parent, selected, P, -1, 1)['data']
        zero = tuple(value-1 for value in local['zero_orders'])
        infinity = tuple(value+2 for value in local['infinity_orders'])
        witness = synthesize(prime, parent, selected, P, -1, 1, zero, infinity)
        assert witness is not None
        kernel = witness['kernel']
        tensor, clearing, gaussian = polynomial_seed(prime, parent, 1, 0, len(selected), P, kernel,
                                                     selected_indices=selected)
        dense = {powers: tuple(value.get(i, 0) for i in range(max(value)+1))
                 for powers, value in tensor.items()}
        actual = physical_cut_profiles(dense, prime, parent-1, P, J)
        profile = profile_data(prime, parent, 1, 0, len(selected), P, J, kernel, selected_indices=selected)
        expected = (tuple(sorted(a+b+2*clearing for a in profile['zero_active']
                                for b in expand_spectators(profile['zero_spectators'], prime))),
                    tuple(sorted(a+b-2*clearing for a in profile['infinity_active']
                                for b in expand_spectators(profile['infinity_spectators'], prime))))
        assert actual == expected
        physical.append({'prime': prime, 'parent': parent, 'selected': selected, 'ramification': P,
                         'aggregate_charge': J, 'matrix_order': prime**(2*(parent-1)),
                         'endpoint_coefficients': len(set(tensor)|set(gaussian)),
                         'coefficient_variables': witness['variables']})
    bad = list(data['zero_orders'])
    bad[0] = data['known_below']+system['lower']+1
    try:
        constraints(system, bad, data['infinity_orders'])
    except ArithmeticError:
        pass
    else:
        raise AssertionError('An insufficient-precision target was accepted.')
    result = {'status': 'Exact finite-frame affine synthesis controls; no pressure improvement claim.',
              'literal_transformed_entry_equalities': coefficient_equalities,
              'complete_Smith_prefix_controls': prefix_controls,
              'complete_affine_assignments_per_query': len(assignments),
              'exhaustive_affine_target_queries': len(queries), 'feasible_exhaustive_queries': feasible,
              'nontrivial_infinity_flag_witnesses': flags,
              'original_physical_controls': physical, 'insufficient_precision_rejection': True}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
