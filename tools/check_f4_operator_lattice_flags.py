#!/usr/bin/env python3
"""Literal constant-kernel tests of F4 operator flags and scalar Pascal sources."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import random

from check_weyl_operator_lattice_synthesis import oracle
from coupled_identity_operator_filters import states
from coupled_identity_weyl_profiles import weyl_data
from f4_identity_operator_filters import normalized_kernel
from f4_weyl_operator_lattice_synthesis import synthesize
from scalar_pascal_operator_basis import operator_basis as scalar_basis
from weyl_flag_filter_synthesis import full_flags
from weyl_lattice_filter_synthesis import polynomial_multiply
from weyl_operator_lattice_synthesis import operator_basis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    scalar_cases, entries = [], 0
    for n, P, r, geometry in ((2, 1, 2, 'mixed'), (4, 1, 2, 'initial'),
                               (8, 3, 2, 'terminal'), (4, 1, 3, 'initial'),
                               (8, 1, 3, 'terminal'), (32, 1, 3, 'split:1:2'),
                               (32, 3, 3, 'split:2:1')):
        compact, general = scalar_basis(n, P, r, geometry), operator_basis(n, P, r, geometry)
        assert compact['gaussian_orders'] == general['gaussian_orders']
        data = weyl_data(2, n, r, geometry)
        numerator = data['inverse_parameter_numerator']
        assert all(not numerator[i][j] for i in range(r) for j in range(r) if i != j)
        assert all(numerator[i][i] == numerator[0][0] for i in range(r))
        unit = 1
        factor = sum(1 << (P*degree) for degree, value in numerator[0][0].items() if value)
        for _ in range(r):
            unit = polynomial_multiply(unit, factor)
        assert unit & 1
        for index, coefficients in compact['coefficients'].items():
            assert {entry: polynomial_multiply(unit, value) for entry, value in coefficients.items()} == general['coefficients'][index]
            entries += len(coefficients)
        scalar_cases.append({'n': n, 'ramification': P, 'modes': r, 'geometry': geometry,
                             'common_polynomial_unit_degree': unit.bit_length()-1})
    alphabet = states(2, 2)
    indices = tuple((row, column) for row, column in product(alphabet, repeat=2)
                    if sum(row) == sum(column))
    kernels = []
    for coefficients in product(range(4), repeat=len(indices)):
        kernel = {index: {0: value} for index, value in zip(indices, coefficients) if value}
        try:
            kernel, inverse = normalized_kernel(2, kernel, (0, 1))
        except ValueError:
            continue
        kernels.append(kernel)
    assert len(kernels) == 324
    data = operator_basis(2, 1, 2, 'mixed')
    rng = random.Random(160412)
    targets = {data['gaussian_orders'], (0, 0, 0, 0), (5, 5, 5, 5)}
    targets.update(tuple(rng.randrange(-1, value+2) for value in data['gaussian_orders']) for _ in range(25))
    infinity_targets = (data['infinity_orders'], (0, -3, -3, -4), (0, 1, -3, -4))
    literal = oracle()
    decisions, feasible, witnesses = 0, 0, 0
    for flag in full_flags(2):
        basis = [[int(i == j) for j in range(4)] for i in range(4)]
        for i in range(2):
            for j in range(2):
                basis[i+1][j+1] = flag[i][j]
        basis = tuple(tuple(row) for row in basis)
        evaluated = [literal(kernel, basis) for kernel in kernels]
        for zero, infinity in product(sorted(targets), infinity_targets):
            expected = any(all(a >= b for a, b in zip(zero_actual, zero)) and
                           all(a >= b for a, b in zip(infinity_actual, infinity))
                           for zero_actual, infinity_actual in evaluated)
            actual = synthesize(2, 1, 2, 'mixed', 0, zero, infinity, basis,
                                field_size=4, upper_power=0)
            extended = synthesize(2, 1, 2, 'mixed', 0, zero, infinity, basis,
                                  construct=False, field_size=16, upper_power=0)
            assert (actual is not None) == bool(extended) == expected
            decisions += 1
            if actual is not None:
                feasible += 1
                zero_actual, infinity_actual = literal(actual['kernel'], basis)
                assert all(a >= b for a, b in zip(zero_actual, zero))
                assert all(a >= b for a, b in zip(infinity_actual, infinity))
                witnesses += 1
    result = {'status': 'Exact scalar-source and nonbinary operator-flag controls; no pressure claim.',
              'scalar_Pascal_cases': scalar_cases, 'exact_common_unit_source_entries': entries,
              'complete_F4_constant_tables': 4**len(indices),
              'admissible_filtered_endpoint_kernels': len(kernels),
              'complete_F4_two_dimensional_flags': len(full_flags(2)),
              'literal_constant_matrix_calculations': len(full_flags(2))*len(kernels),
              'F4_and_F16_feasibility_comparisons': decisions,
              'feasible_decisions': feasible, 'independent_literal_witness_checks': witnesses}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'scalar_Pascal_cases'}))


if __name__ == '__main__':
    main()
