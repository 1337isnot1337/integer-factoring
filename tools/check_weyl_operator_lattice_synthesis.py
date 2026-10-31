#!/usr/bin/env python3
"""Compare operator-lattice feasibility with every constant F4 two-mode kernel."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import random

from coupled_identity_operator_filters import identity_kernel, states
from coupled_identity_weyl_profiles import reduced_components, reduced_matrix, weyl_data
from f4_identity_operator_filters import normalized_kernel
from f4_repaired_identity_filters import MUL, add, multiply
from weyl_lattice_filter_synthesis import polynomial_degrees, regular_smith_basis
from weyl_operator_lattice_synthesis import operator_basis, synthesize


def matrix_multiply(first, second):
    result = [[{} for _ in second[0]] for _ in first]
    for i, row in enumerate(first):
        for k, value in enumerate(row):
            for j, other in enumerate(second[k]):
                result[i][j] = add(result[i][j], multiply(value, other))
    return result


def oracle():
    """Literal weighted Weyl matrices, assembled without the synthesis module."""
    n, P, modes, geometry = 2, 1, 2, 'mixed'
    alphabet = states(2, modes)
    data = weyl_data(2, n, modes, geometry)
    components = reduced_components(2, n, modes, geometry)
    gaussian = reduced_matrix(2, n, P, modes, identity_kernel(2, modes), geometry)
    transforms = []
    for occupation in range(3):
        indices = [i for i, state in enumerate(alphabet) if sum(state) == occupation]
        block = [[sum(1 << degree for degree, value in gaussian[i][j].items() if value)
                  for j in indices] for i in indices]
        orders, left, right, diagonal = regular_smith_basis(block)
        transforms.append((indices,
                           [[{degree: 1 for degree in polynomial_degrees(value)} for value in row] for row in left],
                           [[{degree: 1 for degree in polynomial_degrees(value)} for value in row] for row in right]))
    kappa = tuple(sum((2*j+1-2*n)*digit for j, digit in zip(data['selected'], state))
                  for state in alphabet)

    def evaluate(kernel, infinity_basis):
        matrix = [[{} for _ in alphabet] for _ in alphabet]
        for index, scalar in kernel.items():
            for (row, column), value in components[index].items():
                weight = sum(j*digit for j, digit in zip(data['selected'], alphabet[row]))
                value = {degree+weight: coefficient for degree, coefficient in value.items()}
                matrix[row][column] = add(matrix[row][column], multiply(scalar, value))
        zero = []
        for indices, left, right in transforms:
            block = [[matrix[i][j] for j in indices] for i in indices]
            block = matrix_multiply(matrix_multiply(left, block), right)
            zero.extend(min((degree for value in row for degree in value), default=10**6) for row in block)
        infinity = []
        for row in infinity_basis:
            result = []
            for j, column in enumerate(alphabet):
                polynomial = {}
                for i, coefficient in enumerate(row):
                    if coefficient:
                        polynomial = add(polynomial, kernel.get((alphabet[i], column), {}), coefficient)
                result.extend(kappa[j]-2*degree for degree in polynomial)
            infinity.append(min(result, default=10**6))
        return tuple(zero), tuple(infinity)
    return evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
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
    direct = oracle()
    data = operator_basis(2, 1, 2, 'mixed')
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    bases = [identity]
    for i, j in ((1, 2), (2, 1)):
        basis = [list(row) for row in identity]
        basis[i][j] ^= 1
        bases.append(tuple(tuple(row) for row in basis))
    rng = random.Random(281077)
    targets = {data['gaussian_orders'], (0, 0, 0, 0), (5, 5, 5, 5)}
    targets.update(tuple(rng.randrange(-1, bound+2) for bound in data['gaussian_orders'])
                   for _ in range(40))
    infinity_targets = (data['infinity_orders'], (0, -3, -3, -4), (0, 1, -3, -4))
    decisions, feasible, witnesses = 0, 0, 0
    for basis in bases:
        evaluated = [(kernel, direct(kernel, basis)) for kernel in kernels]
        for zero_target, infinity_target in product(sorted(targets), infinity_targets):
            expected = any(all(a >= b for a, b in zip(zero, zero_target)) and
                           all(a >= b for a, b in zip(infinity, infinity_target))
                           for kernel, (zero, infinity) in evaluated)
            actual = synthesize(2, 1, 2, 'mixed', 0, zero_target, infinity_target,
                                basis, upper_power=0)
            assert (actual is not None) == expected
            decisions += 1
            if actual is not None:
                feasible += 1
                zero, infinity = direct(actual['kernel'], basis)
                assert all(a >= b for a, b in zip(zero, zero_target))
                assert all(a >= b for a, b in zip(infinity, infinity_target))
                assert all(max(value) == min(value) == 0 for value in actual['kernel'].values())
                witnesses += 1
    # A Laurent scalar is admissible at both cuts with precisely opposite shifts.
    laurent = 0
    for n, P, modes, geometry in ((2, 1, 2, 'mixed'), (3, 3, 3, 'split:1:2'),
                                  (32, 1, 3, 'split:1:2')):
        data = operator_basis(n, P, modes, geometry)
        for shift in (-2, 0, 2):
            actual = synthesize(n, P, modes, geometry, shift,
                                tuple(value+shift for value in data['gaussian_orders']),
                                tuple(value-2*shift for value in data['infinity_orders']),
                                upper_power=shift)
            assert actual is not None
            assert actual['kernel'] == {index: {shift: 1} for index in identity_kernel(2, modes)}
            laurent += 1
    result = {'status': 'Exact feasibility comparisons against every admissible constant two-mode F4 operator; no pressure claim.',
              'constant_coefficients_enumerated': 4**len(indices),
              'admissible_filtered_endpoint_kernels': len(kernels),
              'literal_transformed_matrices': len(kernels)*len(bases),
              'exact_feasibility_decisions': decisions, 'feasible_decisions': feasible,
              'independent_witness_verifications': witnesses,
              'Laurent_and_large_interval_controls': laurent}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
