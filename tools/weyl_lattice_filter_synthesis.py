#!/usr/bin/env python3
"""Exact polynomial filter synthesis in a fixed reduced Smith basis.

For four binary states, nonzero endpoint coordinate projections suffice
over F4. For eight states the exact criterion is that two binary endpoint
image vectors have supports covering every coordinate. A physical-Y
diagonal map then restores the endpoint at no profile or degree cost.

Prescribed row valuations in a regular Gaussian Smith basis certify every
Smith prefix. They allow independent targets in the reduced occupation
blocks, instead of one common bound on the original active matrix.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from coupled_identity_weyl_profiles import reduced_components, weyl_data
from f4_repaired_identity_filters import MUL, add, evaluate_one, multiply
from prime_state_filter_profiles import validate


def polynomial_order(value):
    if not value:
        raise ValueError('The zero polynomial has no finite order.')
    return (value & -value).bit_length() - 1


def polynomial_multiply(first, second):
    if first.bit_count() > second.bit_count():
        first, second = second, first
    result = 0
    while first:
        low = first & -first
        result ^= second << (low.bit_length() - 1)
        first ^= low
    return result


def polynomial_degrees(value):
    while value:
        low = value & -value
        yield low.bit_length() - 1
        value ^= low


def matrix_multiply(first, second):
    return [[_xor(polynomial_multiply(first[i][k], second[k][j]) for k in range(len(second)))
             for j in range(len(second[0]))] for i in range(len(first))]


def _xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def regular_smith_basis(source):
    """Exact polynomial U,V over F2 with unit determinants and U source V diagonal."""
    size = len(source)
    matrix = [list(row) for row in source]
    left = [[int(i == j) for j in range(size)] for i in range(size)]
    right = [[int(i == j) for j in range(size)] for i in range(size)]
    orders = []
    for k in range(size):
        candidates = [(polynomial_order(matrix[i][j]), i, j)
                      for i in range(k, size) for j in range(k, size) if matrix[i][j]]
        if not candidates:
            raise ArithmeticError('The promised Gaussian block is singular.')
        order, row, column = min(candidates)
        matrix[k], matrix[row] = matrix[row], matrix[k]
        left[k], left[row] = left[row], left[k]
        for values in (matrix, right):
            for line in values:
                line[k], line[column] = line[column], line[k]
        unit = matrix[k][k] >> order
        assert unit & 1
        for i in range(k + 1, size):
            if not matrix[i][k]:
                continue
            assert polynomial_order(matrix[i][k]) >= order
            quotient = matrix[i][k] >> order
            for values in (matrix, left):
                values[i] = [polynomial_multiply(unit, a) ^ polynomial_multiply(quotient, b)
                             for a, b in zip(values[i], values[k])]
            assert matrix[i][k] == 0
        for j in range(k + 1, size):
            if not matrix[k][j]:
                continue
            assert polynomial_order(matrix[k][j]) >= order
            quotient = matrix[k][j] >> order
            for values in (matrix, right):
                for line in values:
                    line[j] = polynomial_multiply(unit, line[j]) ^ polynomial_multiply(quotient, line[k])
            assert matrix[k][j] == 0
        orders.append(order)
    assert all(not matrix[i][j] for i in range(size) for j in range(size) if i != j)
    assert orders == sorted(orders)
    assert matrix_multiply(matrix_multiply(left, source), right) == matrix
    return tuple(orders), left, right, matrix


@lru_cache(maxsize=16)
def lattice_basis(n, ramification=1, geometry='initial', modes=2):
    """Complete diagonal-filter coefficients in the Gaussian Smith basis."""
    validate(2, n, ramification, (0, 0))
    if modes not in (2, 3):
        raise ValueError('This prototype supports two or three binary identity directions.')
    states = tuple(product(range(2), repeat=modes))
    data = weyl_data(2, n, modes, geometry)
    components = reduced_components(2, n, modes, geometry)
    full = []
    for state in states:
        matrix = [[0 for _ in states] for _ in states]
        for (row, column), polynomial in components[state, state].items():
            weight = ramification * sum(j * digit for j, digit in zip(data['selected'], states[row]))
            matrix[row][column] = sum(1 << (ramification * degree + weight)
                                      for degree, coefficient in polynomial.items() if coefficient)
        full.append(matrix)
    transformed = [dict() for _ in states]
    gaussian_orders, occupations, blocks = [], [], []
    row_offset = 0
    for occupation in range(modes + 1):
        indices = [i for i, state in enumerate(states) if sum(state) == occupation]
        gaussian = [[_xor(matrix[i][j] for matrix in full) for j in indices] for i in indices]
        orders, left, right, diagonal = regular_smith_basis(gaussian)
        blocks.append({'occupation': occupation, 'indices': indices, 'gaussian': gaussian,
                       'left': left, 'right': right, 'diagonal': diagonal,
                       'orders': orders, 'row_offset': row_offset})
        for state_index, matrix in enumerate(full):
            value = matrix_multiply(matrix_multiply(left, [[matrix[i][j] for j in indices]
                                                          for i in indices]), right)
            for i, row in enumerate(value):
                for j, coefficient in enumerate(row):
                    if coefficient:
                        transformed[state_index][row_offset + i, j] = coefficient
        gaussian_orders.extend(orders)
        occupations.extend([occupation] * len(indices))
        row_offset += len(indices)
    return {'states': states, 'selected': data['selected'], 'blocks': blocks,
            'coefficients': tuple(transformed), 'gaussian_orders': tuple(gaussian_orders),
            'row_occupations': tuple(occupations)}


@lru_cache(maxsize=16)
def coefficient_tables(n, ramification, geometry, degrees, modes=2):
    data = lattice_basis(n, ramification, geometry, modes)
    size = len(data['states'])
    if len(degrees) != size or any(not isinstance(value, int) or value < 0 for value in degrees):
        raise ValueError('Use one nonnegative degree bound per occupation state.')
    variables = tuple((state, degree) for state, bound in enumerate(degrees)
                      for degree in range(bound + 1))
    equations = {}
    endpoint_masks = [0] * size
    for column, (state, degree) in enumerate(variables):
        bit = 1 << column
        endpoint_masks[state] |= bit
        for (row, col), coefficient in data['coefficients'][state].items():
            for power in polynomial_degrees(coefficient):
                key = row, col, power + degree
                equations[key] = equations.get(key, 0) ^ bit
    rows = tuple(tuple(sorted((power, value) for (row, col, power), value in equations.items()
                              if row == index and value)) for index in range(size))
    return variables, tuple(endpoint_masks), rows


def homogeneous_system(n, ramification, geometry, degrees, row_targets, modes=2):
    if len(row_targets) != 2 ** modes or any(not isinstance(value, int) for value in row_targets):
        raise ValueError('Use one integer target per reduced row.')
    variables, endpoint_masks, rows = coefficient_tables(n, ramification, geometry, tuple(degrees), modes)
    pivots = {}
    constraints = 0
    for target, row in zip(row_targets, rows):
        for power, equation in row:
            if power >= target:
                break
            constraints += 1
            while equation:
                pivot = equation.bit_length() - 1
                if pivot in pivots:
                    equation ^= pivots[pivot]
                else:
                    pivots[pivot] = equation
                    break
    return variables, pivots, endpoint_masks, constraints


def reduce_row(row, pivots):
    for pivot in sorted(pivots, reverse=True):
        if row & (1 << pivot):
            row ^= pivots[pivot]
    return row


def endpoint_image(variables, pivots, endpoint_masks):
    image = {}
    for free in range(len(variables)):
        if free in pivots:
            continue
        lift = 1 << free
        for pivot, equation in sorted(pivots.items()):
            if (lift & equation).bit_count() % 2:
                lift |= 1 << pivot
        endpoint = sum(((lift & mask).bit_count() % 2) << j for j, mask in enumerate(endpoint_masks))
        while endpoint:
            pivot = endpoint.bit_length() - 1
            if pivot in image:
                endpoint ^= image[pivot][0]
                lift ^= image[pivot][1]
            else:
                image[pivot] = endpoint, lift
                break
    return tuple(image.values())


def binary_image_basis(remainders):
    """Image of the endpoint map from its coefficient functionals on free variables."""
    columns = 0
    for row in remainders:
        columns |= row
    basis = {}
    for column in polynomial_degrees(columns):
        value = sum(((row >> column) & 1) << j for j, row in enumerate(remainders))
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return tuple(basis.values())


def ones_endpoint_possible(remainders):
    pivots = {}
    for row in remainders:
        rhs = 1
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                old, value = pivots[pivot]
                row ^= old
                rhs ^= value
            else:
                pivots[pivot] = row, rhs
                break
        else:
            if rhs:
                return False
    return True


def covering_pair(basis, size):
    """Two words of a binary image cover all coordinates iff its F4 extension has full support."""
    values = [0]
    for vector in basis:
        values += [value ^ vector for value in values]
    full = (1 << size) - 1
    if full in values:
        return full, 0
    supersets = [None] * (1 << size)
    for value in values:
        supersets[value] = value
    for bit in range(size):
        for mask in range(1 << size):
            if not mask & (1 << bit) and supersets[mask] is None:
                supersets[mask] = supersets[mask | (1 << bit)]
    return next(((value, supersets[full ^ value]) for value in values
                 if supersets[full ^ value] is not None), None)


def nonzero_endpoint_lift(image, size=4):
    """For four coordinates over F4, grow support with at most three forbidden scalars."""
    if size != 4:
        pair = covering_pair(tuple(vector for vector, lift in image), size)
        if pair is None:
            return None
        lifts = {0: 0}
        for vector, binary_lift in image:
            lifts.update({value ^ vector: lift ^ binary_lift for value, lift in tuple(lifts.items())})
        first, second = pair
        v, w = lifts[first], lifts[second]
        endpoint = [((first >> j) & 1) | (((second >> j) & 1) << 1) for j in range(size)]
        lift = {column: ((v >> column) & 1) | (((w >> column) & 1) << 1)
                for column in polynomial_degrees(v | w)}
        return endpoint, lift
    endpoint = [0] * size
    lift = {}
    for vector, binary_lift in image:
        support = tuple(j for j in range(4) if vector & (1 << j))
        if len(support) == 4:
            return [1] * 4, {column: 1 for column in polynomial_degrees(binary_lift)}
        union = tuple(j for j in range(4) if endpoint[j] or j in support)
        if all(endpoint[j] for j in union):
            continue
        scalar = next(value for value in range(4)
                      if all(endpoint[j] ^ (value if j in support else 0) for j in union))
        for j in support:
            endpoint[j] ^= scalar
        for column in polynomial_degrees(binary_lift):
            lift[column] = lift.get(column, 0) ^ scalar
            if not lift[column]:
                del lift[column]
        if all(endpoint):
            return endpoint, lift
    raise ArithmeticError('The full-coordinate endpoint image did not produce a nonzero F4 endpoint.')


def verify_witness(n, ramification, geometry, degrees, row_targets, alpha, modes=2):
    data = lattice_basis(n, ramification, geometry, modes)
    result = {}
    for state, basis in zip(data['states'], data['coefficients']):
        for entry, coefficient in basis.items():
            polynomial = {degree: 1 for degree in polynomial_degrees(coefficient)}
            result[entry] = add(result.get(entry, {}), multiply(polynomial, alpha[state]))
    assert all(evaluate_one(alpha[state]) and max(alpha[state]) <= bound
               for state, bound in zip(data['states'], degrees))
    assert all(not value or min(value) >= row_targets[row]
               for (row, column), value in result.items())
    return {'checked_transformed_entries': len(result),
            'nonzero_transformed_entries': sum(bool(value) for value in result.values())}


def synthesize(n, ramification, geometry, degrees, row_targets, construct=True, modes=2):
    degrees, row_targets = tuple(degrees), tuple(row_targets)
    variables, pivots, endpoint_masks, constraints = homogeneous_system(
        n, ramification, geometry, degrees, row_targets, modes)
    remainders = tuple(reduce_row(mask, pivots) for mask in endpoint_masks)
    if not all(remainders):
        return None
    if len(endpoint_masks) > 4 and not ones_endpoint_possible(remainders):
        if covering_pair(binary_image_basis(remainders), len(endpoint_masks)) is None:
            return None
    if not construct:
        return True
    image = endpoint_image(variables, pivots, endpoint_masks)
    endpoint, lift = nonzero_endpoint_lift(image, len(endpoint_masks))
    data = lattice_basis(n, ramification, geometry, modes)
    alpha = {state: {} for state in data['states']}
    for column, coefficient in lift.items():
        state, degree = variables[column]
        alpha[data['states'][state]][degree] = coefficient
    result = {'n': n, 'ramification': ramification, 'geometry': geometry,
              'coefficient_field': 'F4', 'modes': modes,
              'degree_bounds': degrees, 'row_targets': row_targets,
              'gaussian_row_orders': data['gaussian_orders'],
              'coefficient_variables': len(variables), 'coefficient_constraints': constraints,
              'homogeneous_rank': len(pivots), 'endpoint_image_dimension': len(image),
              'endpoint_values': endpoint, 'alpha': alpha,
              'actual_degrees': tuple(max(alpha[state]) for state in data['states'])}
    result['witness_verification'] = verify_witness(
        n, ramification, geometry, degrees, row_targets, alpha, modes)
    return result


def bound_blocks(n, ramification, geometry, row_targets, modes=2):
    data = lattice_basis(n, ramification, geometry, modes)
    result = defaultdict(list)
    for auxiliary in data['states']:
        shift = -2 * ramification * sum(j * digit for j, digit in zip(data['selected'], auxiliary))
        for occupation, value in zip(data['row_occupations'], row_targets):
            result[occupation - sum(auxiliary)].append(2 * value + shift)
    return tuple({'charge': charge, 'orders': tuple(sorted(orders))}
                 for charge, orders in sorted(result.items()))


def json_ready(result):
    if result is None:
        return None
    result = dict(result)
    result['alpha'] = [{'state': state, 'polynomial': sorted(value.items())}
                       for state, value in result['alpha'].items()]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=8)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--geometry', default='initial',
                        help='A proved geometry, including split:a:b.')
    parser.add_argument('--degrees', help='One degree cap per state; default is zero.')
    parser.add_argument('--targets', help='One row target per state; default is the Gaussian orders.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    degrees = tuple(map(int, args.degrees.split(','))) if args.degrees else (0,) * (2 ** args.modes)
    targets = (tuple(map(int, args.targets.split(','))) if args.targets else
               lattice_basis(args.n, args.ramification, args.geometry, args.modes)['gaussian_orders'])
    result = synthesize(args.n, args.ramification, args.geometry, degrees, targets, modes=args.modes)
    encoded = json.dumps({'feasible': result is not None, 'witness': json_ready(result)}, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
