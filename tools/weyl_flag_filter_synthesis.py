#!/usr/bin/env python3
"""Diagonal polynomial synthesis after constant F4 row flags.

The reference Gaussian Smith bases come from weyl_lattice_filter_synthesis.
Each occupation block may be multiplied on the left by an invertible
constant F4 matrix.  Row-divisibility targets still bound every Smith
prefix.  A fixed F16 coefficient extension supplies nonzero endpoints
whenever no endpoint coordinate is identically zero on the solution space.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product

import f4_lattice_linear as linear
from weyl_lattice_filter_synthesis import lattice_basis, polynomial_degrees


def dot_vector(first, second):
    value = 0
    for a, b in zip(first, second):
        value ^= linear.MUL[a][b]
    return value


def matrix_rank(matrix):
    return len(linear.eliminate([linear.encode(row) for row in matrix]))


@lru_cache(maxsize=3)
def full_flags(size):
    """Represent all full flags in F4^size, for size at most three.

    Rows are ordered low, middle, high.  The last row spans the final
    line; the last two rows span its plane.  Flags, rather than all bases,
    suffice when target row orders are nondecreasing in this order.
    """
    if size == 1:
        return (((1,),),)
    points = tuple(vector for vector in product(range(4), repeat=size)
                   if any(vector) and next(value for value in vector if value) == 1)
    result = []
    if size == 2:
        for high in points:
            low = next(value for value in points if matrix_rank((value, high)) == 2)
            result.append((low, high))
    elif size == 3:
        for high in points:
            for normal in points:
                if dot_vector(high, normal):
                    continue
                middle = next(value for value in points
                              if not dot_vector(value, normal) and matrix_rank((value, high)) == 2)
                low = next(value for value in points if dot_vector(value, normal))
                result.append((low, middle, high))
    else:
        raise ValueError('This prototype enumerates flags only through dimension three.')
    assert all(matrix_rank(matrix) == size for matrix in result)
    assert len(result) == {2: 5, 3: 105}[size]
    return tuple(result)


def identity_flag(size):
    matrix = tuple(tuple(int(i == j) for j in range(size)) for i in range(size))
    return full_flags(size).index(matrix)


@lru_cache(maxsize=48)
def flag_basis(n, ramification=1, geometry='mixed', modes=2, flags=None):
    base = lattice_basis(n, ramification, geometry, modes)
    if flags is None:
        flags = tuple(identity_flag(len(block['indices'])) for block in base['blocks'])
    if len(flags) != len(base['blocks']):
        raise ValueError('Give one flag index per occupation block.')
    transformed = [dict() for _ in base['states']]
    for block, flag_index in zip(base['blocks'], flags):
        size, offset = len(block['indices']), block['row_offset']
        flag = full_flags(size)[flag_index]
        for state_index, coefficients in enumerate(base['coefficients']):
            for i in range(size):
                for j in range(size):
                    value = 0, 0
                    for k in range(size):
                        term = coefficients.get((offset + k, j), 0)
                        if term and flag[i][k]:
                            value = linear.add(value, linear.scale((term, 0), flag[i][k]))
                    if value[0] | value[1]:
                        transformed[state_index][offset + i, j] = value
    minimum = tuple(min(min(polynomial_degrees(low | high))
                        for coefficients in transformed
                        for (row, column), (low, high) in coefficients.items() if row == index)
                    for index in range(len(base['states'])))
    gaussian_rows = []
    for block, flag_index in zip(base['blocks'], flags):
        flag = full_flags(len(block['indices']))[flag_index]
        gaussian_rows.extend(min(order for order, coefficient in zip(block['orders'], row) if coefficient)
                             for row in flag)
    return {'base': base, 'flags': tuple(flags), 'coefficients': tuple(transformed),
            'minimum_rows': minimum, 'gaussian_rows': tuple(gaussian_rows)}


@lru_cache(maxsize=24)
def coefficient_tables(n, ramification, geometry, modes, flags, degrees):
    data = flag_basis(n, ramification, geometry, modes, flags)
    size = len(data['base']['states'])
    if len(degrees) != size or min(degrees) < 0:
        raise ValueError('Give one nonnegative degree bound per state.')
    variables = tuple((state, degree) for state, bound in enumerate(degrees)
                      for degree in range(bound + 1))
    equations = {}
    masks = [0] * size
    for column, (state, degree) in enumerate(variables):
        bit = 1 << column
        masks[state] |= bit
        for (row, col), (low, high) in data['coefficients'][state].items():
            for power in polynomial_degrees(low | high):
                value = ((bit if (low >> power) & 1 else 0),
                         (bit if (high >> power) & 1 else 0))
                key = row, col, power + degree
                equations[key] = linear.add(equations.get(key, (0, 0)), value)
    rows = tuple(tuple(sorted((power, value) for (row, col, power), value in equations.items()
                              if row == index and (value[0] | value[1]))) for index in range(size))
    return variables, tuple(masks), rows


def homogeneous_system(n, ramification, geometry, degrees, targets, flags, modes=2):
    variables, masks, rows = coefficient_tables(n, ramification, geometry, modes, tuple(flags), tuple(degrees))
    equations = [equation for target, row in zip(targets, rows)
                 for power, equation in row if power < target]
    pivots = linear.eliminate(equations)
    remainders = linear.endpoint_remainders(masks, pivots)
    return variables, masks, equations, pivots, remainders


def synthesize(n, ramification, geometry, degrees, targets, flags, modes=2,
               construct=True, field_size=16):
    data = flag_basis(n, ramification, geometry, modes, tuple(flags))
    variables, masks, equations, pivots, remainders = homogeneous_system(
        n, ramification, geometry, tuple(degrees), tuple(targets), tuple(flags), modes)
    if any(not (low | high) for low, high in remainders):
        return None
    if not construct and field_size >= len(remainders):
        return True
    lift = linear.full_support_lift(remainders, pivots, field_size)
    if lift is None:
        return None
    if not construct:
        return True
    endpoints, assignment, dimension = lift
    assert linear.verify_assignment(equations, assignment)
    actual_endpoint = [0] * len(masks)
    alpha = {state: {} for state in data['base']['states']}
    for column, coefficient in assignment.items():
        state, degree = variables[column]
        alpha[data['base']['states'][state]][degree] = coefficient
        actual_endpoint[state] ^= coefficient
    assert tuple(actual_endpoint) == endpoints and all(endpoints)
    return {'n': n, 'ramification': ramification, 'geometry': geometry, 'modes': modes,
            'coefficient_field': f'F{field_size}', 'flags': tuple(flags),
            'degree_bounds': tuple(degrees), 'row_targets': tuple(targets),
            'gaussian_row_orders': data['base']['gaussian_orders'],
            'endpoint_values': endpoints, 'endpoint_image_dimension': dimension,
            'coefficient_variables': len(variables), 'homogeneous_rank': len(pivots),
            'coefficient_equations': len(equations), 'alpha': alpha,
            'actual_degrees': tuple(max(alpha[state]) for state in data['base']['states'])}
