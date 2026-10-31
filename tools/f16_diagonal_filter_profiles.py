#!/usr/bin/env python3
"""Complete F16 polynomial profiles for two/three binary identity directions.

All reduced occupation blocks have size at most three.  Exact valuations
of every minor therefore give their Smith orders without a series cutoff.
Polynomial coefficients use four integer bit planes and the fixed field
F16=F4[u]/(u^2+u+omega).
"""
from __future__ import annotations

from itertools import combinations, permutations

from f4_lattice_linear import f16_multiply
from weyl_lattice_filter_synthesis import bound_blocks, lattice_basis, polynomial_multiply


ZERO = (0, 0, 0, 0)
ONE = (1, 0, 0, 0)
MUL = tuple(tuple(f16_multiply(a, b) for b in range(16)) for a in range(16))


def encode(polynomial):
    return tuple(sum(((coefficient >> bit) & 1) << degree
                     for degree, coefficient in polynomial.items()) for bit in range(4))


def add(first, second):
    return tuple(a ^ b for a, b in zip(first, second))


def multiply(first, second):
    result = [0] * 4
    for i, a in enumerate(first):
        if not a:
            continue
        for j, b in enumerate(second):
            if not b:
                continue
            product = polynomial_multiply(a, b)
            coefficient = MUL[1 << i][1 << j]
            for bit in range(4):
                if (coefficient >> bit) & 1:
                    result[bit] ^= product
    return tuple(result)


def order(value):
    bits = value[0] | value[1] | value[2] | value[3]
    return (bits & -bits).bit_length() - 1 if bits else None


def determinant(matrix):
    result = ZERO
    for permutation in permutations(range(len(matrix))):
        term = ONE
        for i, j in enumerate(permutation):
            term = multiply(term, matrix[i][j])
            if term == ZERO:
                break
        result = add(result, term)
    return result


def smith_orders(matrix):
    size = len(matrix)
    prefixes = [0]
    for k in range(1, size + 1):
        values = []
        for rows in combinations(range(size), k):
            for columns in combinations(range(size), k):
                value = order(determinant([[matrix[i][j] for j in columns] for i in rows]))
                if value is not None:
                    values.append(value)
        if not values:
            raise ArithmeticError('A repaired nonzero-endpoint filter must have full reduced rank.')
        prefixes.append(min(values))
    result = tuple(b - a for a, b in zip(prefixes, prefixes[1:]))
    assert result == tuple(sorted(result))
    return result


def reduced_orders(n, ramification, modes, alpha, geometry='initial'):
    base = lattice_basis(n, ramification, geometry, modes)
    polynomials = tuple(encode(alpha[state]) for state in base['states'])
    result = []
    for block in base['blocks']:
        size, offset = len(block['indices']), block['row_offset']
        matrix = [[ZERO for _ in range(size)] for _ in range(size)]
        for coefficients, polynomial in zip(base['coefficients'], polynomials):
            for i in range(size):
                for j in range(size):
                    coefficient = coefficients.get((offset + i, j), 0)
                    if coefficient:
                        matrix[i][j] = add(matrix[i][j], multiply((coefficient, 0, 0, 0), polynomial))
        result.extend(smith_orders(matrix))
    return tuple(result)


def zero_profiles(n, ramification, modes, alpha, geometry='initial'):
    values = reduced_orders(n, ramification, modes, alpha, geometry)
    return bound_blocks(n, ramification, geometry, values, modes)


def verify_flag_witness(witness):
    from weyl_flag_filter_synthesis import flag_basis
    n, P, modes, geometry = (witness[key] for key in ('n', 'ramification', 'modes', 'geometry'))
    data = flag_basis(n, P, geometry, modes, witness['flags'])
    matrix = {}
    for state, coefficients in zip(data['base']['states'], data['coefficients']):
        polynomial = encode(witness['alpha'][state])
        for entry, coefficient in coefficients.items():
            matrix[entry] = add(matrix.get(entry, ZERO), multiply(coefficient + (0, 0), polynomial))
    assert all(order(value) is None or order(value) >= witness['row_targets'][row]
               for (row, column), value in matrix.items())
    orders = reduced_orders(n, P, modes, witness['alpha'], geometry)
    for block in data['base']['blocks']:
        start, size = block['row_offset'], len(block['indices'])
        targets = sorted(witness['row_targets'][start:start + size])
        assert all(sum(orders[start:start + k]) >= sum(targets[:k]) for k in range(1, size + 1))
    return {'transformed_entries': len(matrix), 'exact_reduced_orders': orders}
