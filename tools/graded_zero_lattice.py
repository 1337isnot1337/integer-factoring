#!/usr/bin/env python3
"""An explicit graded lattice change for the retained five-letter crossing.

The resulting polynomial kernel has a stronger zero profile and a worse
opposite profile. No improved factoring bound follows from this constructor.
"""
from __future__ import annotations

from functools import reduce

from joint_energy_conjugations import IDENTITY, N, binary_product, model
from tensor_pascal_box_crossings import matrix_product
from weighted_pascal_flags import binary_inverse


def coefficient(data, degree):
    return sum(((data['operator'][r][c][0] >> (degree-data['energies'][r])) & 1)
               << (N*r+c) for r in range(N) for c in range(N)
               if degree >= data['energies'][r])


def ramify(polynomial, denominator):
    return sum(1 << (denominator*i) for i in range(polynomial.bit_length())
               if polynomial >> i & 1)


def construct():
    data = model()
    leading, following = coefficient(data, 8), coefficient(data, 9)
    assert binary_product(leading, leading) == leading
    complement = IDENTITY ^ leading
    normal = binary_product(binary_product(complement, following), complement)
    pivot_r, pivot_c = next((r, c) for r in range(N) for c in range(N)
                            if normal >> (N*r+c) & 1)
    vector = [(normal >> (N*r+pivot_c)) & 1 for r in range(N)]
    covector = [(normal >> (N*pivot_r+c)) & 1 for c in range(N)]
    assert normal == sum((vector[r] & covector[c]) << (N*r+c)
                         for r in range(N) for c in range(N))
    assert binary_product(normal, normal) == 0
    basis = [[int(i == j) for j in range(N)] for i in range(N)]
    potentials, changed_blocks = [0]*N, []
    for block in data['blocks']:
        if not any(vector[i] for i in block):
            continue
        support = [i for i in block if vector[i]]
        assert len(support) == 2 and all(vector[i] == covector[i] for i in block)
        columns = [vector, [int(i == support[0]) for i in range(N)]]
        columns += [[int(i == j) for i in range(N)] for j in block if j not in support]
        assert len(columns) == len(block)
        for j, column in zip(block, columns):
            for i in block:
                basis[i][j] = column[i]
        potentials[block[0]], potentials[block[1]] = 1, -1
        changed_blocks.append({'energy': data['energies'][block[0]], 'indices': block,
                               'normal_vector_support': support})
    inverse = binary_inverse(basis)
    assert matrix_product(inverse, basis) == [[int(i == j) for j in range(N)] for i in range(N)]
    assert all(not basis[i][j] or data['energies'][i] == data['energies'][j]
               for i in range(N) for j in range(N))
    original = [[value[0] for value in row] for row in data['operator']]
    changed = matrix_product(matrix_product(inverse, original), basis)
    denominator, clearing = 2, 2
    matrix = [[ramify(value, denominator) << (clearing+potentials[r]-potentials[c])
               for c, value in enumerate(row)] for r, row in enumerate(changed)]
    values = []
    for index in range(625):
        x, y, u, v = ((index // 5**i) % 5 for i in range(4))
        value = matrix[x+5*y][v+5*u]
        assert value.bit_count() % 2 == int(x == v and y == u)
        values.append((value, 0))
    return {'values': values, 'basis_columns': basis, 'inverse_basis': inverse,
            'ramified_joint_potentials': potentials, 'ramification': denominator,
            'common_monomial_clearing': clearing, 'changed_energy_blocks': changed_blocks,
            'normal_coefficient_rank': 1, 'normal_coefficient_square_zero': True,
            'polynomial_degree': max(value[0].bit_length()-1 for value in values)}


def positive_dress(values, generator, power):
    """Polynomial units preserve the whole zero profile at the fixed heights."""
    from cube_f4_bivariate_smith import F4_MUL
    from f4_pascal_flags import matrix_product as field_product

    if type(power) is not int or power < 0 or len(generator) != N:
        raise ValueError('Use a nonnegative power and one generator coordinate per joint state.')
    if any(len(row) != N or any(type(value) is not int or value not in range(4)
                               for value in row) for row in generator):
        raise ValueError('Use a square F4 generator with entries 0,1,2,3.')
    energies = model()['energies']
    assert all(not value or energies[i] == energies[j]
               for i, row in enumerate(generator) for j, value in enumerate(row))
    assert all(not reduce(int.__xor__,
               (F4_MUL[generator[i][k]][generator[k][j]] for k in range(N)), 0)
               for i in range(N) for j in range(N))
    matrix = [[(0, 0)] * N for _ in range(N)]
    for index, value in enumerate(values):
        x, y, u, v = ((index // 5**i) % 5 for i in range(4))
        matrix[x+5*y][v+5*u] = tuple(value)
    field_generator = [[(value & 1, value >> 1) for value in row] for row in generator]
    left, right = field_product(field_generator, matrix), field_product(matrix, field_generator)
    sandwich = field_product(left, field_generator)
    output = []
    for index in range(625):
        x, y, u, v = ((index // 5**i) % 5 for i in range(4))
        r, c = x+5*y, v+5*u
        entry = tuple(matrix[r][c][p] ^ ((left[r][c][p] ^ right[r][c][p]) << power)
                      ^ (sandwich[r][c][p] << (2*power)) for p in range(2))
        endpoint = entry[0].bit_count() % 2 | ((entry[1].bit_count() % 2) << 1)
        assert endpoint == int(x == v and y == u)
        output.append(entry)
    return output
