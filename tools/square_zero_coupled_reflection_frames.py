#!/usr/bin/env python3
"""Complete coupled physical frames for arbitrary square-zero F4 matrices.

The coefficient action uses every matrix entry and every interaction among
the four physical legs. Polynomiality is checked only after one common
Laurent shift. The first-jet lemma certifies pole one for the checked
reflection seeds; larger poles are accepted only by the exact full test.
"""
from __future__ import annotations

from physical_split_crossings import F4_PRODUCT
from reflection_laurent_physical_frames import (
    add, bell_normalize, bell_restore, degree, scale, valuation,
)


def product(first, second):
    rows, inner, columns = len(first), len(second), len(second[0])
    if any(len(row) != inner for row in first):
        raise ValueError('Incompatible finite-field matrix dimensions.')
    result = [[0] * columns for _ in range(rows)]
    for i in range(rows):
        for k in range(inner):
            if first[i][k]:
                for j in range(columns):
                    result[i][j] ^= F4_PRODUCT[first[i][k]][second[k][j]]
    return result


def transpose(matrix):
    return list(map(list, zip(*matrix)))


def square_zero_matrix(matrix):
    """Return whether a nonzero square F4 matrix satisfies G squared = zero."""
    q = len(matrix)
    if not q or any(len(row) != q or any(value not in range(4) for value in row) for row in matrix):
        return False
    return any(value for row in matrix for value in row) and not any(value for row in product(matrix, matrix) for value in row)


def rank(matrix):
    work = [list(row) for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        selected = next((i for i in range(pivot_row, len(work)) if work[i][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        inverse = (0, 1, 3, 2)[work[pivot_row][column]]
        work[pivot_row] = [F4_PRODUCT[inverse][value] for value in work[pivot_row]]
        for i in range(pivot_row + 1, len(work)):
            coefficient = work[i][column]
            if coefficient:
                work[i] = [a ^ F4_PRODUCT[coefficient][b] for a, b in zip(work[i], work[pivot_row])]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def matching_generator(permutation):
    """G=I+P for an involutive letter permutation, with optional fixed points."""
    q = len(permutation)
    if sorted(permutation) != list(range(q)) or any(permutation[permutation[i]] != i for i in range(q)):
        raise ValueError('Use an involutive permutation of the physical alphabet.')
    matrix = [[int(i == j) ^ int(i == permutation[j]) for j in range(q)] for i in range(q)]
    if not square_zero_matrix(matrix):
        raise ValueError('The identity permutation gives a trivial generator.')
    return matrix


def apply_matrix(values, q, leg, matrix):
    columns = [[(i, matrix[i][j]) for i in range(q) if matrix[i][j]] for j in range(q)]
    stride = q ** leg
    output = {}
    for index, value in values.items():
        source = index // stride % q
        for target, coefficient in columns[source]:
            destination = index + (target - source) * stride
            changed = add(output.get(destination, (0, 0)), scale(value, coefficient))
            if changed != (0, 0):
                output[destination] = changed
            elif destination in output:
                del output[destination]
    return output


def matrix_terms(bell_seed, generator):
    """The five coefficients in epsilon for A_x A^-T_y A_u A^-T_v."""
    q = bell_seed['q_x']
    if bell_seed['q_u'] != q or len(generator) != q or not square_zero_matrix(generator):
        raise ValueError('Use a nonzero square-zero generator on the complete common alphabet.')
    transposed = transpose(generator)
    terms = [{i: value for i, value in enumerate(bell_seed['values']) if value != (0, 0)}]
    terms += [{} for _ in range(4)]
    for leg, matrix in enumerate((generator, transposed, generator, transposed)):
        for order in range(leg, -1, -1):
            changed = apply_matrix(terms[order], q, leg, matrix)
            for index, value in changed.items():
                updated = add(terms[order + 1].get(index, (0, 0)), value)
                if updated != (0, 0):
                    terms[order + 1][index] = updated
                elif index in terms[order + 1]:
                    del terms[order + 1][index]
    return terms


def from_terms(seed, terms, pole=1):
    if not isinstance(pole, int) or pole < 1:
        raise ValueError('Use a positive integer pole order.')
    offset, mask = 4 * pole, (1 << (4 * pole)) - 1
    values = [(0, 0)] * len(seed['values'])
    for order, coefficients in enumerate(terms):
        for index, (a, b) in coefficients.items():
            before = values[index]
            shift = offset - order * pole
            values[index] = before[0] ^ (a << shift), before[1] ^ (b << shift)
    if any((a | b) & mask for a, b in values):
        raise ArithmeticError('The complete square-zero coupled frame retains an actual pole.')
    values = [(a >> offset, b >> offset) for a, b in values]
    q = seed['q_x']
    for index, (a, b) in enumerate(values):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        if (a.bit_count() % 2, b.bit_count() % 2) != (int(x == v and y == u), 0):
            raise ArithmeticError('The complete square-zero frame changed the strict endpoint.')
    actual_degree = max(degree(value) for value in values)
    if actual_degree != seed['degree']:
        raise ArithmeticError('A negative-power frame unexpectedly changed the leading coefficient degree.')
    return dict(seed, values=values, degree=actual_degree,
                monomial_content=min(valuation(value) for value in values if value != (0, 0)))


def construct(seed, generator, coordinates='bell', restore=True, pole=1, bell_data=None):
    """Build a complete polynomial candidate from an exact F4 generator.

    coordinates='original' specifies G on the original x leg. The Bell
    generator is Q^-1 G Q. With restore=True, returned values use the
    original physical coordinates. Precomputed (bell_seed, factors) may be
    passed through bell_data for a search that retains one source seed.
    """
    if not square_zero_matrix(generator):
        raise ValueError('The supplied generator is not a nonzero square-zero F4 matrix.')
    bell, factors = bell_normalize(seed) if bell_data is None else bell_data
    if coordinates == 'bell':
        bell_generator = [list(row) for row in generator]
    elif coordinates == 'original':
        bell_generator = product(product(factors['inverse_factor'], generator), factors['first_factor'])
    else:
        raise ValueError('Choose bell or original generator coordinates.')
    terms = matrix_terms(bell, bell_generator)
    result = from_terms(bell, terms, pole)
    content = result['monomial_content']
    if restore:
        result = bell_restore(result, factors)
        result['monomial_content'] = content
    result['square_zero_frame'] = {'generator': generator, 'generator_coordinates': coordinates,
                                   'bell_generator': bell_generator, 'rank': rank(generator),
                                   'pole': pole, 'restored_coordinates': restore}
    return result


def first_jet_identity(seed, bell_data=None):
    """Verify M=I+w²(ad_H+cI)+w⁴R in the third flattening exactly."""
    bell, factors = bell_normalize(seed) if bell_data is None else bell_data
    q = bell['q_x']

    def coefficient(x, y, u, v, exponent):
        a, b = bell['values'][x + q * (y + q * (u + q * v))]
        return ((a >> exponent) & 1) + 2 * ((b >> exponent) & 1)

    scalar = coefficient(0, 0, 0, 0, 2)
    operator = [[coefficient(x, y, 0, 0, 2) ^ (scalar if x == y else 0)
                 for y in range(q)] for x in range(q)]
    for x in range(q):
        for y in range(q):
            for u in range(q):
                for v in range(q):
                    if coefficient(x, y, u, v, 0) != int(x == y and u == v):
                        raise ArithmeticError('The constant third-cut coefficient is not identity.')
                    if coefficient(x, y, u, v, 1) or coefficient(x, y, u, v, 3):
                        raise ArithmeticError('The first or third coefficient is not zero.')
                    expected = ((operator[x][y] if u == v else 0)
                                ^ (operator[u][v] if x == y else 0)
                                ^ (scalar if x == y and u == v else 0))
                    if coefficient(x, y, u, v, 2) != expected:
                        raise ArithmeticError('The second third-cut coefficient is not ad_H+cI.')
    return {'alphabet': q, 'degree': seed['degree'], 'operator': operator, 'scalar': scalar,
            'exact_coefficient_entries': 4 * q ** 4, 'bell_factors': factors}
