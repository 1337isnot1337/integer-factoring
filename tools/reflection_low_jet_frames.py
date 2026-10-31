#!/usr/bin/env python3
"""Uniform q-by-q reflection jets for complete coupled-pole admissibility.

Polynomial entries are binary coefficient bit masks, reduced modulo
w**cutoff.  The physical alphabet has q=2**(M/2) letters.  Only selected
half-size minors are constructed; the full q**4 crossing is not needed.
"""
from __future__ import annotations

from reflection_dual_rail_gaussian import reflection_letters, validate
from square_zero_coupled_reflection_frames import product, square_zero_matrix


def binary_inverse(matrix):
    size = len(matrix)
    rows = [sum(value << column for column, value in enumerate(row)) | (1 << (size + i))
            for i, row in enumerate(matrix)]
    for column in range(size):
        selected = next((i for i in range(column, size) if rows[i] >> column & 1), None)
        if selected is None:
            raise ArithmeticError('The constant selected-minor matrix is singular.')
        rows[column], rows[selected] = rows[selected], rows[column]
        for i in range(size):
            if i != column and rows[i] >> column & 1:
                rows[i] ^= rows[column]
    return [[rows[i] >> (size + j) & 1 for j in range(size)] for i in range(size)]


def left_binary(binary, polynomials):
    size = len(binary)
    output = [[0] * size for _ in range(size)]
    for i, row in enumerate(binary):
        for k, value in enumerate(row):
            if value:
                output[i] = [a ^ b for a, b in zip(output[i], polynomials[k])]
    return output


def right_binary(polynomials, binary):
    # Transposition changes a right product into a left product.
    transposed = left_binary(list(map(list, zip(*binary))),
                             list(map(list, zip(*polynomials))))
    return list(map(list, zip(*transposed)))


def selected_power_minor(rows, columns, cutoff):
    """det(w**(2*(row&column))) by subset DP, in characteristic two."""
    count = len(rows)
    limit = (1 << cutoff) - 1
    values = [0] * (1 << count)
    values[0] = 1
    for subset in range(len(values) - 1):
        current = values[subset]
        if not current:
            continue
        row = rows[subset.bit_count()]
        remaining = ((1 << count) - 1) ^ subset
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            column = columns[bit.bit_length() - 1]
            exponent = 2 * (row & column)
            if exponent < cutoff:
                values[subset | bit] ^= (current << exponent) & limit
    return values[-1]


def reflection_jets(interval, cutoff=None):
    """Construct B=Q(0)^-1 Q and its inverse modulo w**(2M).

    Literal x,y,u,v indexing has Q=rho((1+w)H_M), R=rho((1+w)G_M).
    Since H=J G J, the two selected-minor matrices differ by complementing
    both physical bit labels.  Each minor gains (1+w)**(M/2)=1+w**(M/2).
    """
    validate(interval)
    cutoff = 2 * interval if cutoff is None else cutoff
    if not isinstance(cutoff, int) or not 1 <= cutoff <= 2 * interval:
        raise ValueError('Use an integer cutoff between one and 2M.')
    half = interval // 2
    letters = [[i for i in range(interval) if mask >> i & 1]
               for mask in reflection_letters(interval)]
    size = len(letters)
    limit = (1 << cutoff) - 1
    r_g = [[0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i, size):
            value = selected_power_minor(letters[i], letters[j], cutoff)
            value = (value ^ (value << half)) & limit
            r_g[i][j] = r_g[j][i] = value
    opposite = size - 1
    q_h = [[r_g[i ^ opposite][j ^ opposite] for j in range(size)] for i in range(size)]
    q0 = [[value & 1 for value in row] for row in q_h]
    q0_inverse = binary_inverse(q0)
    if q0_inverse != [[value & 1 for value in row] for row in r_g]:
        raise ArithmeticError('The two constant physical block factors are not inverses.')
    b = left_binary(q0_inverse, q_h)
    inverse_b = right_binary(r_g, q0)
    if [[value & 1 for value in row] for row in b] != [[int(i == j) for j in range(size)] for i in range(size)]:
        raise ArithmeticError('The normalized physical matrix has the wrong constant jet.')
    return {'interval': interval, 'alphabet': size, 'cutoff': cutoff,
            'Q_H': q_h, 'R_G': r_g, 'Q0': q0, 'Q0_inverse': q0_inverse,
            'B': b, 'inverse_B': inverse_b}


def matrix_coefficient(polynomials, exponent):
    if exponent < 0:
        return [[0] * len(polynomials) for _ in polynomials]
    return [[value >> exponent & 1 for value in row] for row in polynomials]


def add_matrices(first, second):
    return [[a ^ b for a, b in zip(row_a, row_b)] for row_a, row_b in zip(first, second)]


def pole_conditions(jets, generator, maximum_pole=None):
    """Necessary and sufficient matrix-jet tests for every k<=M/2.

    Returns each nonzero negative coefficient of
    B+w^-k(GB+BG)+w^-2k GBG.  The theorem identifies precisely these
    obstructions with polynomiality of the original four-leg crossing.
    """
    interval, size = jets['interval'], jets['alphabet']
    maximum_pole = interval // 2 if maximum_pole is None else maximum_pole
    if not isinstance(maximum_pole, int) or not 1 <= maximum_pole <= interval // 2:
        raise ValueError('The uniform criterion applies to 1<=k<=M/2.')
    if jets['cutoff'] < 2 * maximum_pole:
        raise ValueError('The supplied physical jet cutoff is too short.')
    if (len(generator) != size
            or any(len(row) != size or any(value not in range(4) for value in row) for row in generator)):
        raise ValueError('Use an F4 matrix on the whole physical alphabet.')
    if not any(value for row in generator for value in row):
        return [{'pole': pole, 'polynomial': True, 'negative_coefficients': {}}
                for pole in range(1, maximum_pole + 1)]
    if not square_zero_matrix(generator):
        raise ValueError('The full physical generator must satisfy G squared = zero.')
    zero = [[0] * size for _ in range(size)]
    commutators, quadratics = [], []
    for exponent in range(2 * maximum_pole):
        coefficient = matrix_coefficient(jets['B'], exponent)
        left = product(generator, coefficient)
        right = product(coefficient, generator)
        commutators.append(add_matrices(left, right))
        quadratics.append(product(left, generator))
    results = []
    for pole in range(1, maximum_pole + 1):
        obstructions = {}
        for exponent in range(-2 * pole, 0):
            linear = commutators[exponent + pole] if exponent + pole >= 0 else zero
            quadratic = quadratics[exponent + 2 * pole]
            value = add_matrices(linear, quadratic)
            if any(entry for row in value for entry in row):
                obstructions[exponent] = value
        results.append({'pole': pole, 'polynomial': not obstructions,
                        'negative_coefficients': obstructions})
    return results
