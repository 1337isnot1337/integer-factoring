#!/usr/bin/env python3
"""Exact physical split maps and exploratory pressure for F4 crossings.

Each polynomial is represented by two binary bit planes, with omega^2=omega+1.
The endpoint is checked before either full-rank Smith oracle is constructed.
The two physical alphabets and their charged heights can be different.
"""

from itertools import product

from native_matrix_smith_f4 import MatrixOracle


F4_PRODUCT = (
    (0, 0, 0, 0),
    (0, 1, 2, 3),
    (0, 2, 3, 1),
    (0, 3, 1, 2),
)


def field_dot(left, right):
    result = 0
    for a, b in zip(left, right):
        result ^= F4_PRODUCT[a][b]
    return result


def scalar_multiply(value, scalar):
    low, high = value
    return ((0, 0), (low, high), (high, low ^ high), (low ^ high, low))[scalar]


def linear_combination(terms):
    low = high = 0
    for value, scalar in terms:
        a, b = scalar_multiply(value, scalar)
        low ^= a
        high ^= b
    return low, high


def check_split(left, right):
    """right is stored by its columns, so both maps are lists of vectors."""
    if not left or len(left) != len(right):
        raise ValueError("Supply nonempty left rows and right columns of equal count.")
    ambient = len(left[0])
    if any(len(row) != ambient for row in left + right):
        raise ValueError("Every split vector must have the same ambient dimension.")
    if any(value not in range(4) for row in left + right for value in row):
        raise ValueError("Field coefficients must be encoded by 0,1,2,3.")
    for i, row in enumerate(left):
        for j, column in enumerate(right):
            if field_dot(row, column) != int(i == j):
                raise ValueError("The physical maps do not satisfy L R = I.")
    return ambient


def transform_leg(values, dimensions, axis, matrix):
    new_dimensions = list(dimensions)
    new_dimensions[axis] = len(matrix)
    result = []
    for reversed_coordinates in product(*(range(d) for d in new_dimensions[::-1])):
        coordinates = list(reversed_coordinates[::-1])
        output = coordinates[axis]
        terms = []
        for source, scalar in enumerate(matrix[output]):
            if not scalar:
                continue
            coordinates[axis] = source
            index, stride = 0, 1
            for coordinate, dimension in zip(coordinates, dimensions):
                index += coordinate * stride
                stride *= dimension
            terms.append((values[index], scalar))
        result.append(linear_combination(terms))
    return result, new_dimensions


def retract(seed, left_x, right_x, left_y=None, right_y=None):
    """K' = (Lx,Ly,Ry^T,Rx^T) K, with all endpoint coefficients checked."""
    if left_y is None:
        left_y, right_y = left_x, right_x
    if right_y is None:
        raise ValueError("Both maps of the second split are required.")
    ambient_x = check_split(left_x, right_x)
    ambient_y = check_split(left_y, right_y)
    expected = [ambient_x, ambient_y, ambient_y, ambient_x]
    if list(seed['dimensions']) != expected:
        raise ValueError("The ambient crossing dimensions do not match the split maps.")
    values, dimensions = list(seed['values']), expected
    for axis, matrix in enumerate((left_x, left_y, right_y, right_x)):
        values, dimensions = transform_leg(values, dimensions, axis, matrix)
    qx, qy = len(left_x), len(left_y)
    for index, value in enumerate(values):
        x = index % qx
        y = index // qx % qy
        u = index // (qx * qy) % qy
        v = index // (qx * qy * qy)
        endpoint = (value[0].bit_count() % 2, value[1].bit_count() % 2)
        if endpoint != (int(x == v and y == u), 0):
            raise ArithmeticError("A complete retracted crossing endpoint entry failed.")
    return {
        'values': values, 'dimensions': dimensions, 'q_x': qx, 'q_u': qy,
        'degree': max((part.bit_length() - 1 for value in values for part in value), default=0),
        'denominator_degree': 0, 'strict_endpoint': True,
    }


def flatten(seed, cut):
    qx, qy = seed['q_x'], seed['q_u']
    rows = list(product(range(qx), range(qy)))
    columns = list(product(range(qy), range(qx)))
    matrix = []
    for a, b in rows:
        row = []
        for c, d in columns:
            x, y, u, v = (a, b, c, d) if cut == 0 else (a, c, b, d)
            row.append(seed['values'][x + qx * y + qx * qy * u + qx * qy * qy * v])
        matrix.append(row)
    return matrix, rows, columns


class SplitProfile:
    def __init__(self, seed):
        self.seed = seed
        self.matrices = [flatten(seed, cut) for cut in (0, 1)]
        order = seed['q_x'] * seed['q_u']
        self.oracles = [MatrixOracle(matrix, order) for matrix, _, _ in self.matrices]

    def profiles(self, energy_x, energy_y=None):
        if energy_y is None:
            energy_y = energy_x
        if len(energy_x) != self.seed['q_x'] or len(energy_y) != self.seed['q_u']:
            raise ValueError("One charged energy is required for each physical symbol.")
        result = []
        for cut, (_, rows, columns) in enumerate(self.matrices):
            row_weights = [(energy_x[a] + energy_y[b]) / 2 if cut == 0
                           else energy_y[b] / 2 for a, b in rows]
            column_weights = [0 if cut == 0 else energy_x[b] / 2 for _, b in columns]
            result.append(self.oracles[cut].profile(row_weights, column_weights,
                                                    infinity=bool(cut))['orders'])
        return result

    def pressure(self, energy_x, energy_y=None, profiles=None):
        """Numerical scalar root and scale search; this is not a strict certificate."""
        import numpy as np
        from scipy.optimize import brentq, minimize_scalar
        from scipy.special import logsumexp

        if energy_y is None:
            energy_y = energy_x
        if profiles is None:
            profiles = self.profiles(energy_x, energy_y)
        nu, mu = (np.asarray(profile, float) for profile in profiles)
        ex, ey = np.asarray(energy_x, float), np.asarray(energy_y, float)

        def rho(log_scale):
            constant = (logsumexp(-log_scale * nu) + logsumexp(-log_scale * mu)
                        + log_scale * (min(ex) + min(ey)))
            gx, gy = log_scale * (ex - min(ex)), log_scale * (ey - min(ey))
            if constant < 1e-12:
                return 1.0

            def difference(rate):
                return constant - rate * (logsumexp(-gx / rate) + logsumexp(-gy / rate))

            high = 2.0
            while difference(high) > 0:
                high *= 2
            return brentq(difference, 1e-6, high, xtol=1e-10) / 2

        optimum = minimize_scalar(rho, bounds=(1e-5, 2), method='bounded',
                                  options={'xatol': 1e-8})
        return {'rho': float(optimum.fun), 'log_scale': float(optimum.x),
                'profiles': profiles}


def identity_split(dimension):
    return [[int(i == j) for j in range(dimension)] for i in range(dimension)]


def pair_split(dimension, pair, left_pair, right_pair):
    a, b = pair
    alphabet = [i for i in range(dimension) if i != b]
    maps = []
    for coefficients in (left_pair, right_pair):
        matrix = []
        for symbol in alphabet:
            row = [0] * dimension
            if symbol == a:
                row[a], row[b] = coefficients
            else:
                row[symbol] = 1
            matrix.append(row)
        maps.append(matrix)
    check_split(*maps)
    return alphabet, *maps


def hyperplane_split(functional, complement):
    if field_dot(functional, complement) != 1:
        raise ValueError("The functional must take value one on its complement.")
    pivot = next(i for i, value in enumerate(functional) if value)
    if functional[pivot] != 1:
        raise ValueError("Normalize the first nonzero functional coefficient to one.")
    dimension = len(functional)
    alphabet = [i for i in range(dimension) if i != pivot]
    right = [[int(i == c) ^ (functional[c] if i == pivot else 0)
              for i in range(dimension)] for c in alphabet]
    left = [[int(i == c) ^ F4_PRODUCT[complement[c]][functional[i]]
             for i in range(dimension)] for c in alphabet]
    check_split(left, right)
    return alphabet, left, right
