#!/usr/bin/env python3
"""Three active modes for even defects on the hottest reflection-logical mode.

Input factor: {physical_mask: {w_exponent: F4_coefficient}}, in x,y,u,v
bit order. The scalar coefficient equals one at w=1; every nonempty
coefficient vanishes there. All Laurent coefficients and spectator
modes are retained exactly. Floating pressure optimization is separate.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

from check_rectangular_f4_repair import fmul
from check_rectangular_defect_support import Rat, Z, ZERO, ONE
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import validate, profile_factors, physical_bit_energies


EVEN_MASKS = (0, 3, 5, 6, 9, 10, 12, 15)


def normalize_factor(factor):
    result = {}
    for mask, polynomial in factor.items():
        mask = int(mask)
        if mask not in EVEN_MASKS:
            raise ValueError('Use even exterior monomials in the four selected physical variables.')
        values = {}
        for exponent, coefficient in polynomial.items():
            exponent = int(exponent)
            if not isinstance(coefficient, int) or coefficient not in (0, 1, 2, 3):
                raise ValueError('Every coefficient must be an F4 element encoded as 0,1,2,3.')
            if coefficient:
                values[exponent] = values.get(exponent, 0) ^ coefficient
        values = {exponent: coefficient for exponent, coefficient in values.items() if coefficient}
        if values:
            result[mask] = values
    result.setdefault(0, {0: 1})
    for mask, polynomial in result.items():
        endpoint = 0
        for coefficient in polynomial.values():
            endpoint ^= coefficient
        if endpoint != int(mask == 0):
            raise ValueError('The exterior factor must equal one coefficientwise at w=1.')
    return result


def add(first, second):
    result = dict(first)
    for degree, coefficient in second.items():
        value = result.get(degree, 0) ^ coefficient
        if value:
            result[degree] = value
        else:
            result.pop(degree, None)
    return result


def multiply(first, second):
    result = {}
    for a, x in first.items():
        for b, y in second.items():
            degree = a + b
            value = result.get(degree, 0) ^ fmul(x, y)
            if value:
                result[degree] = value
            else:
                result.pop(degree, None)
    return result


def rat_coefficients(value):
    """The three-mode images have only monomial denominators."""
    if not value:
        return {}
    if value.d & (value.d - 1):
        raise ArithmeticError('An unexpected unit denominator entered the explicit active images.')
    shift = 1 - value.d.bit_length()
    source, result = value.n, {}
    while source:
        bit = source & -source
        source ^= bit
        result[bit.bit_length() - 1 + shift] = 1
    return result


def active_images(interval, ramification=1):
    validate(interval)
    if not isinstance(ramification, int) or ramification < 1:
        raise ValueError('Use a positive integer substitution degree.')
    p, m = ramification, interval // 2
    c, s = ONE + Z ** (2 * p), Z ** (4 * p)
    geometric = sum((s ** j for j in range(m - 1)), ZERO)
    small, large = Z ** (2 * p), Z ** ((2 * interval + 2) * p)
    empty = (ZERO, ZERO, ZERO)
    return {
        'left_first': ((ZERO, ONE, ZERO), empty),
        'left_second': ((ONE, c, c * s), empty),
        'right_first': (empty, (ONE / small, c / large, c * geometric / large)),
        'right_second': (empty, (ZERO, ONE / large, ZERO)),
    }


def act(matrix, image):
    """Left action of one transported CAR operator on a complete matrix."""
    size = len(matrix)
    result = [[{} for _ in range(size)] for _ in range(size)]
    for creation, coefficients in ((True, image[0]), (False, image[1])):
        for mode, coefficient in enumerate(coefficients):
            polynomial = rat_coefficients(coefficient)
            if not polynomial:
                continue
            bit = 1 << mode
            for source in range(size):
                if bool(source & bit) == creation:
                    continue
                target = source ^ bit
                for column, value in enumerate(matrix[source]):
                    if value:
                        result[target][column] = add(result[target][column], multiply(polynomial, value))
    return result


@lru_cache(maxsize=64)
def components(interval, ramification=1, place='zero'):
    if place not in ('zero', 'infinity'):
        raise ValueError('Choose zero or infinity.')
    images = active_images(interval, ramification)
    orders = (2 * ramification, (2 * interval + 2) * ramification,
              (2 * interval + 6) * ramification)
    size = 8
    initial = [[{} for _ in range(size)] for _ in range(size)]
    for state in range(size):
        initial[state][state] = {sum(order for i, order in enumerate(orders) if state >> i & 1): 1}
    left = ('x', 'y') if place == 'zero' else ('u', 'x')
    right = ('u', 'v') if place == 'zero' else ('v', 'y')
    names = {left[0]: 'left_first', left[1]: 'left_second',
             right[0]: 'right_first', right[1]: 'right_second'}
    answer = {}
    for mask in EVEN_MASKS:
        matrix = initial
        for letter in right + left:
            if mask >> 'xyuv'.index(letter) & 1:
                matrix = act(matrix, images[names[letter]])
        answer[mask] = matrix
    return answer


def compile_interval(interval, place='zero', ramification=1):
    """All eight exact gauged monomial matrices, before scalar polynomials.

    Every matrix entry is a complete {integer_z_power: F4_coefficient}
    Laurent polynomial. Multiplication by w^e means z^(+e) at zero
    and z^(-e) at infinity. No series truncation or unit denominator
    occurs in these compiled matrices.
    """
    base = components(interval, ramification, place)
    gauged = (0, 1) if place == 'zero' else (2, 3)
    output = {}
    for mask, matrix in base.items():
        height = 2 * ramification * sum((mask >> leg) & 1 for leg in gauged)
        output[mask] = [[{degree + height: value for degree, value in polynomial.items()}
                         for polynomial in row] for row in matrix]
    m = interval // 2
    indices = (0, m, m + 1)
    all_modes = tuple(ramification * order for order in profile_factors(interval))
    return {'M': interval, 'place': place, 'ramification': ramification,
            'components': output, 'active_indices': indices,
            'active_orders': tuple(all_modes[i] for i in indices),
            'spectator_orders': tuple(order for i, order in enumerate(all_modes) if i not in indices),
            'scalar_order': (1 if place == 'zero' else -1) * ramification * interval * interval // 2,
            'exact_laurent_matrix': True, 'unit_denominator': 1,
            'physical_monomial_bit_order': 'x,y,u,v', 'logical_core_index': m - 1}


def active_matrix(interval, ramification, place, factor):
    factor = normalize_factor(factor)
    base = components(interval, ramification, place)
    gauged = (0, 1) if place == 'zero' else (2, 3)
    sign = 1 if place == 'zero' else -1
    result = [[{} for _ in range(8)] for _ in range(8)]
    for mask, scalar in factor.items():
        height = 2 * ramification * sum((mask >> leg) & 1 for leg in gauged)
        polynomial = {sign * degree + height: coefficient for degree, coefficient in scalar.items()}
        for i in range(8):
            for j in range(8):
                if base[mask][i][j]:
                    result[i][j] = add(result[i][j], multiply(polynomial, base[mask][i][j]))
    return result


def smith_profile(matrix):
    shifts = [[min(value) if value else 0 for value in row] for row in matrix]
    packed = [[tuple(sum(1 << (degree - shifts[i][j]) for degree, value in entry.items() if value >> bit & 1)
                     for bit in range(2)) for j, entry in enumerate(row)] for i, row in enumerate(matrix)]
    bits = max(1, max(value.bit_length() for row in packed for pair in row for value in pair))
    precision = max(128, 1 << (bits - 1).bit_length())
    if precision > 8192:
        raise ArithmeticError('The complete Laurent source exceeds the native 8192-coefficient limit.')
    while precision <= 8192:
        try:
            result = MatrixOracle(packed, len(matrix), shifts, precision).profile()
            orders = result['orders']
            if any(order != round(order) for order in orders):
                raise ArithmeticError('An integral profile returned a nonintegral order.')
            return tuple(map(int, orders)), {'exact_laurent_matrix': True,
                'source_maximum_width': bits, 'profile_precision': precision,
                'minimum_entry_order': min(min(value) for row in matrix for value in row if value),
                'maximum_entry_degree': max(max(value) for row in matrix for value in row if value)}
        except ArithmeticError:
            precision *= 2
    raise ArithmeticError('The native window did not resolve every promised active pivot.')


def profile_data(interval, ramification, factor):
    factor = normalize_factor(factor)
    validate(interval)
    active_one = (2 * ramification, (2 * interval + 2) * ramification,
                  (2 * interval + 6) * ramification)
    all_modes = [ramification * order for order in profile_factors(interval)]
    spectators = list(all_modes)
    for order in active_one:
        spectators.remove(order)
    clearing = max(0, -min(degree for value in factor.values() for degree in value))
    scalar = ramification * interval * interval // 2 + clearing
    result = {'interval': interval, 'ramification': ramification,
              'logical_core_index': interval // 2 - 1, 'active_one_particle_orders': active_one,
              'spectator_orders': tuple(spectators), 'common_laurent_clearing': clearing,
              'scalar_orders': (scalar, -scalar), 'active_matrix_order': 8,
              'polynomial_degree_bound': 2 * ramification * interval * interval
                    + max(degree for value in factor.values() for degree in value) + clearing,
              'factor': factor, 'profiles': [], 'metadata': []}
    for place in ('zero', 'infinity'):
        orders, metadata = smith_profile(active_matrix(interval, ramification, place, factor))
        result['profiles'].append(orders)
        result['metadata'].append(metadata)
    return result


def pressure(data, log_s):
    if not math.isfinite(log_s) or log_s <= 0:
        raise ValueError('Use a positive finite logarithmic scale.')
    log_b = sum(float(logsumexp(-log_s * np.asarray(profile))) for profile in data['profiles'])
    log_b += 2 * float(np.logaddexp(0, -log_s * np.asarray(data['spectator_orders'])).sum())
    energies = data['ramification'] * np.asarray(physical_bit_energies(data['interval']), dtype=float)
    def equation(rate):
        return log_b - 2 * rate * float(np.logaddexp(0, -log_s * energies / rate).sum())
    lower, upper = 1e-9, 4.0
    if equation(lower) <= 0:
        raise ArithmeticError('A nonpositive pressure root requires separate review.')
    while equation(upper) > 0:
        upper *= 2
    return brentq(equation, lower, upper, xtol=1e-13) / 2


def optimize_scale(data):
    result = minimize_scalar(lambda x: pressure(data, math.exp(x)),
                             bounds=(math.log(1e-5 / data['ramification']),
                                     math.log(10 / data['ramification'])), method='bounded',
                             options={'xatol': 1e-11})
    return {'rho': float(result.fun), 'log_s': math.exp(result.x),
            'evaluations': int(result.nfev), 'success': bool(result.success)}
