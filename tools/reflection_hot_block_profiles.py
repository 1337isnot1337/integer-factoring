#!/usr/bin/env python3
"""At most 3k active modes for a power-of-two hot reflection-logical block.

The practical core=2 case has six modes and 64-by-64 matrices at both
places, independent of the ambient interval. All coefficients below are
complete Laurent polynomials over F4, with no series truncation.
"""
from __future__ import annotations

from functools import lru_cache

from check_rectangular_defect_support import Z, ZERO, ONE
from reflection_dual_rail_gaussian import validate, profile_factors
from reflection_hot_core_profiles import add, multiply, act, smith_profile, pressure, optimize_scale


def validate_core(interval, core):
    validate(interval)
    if not isinstance(core, int) or core < 1 or core & (core - 1) or core > interval // 2:
        raise ValueError('Use a power-of-two number of hot logical modes, at most M/2.')


def even_masks(core):
    return tuple(mask for mask in range(1 << (4 * core)) if mask.bit_count() % 2 == 0)


def normalize_factor(factor, core=2):
    result = {}
    for raw_mask, polynomial in factor.items():
        mask = int(raw_mask)
        if mask < 0 or mask >= 1 << (4 * core) or mask.bit_count() % 2:
            raise ValueError('Use even exterior monomials on the four physical hot blocks.')
        values = {}
        for raw_exponent, coefficient in polynomial.items():
            exponent = int(raw_exponent)
            if not isinstance(coefficient, int) or coefficient not in (0, 1, 2, 3):
                raise ValueError('Use F4 coefficients encoded as 0,1,2,3.')
            if coefficient:
                values[exponent] = values.get(exponent, 0) ^ coefficient
        values = {degree: value for degree, value in values.items() if value}
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


def active_images(interval, ramification=1, core=2):
    validate_core(interval, core)
    if not isinstance(ramification, int) or ramification < 1:
        raise ValueError('Use a positive integer substitution degree.')
    p, m, k = ramification, interval // 2, core
    groups = m // k
    retained = min(m, 2 * k)
    total = k + retained
    c, s = ONE + Z ** (2 * p), Z ** (4 * p)
    large = Z ** ((2 * interval + 2) * p)
    geometric = sum(((s ** k) ** h for h in range(groups - 1)), ZERO)
    images = {}
    def vector():
        return [ZERO] * total
    for j in range(k):
        creation, annihilation = vector(), vector()
        creation[k + j] = ONE
        images['left_first', j] = (creation, annihilation)
        creation, annihilation = vector(), vector()
        creation[j] = ONE
        for t in range(k):
            coefficient = c * s ** (t - (j & t))
            creation[k + t] = coefficient
            if groups >= 2:
                creation[2 * k + t] = coefficient * s ** k
        images['left_second', j] = (creation, annihilation)
        creation, annihilation = vector(), vector()
        annihilation[j] = Z ** (-(2 + 4 * j) * p)
        for t in range(k):
            coefficient = c * s ** (-(j & t)) / large
            annihilation[k + t] = coefficient
            if groups >= 2:
                annihilation[2 * k + t] = coefficient * geometric
        images['right_first', j] = (creation, annihilation)
        creation, annihilation = vector(), vector()
        annihilation[k + j] = ONE / (large * s ** j)
        images['right_second', j] = (creation, annihilation)
    indices = tuple(range(k)) + tuple(m + i for i in range(retained))
    orders = tuple(ramification * profile_factors(interval)[i] for i in indices)
    return images, indices, orders


@lru_cache(maxsize=2048)
def component(interval, ramification, place, core, mask):
    if place not in ('zero', 'infinity'):
        raise ValueError('Choose zero or infinity.')
    if mask not in even_masks(core):
        raise ValueError('Use an even physical monomial mask.')
    images, indices, orders = active_images(interval, ramification, core)
    size = 1 << len(indices)
    if size > 1024:
        raise ValueError('The dense native profile path currently supports at most ten active modes.')
    matrix = [[{} for _ in range(size)] for _ in range(size)]
    for state in range(size):
        matrix[state][state] = {sum(order for i, order in enumerate(orders) if state >> i & 1): 1}
    left = ('x', 'y') if place == 'zero' else ('u', 'x')
    right = ('u', 'v') if place == 'zero' else ('v', 'y')
    names = {left[0]: 'left_first', left[1]: 'left_second',
             right[0]: 'right_first', right[1]: 'right_second'}
    for letter in right + left:
        for j in range(core):
            if mask >> ('xyuv'.index(letter) * core + j) & 1:
                matrix = act(matrix, images[names[letter], j])
    gauged = (0, 1) if place == 'zero' else (2, 3)
    height = ramification * sum((2 + 4 * j) for leg in gauged for j in range(core)
                               if mask >> (leg * core + j) & 1)
    return [[{degree + height: value for degree, value in polynomial.items()}
             for polynomial in row] for row in matrix]


def compile_interval(interval, place='zero', ramification=1, core=2, masks=None):
    validate_core(interval, core)
    _, indices, orders = active_images(interval, ramification, core)
    all_modes = tuple(ramification * order for order in profile_factors(interval))
    selected = even_masks(core) if masks is None else tuple(masks)
    return {'M': interval, 'place': place, 'ramification': ramification, 'core': core,
            'components': {mask: component(interval, ramification, place, core, mask) for mask in selected},
            'active_indices': indices, 'active_orders': orders,
            'spectator_orders': tuple(order for i, order in enumerate(all_modes) if i not in indices),
            'scalar_order': (1 if place == 'zero' else -1) * ramification * interval * interval // 2,
            'exact_laurent_matrix': True, 'unit_denominator': 1,
            'physical_monomial_bit_order': 'x hot block, y hot block, u hot block, v hot block',
            'logical_core_indices': tuple(interval // 2 - 1 - j for j in range(core))}


def active_matrix(interval, ramification, place, factor, core=2):
    factor = normalize_factor(factor, core)
    data = compile_interval(interval, place, ramification, core, factor)
    size = 1 << len(data['active_indices'])
    result = [[{} for _ in range(size)] for _ in range(size)]
    sign = 1 if place == 'zero' else -1
    for mask, scalar in factor.items():
        polynomial = {sign * degree: value for degree, value in scalar.items()}
        matrix = data['components'][mask]
        for i in range(size):
            for j in range(size):
                if matrix[i][j]:
                    result[i][j] = add(result[i][j], multiply(polynomial, matrix[i][j]))
    return result


def profile_data(interval, ramification, factor, core=2):
    factor = normalize_factor(factor, core)
    data = compile_interval(interval, 'zero', ramification, core, ())
    clearing = max(0, -min(degree for value in factor.values() for degree in value))
    scalar = ramification * interval * interval // 2 + clearing
    result = {'interval': interval, 'ramification': ramification, 'core': core,
              'logical_core_indices': data['logical_core_indices'],
              'active_one_particle_orders': data['active_orders'],
              'spectator_orders': data['spectator_orders'],
              'common_laurent_clearing': clearing, 'scalar_orders': (scalar, -scalar),
              'active_matrix_order': 1 << len(data['active_indices']),
              'polynomial_degree_bound': 2 * ramification * interval * interval
                    + max(degree for value in factor.values() for degree in value) + clearing,
              'factor': factor, 'profiles': [], 'metadata': []}
    for place in ('zero', 'infinity'):
        orders, metadata = smith_profile(active_matrix(interval, ramification, place, factor, core))
        result['profiles'].append(orders)
        result['metadata'].append(metadata)
    return result
