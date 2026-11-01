#!/usr/bin/env python3
"""Convolution-free exact rational bounds for finite typed adaptive flags.

Every local estimate is a rational coefficient bound, rounded upward only.
The input local Smith profiles retain the independent proof obligations of
finite_typed_adaptive_certificate. No floating rank estimate is accepted.
"""
from __future__ import annotations

from fractions import Fraction
from heapq import heappop, heappush

from finite_typed_adaptive_certificate import amplification_bounds
from weighted_adaptive_mixed_place import type_class_size


def rational_jet_costs(profile, copies, total_degree, parameter, cut_dimension):
    """Bound [x^Q]P(x)^m/(1-x)^2 at one exact rational 0<parameter<1."""
    parameter = Fraction(parameter)
    if (not isinstance(copies, int) or copies < 1 or not isinstance(total_degree, int)
            or total_degree < 0 or not isinstance(cut_dimension, int) or cut_dimension < 1
            or not 0 < parameter < 1):
        raise ValueError('Use positive length/dimension, nonnegative degree, and a rational parameter in (0,1).')
    shift = copies * Fraction(profile.mean_shift)
    if shift.denominator != 1:
        raise ValueError('The tensor length must clear the exact type mean.')
    minimum = min(profile.orders)
    offset = int(shift) - copies * minimum
    first = max(0, (-offset + profile.ramification - 1) // profile.ramification)
    if first > total_degree:
        return [0] * (total_degree + 1)
    polynomial = sum((parameter ** (order - minimum) for order in profile.orders), Fraction(0))
    power = polynomial ** copies / (1 - parameter) ** 2
    index = profile.ramification * first + offset
    value = power / parameter ** index
    numerator, denominator = value.numerator, value.denominator
    step_numerator = parameter.denominator ** profile.ramification
    step_denominator = parameter.numerator ** profile.ramification
    output = [0] * first
    for j in range(first, total_degree + 1):
        if numerator >= cut_dimension * denominator:
            output.extend([cut_dimension] * (total_degree - j + 1))
            break
        output.append((numerator + denominator - 1) // denominator)
        numerator *= step_numerator
        denominator *= step_denominator
    if len(output) != total_degree + 1 or output != sorted(output):
        raise ArithmeticError('The exact rational local sequence is not monotone or complete.')
    return output


def rational_catalogue_flag(catalogues, copies, total_degree, cut_dimension):
    """Each place has (GaugeProfile, positive rational parameter list) pairs."""
    if not catalogues or any(not place for place in catalogues):
        raise ValueError('Both sides require nonempty local gauge catalogues.')
    place_costs, place_choices = [], []
    for place in catalogues:
        choices, alternatives = [], []
        for gauge_index, (profile, parameters) in enumerate(place):
            if not parameters:
                raise ValueError('Every gauge requires a nonempty rational parameter list.')
            for parameter_index, parameter in enumerate(parameters):
                alternatives.append(rational_jet_costs(profile, copies, total_degree, parameter, cut_dimension))
                choices.append((gauge_index, parameter_index))
        selected = [min((values[j], i) for i, values in enumerate(alternatives))
                    for j in range(total_degree + 1)]
        place_costs.append([pair[0] for pair in selected])
        place_choices.append([choices[pair[1]] for pair in selected])
    queue = []
    for place, values in enumerate(place_costs):
        heappush(queue, (values[0], place, 0))
    prefix, flag = 0, []
    for _ in range(total_degree + 1):
        cost, place, order = heappop(queue)
        gauge, parameter = place_choices[place][order]
        prefix = min(cut_dimension, prefix + cost)
        flag.append({'place_index': place, 'order': order, 'gauge_index': gauge,
                     'parameter_index': parameter, 'cost': cost, 'prefix_cost': prefix})
        if order < total_degree:
            heappush(queue, (place_costs[place][order + 1], place, order + 1))
    return flag


def finite_rational_certificate(row_catalogues, column_catalogues, copies, degree,
                                xcounts, ycounts, exponent=(140, 99), precision=64):
    if (not isinstance(copies, int) or copies < 1 or not isinstance(degree, int) or degree < 0
            or sum(xcounts) != copies or sum(ycounts) != copies):
        raise ValueError('Use positive length, nonnegative degree, and common-length types.')
    dx, dy = type_class_size(xcounts), type_class_size(ycounts)
    if min(dx, dy) < 2:
        raise ValueError('Both finite external label sets must contain at least two words.')
    channels = len(xcounts) * len(ycounts)
    if any(len(profile.orders) != channels
           for side in (row_catalogues, column_catalogues) for place in side for profile, _ in place):
        raise ValueError('The local profiles must contain every generic cut channel.')
    total_degree, dimension = copies * degree, dx * dy
    row = rational_catalogue_flag(row_catalogues, copies, total_degree, dimension)
    column = rational_catalogue_flag(column_catalogues, copies, total_degree, dimension)
    pairs = [(row[k]['prefix_cost'], column[total_degree - k]['prefix_cost'])
             for k in range(total_degree + 1)]
    gate = amplification_bounds(pairs, dimension, exponent, precision)
    return {'copies': copies, 'seed_degree': degree, 'tensor_degree': total_degree,
            'type_counts_x': list(xcounts), 'type_counts_y': list(ycounts),
            'external_label_counts': [dx, dy], 'cut_dimension': dimension,
            'rank_pairs': pairs, 'amplification': gate, 'row_flag': row, 'column_flag': column}
