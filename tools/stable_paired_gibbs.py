#!/usr/bin/env python3
"""Numerically stable paired Gibbs roots and envelope gradients.

Native selected prefixes remain floating-weight search outputs. Their
orders are reconstructed as exact affine forms in the binary-exact input
parameters, and necessary consistency checks are enforced. This does not
prove that all selected minors attain the global Smith profile.

No negative leading gap is clipped. Tiny positive partition tails and
entropy factors are retained logarithmically. The compatible fallback
adapter leaves the historical ordinary evaluator unchanged.
"""
from __future__ import annotations

import math


LOG_TWO = math.log(2.)


def log_sum(values):
    values = tuple(values)
    if not values:
        return -math.inf
    maximum = max(values)
    if maximum == -math.inf:
        return maximum
    return maximum+math.log(math.fsum(math.exp(value-maximum) for value in values))


def log_softplus(value):
    """log(log(1+exp(value))), without underflowing its small positive value."""
    if value == -math.inf:
        return value
    if value < -36:
        # The relative correction to exp(value) is below double precision.
        return value
    if value > 36:
        return math.log(value+math.log1p(math.exp(-value)))
    return math.log(math.log1p(math.exp(value)))


def partition_shape(gaps, rate=1.):
    """Return log(log Z), log Z, and all log probabilities for min-shifted data."""
    assert gaps and min(gaps) == 0 and rate > 0
    first = gaps.index(0)
    exponents = [-value/rate for value in gaps]
    log_tail = log_sum(value for index, value in enumerate(exponents) if index != first)
    log_log_z = log_softplus(log_tail)
    log_z = math.exp(log_log_z) if log_log_z > -746 else 0.
    return log_log_z, log_z, tuple(value-log_z for value in exponents)


def subtract_forms(first, second):
    result = dict(first)
    for index, value in second:
        result[index] = result.get(index, 0)-value
        if not result[index]:
            del result[index]
    return tuple(sorted(result.items()))


def add_forms(forms):
    result = {}
    for factor, form in forms:
        for index, value in form:
            result[index] = result.get(index, 0)+factor*value
            if not result[index]:
                del result[index]
    return tuple(sorted(result.items()))


def exact_affine_data(oracle, parameters, cut):
    if oracle.denominator_degree:
        raise ValueError('This stable endpoint-bound check currently requires a polynomial strict crossing.')
    params = tuple(map(float, parameters))
    prefixes = oracle.prefixes(params, cut)
    ratios = [value.as_integer_ratio() for value in params]
    denominator = max(second for _, second in ratios)
    units = [first*(denominator//second) for first, second in ratios]
    n, a, b = oracle.n, oracle.A, oracle.B
    order_forms, order_units = [], []
    for side in (0, 1):
        selected = cut if side == 0 else 1-cut
        previous = 0
        forms, values = [], []
        for index in range(n):
            position = side*n+index
            row, column = prefixes['rows'][position], prefixes['cols'][position]
            x, first = divmod(row, b)
            second, v = divmod(column, a)
            labels = (x, first, second, v) if selected == 0 else (x, second, first, v)
            base = prefixes['base_prefix'][position]-previous
            previous = prefixes['base_prefix'][position]
            form = [(len(params)-1, base)] if base else []
            form += [(side*oracle.stride+oracle.offsets[leg]+label-1, 1)
                     for leg, label in enumerate(labels) if label]
            form = tuple(sorted(form))
            value = sum(coefficient*units[at] for at, coefficient in form)
            reconstructed = value/denominator
            observed = prefixes['orders'][position]
            if abs(reconstructed-observed) > 2e-8*max(1., abs(reconstructed), abs(observed)):
                raise ArithmeticError('A native selected order disagrees with its exact affine reconstruction.')
            forms.append(form)
            values.append(value)
        if values != sorted(values):
            raise ArithmeticError('Exactly reconstructed selected orders are not nondecreasing.')
        order_forms.append(tuple(forms))
        order_units.append(tuple(values))
    energy_forms, energy_units = [], []
    for size, legs in ((a, (0, 3)), (b, (1, 2))):
        forms = [tuple((side*oracle.stride+oracle.offsets[leg]+label-1, 1)
                       for side in (0, 1) for leg in legs) if label else ()
                 for label in range(size)]
        values = [sum(coefficient*units[at] for at, coefficient in form) for form in forms]
        energy_forms.append(tuple(forms))
        energy_units.append(tuple(values))
    minima = [values.index(min(values)) for values in energy_units]
    delta = sum(min(values) for values in energy_units)-sum(values[0] for values in order_units)
    if delta < 0:
        raise ArithmeticError('Exactly reconstructed selected minima violate the polynomial endpoint bound.')
    delta_form = add_forms([(1, energy_forms[side][minima[side]]) for side in (0, 1)]+
                          [(-1, order_forms[side][0]) for side in (0, 1)])
    assert delta == sum(coefficient*units[at] for at, coefficient in delta_form)
    return {'denominator': denominator, 'parameter_units': tuple(units),
            'order_forms': tuple(order_forms), 'order_units': tuple(order_units),
            'energy_forms': tuple(energy_forms), 'energy_units': tuple(energy_units),
            'energy_minimum_indices': tuple(minima), 'delta_units': delta, 'delta_form': delta_form,
            'native_selected_prefixes': prefixes}


def stable_paired_result(oracle, parameters, cut=1):
    """Return a positive root, log types, and a stable gradient when representable."""
    from scipy.optimize import brentq
    data = exact_affine_data(oracle, parameters, cut)
    denominator = data['denominator']
    order_gaps = [tuple((value-values[0])/denominator for value in values)
                  for values in data['order_units']]
    energy_gaps = [tuple((value-min(values))/denominator for value in values)
                   for values in data['energy_units']]
    profile_shapes = [partition_shape(gaps) for gaps in order_gaps]
    delta = data['delta_units']
    log_delta = math.log(delta)-math.log(denominator) if delta else -math.inf
    log_constant = log_sum([log_delta]+[shape[0] for shape in profile_shapes])
    if not math.isfinite(log_constant):
        raise ArithmeticError('The positive leading-gap partition remainder is not numerically representable in log space.')

    def equation(log_rate):
        rate = math.exp(log_rate)
        shapes = [partition_shape(gaps, rate) for gaps in energy_gaps]
        return log_constant-log_rate-log_sum(shape[0] for shape in shapes)

    low, high = math.log(1e-10), math.log(8.)
    while equation(low) <= 0:
        low -= 10
        if low < -700:
            raise ArithmeticError('The positive Gibbs root lies below the safe floating log-rate range.')
    while equation(high) > 0:
        high += LOG_TWO
        if high > math.log(1e8):
            raise ArithmeticError('The positive Gibbs root exceeds the configured numerical range.')
    log_rate = brentq(equation, low, high, xtol=2e-14, rtol=2e-14)
    rate = math.exp(log_rate)
    type_shapes = [partition_shape(gaps, rate) for gaps in energy_gaps]
    log_entropies = []
    for gaps, shape in zip(energy_gaps, type_shapes):
        terms = [shape[0]]
        for gap, log_probability in zip(gaps, shape[2]):
            if gap:
                terms.append(log_probability+math.log(gap)-log_rate)
        log_entropies.append(log_sum(terms))
    log_entropy = log_sum(log_entropies)
    if not math.isfinite(log_entropy):
        raise ArithmeticError('The positive physical entropy is not numerically representable in log space.')

    contributions = [[] for _ in parameters]
    gradient_status = 'finite_log_scaled_envelope_gradient'
    try:
        if data['delta_form']:
            factor = math.exp(-LOG_TWO-log_entropy)
            for index, coefficient in data['delta_form']:
                contributions[index].append(coefficient*factor)
        for side in (0, 1):
            first = data['energy_minimum_indices'][side]
            for form, log_probability in zip(data['energy_forms'][side], type_shapes[side][2]):
                difference = subtract_forms(form, data['energy_forms'][side][first])
                if difference:
                    factor = math.exp(log_probability-LOG_TWO-log_entropy)
                    for index, coefficient in difference:
                        contributions[index].append(coefficient*factor)
            for form, log_probability in zip(data['order_forms'][side], profile_shapes[side][2]):
                difference = subtract_forms(form, data['order_forms'][side][0])
                if difference:
                    factor = math.exp(log_probability-LOG_TWO-log_entropy)
                    for index, coefficient in difference:
                        contributions[index].append(-coefficient*factor)
        gradient = [math.fsum(values) for values in contributions]
        if not all(math.isfinite(value) for value in gradient):
            raise OverflowError('The stable gradient exceeded the finite double range.')
    except OverflowError:
        gradient = None
        gradient_status = 'positive_root_retained_gradient_outside_float_range'
    return {'rho': rate/2, 'gradient': gradient, 'gradient_status': gradient_status,
            'type_x': [math.exp(value) for value in type_shapes[0][2]],
            'type_y': [math.exp(value) for value in type_shapes[1][2]],
            'log_type_x': list(type_shapes[0][2]), 'log_type_y': list(type_shapes[1][2]),
            'log_paired_entropy': log_entropy, 'log_zero_rate_constant': log_constant,
            'delta_numerator': str(delta), 'delta_denominator': str(denominator),
            'exact_affine_consistency': {'orders_nondecreasing': True, 'delta_nonnegative': True,
                                         'endpoint_minimum_bound': True},
            'orders_times_log_s': [[value/denominator for value in row] for row in data['order_units']],
            'physical_energies_times_log_s': [[value/denominator for value in row] for row in data['energy_units']],
            'scope': 'Stable floating evaluation of selected native profiles with exact affine necessary checks; not a proof of global profile attainment or a certified pressure sign.'}


def stable_paired_value_gradient(oracle, parameters, cut=1):
    """Drop-in four-item result for the historical unconstrained objective."""
    import numpy as np
    result = stable_paired_result(oracle, parameters, cut)
    if result['gradient'] is None:
        raise ArithmeticError('The Gibbs root is positive but its envelope gradient exceeds the floating range.')
    return result['rho'], np.asarray(result['gradient']), result['type_x'], result['type_y']


def paired_value_gradient_with_fallback(oracle, parameters, cut=1, entropy_switch=1e-7):
    """Replay root guards and small-entropy results through the stable evaluator."""
    from search_rectangular_f4_pressure import paired_value_gradient
    try:
        result = paired_value_gradient(oracle, parameters, cut)
        entropy = math.fsum(-value*math.log(value) for value in result[2]+result[3] if value > 0)
        if not math.isfinite(entropy) or entropy < entropy_switch:
            return stable_paired_value_gradient(oracle, parameters, cut)
        return result
    except ArithmeticError as error:
        if str(error) not in ('Nonpositive paired Gibbs root needs exact investigation.',
                              'Type entropy is too small for stable optimization.'):
            raise
        return stable_paired_value_gradient(oracle, parameters, cut)
