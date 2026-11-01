#!/usr/bin/env python3
"""Implicit residue-count ranges and charged inverse precision on odd domains."""

from pole_cauchy_analytic_kernels import compile_regular_odd_kernel
from pole_cauchy_analytic_kernels import inverse_polynomial_kernel_contraction
from pole_cauchy_polynomial_solve import _validate


def residue_count_range(domain, level):
    """Min/max counts among all odd classes modulo 2^level, without listing them."""
    if not isinstance(level, int) or level < 1:
        raise ValueError('Use a positive integer residue level.')
    if not domain:
        raise ValueError('Use a nonempty odd-progression domain.')
    modulus, baseline, events, previous = 1 << (level - 1), 0, {}, 0

    def add_interval(start, end):
        if start < end:
            events[start] = events.get(start, 0) + 1
            events[end] = events.get(end, 0) - 1

    for first, length in domain:
        if any(not isinstance(x, int) for x in (first, length)) or first <= previous or first % 2 == 0 or length < 1:
            raise ValueError('Use ordered disjoint positive odd progressions.')
        previous = first + 2 * (length - 1)
        quotient, remainder = divmod(length, modulus)
        baseline += quotient
        start = ((first - 1) // 2) % modulus
        end = start + remainder
        add_interval(start, min(end, modulus))
        if end > modulus:
            add_interval(0, end - modulus)
    events.setdefault(0, 0)
    events.setdefault(modulus, 0)
    points = sorted(events)
    current, minimum, maximum = baseline, None, None
    for index, position in enumerate(points[:-1]):
        current += events[position]
        if position < points[index + 1]:
            minimum = current if minimum is None else min(minimum, current)
            maximum = current if maximum is None else max(maximum, current)
    return {'minimum': minimum, 'maximum': maximum, 'event_positions': len(points),
            'residue_classes_enumerated': 0}


def inverse_precision_certificate(n, left, right, k):
    """A precision bound polynomial in k and the explicit progression count."""
    if not right:
        raise ValueError('Use nonempty domains.')
    size = _validate(n, left, right, right[0][0], k, 1)
    levels, loss = [], 0
    for level in range(1, k + 1):
        row = residue_count_range(left, level)
        column = residue_count_range(right, level)
        imbalance = max(row['maximum'] - column['minimum'],
                        column['maximum'] - row['minimum'])
        contribution = max(0, 2 * imbalance - 1)
        loss += contribution
        levels.append({'level': level, 'row_range': [row['minimum'], row['maximum']],
                       'column_range': [column['minimum'], column['maximum']],
                       'imbalance_bound': imbalance, 'loss_contribution': contribution,
                       'event_positions': row['event_positions'] + column['event_positions']})
    balanced = all(row['row_range'][1] - row['row_range'][0] <= 1
                   and row['column_range'][1] - row['column_range'][0] <= 1 for row in levels)
    if balanced:
        expected = max(0, k - 1 - ((size & -size).bit_length() - 1))
        if loss != expected:
            raise ArithmeticError('The balanced-size precision identity failed.')
    representation_bound = (k - 1) * max(0, 2 * len(left) + 2 * len(right) - 3)
    if loss > representation_bound:
        raise ArithmeticError('The progression-count precision bound failed.')
    return {'domain_size': size, 'inverse_loss_bound': loss, 'balanced_through_displacement': balanced,
            'representation_loss_bound': representation_bound,
            'levels': levels, 'enumerated_matrix_nodes': 0}


def inverse_certified_analytic_kernel(n, left, right, k, circuit, absolute_precision,
                                      *, denominator_bits=0, maximum_compilation_precision=None):
    """Contract a regular joint circuit on certified progression unions.

    A user-supplied cap can reject a large certified precision before any
    compiler allocation. The two jet calculation always charges k as well.
    """
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError('Use positive absolute output precision.')
    if not isinstance(denominator_bits, int) or denominator_bits < 0:
        raise ValueError('Use a nonnegative external denominator exponent.')
    certificate = inverse_precision_certificate(n, left, right, k)
    compilation_precision = absolute_precision + denominator_bits + certificate['inverse_loss_bound']
    if maximum_compilation_precision is not None:
        if not isinstance(maximum_compilation_precision, int) or maximum_compilation_precision < 1:
            raise ValueError('The optional precision cap must be a positive integer.')
        if compilation_precision > maximum_compilation_precision:
            raise ValueError('The certified compiler precision exceeds the supplied cap.')
    moment_precision = absolute_precision + k + denominator_bits
    compiled = compile_regular_odd_kernel(circuit, compilation_precision)
    result = inverse_polynomial_kernel_contraction(n, left, right, k, compiled['terms'], moment_precision)
    result['scale_exponent'] -= denominator_bits
    result.update(absolute_precision=absolute_precision, compiler=compiled,
                  precision_certificate=certificate, compilation_precision=compilation_precision,
                  moment_precision=moment_precision, denominator_bits=denominator_bits,
                  enumerated_matrix_nodes=0,
                  scope='Progression unions with an explicit, fully charged residue-range inverse bound.')
    return result
