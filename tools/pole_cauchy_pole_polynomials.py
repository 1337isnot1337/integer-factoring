#!/usr/bin/env python3
"""Batch proper rational weights through two normalized product jets.

Coefficient index i denotes pole order i+1. Coefficients and poles are
integers or integer coordinate pairs in Z[w], w^2+w+1=0. Every numerical
degree and coefficient bit length is charged; matrix nodes are not listed.
"""

from interval_unit_jets import QuadraticJets
from pole_cauchy_domain_precision import inverse_precision_certificate
from pole_cauchy_quadratic_selectors import (
    _pair_product, _pair_valuation, _pole, affine_quadratic_strata,
    affine_quadratic_jet_product,
)
from pole_cauchy_rational_poles import divisibility_count
from pole_cauchy_residue_blocks import selector_rounding_coefficients


def _coefficients(values):
    result = [_pole(value) for value in values]
    while result and result[-1] == (0, 0):
        result.pop()
    if not result:
        raise ValueError('Require a nonzero proper pole polynomial')
    return result


def _integral_weights(coefficients, shift):
    floor = min(_pair_valuation(value)-index*shift
                for index, value in enumerate(coefficients) if value != (0, 0))
    weights = []
    for index, value in enumerate(coefficients):
        exponent = -index*shift-floor
        if exponent >= 0:
            weights.append(tuple(x << exponent for x in value))
        else:
            assert all(x % (1 << -exponent) == 0 for x in value)
            weights.append(tuple(x >> -exponent for x in value))
    return floor, weights


def inverse_pole_polynomial_contraction(n, left, right, k, row_pole, column_pole,
                                        row_coefficients, column_coefficients,
                                        absolute_precision):
    """Contract M^-1 against two proper single-pole rational polynomials.

The weights are sum_i alpha_i/(a-c)^(i+1) and sum_j beta_j/(b-d)^(j+1).
Only two product jets are assembled, regardless of the number of nonzero
coefficients. A separate exact scale records all coefficient divisions.
"""
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError('Require positive absolute precision')
    alpha, beta = _coefficients(row_coefficients), _coefficients(column_coefficients)
    rows, columns = len(alpha), len(beta)
    certificate = inverse_precision_certificate(n, left, right, k)
    c, d, z = _pole(row_pole), _pole(column_pole), (n, 1 << k)
    cd = _pair_product(c, d)
    denominator = tuple(a-b for a, b in zip(cd, z))
    v0 = _pair_valuation(denominator)
    specifications = []
    for first, length in left:
        specifications.append((1, 1, (n-first*d[0], (1 << k)-first*d[1]),
                               (-2*d[0], -2*d[1]), (-first, 0), (-2, 0), length))
        specifications.append((-1, 0, (first-c[0], -c[1]), (2, 0),
                               (-1, 0), (0, 0), length))
    for first, length in right:
        specifications.append((1, 0, (n-first*c[0], (1 << k)-first*c[1]),
                               (-2*c[0], -2*c[1]), (-first, 0), (-2, 0), length))
        specifications.append((-1, 1, (first-d[0], -d[1]), (2, 0),
                               (-1, 0), (0, 0), length))
    prepared, q, maximum = [], 0, 0
    pole_bounds = [0, 0]
    for sign, variable, first, step, first_derivative, step_derivative, length in specifications:
        strata = affine_quadratic_strata(first, step, length)
        largest = max(value for _, _, _, value in strata)
        q += sign*sum(count*value for _, _, count, value in strata)
        maximum = max(maximum, largest)
        if sign == -1:
            pole_bounds[variable] = max(pole_bounds[variable], largest)
        prepared.append((sign, variable, first, step, first_derivative, step_derivative, strata))
    lower = min(0, k+v0-certificate['inverse_loss_bound']-sum(pole_bounds))
    if q < lower:
        raise ArithmeticError('The product valuation violates inverse conditioning')
    shift = 1+max(v0, maximum)
    first_floor, first_weights = _integral_weights(alpha, shift)
    second_floor, second_weights = _integral_weights(beta, shift)
    e = min(0, q)
    scale = first_floor+second_floor+e-k-v0
    precision = absolute_precision-scale
    report = {'absolute_precision': absolute_precision, 'scale_exponent': scale,
              'working_precision': max(0, precision),
              'row_pole': list(c), 'column_pole': list(d),
              'row_degree': rows, 'column_degree': columns,
              'nonzero_row_coefficients': sum(x != (0, 0) for x in alpha),
              'nonzero_column_coefficients': sum(x != (0, 0) for x in beta),
              'row_weight_scale': first_floor, 'column_weight_scale': second_floor,
              'jet_shift_exponent': shift, 'product_valuation': q,
              'product_valuation_lower_bound': lower, 'denominator_valuation': v0,
              'row_denominator_bound': pole_bounds[0], 'column_denominator_bound': pole_bounds[1],
              'affine_valuation_strata': sum(len(item[-1]) for item in prepared),
              'precision_certificate': certificate, 'enumerated_matrix_nodes': 0,
              'scope': 'Explicit proper rational coefficient arrays; no product-congruence or rank oracle.'}
    if precision <= 0:
        return {**report, 'residue': [0, 0], 'product_unit_precision': 0,
                'assembled_product_jets': 0, 'denominator_coefficient_cells': 0,
                'vanishes_at_requested_precision': True}
    ring = QuadraticJets(precision, 0)
    weight_a = [ring.constant(*value) for value in first_weights]
    weight_b = [ring.constant(*value) for value in second_weights]
    product_precision = max(0, precision-(q-e))
    cache = {}
    if product_precision:
        jet_rings = (QuadraticJets(product_precision, rows-1),
                     QuadraticJets(product_precision, columns-1))
        jets = [jet_ring.constant(1) for jet_ring in jet_rings]
        for sign, variable, first, step, first_derivative, step_derivative, strata in prepared:
            jet_ring = jet_rings[variable]
            factor = affine_quadratic_jet_product(jet_ring, first, step, first_derivative,
                                                  step_derivative, strata, shift, cache)
            jets[variable] = jet_ring.multiply(jets[variable],
                                               factor if sign == 1 else jet_ring.inverse(factor))
        jets[0] = jet_rings[0].multiply(jets[0], jet_rings[0].inverse(
            jet_rings[0].power(jet_rings[0].constant(*z), certificate['domain_size'])))
    inverse_constant = ring.inverse(ring.constant(*(x >> v0 for x in denominator)))
    linear_scale = 1 << (shift-v0)
    x_term, y_term = (ring.constant(*(x*linear_scale for x in value)) for value in (d, c))
    xy_term = ring.constant(1 << (2*shift-v0))
    table = [[ring.constant(0) for _ in range(columns)] for _ in range(rows)]
    table[0][0] = inverse_constant
    for i in range(rows):
        for j in range(columns):
            if not (i or j):
                continue
            value = ring.constant(0)
            if i:
                value = ring.add(value, ring.multiply(x_term, table[i-1][j]))
            if j:
                value = ring.add(value, ring.multiply(y_term, table[i][j-1]))
            if i and j:
                value = ring.add(value, ring.multiply(xy_term, table[i-1][j-1]))
            table[i][j] = ring.scale(ring.multiply(inverse_constant, value), -1)
    def contract(a, b):
        value = ring.constant(0)
        for i in range(rows):
            for j in range(columns):
                value = ring.add(value, ring.multiply(table[i][j], ring.multiply(a[i], b[j])))
        return value
    answer = ring.scale(contract(weight_a, weight_b), 1 << -e)
    if product_precision:
        correlated = []
        for weights, jet in zip((weight_a, weight_b), jets):
            vector = []
            for first in range(len(weights)):
                value = ring.constant(0)
                for last in range(first, len(weights)):
                    value = ring.add(value, ring.multiply(weights[last], ring.constant(*jet[last-first])))
                vector.append(value)
            correlated.append(vector)
        answer = ring.add(answer, ring.scale(contract(*correlated), -(1 << (q-e))))
    answer = ring.multiply(answer, ring.multiply(ring.constant(*z), ring.inverse(ring.constant(0, 1))))
    return {**report, 'residue': list(answer[0]), 'product_unit_precision': product_precision,
            'assembled_product_jets': 2 if product_precision else 0,
            'denominator_coefficient_cells': rows*columns,
            'coefficient_correlation_terms': (rows*(rows+1)+columns*(columns+1))//2,
            'distinct_elementary_intervals': len(cache),
            'vanishes_at_requested_precision': False}


def inverse_residue_block_batched(n, left, right, k, row_residue, row_bits,
                                  column_residue, column_bits, precision):
    """Full-precision normalized inverse block, using one pair of product jets."""
    if any(not isinstance(bits, int) or bits < 0 for bits in (row_bits, column_bits)):
        raise ValueError('Require nonnegative residue precisions')
    if any(not isinstance(value, int) for value in (row_residue, column_residue)):
        raise ValueError('Require integer residues')
    coefficients = selector_rounding_coefficients(precision)
    certificate = inverse_precision_certificate(n, left, right, k)
    loss = certificate['inverse_loss_bound']
    row_residue %= 1 << row_bits
    column_residue %= 1 << column_bits
    row_weights, column_weights = [0]*(6*precision-3), [0]*(6*precision-3)
    for degree, coefficient in coefficients.items():
        row_weights[3*degree-1] = coefficient << (loss+3*degree*row_bits)
        column_weights[3*degree-1] = coefficient << (3*degree*column_bits)
    result = inverse_pole_polynomial_contraction(
        n, left, right, k, (row_residue, -(1 << row_bits)),
        (column_residue, -(1 << column_bits)), row_weights, column_weights, precision)
    scale = result['scale_exponent']
    if scale >= precision:
        residue = [0, 0]
    elif scale >= 0:
        residue = [(value << scale) % (1 << precision) for value in result['residue']]
    else:
        assert all(value % (1 << -scale) == 0 for value in result['residue'])
        assert result['working_precision']+scale >= precision
        residue = [(value >> -scale) % (1 << precision) for value in result['residue']]
    return {'residue': residue, 'precision': precision,
            'inverse_normalization_exponent': loss,
            'row_residue': row_residue, 'row_bits': row_bits,
            'column_residue': column_residue, 'column_bits': column_bits,
            'selected_rows': sum(divisibility_count(first-row_residue, 2, length, row_bits)
                                 for first, length in left),
            'selected_columns': sum(divisibility_count(first-column_residue, 2, length, column_bits)
                                    for first, length in right),
            'selector_degree': 2*precision-1, 'maximum_pole_order': 6*precision-3,
            'batched_contraction': result, 'enumerated_matrix_nodes': 0,
            'scope': 'Full-precision normalized inverse block; no rank or factor-count inference.'}
