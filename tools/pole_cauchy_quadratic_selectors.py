#!/usr/bin/env python3
"""Exact quadratic-pole jets and residue-block sums of normalized inverses.

All poles have integer coordinates in the basis 1,w, with w^2+w+1=0.
The product factors must be nonzero at the chosen poles. Pole orders
are numerical parameters; no matrix nodes are enumerated.
"""

from interval_unit_jets import QuadraticJets
from modular_interval_jets import interval_elementaries_mod
from pole_cauchy_domain_precision import inverse_precision_certificate
from pole_cauchy_rational_poles import divisibility_count


def _valuation(value):
    assert value
    return (value & -value).bit_length() - 1


def _pair_valuation(pair):
    nonzero = [coordinate for coordinate in pair if coordinate]
    if not nonzero:
        raise ValueError("A quadratic denominator or product factor is zero")
    return min(_valuation(value) for value in nonzero)


def _pair_product(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def _pole(value):
    if isinstance(value, int):
        return value, 0
    if isinstance(value, (tuple, list)) and len(value) == 2 and all(isinstance(x, int) for x in value):
        return tuple(value)
    raise ValueError("A quadratic pole must have two integer coordinates")


def affine_quadratic_strata(first, step, length):
    """(first index, stride, count, valuation), with unit even-step factors.

At every parity split at most one branch has both coordinates even.
Only that branch continues. A singleton is handled directly, including
an isolated very high valuation or a forbidden zero.
"""
    if not isinstance(length, int) or length < 0:
        raise ValueError("Require a nonnegative integer length")
    if not length:
        return []
    first, step = tuple(first), tuple(step)
    offset, stride, value = 0, 1, 0
    result = []
    while length:
        if length == 1:
            result.append((offset, 0, 1, value + _pair_valuation(first)))
            break
        nonzero = [x for x in first + step if x]
        if not nonzero:
            raise ValueError("An affine quadratic product contains zero")
        common = min(_valuation(x) for x in nonzero)
        value += common
        first = tuple(x >> common for x in first)
        step = tuple(x >> common for x in step)
        if all(x % 2 == 0 for x in step):
            assert any(x % 2 for x in first)
            result.append((offset, stride, length, value))
            break
        following = None
        for parity in (0, 1):
            count = (length + 1 - parity) // 2
            if not count:
                continue
            branch = tuple(a + parity * b for a, b in zip(first, step))
            if any(x % 2 for x in branch):
                result.append((offset + parity * stride, 2 * stride, count, value))
            else:
                assert following is None
                following = (tuple(x // 2 for x in branch), count,
                             offset + parity * stride, 2 * stride)
        if following is None:
            break
        first, length, offset, stride = following
        value += 1
    return result


def _unit_product(ring, unit, step, length, cache):
    nonzero = [(index, coordinate) for index, pair in enumerate(step)
               for coordinate in pair if coordinate]
    if not nonzero:
        return ring.power(unit, length)
    assert all(coordinate % 2 == 0 for _, coordinate in nonzero)
    degree = min((ring.precision - 1) // min(_valuation(x) for _, x in nonzero), length)
    formal_order = min(index for index, _ in nonzero)
    if formal_order:
        degree = min(degree, (ring.size - 1) // formal_order)
    key = (length, degree, ring.precision)
    if key not in cache:
        cache[key] = interval_elementaries_mod(0, length, degree, ring.precision)
    elementary = cache[key]
    ratio = ring.multiply(step, ring.inverse(unit))
    series, power = ring.constant(1), ring.constant(1)
    for j in range(1, degree + 1):
        power = ring.multiply(power, ratio)
        series = ring.add(series, ring.scale(power, elementary[j]))
    return ring.multiply(ring.power(unit, length), series)


def affine_quadratic_jet_product(ring, first, step, first_derivative, step_derivative,
                                 strata, shift_exponent, cache=None):
    """Unit jet of prod_i(first+i*step+2^S X*(first'+i*step'))."""
    cache = {} if cache is None else cache
    result = ring.constant(1)
    for offset, stride, count, value in strata:
        unit = [((a + offset * b) >> value) % ring.modulus for a, b in zip(first, step)]
        increment = [(stride * b >> value) % ring.modulus for b in step]
        unit_jet, step_jet = list(ring.constant(*unit)), list(ring.constant(*increment))
        if ring.size > 1:
            assert shift_exponent > value
            scale = 1 << (shift_exponent - value)
            unit_jet[1] = tuple((a + offset * b) * scale % ring.modulus
                                for a, b in zip(first_derivative, step_derivative))
            step_jet[1] = tuple(stride * b * scale % ring.modulus for b in step_derivative)
        result = ring.multiply(result, _unit_product(ring, tuple(unit_jet), tuple(step_jet), count, cache))
    return result


def inverse_quadratic_pole_contraction(n, left, right, k, row_pole, column_pole,
                                       row_order, column_order, absolute_precision):
    """Full inverse contraction with repeated, integral quadratic poles."""
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError("Require positive absolute precision")
    if any(not isinstance(order, int) or order < 1 for order in (row_order, column_order)):
        raise ValueError("Require positive numerical pole orders")
    certificate = inverse_precision_certificate(n, left, right, k)
    c, d = _pole(row_pole), _pole(column_pole)
    z = (n, 1 << k)
    cd = _pair_product(c, d)
    denominator = (cd[0] - z[0], cd[1] - z[1])
    denominator_valuation = _pair_valuation(denominator)
    specifications = []
    for first, length in left:
        specifications.append((1, 1, (n - first * d[0], (1 << k) - first * d[1]),
                               (-2 * d[0], -2 * d[1]), (-first, 0), (-2, 0), length))
        specifications.append((-1, 0, (first - c[0], -c[1]), (2, 0), (-1, 0), (0, 0), length))
    for first, length in right:
        specifications.append((1, 0, (n - first * c[0], (1 << k) - first * c[1]),
                               (-2 * c[0], -2 * c[1]), (-first, 0), (-2, 0), length))
        specifications.append((-1, 1, (first - d[0], -d[1]), (2, 0), (-1, 0), (0, 0), length))
    prepared, product_valuation, maximum_factor_valuation = [], 0, 0
    pole_bounds = [0, 0]
    for sign, variable, first, step, first_derivative, step_derivative, length in specifications:
        strata = affine_quadratic_strata(first, step, length)
        maximum = max(value for _, _, _, value in strata)
        product_valuation += sign * sum(count * value for _, _, count, value in strata)
        maximum_factor_valuation = max(maximum_factor_valuation, maximum)
        if sign == -1:
            pole_bounds[variable] = max(pole_bounds[variable], maximum)
        prepared.append((sign, variable, first, step, first_derivative, step_derivative, strata))
    certified_lower = min(0, k + denominator_valuation - certificate["inverse_loss_bound"]
                          - sum(pole_bounds))
    if product_valuation < certified_lower:
        raise ArithmeticError("The quadratic product valuation violates inverse conditioning")
    shift_exponent = 1 + max(denominator_valuation, maximum_factor_valuation)
    derivative_degree = row_order + column_order - 2
    coefficient_division = derivative_degree * shift_exponent
    target_precision = absolute_precision + coefficient_division + k + denominator_valuation
    scale_floor = min(0, product_valuation)
    working_precision = target_precision - scale_floor
    product_precision = max(0, target_precision - product_valuation)
    ring = QuadraticJets(working_precision, 0)
    cache = {}
    if product_precision:
        jet_rings = (QuadraticJets(product_precision, row_order - 1),
                     QuadraticJets(product_precision, column_order - 1))
        jets = [jet_ring.constant(1) for jet_ring in jet_rings]
        for sign, variable, first, step, first_derivative, step_derivative, strata in prepared:
            jet_ring = jet_rings[variable]
            factor = affine_quadratic_jet_product(jet_ring, first, step, first_derivative,
                                                  step_derivative, strata, shift_exponent, cache)
            if sign == -1:
                factor = jet_ring.inverse(factor)
            jets[variable] = jet_ring.multiply(jets[variable], factor)
        jet_ring = jet_rings[0]
        jets[0] = jet_ring.multiply(jets[0], jet_ring.inverse(
            jet_ring.power(jet_ring.constant(*z), certificate["domain_size"])))
    unit_denominator = ring.constant(*(x >> denominator_valuation for x in denominator))
    inverse_constant = ring.inverse(unit_denominator)
    linear_scale = 1 << (shift_exponent - denominator_valuation)
    x_coefficient = ring.constant(*(x * linear_scale for x in d))
    y_coefficient = ring.constant(*(x * linear_scale for x in c))
    xy_coefficient = ring.constant(1 << (2 * shift_exponent - denominator_valuation))
    inverse_coefficients = [[ring.constant(0) for _ in range(column_order)] for _ in range(row_order)]
    inverse_coefficients[0][0] = inverse_constant
    for i in range(row_order):
        for j in range(column_order):
            if not (i or j):
                continue
            total = ring.constant(0)
            if i:
                total = ring.add(total, ring.multiply(x_coefficient, inverse_coefficients[i - 1][j]))
            if j:
                total = ring.add(total, ring.multiply(y_coefficient, inverse_coefficients[i][j - 1]))
            if i and j:
                total = ring.add(total, ring.multiply(xy_coefficient, inverse_coefficients[i - 1][j - 1]))
            inverse_coefficients[i][j] = ring.scale(ring.multiply(inverse_constant, total), -1)
    result = ring.scale(inverse_coefficients[-1][-1], 1 << (-scale_floor))
    if product_precision:
        convolution = ring.constant(0)
        for i in range(row_order):
            for j in range(column_order):
                term = ring.multiply(ring.constant(*jets[0][i]), ring.constant(*jets[1][j]))
                term = ring.multiply(term, inverse_coefficients[row_order - 1 - i][column_order - 1 - j])
                convolution = ring.add(convolution, term)
        result = ring.add(result, ring.scale(convolution, -(1 << (product_valuation - scale_floor))))
    result = ring.multiply(result, ring.multiply(ring.constant(*z), ring.inverse(ring.constant(0, 1))))
    return {"residue": list(result[0]),
            "scale_exponent": scale_floor - k - denominator_valuation - coefficient_division,
            "working_precision": working_precision, "absolute_precision": absolute_precision,
            "row_pole": list(c), "column_pole": list(d),
            "row_order": row_order, "column_order": column_order,
            "row_denominator_bound": pole_bounds[0], "column_denominator_bound": pole_bounds[1],
            "product_valuation": product_valuation, "product_valuation_lower_bound": certified_lower,
            "product_unit_precision": product_precision,
            "denominator_valuation": denominator_valuation,
            "maximum_affine_factor_valuation": maximum_factor_valuation,
            "jet_shift_exponent": shift_exponent, "coefficient_division_valuation": coefficient_division,
            "affine_valuation_strata": sum(len(item[-1]) for item in prepared),
            "distinct_elementary_intervals": len(cache),
            "precision_certificate": certificate, "enumerated_matrix_nodes": 0,
            "scope": "Integral quadratic poles with nonzero product factors; numerical pole orders are charged."}


def inverse_residue_block_parity(n, left, right, k, row_residue, row_bits,
                                 column_residue, column_bits):
    """Sum of 2^L M^-1 entries in two residue classes, modulo 2."""
    if any(not isinstance(bits, int) or bits < 0 for bits in (row_bits, column_bits)):
        raise ValueError("Require nonnegative residue precisions")
    if any(not isinstance(residue, int) for residue in (row_residue, column_residue)):
        raise ValueError("Require integer residues")
    row_residue %= 1 << row_bits
    column_residue %= 1 << column_bits
    c, d = (row_residue, -(1 << row_bits)), (column_residue, -(1 << column_bits))
    contraction = inverse_quadratic_pole_contraction(n, left, right, k, c, d, 3, 3, 1)
    loss = contraction["precision_certificate"]["inverse_loss_bound"]
    normalization = loss + 3 * (row_bits + column_bits)
    scale = contraction["scale_exponent"] + normalization
    if scale > 0:
        parity = [0, 0]
    elif scale == 0:
        parity = [coordinate % 2 for coordinate in contraction["residue"]]
    else:
        assert all(coordinate % (1 << (-scale)) == 0 for coordinate in contraction["residue"])
        parity = [(coordinate >> (-scale)) % 2 for coordinate in contraction["residue"]]
    return {"parity": parity, "inverse_normalization_exponent": loss,
            "row_residue": row_residue, "row_bits": row_bits,
            "column_residue": column_residue, "column_bits": column_bits,
            "selected_rows": sum(divisibility_count(first - row_residue, 2, length, row_bits)
                                  for first, length in left),
            "selected_columns": sum(divisibility_count(first - column_residue, 2, length, column_bits)
                                     for first, length in right),
            "quadratic_pole_contraction": contraction,
            "enumerated_matrix_nodes": 0,
            "scope": "Residue-block sum of the normalized inverse, not the factor-pair incidence matrix."}
