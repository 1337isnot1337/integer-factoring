#!/usr/bin/env python3
"""A single entry of M^-1 v for polynomial v, without enumerating domains."""

from __future__ import annotations

from interval_unit_jets import QuadraticJets
from modular_interval_jets import unit_progression_jet
from pole_cauchy_conditioning import affine_product, domain_size


def _validate(n: int, left: tuple, right: tuple, b: int, k: int, precision: int) -> int:
    if n < 1 or n % 2 == 0 or precision < 1:
        raise ValueError("Require positive odd N and positive precision")
    for domain in (left, right):
        if not domain:
            raise ValueError("Domains must be nonempty")
        previous = 0
        for first, length in domain:
            if first <= previous or first % 2 == 0 or length < 1:
                raise ValueError("Require ordered disjoint positive odd progressions")
            previous = first + 2 * (length - 1)
    size = domain_size(left)
    if size != domain_size(right):
        raise ValueError("The matrix must be square")
    if not any(first <= b <= first + 2 * (length - 1) and (b - first) % 2 == 0
               for first, length in right):
        raise ValueError("The requested output column must belong to B")
    max_a = left[-1][0] + 2 * (left[-1][1] - 1)
    max_b = right[-1][0] + 2 * (right[-1][1] - 1)
    bound = max(abs(n - left[0][0] * right[0][0]), abs(n - max_a * max_b))
    if k < max(1, bound.bit_length()):
        raise ValueError("Displacement must separate all nonzero product differences")
    return size


def inverse_constant_vector_entry(n: int, left: tuple, right: tuple,
                                  b: int, k: int, precision: int) -> dict:
    """Normalized entry of M^-1 1: an exact valuation and a unit modulo 2^p."""
    size = _validate(n, left, right, b, k, precision)
    ring = QuadraticJets(precision, 0)
    z = ring.constant(n, 1 << k)
    cross_exponent = 0
    cross_unit = ring.constant(1)
    for first, length in left:
        exponent, unit = affine_product(ring, n - first * b, -2 * b, length, k)
        cross_exponent += exponent
        cross_unit = ring.multiply(cross_unit, unit)
    derivative_exponent = 0
    derivative_unit = ring.constant(1)
    node_product = ring.constant(1)
    for first, length in right:
        exponent, unit = affine_product(ring, b - first, -2, length, None, omit_zero=True)
        derivative_exponent += exponent
        derivative_unit = ring.multiply(derivative_unit, unit)
        exponent, unit = affine_product(ring, first, 2, length, None)
        assert exponent == 0
        node_product = ring.multiply(node_product, unit)
    numerator = ring.multiply(cross_unit, node_product)
    denominator = ring.multiply(derivative_unit, ring.constant(0, b))
    denominator = ring.multiply(denominator, ring.power(z, size - 1))
    unit = ring.multiply(numerator, ring.inverse(denominator))
    if size % 2 == 0:
        unit = ring.scale(unit, -1)
    assert unit[0][0] % 2 or unit[0][1] % 2
    return {"valuation": cross_exponent - derivative_exponent - k,
            "unit": list(unit[0]), "precision": precision}


def normalized_interpolation_series(n: int, left: tuple, right: tuple,
                                    k: int, precision: int, degree: int) -> tuple:
    """Truncate prod_b(1-Y/b) / prod_a(1-aY/z) through the given degree."""
    ring = QuadraticJets(precision, degree)
    y = [(0, 0)] * ring.size
    if degree:
        y[1] = (1, 0)
    y = tuple(y)
    z_inverse = ring.inverse(ring.constant(n, 1 << k))
    y_over_z = ring.multiply(y, z_inverse)
    row_product = ring.constant(1)
    for first, length in left:
        base = ring.add(ring.constant(1), ring.scale(y_over_z, -first))
        step = ring.scale(y_over_z, -2)
        row_product = ring.multiply(row_product, unit_progression_jet(
            ring, (0, length), base, step))
    column_product = ring.constant(1)
    for first, length in right:
        base = ring.add(ring.constant(first), ring.scale(y, -1))
        factor = unit_progression_jet(ring, (0, length), base, ring.constant(2))
        factor = ring.multiply(factor, ring.inverse(ring.constant(*factor[0])))
        column_product = ring.multiply(column_product, factor)
    return ring.multiply(column_product, ring.inverse(row_product))


def inverse_polynomial_vector_entry(n: int, left: tuple, right: tuple,
                                    b: int, k: int, coefficients: tuple[int, ...],
                                    precision: int) -> dict:
    """Return an integral scaled residue of (M^-1 v)[b], v(a)=sum c_d a^d.

    If the returned scale is E and residue is U, the answer is 2^E U
    with error of valuation at least E+precision. U may be even or zero:
    cancellation in the polynomial input does not claim a new exact
    valuation. The numerical node counts are never enumerated.
    """
    if not coefficients or any(not isinstance(c, int) for c in coefficients):
        raise ValueError("Provide a nonempty tuple of integer coefficients")
    constant = inverse_constant_vector_entry(n, left, right, b, k, precision)
    degree = len(coefficients) - 1
    series = normalized_interpolation_series(n, left, right, k, precision, degree)
    ring = QuadraticJets(precision, 0)
    z_over_b = ring.scale(ring.constant(n, 1 << k), pow(b, -1, ring.modulus))
    prefix = ring.constant(0)
    b_power = 1
    multiplier_power = ring.constant(1)
    multiplier = ring.constant(0)
    for index, coefficient in enumerate(coefficients):
        prefix = ring.add(prefix, ring.scale(ring.constant(*series[index]), b_power))
        multiplier = ring.add(multiplier, ring.scale(ring.multiply(multiplier_power, prefix), coefficient))
        multiplier_power = ring.multiply(multiplier_power, z_over_b)
        b_power = b_power * b % ring.modulus
    result = ring.multiply(ring.constant(*constant["unit"]), multiplier)
    return {"scale_exponent": constant["valuation"], "residue": list(result[0]),
            "precision": precision, "polynomial_degree": degree}


def inverse_polynomial_bilinear_form(n: int, left: tuple, right: tuple, k: int,
                                   row_coefficients: tuple[int, ...],
                                   column_coefficients: tuple[int, ...],
                                   precision: int) -> dict:
    """Compute 2^k g(B)^T M^-1 f(A) modulo 2^precision.

    Both complete matrix axes are summed. The explicit scale -k does
    not claim an exact valuation; coefficients may cancel.
    """
    if not row_coefficients or not column_coefficients:
        raise ValueError("Provide nonempty coefficient tuples")
    if any(not isinstance(c, int) for c in row_coefficients + column_coefficients):
        raise ValueError("Polynomial coefficients must be integers")
    size = _validate(n, left, right, right[0][0], k, precision)
    row_degree = len(row_coefficients) - 1
    column_degree = len(column_coefficients) - 1
    origin_series = normalized_interpolation_series(n, left, right, k, precision, row_degree)
    ring = QuadraticJets(precision, column_degree)
    y = [(0, 0)] * ring.size
    if column_degree:
        y[1] = (1, 0)
    y = tuple(y)
    z = ring.constant(n, 1 << k)
    zy = ring.multiply(z, y)
    numerator = ring.constant(1)
    denominator = ring.constant(1)
    for first, length in left:
        factor = unit_progression_jet(ring, (0, length),
                                      ring.add(ring.constant(first), ring.scale(zy, -1)),
                                      ring.constant(2))
        factor = ring.multiply(factor, ring.inverse(ring.constant(*factor[0])))
        numerator = ring.multiply(numerator, factor)
    for first, length in right:
        factor = unit_progression_jet(ring, (0, length),
                                      ring.add(ring.constant(1), ring.scale(y, -first)),
                                      ring.scale(y, -2))
        denominator = ring.multiply(denominator, factor)
    infinity_series = ring.multiply(numerator, ring.inverse(denominator))
    scalar = QuadraticJets(precision, 0)
    z = scalar.constant(n, 1 << k)
    origin_constant = scalar.constant(1)
    for domain in (left, right):
        for first, length in domain:
            origin_constant = scalar.multiply(origin_constant, unit_progression_jet(
                scalar, (0, length), scalar.constant(first), scalar.constant(2)))
    origin_constant = scalar.multiply(origin_constant, scalar.inverse(scalar.power(z, size)))
    origin_series = tuple(scalar.pair_product(origin_constant[0], value) for value in origin_series)
    result = scalar.constant(0)
    factor = scalar.multiply(z, scalar.inverse(scalar.constant(0, 1)))
    for d, row_coefficient in enumerate(row_coefficients):
        for e, column_coefficient in enumerate(column_coefficients):
            value = scalar.constant(int(d == e))
            for j in range(max(0, d - e), d + 1):
                product = scalar.pair_product(origin_series[j], infinity_series[e - d + j])
                value = scalar.add(value, scalar.constant(-product[0], -product[1]))
            value = scalar.scale(scalar.multiply(factor, value), row_coefficient * column_coefficient)
            result = scalar.add(result, value)
        factor = scalar.multiply(factor, z)
    return {"scale_exponent": -k, "residue": list(result[0]), "precision": precision,
            "row_degree": row_degree, "column_degree": column_degree}
