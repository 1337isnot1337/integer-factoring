#!/usr/bin/env python3
"""Exact polynomial realization in one Boolean Cauchy packet.

All arithmetic is rational or modular integer arithmetic. This module does not
claim an improved factoring algorithm or a new asymptotic precision bound.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb


def trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result or [0]


def add(left, right):
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return trim(result)


def scale(polynomial, scalar):
    return trim([scalar * value for value in polynomial])


def multiply(left, right, modulus=None):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    result[i + j] += a * b
    if modulus:
        result = [value % modulus for value in result]
    return trim(result)


def evaluate(polynomial, point):
    result = 0
    for coefficient in reversed(polynomial):
        result = result * point + coefficient
    return result


def shift(polynomial, center, dilation):
    result = [0]
    for coefficient in reversed(polynomial):
        result = add(multiply(result, [center, dilation]), [coefficient])
    return result


def valuation(value):
    value = Fraction(value)
    if not value:
        raise ValueError("The valuation of zero is not finite.")
    numerator = abs(value.numerator)
    denominator = value.denominator
    return ((numerator & -numerator).bit_length()
            - (denominator & -denominator).bit_length())


def content(polynomial):
    nonzero = [valuation(value) for value in polynomial if value]
    if not nonzero:
        raise ValueError("The zero polynomial has no finite content.")
    return min(nonzero)


def validate_packet_size(length):
    if not isinstance(length, int) or length < 2 or length & (length - 1):
        raise ValueError("Packet length must be a power of two at least two.")


def cauchy_polynomials(length):
    """Return W_q(T), with q=y-x, by direct multiplication of its roots."""
    validate_packet_size(length)
    grid = range(-length, length)
    derivative_at_zero = 1
    for root in grid:
        if root:
            derivative_at_zero *= -root
    result = {}
    for omitted in range(-length + 1, length):
        polynomial = [1]
        for root in grid:
            if root != omitted:
                polynomial = multiply(polynomial, [-root, 1])
        result[omitted] = [Fraction(value, derivative_at_zero) for value in polynomial]
    return result


def affine_realization(polynomial, length):
    """Construct the affine row and zero/one-hot column tensors.

    The polynomial must have degree at most length and be 2-adically integral
    at every integer in [-length,length). Returned values have odd denominators.
    """
    validate_packet_size(length)
    polynomial = trim([Fraction(value) for value in polynomial])
    if len(polynomial) > length + 1:
        raise ValueError("The affine realization requires degree at most L.")
    values = [evaluate(polynomial, point) for point in range(-length, length)]
    if any(value.denominator % 2 == 0 for value in values):
        raise ValueError("The polynomial is not integral on the complete grid.")

    differences = values[:length + 1]
    forward = []
    for _ in range(length + 1):
        forward.append(differences[0])
        differences = [right - left for left, right in zip(differences, differences[1:])]
    denominator = comb(2 * length - 1, length)
    column = [Fraction(0)] * length
    for order in range(1, length + 1):
        coefficient = Fraction((-1)**order * comb(2 * length - 1, order), denominator) * forward[order]
        for power in range(length - order + 1):
            column[order - 1 + power] += coefficient * (-1)**power * comb(length - order, power)
    row = [(-1)**(length - 1 - x) * comb(length - 1, x) for x in range(length)]
    constant = values[0]
    zero = 1 - sum(column)
    assert all(Fraction(value).denominator & 1 for value in [constant, zero] + row + column)
    return {"L": length, "constant": constant, "row": row, "column": column, "column_zero": zero}


def reconstructed_affine_polynomial(realization):
    length = realization["L"]
    kernels = cauchy_polynomials(length)
    result = [realization["constant"]]
    for x, a in enumerate(realization["row"]):
        for y, b in enumerate(realization["column"]):
            result = add(result, scale(kernels[y - x], a * b))
    return trim(result)


def direct_boolean_affine_polynomial(realization):
    """Literal finite sum over Boolean row states and supported column states."""
    length = realization["L"]
    kernels = cauchy_polynomials(length)
    result = [0]
    supported_columns = [(0, realization["column_zero"])]
    supported_columns += [(1 << y, value) for y, value in enumerate(realization["column"])]
    for column_state, column_weight in supported_columns:
        if not column_weight:
            continue
        images = []
        for x in range(length):
            image = [0]
            for y in range(length):
                if column_state >> y & 1:
                    image = add(image, kernels[y - x])
            images.append(image)
        for row_state in range(1 << length):
            weight = realization["constant"] + sum(
                realization["row"][x] for x in range(length) if row_state >> x & 1
            )
            term = [weight * column_weight]
            for x, image in enumerate(images):
                factor = image if row_state >> x & 1 else add([1], scale(image, -1))
                term = multiply(term, factor)
            result = add(result, term)
    return trim(result)


def odd_root_polynomial(length):
    validate_packet_size(length)
    result = [Fraction(1)]
    for root in range(-length, length):
        if root & 1:
            result = multiply(result, [1, Fraction(-1, root)])
    return result


def modular_packet_polynomial(length, center, precision, row_state, column_tensor):
    """Compute every coefficient of a shifted packet using integer residues.

    The row tensor is the indicator of row_state. The column tensor is given
    explicitly on all Boolean states. Rational kernel entries are independently
    constructed globally, shifted, and then reduced using their odd denominators.
    """
    if len(column_tensor) != 1 << length:
        raise ValueError("Incorrect column tensor length.")
    modulus = 1 << precision
    kernels = {}
    for difference, polynomial in cauchy_polynomials(length).items():
        shifted = shift(polynomial, center, 2 * length)
        assert all(Fraction(value).denominator & 1 for value in shifted)
        kernels[difference] = [
            Fraction(value).numerator * pow(Fraction(value).denominator, -1, modulus) % modulus
            for value in shifted
        ]
    degree = length * (2 * length - 1)
    result = [0] * (degree + 1)
    for state, coefficient in enumerate(column_tensor):
        coefficient %= modulus
        if not coefficient:
            continue
        term = [coefficient]
        for x in range(length):
            image = [0] * (2 * length)
            for y in range(length):
                if state >> y & 1:
                    for index, value in enumerate(kernels[y - x]):
                        image[index] += value
            if not (row_state >> x & 1):
                image = [-value for value in image]
                image[0] += 1
            term = multiply(term, [value % modulus for value in image], modulus)
        for index, value in enumerate(term):
            result[index] = (result[index] + value) % modulus
    return trim(result)
