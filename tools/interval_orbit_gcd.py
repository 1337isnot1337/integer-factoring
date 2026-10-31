#!/usr/bin/env python3
"""Dense rational reference for interval orbit products; not a compressed algorithm."""

from math import gcd


def multiply_linear(polynomial, constant, linear):
    result = [0] * (len(polynomial) + 1)
    for i, coefficient in enumerate(polynomial):
        result[i] += constant * coefficient
        result[i + 1] += linear * coefficient
    return tuple(result)


def interval_product(interval, number=None):
    """Product (X-a), or (number-a*X), over a positive integer interval."""
    lower, upper = interval
    assert 1 <= lower <= upper
    result = (1,)
    for a in range(lower, upper):
        result = multiply_linear(result, -a if number is None else number,
                                  1 if number is None else -a)
    return result


def primitive(polynomial):
    values = list(polynomial)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    content = gcd(*values)
    if not content:
        return (0,)
    if values[-1] < 0:
        content = -content
    return tuple(a // content for a in values)


def remainder(left, right):
    """Primitive pseudo-remainder, preserving the gcd over Q."""
    assert right != (0,)
    if len(right) == 1:
        return (0,)
    while len(left) >= len(right) and left != (0,):
        shift = len(left) - len(right)
        leading = left[-1]
        reduced = [right[-1] * a for a in left]
        for i, a in enumerate(right):
            reduced[i + shift] -= leading * a
        left = primitive(reduced)
    return left


def polynomial_gcd(left, right):
    left, right = primitive(left), primitive(right)
    while right != (0,):
        left, right = right, remainder(left, right)
    return primitive(left)


def interval_factor_gcd(number, x_interval, y_interval):
    """Monic product (X-d) over factors d in the queried rectangle."""
    assert number >= 1
    result = polynomial_gcd(interval_product(x_interval),
                             interval_product(y_interval, number))
    assert result[-1] == 1
    return result
