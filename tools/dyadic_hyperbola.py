#!/usr/bin/env python3
"""Exact dyadic inverse fibres and interval counts; no faster factoring claim."""

from math import isqrt


def floor_sum(length, modulus, slope, offset):
    """Sum floor((slope*i+offset)/modulus), using Euclidean descent."""
    total = 0
    while True:
        quotient, slope = divmod(slope, modulus)
        total += quotient * length * (length - 1) // 2
        quotient, offset = divmod(offset, modulus)
        total += quotient * length
        height = slope * length + offset
        if height < modulus:
            return total
        length, offset = divmod(height, modulus)
        modulus, slope = slope, modulus


def affine_interval_count(modulus, slope, offset, start, stop, lower, upper):
    """Count start<=t<stop with lower<=(slope*t+offset)%modulus<upper."""
    assert 0 <= start <= stop <= modulus
    assert 0 <= lower <= upper <= modulus
    length = stop - start
    offset += slope * start

    def below(bound):
        return (length + floor_sum(length, modulus, slope, offset)
                - floor_sum(length, modulus, slope, offset + modulus - bound))

    return below(upper) - below(lower)


def inverse_fibre(number, bits, low_bits, root):
    """Return b,c,h with y=b+2^s*((c+h*t)%2^(m-s)) for x=root+2^s*t."""
    assert bits >= 1 and (bits + 1) // 2 <= low_bits <= bits
    base, modulus = 1 << low_bits, 1 << bits
    assert number % 2 and 0 < root < base and root % 2
    inverse = pow(root, -1, modulus)
    value = number * inverse % modulus
    high_modulus = 1 << (bits - low_bits)
    return (value % base, value // base,
            -number * inverse * inverse % high_modulus)


def inverse_polynomial(number, bits, low_bits, root):
    """An exact degree floor((bits-1)/low_bits) high-bit fibre polynomial."""
    assert bits >= low_bits >= 1 and number % 2
    base, modulus = 1 << low_bits, 1 << bits
    assert 0 < root < base and root % 2
    inverse = pow(root, -1, modulus)
    value = number * inverse % modulus
    high_modulus = modulus // base
    coefficients = [value // base]
    coefficient = -number * inverse * inverse % high_modulus
    for _ in range((bits - 1) // low_bits):
        coefficients.append(coefficient)
        coefficient = -coefficient * base * inverse % high_modulus
    return value % base, tuple(coefficients)


def residue_interval(lower, upper, root, base, high_modulus):
    """Intersect [lower,upper) with root+base*[0,high_modulus)."""
    start = max(0, min(high_modulus, -((root - lower) // base)))
    stop = max(0, min(high_modulus, -((root - upper) // base)))
    return start, stop


def rectangle_count(number, bits, x_interval, y_interval, low_bits=None):
    """Exact modular-hyperbola count and the number of low roots visited."""
    assert number % 2 and bits >= 1
    low_bits = (bits + 1) // 2 if low_bits is None else low_bits
    assert (bits + 1) // 2 <= low_bits <= bits
    modulus, base = 1 << bits, 1 << low_bits
    high_modulus = modulus // base
    x0, x1 = x_interval
    y0, y1 = y_interval
    assert 0 <= x0 <= x1 <= modulus and 0 <= y0 <= y1 <= modulus
    if x0 == x1 or y0 == y1:
        return 0, 0
    if x1 - x0 >= base:
        roots = range(1, base, 2)
    else:
        roots = (x % base for x in range(x0 + (1 - x0 % 2), x1, 2))
    count = visited = 0
    for root in roots:
        visited += 1
        b, c, slope = inverse_fibre(number, bits, low_bits, root)
        start, stop = residue_interval(x0, x1, root, base, high_modulus)
        lower, upper = residue_interval(y0, y1, b, base, high_modulus)
        count += affine_interval_count(high_modulus, slope, c, start, stop,
                                       lower, upper)
    return count, visited


def mobius_section(matrix, digit, bits):
    """Section of (a*t+b)/(c*t+1) after one low input/output bit."""
    a, b, c = matrix
    assert bits >= 1 and c % 2 == 0 and a % 2
    output_digit = (a * digit + b) % 2
    denominator = c * digit + 1
    modulus = 1 << (bits - 1)
    scale = pow(denominator, -1, modulus)
    next_matrix = ((a - output_digit * c) * scale % modulus,
                   ((a * digit + b - output_digit * denominator) // 2)
                   * scale % modulus,
                   2 * c * scale % modulus)
    return output_digit, next_matrix


def mobius_interval_count(matrix, bits, x_interval, y_interval, stats=None):
    """Exact recursive interval count, stopping at affine or full fibres."""
    modulus = 1 << bits
    a, b, c = (v % modulus for v in matrix)
    if stats is not None:
        stats['nodes'] = stats.get('nodes', 0) + 1
    x0, x1 = x_interval
    y0, y1 = y_interval
    assert 0 <= x0 <= x1 <= modulus and 0 <= y0 <= y1 <= modulus
    if x0 == x1 or y0 == y1:
        return 0
    if x_interval == (0, modulus):
        return y1 - y0
    if y_interval == (0, modulus):
        return x1 - x0
    assert a % 2 and c % 2 == 0
    if c == 0:
        if stats is not None:
            stats['affine_leaves'] = stats.get('affine_leaves', 0) + 1
        return affine_interval_count(modulus, a, b, x0, x1, y0, y1)
    total = 0
    for digit in (0, 1):
        output_digit, section = mobius_section((a, b, c), digit, bits)
        x_child = residue_interval(x0, x1, digit, 2, modulus // 2)
        y_child = residue_interval(y0, y1, output_digit, 2, modulus // 2)
        total += mobius_interval_count(section, bits - 1, x_child, y_child, stats)
    return total


def recursive_rectangle_count(number, bits, x_interval, y_interval, stats=None):
    """Odd inverse graph as a determinant-preserving Mobius recursion."""
    assert bits >= 1 and number % 2
    modulus = 1 << (bits - 1)
    x_child = residue_interval(*x_interval, 1, 2, modulus)
    y_child = residue_interval(*y_interval, 1, 2, modulus)
    return mobius_interval_count((-1, (number - 1) // 2, 2), bits - 1,
                                 x_child, y_child, stats)


def divisor_count(number, lower=2, upper=None):
    """Count divisors in [lower,upper) up to sqrt(number), for odd number>=3."""
    assert number >= 3 and number % 2
    limit = isqrt(number) + 1
    upper = limit if upper is None else min(upper, limit)
    lower = max(lower, 2)
    total = 0
    power = 2
    while power < upper:
        x0, x1 = max(lower, power), min(upper, 2 * power)
        if x0 < x1:
            y_interval = (number // (2 * power) + 1, number // power + 1)
            total += rectangle_count(number, number.bit_length(),
                                     (x0, x1), y_interval)[0]
        power *= 2
    return total


def smallest_factor(number):
    """Demonstrate exact recovery. Worst-case cost remains 2^(n/2)*poly(n)."""
    assert number >= 2
    if number % 2 == 0:
        return 2
    if not divisor_count(number):
        return number
    lower, upper = 2, isqrt(number)
    while lower < upper:
        middle = (lower + upper) // 2
        if divisor_count(number, lower, middle + 1):
            upper = middle
        else:
            lower = middle + 1
    return lower
