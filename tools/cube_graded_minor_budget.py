#!/usr/bin/env python3
"""Exact shared-minor envelopes for graded eight-leg cube tensors.

This constructs a necessary profile budget, not an actual crossing tensor.
Alphabet multiplicities and matrix dimensions stay compressed as integers.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256


def convolution(first, second):
    result = Counter()
    for a, x in first.items():
        for b, y in second.items():
            result[a + b] += x * y
    return result


def power(polynomial, exponent):
    result = Counter({0: 1})
    for _ in range(exponent):
        result = convolution(result, polynomial)
    return result


def raw_knots(energies):
    if not energies or energies.get(0, 0) < 1:
        raise ValueError('The alphabet must contain a zero-energy letter.')
    if any(not isinstance(e, int) or e < 0 or not isinstance(n, int) or n <= 0
           for e, n in energies.items()):
        raise ValueError('Use nonnegative integer energies and positive integer multiplicities.')
    if sum(energies.values()) < 2:
        raise ValueError('The physical alphabet must contain at least two letters.')
    squares = convolution(energies, energies)
    fourths = convolution(squares, squares)
    records = []
    for wires, factor, counts in ((1, Fraction(2, 3), energies),
                                  (2, Fraction(1), squares),
                                  (4, Fraction(2), fourths)):
        size = total = 0
        knots = [(0, Fraction(0))]
        for energy, count in sorted(counts.items()):
            size += count
            total += count * energy
            knots.append((size, factor * total))
        records.append({'wires': wires, 'knots': knots})
    return records


def envelope(energies):
    records = raw_knots(energies)
    bounds = {}
    for record in records:
        for size, value in record['knots']:
            if size not in bounds or value < bounds[size]:
                bounds[size] = value
    hull = []
    for size, value in sorted(bounds.items()):
        while len(hull) >= 2:
            first, second = hull[-2:]
            if ((second[1] - first[1]) * (size - second[0])
                    >= (value - second[1]) * (second[0] - first[0])):
                hull.pop()
            else:
                break
        hull.append((size, value))
    blocks = [(second[0] - first[0],
               (second[1] - first[1]) / (second[0] - first[0]))
              for first, second in zip(hull, hull[1:])]
    assert sum(count for count, _ in blocks) == sum(energies.values()) ** 4
    assert all(a[1] <= b[1] for a, b in zip(blocks, blocks[1:]))
    return {'alphabet': sum(energies.values()), 'raw': records,
            'hull': hull, 'blocks': blocks}


def cumulative(blocks, count):
    value = Fraction(0)
    for multiplicity, order in blocks:
        used = min(count, multiplicity)
        value += used * order
        count -= used
        if not count:
            return value
    if count:
        raise ValueError('The requested prefix exceeds the profile dimension.')
    return value


def integer_relaxation(result):
    """A feasible integer abstract profile; not a universal lower envelope."""
    profile = Counter()
    for first, second in zip(result['hull'], result['hull'][1:]):
        count = second[0] - first[0]
        total = second[1] - first[1]
        if total.denominator != 1:
            raise ValueError('Integral block totals are required for this rounding rule.')
        order, upper = divmod(int(total), count)
        profile[order] += count - upper
        profile[order + 1] += upper
    profile = Counter({order: count for order, count in profile.items() if count})
    blocks = [(count, Fraction(order)) for order, count in sorted(profile.items())]
    for record in result['raw']:
        for count, bound in record['knots']:
            assert cumulative(blocks, count) <= bound
    return profile


def positive_interval_certificate(polynomial):
    polynomial = {degree: value for degree, value in polynomial.items() if value}
    lead = min(polynomial)
    degree = max(polynomial) - lead
    quotient = {power - lead: value for power, value in polynomial.items()}
    transformed = [0] * (degree + 1)
    for i, value in quotient.items():
        binomial = 1
        for k in range(i, degree + 1):
            transformed[k] += value * binomial
            if k < degree:
                binomial = binomial * (degree - k) // (k + 1 - i)
    if not all(value > 0 for value in transformed):
        raise ArithmeticError('This polynomial does not have the asserted positive coefficient certificate.')
    digest = sha256()
    for value in transformed:
        data = value.to_bytes((value.bit_length() + 7) // 8, 'big')
        digest.update(len(data).to_bytes(8, 'big'))
        digest.update(data)
    return {'leading_degree': lead, 'quotient_degree': degree,
            'positive_coefficients': len(transformed),
            'minimum_positive_coefficient': min(transformed),
            'coefficient_encoding': 'eight-byte length followed by unsigned big-endian integer',
            'positive_coefficients_sha256': digest.hexdigest()}


def integer_hash(value):
    return sha256(value.to_bytes((value.bit_length() + 7) // 8, 'big')).hexdigest()


def rational_blocks(result):
    return [{'multiplicity': count, 'order': [order.numerator, order.denominator]}
            for count, order in result['blocks']]
