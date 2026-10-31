#!/usr/bin/env python3
"""Exact valuation-separated jets at scaled integer centres."""

from math import gcd

from interval_unit_jets import QuadraticJets, unit_progression_jet


def valuation_two(value):
    assert value
    value = abs(value)
    return (value & -value).bit_length() - 1


def nonzero_affine_jet(ring, interval, a, b, c, d, scale_bits):
    """Normalize product (a+b*t+2^scale_bits*(c+d*t)*Z), no constant zeros."""
    lower, upper = interval
    if lower == upper:
        return 0, ring.constant(1)
    common = valuation_two(gcd(a, b))
    a, b = a >> common, b >> common
    scale_bits -= common
    z = list(ring.constant())
    if ring.size > 1:
        z[1] = (1, 0)
    z = tuple(z)
    assert scale_bits >= 1
    if b % 2 == 0:
        assert a % 2
        unit = ring.add(ring.constant(a), ring.scale(z, (1 << scale_bits) * c))
        step = ring.add(ring.constant(b), ring.scale(z, (1 << scale_bits) * d))
        return common * (upper - lower), unit_progression_jet(ring, interval, unit, step)
    bound = max(abs(a + b * lower), abs(a + b * (upper - 1))).bit_length()
    total_valuation, product = common * (upper - lower), ring.constant(1)
    covered = 0
    for v in range(bound):
        stride = 1 << (v + 1)
        root = ((1 << v) - a) * pow(b, -1, stride) % stride
        start = -((root - lower) // stride)
        stop = -((root - upper) // stride)
        if start == stop:
            continue
        assert scale_bits > v
        odd_base = (a + b * root) >> v
        unit = ring.add(ring.constant(odd_base),
                         ring.scale(z, (1 << (scale_bits - v)) * (c + d * root)))
        step = ring.add(ring.constant(2 * b), ring.scale(z, (1 << (scale_bits + 1)) * d))
        product = ring.multiply(product, unit_progression_jet(ring, (start, stop), unit, step))
        total_valuation += v * (stop - start)
        covered += stop - start
    assert covered == upper - lower
    return total_valuation, product


def normalized_product_jet(interval, center, precision, degree, number=None):
    """Return F(center+2^k*Z)=2^valuation*Z^zero_order*unit, with unit(0) odd."""
    lower, upper = interval
    assert 1 <= lower <= upper and (number is None or number >= 1)
    ring = QuadraticJets(precision, degree)
    a, b, c, d = (center, -1, 1, 0) if number is None else (number, -center, 0, -1)
    largest = max(abs(a + b * lower), abs(a + b * (upper - 1)))
    scale_bits = max(1, largest.bit_length())
    total_valuation, zero_order, product = 0, 0, ring.constant(1)
    pieces = [interval]
    if b and (-a) % b == 0 and lower <= -a // b < upper:
        root = -a // b
        coefficient = c + d * root
        assert coefficient
        extra = valuation_two(coefficient)
        total_valuation = scale_bits + extra
        zero_order = 1
        product = ring.constant(coefficient >> extra)
        pieces = [(lower, root), (root + 1, upper)]
    for piece in pieces:
        value, unit = nonzero_affine_jet(ring, piece, a, b, c, d, scale_bits)
        total_valuation += value
        product = ring.multiply(product, unit)
    assert product[0][0] % 2 and all(b == 0 for _, b in product)
    return {'interval': list(interval), 'center': center, 'number': number,
            'precision': precision, 'jet_degree': degree, 'scale_bits': scale_bits,
            'valuation': total_valuation, 'zero_order': zero_order,
            'unit': [a for a, _ in product]}
