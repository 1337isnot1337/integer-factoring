#!/usr/bin/env python3
"""Whole-matrix inverse contractions with separated rational pole weights.

No domain is enumerated. Exact powers of two are retained as binary
exponents, and the returned residue has the requested absolute precision.
"""

from __future__ import annotations

from math import gcd

from interval_unit_jets import QuadraticJets
from modular_interval_jets import unit_progression_jet
from pole_cauchy_conditioning import affine_product, domain_size, residue_count


def _v2(value: int) -> int:
    if not value:
        raise ValueError("A zero integer has no finite valuation")
    value = abs(value)
    return (value & -value).bit_length() - 1


def _power_two(ring, exponent: int) -> int:
    return 0 if exponent >= ring.precision else 1 << exponent


def _unit_progression(ring, first: int, step: int, length: int,
                      w_coefficient: int) -> tuple:
    """Multiply a unit progression, splitting an odd step into two classes."""
    if not length:
        return ring.constant(1)
    if step % 2 == 0:
        return unit_progression_jet(ring, (0, length),
                                    ring.constant(first, w_coefficient),
                                    ring.constant(step))
    result = ring.constant(1)
    for parity in (0, 1):
        count = (length + 1 - parity) // 2
        if count:
            result = ring.multiply(result, unit_progression_jet(
                ring, (0, count), ring.constant(first + step * parity, w_coefficient),
                ring.constant(2 * step)))
    return result


def normalized_shifted_affine_product(ring, first: int, step: int,
                                      length: int, shift_bits: int) -> tuple:
    """Valuation and unit of product(first+step*i+2^shift_bits*w).

    Unlike the separated-pole primitive, this allows nonzero rational
    parts divisible by 2^shift_bits. Such a tail is a unit progression
    after division by 2^shift_bits, since w and 1+w are both units.
    """
    if step == 0 or length < 0 or shift_bits < 0:
        raise ValueError("Require nonzero step and nonnegative length and shift")
    if not length:
        return 0, ring.constant(1)
    common = _v2(gcd(first, step))
    if common >= shift_bits:
        result = _unit_progression(ring, first >> shift_bits, step >> shift_bits,
                                   length, 1)
        return shift_bits * length, result
    base, increment = first >> common, step >> common
    remaining = shift_bits - common
    if increment % 2 == 0:
        assert base % 2
        result = _unit_progression(ring, base, increment, length,
                                   _power_two(ring, remaining))
        return common * length, result

    exponent = common * length
    result = ring.constant(1)
    covered = 0
    for valuation in range(remaining):
        stride = 1 << (valuation + 1)
        root = ((1 << valuation) - base) * pow(increment, -1, stride) % stride
        count = 0 if root >= length else 1 + (length - 1 - root) // stride
        if count:
            unit = _unit_progression(ring, (base + increment * root) >> valuation,
                                     2 * increment, count,
                                     _power_two(ring, remaining - valuation))
            result = ring.multiply(result, unit)
            exponent += valuation * count
            covered += count
    stride = 1 << remaining
    root = -base * pow(increment, -1, stride) % stride
    count = 0 if root >= length else 1 + (length - 1 - root) // stride
    if count:
        unit = _unit_progression(ring, (base + increment * root) >> remaining,
                                 increment, count, 1)
        result = ring.multiply(result, unit)
        exponent += remaining * count
        covered += count
    assert covered == length
    assert result[0][0] % 2 or result[0][1] % 2
    return exponent, result


def _validate(n: int, left: tuple, right: tuple, xi: int, eta: int,
              k: int, absolute_precision: int) -> tuple:
    if n < 1 or n % 2 == 0 or xi < 1 or eta < 1 or k < 1 or absolute_precision < 1:
        raise ValueError("Require positive odd N, positive poles, displacement, and precision")
    horizon = 0
    for domain, pole in ((left, xi), (right, eta)):
        if not domain:
            raise ValueError("Domains must be nonempty")
        previous = 0
        for first, length in domain:
            if first <= previous or first % 2 == 0 or length < 1:
                raise ValueError("Require ordered disjoint positive odd progressions")
            last = first + 2 * (length - 1)
            if first <= pole <= last and (pole - first) % 2 == 0:
                raise ValueError("A rational pole must lie outside its node domain")
            if pole % 2:
                horizon = max(horizon, abs(pole - first).bit_length(),
                              abs(pole - last).bit_length())
            previous = last
    size = domain_size(left)
    if size != domain_size(right):
        raise ValueError("Domains must have equal positive size")
    if k < horizon:
        raise ValueError("Displacement must cover the valuations of integer pole distances")
    return size, horizon


def _ratio_valuation(n: int, left: tuple, right: tuple,
                      xi: int, eta: int, k: int) -> int:
    result = 0
    for level in range(1, k + 1):
        modulus = 1 << level
        if eta % 2:
            result += residue_count(left, n * pow(eta, -1, modulus), level)
            result -= residue_count(right, eta, level)
        if xi % 2:
            result += residue_count(right, n * pow(xi, -1, modulus), level)
            result -= residue_count(left, xi, level)
    return result


def inverse_resolvent_bilinear_form(n: int, left: tuple, right: tuple,
                                   xi: int, eta: int, k: int,
                                   absolute_precision: int) -> dict:
    """Return sum(M^-1[b,a]/((xi-a)*(eta-b))) with absolute 2-adic precision.

    If the output scale is s and residue is U modulo 2^P, the true
    scalar differs from 2^s U by valuation at least absolute_precision.
    The scale need not be the scalar's exact valuation.
    """
    size, horizon = _validate(n, left, right, xi, eta, k, absolute_precision)
    ratio_valuation = _ratio_valuation(n, left, right, xi, eta, k)
    parts = len(left) + len(right)
    assert abs(ratio_valuation) <= parts * (k - 1)
    rational_part = xi * eta - n
    displacement_valuation = min(k, _v2(rational_part)) if rational_part else k
    scale = min(0, ratio_valuation) - k - displacement_valuation
    precision = absolute_precision - scale
    ring = QuadraticJets(precision, 0)
    z = ring.constant(n, _power_two(ring, k))
    unit = ring.constant(1)
    checked_valuation = 0

    def include(exponent_unit: tuple, sign: int) -> None:
        nonlocal checked_valuation, unit
        exponent, factor = exponent_unit
        checked_valuation += sign * exponent
        unit = ring.multiply(unit, factor if sign == 1 else ring.inverse(factor))

    for first, length in left:
        include(normalized_shifted_affine_product(ring, n - eta * first,
                                                 -2 * eta, length, k), 1)
        include(affine_product(ring, xi - first, -2, length, None), -1)
    for first, length in right:
        include(normalized_shifted_affine_product(ring, n - xi * first,
                                                 -2 * xi, length, k), 1)
        include(affine_product(ring, eta - first, -2, length, None), -1)
    assert checked_valuation == ratio_valuation
    unit = ring.multiply(unit, ring.inverse(ring.power(z, size)))
    if ratio_valuation >= 0:
        difference = ring.add(ring.constant(1),
                              ring.scale(unit, -_power_two(ring, ratio_valuation)))
    else:
        difference = ring.add(ring.constant(_power_two(ring, -ratio_valuation)),
                              ring.scale(unit, -1))
    displaced_unit = ring.constant(rational_part >> displacement_valuation,
                                   -_power_two(ring, k - displacement_valuation))
    denominator = ring.multiply(ring.constant(0, 1), displaced_unit)
    result = ring.multiply(difference, ring.multiply(z, ring.inverse(denominator)))
    return {"scale_exponent": scale, "residue": list(result[0]),
            "precision": precision, "absolute_precision": absolute_precision,
            "ratio_valuation": ratio_valuation,
            "query_displacement_valuation": displacement_valuation,
            "integer_distance_horizon": horizon}
