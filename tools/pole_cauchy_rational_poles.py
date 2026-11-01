#!/usr/bin/env python3
"""Implicit inverse contractions for two arbitrary external integer poles.

Even denominator values are permitted. Exact valuations and normalized
affine products replace a divergent formal-series substitution.
"""

from __future__ import annotations

from math import gcd

from interval_unit_jets import QuadraticJets
from modular_interval_jets import interval_elementaries_mod, unit_progression_jet
from pole_cauchy_conditioning import valuation_two
from pole_cauchy_domain_precision import inverse_precision_certificate


def divisibility_count(first: int, step: int, length: int, level: int) -> int:
    modulus = 1 << level
    common = gcd(step, modulus)
    if first % common:
        return 0
    reduced_modulus = modulus // common
    residue = (0 if reduced_modulus == 1 else
               (-first // common) * pow(step // common, -1, reduced_modulus) % reduced_modulus)
    return 0 if residue >= length else 1 + (length - 1 - residue) // reduced_modulus


def affine_valuation(first: int, step: int, length: int, shift_bits: int | None) -> int:
    """Sum min(v2(first+step*i),k), or the exact integer-product valuation."""
    if length < 0:
        raise ValueError("Require a nonnegative length")
    if not length:
        return 0
    if shift_bits is None:
        if step == 0 and first == 0:
            raise ValueError("An integer product contains zero")
        if step and (-first) % step == 0 and 0 <= -first // step < length:
            raise ValueError("An integer product contains zero")
        cutoff = max(abs(first), abs(first + step * (length - 1))).bit_length()
    else:
        if not isinstance(shift_bits, int) or shift_bits < 0:
            raise ValueError("Require a nonnegative displacement exponent")
        cutoff = shift_bits
    return sum(divisibility_count(first, step, length, level) for level in range(1, cutoff + 1))


def affine_quadratic_product(ring, first: int, step: int, length: int,
                            shift_bits: int | None, *, elementary_cache=None) -> tuple[int, tuple]:
    """Valuation and unit of prod(first+step*i+2^k*w), without listing i.

When k is absent the factors are nonzero integers. At integer valuation
>=k, the quadratic factors have valuation exactly k; split the remaining
unit progression by parity so its step is even.
"""
    if length < 0:
        raise ValueError("Require a nonnegative length")
    if not length:
        return 0, ring.constant(1)
    if step == 0:
        if first == 0:
            if shift_bits is None:
                raise ValueError("An integer product contains zero")
            value, unit = shift_bits, ring.constant(0, 1)
        else:
            value = valuation_two(first)
            if shift_bits is not None:
                value = min(value, shift_bits)
            unit = ring.constant(first >> value,
                                 0 if shift_bits is None else 1 << (shift_bits - value))
        return value * length, ring.power(unit, length)
    if shift_bits is None and (-first) % step == 0 and 0 <= -first // step < length:
        raise ValueError("An integer product contains zero")
    common = valuation_two(gcd(first, step))
    if shift_bits is not None:
        common = min(common, shift_bits)
    base, increment = first >> common, step >> common
    clipped = None if shift_bits is None else shift_bits - common
    exponent, product = common * length, ring.constant(1)

    def include_unit(start, increment, count, imaginary):
        nonlocal product
        if count:
            assert increment % 2 == 0
            unit = ring.constant(start, imaginary)
            if elementary_cache is None:
                factor = unit_progression_jet(ring, (0, count), unit, ring.constant(increment))
            else:
                # Scalar specialization of the same short-jet product. Reuse
                # elementary coefficients when several valuation classes
                # have the same interval length and precision.
                degree = (0 if increment == 0 else
                          min((ring.precision - 1) // valuation_two(increment), count))
                key = (count, degree, ring.precision)
                if key not in elementary_cache:
                    elementary_cache[key] = interval_elementaries_mod(0, count, degree, ring.precision)
                elementary = elementary_cache[key]
                ratio = ring.multiply(ring.constant(increment), ring.inverse(unit))
                series, power = ring.constant(1), ring.constant(1)
                for j in range(1, degree + 1):
                    power = ring.multiply(power, ratio)
                    series = ring.add(series, ring.scale(power, elementary[j]))
                factor = ring.multiply(ring.power(unit, count), series)
            product = ring.multiply(product, factor)

    def include_final_ball(start, increment, count):
        if increment % 2 == 0:
            include_unit(start, increment, count, 1)
        else:
            for parity in (0, 1):
                include_unit(start + increment * parity, 2 * increment,
                             (count + 1 - parity) // 2, 1)

    if clipped == 0:
        include_final_ball(base, increment, length)
        return exponent, product
    if increment % 2 == 0:
        assert base % 2
        include_unit(base, increment, length, 0 if clipped is None else 1 << clipped)
        return exponent, product
    bound = (max(abs(base), abs(base + increment * (length - 1))).bit_length()
             if clipped is None else clipped)
    covered = 0
    for value in range(bound):
        stride = 1 << (value + 1)
        residue = ((1 << value) - base) * pow(increment, -1, stride) % stride
        count = 0 if residue >= length else 1 + (length - 1 - residue) // stride
        include_unit((base + increment * residue) >> value, 2 * increment, count,
                     0 if clipped is None else 1 << (clipped - value))
        exponent += value * count
        covered += count
    if clipped is not None:
        stride = 1 << clipped
        residue = -base * pow(increment, -1, stride) % stride
        count = 0 if residue >= length else 1 + (length - 1 - residue) // stride
        include_final_ball((base + increment * residue) >> clipped, increment, count)
        exponent += clipped * count
        covered += count
    assert covered == length
    return exponent, product


def pole_denominator_bound(domain, pole):
    """The exact maximum v2(a-pole), with pole membership rejected."""
    if not isinstance(pole, int):
        raise ValueError("The external pole must be an integer")
    for first, length in domain:
        if first <= pole <= first + 2 * (length - 1) and (pole - first) % 2 == 0:
            raise ValueError("An external pole lies on its matrix domain")
    maximum_distance = max(max(abs(first - pole), abs(first + 2 * (length - 1) - pole))
                           for first, length in domain)
    low, high = 0, maximum_distance.bit_length() + 1
    while high - low > 1:
        level = (low + high) // 2
        if any(divisibility_count(first - pole, 2, length, level) for first, length in domain):
            low = level
        else:
            high = level
    return low


def inverse_two_pole_contraction(n, left, right, k, row_pole, column_pole,
                                 absolute_precision):
    """sum_ab (M^-1)[b,a]/((a-row_pole)(b-column_pole)) to absolute precision."""
    if not isinstance(absolute_precision, int) or absolute_precision < 1:
        raise ValueError("Require positive absolute precision")
    certificate = inverse_precision_certificate(n, left, right, k)
    size = certificate["domain_size"]
    c, d = row_pole, column_pole
    u, v = pole_denominator_bound(left, c), pole_denominator_bound(right, d)
    difference = c * d - n
    displacement_valuation = k if difference == 0 else min(valuation_two(difference), k)
    specifications = []
    for first, length in left:
        specifications.append((1, n - d * first, -2 * d, length, k))
        specifications.append((-1, first - c, 2, length, None))
    for first, length in right:
        specifications.append((1, n - c * first, -2 * c, length, k))
        specifications.append((-1, first - d, 2, length, None))
    product_valuation = sum(sign * affine_valuation(first, step, length, shift)
                            for sign, first, step, length, shift in specifications)
    certified_lower = min(0, k + displacement_valuation
                          - certificate["inverse_loss_bound"] - u - v)
    if product_valuation < certified_lower:
        raise ArithmeticError("The product valuation violates the inverse conditioning bound")
    target_precision = absolute_precision + k + displacement_valuation
    scale_floor = min(0, product_valuation)
    working_precision = target_precision - scale_floor
    ring = QuadraticJets(working_precision, 0)
    bracket = ring.constant(1 << (-scale_floor))
    product_precision = max(0, target_precision - product_valuation)
    elementary_cache = {}
    if product_precision:
        product_ring = QuadraticJets(product_precision, 0)
        unit = product_ring.constant(1)
        checked_valuation = 0
        for sign, first, step, length, shift in specifications:
            value, factor = affine_quadratic_product(product_ring, first, step, length, shift,
                                                      elementary_cache=elementary_cache)
            checked_valuation += sign * value
            unit = product_ring.multiply(unit, factor if sign == 1 else product_ring.inverse(factor))
        assert checked_valuation == product_valuation
        unit = product_ring.multiply(unit, product_ring.inverse(
            product_ring.power(product_ring.constant(n, 1 << k), size)))
        bracket = ring.add(bracket, ring.scale(ring.constant(*unit[0]),
                                              -(1 << (product_valuation - scale_floor))))
    denominator_unit = ring.constant(difference >> displacement_valuation,
                                     -(1 << (k - displacement_valuation)))
    prefactor = ring.multiply(ring.constant(n, 1 << k),
                              ring.inverse(ring.multiply(ring.constant(0, 1), denominator_unit)))
    result = ring.multiply(prefactor, bracket)
    return {"residue": list(result[0]), "scale_exponent": scale_floor - k - displacement_valuation,
            "working_precision": working_precision, "absolute_precision": absolute_precision,
            "product_unit_precision": product_precision, "product_valuation": product_valuation,
            "product_valuation_lower_bound": certified_lower,
            "distinct_elementary_intervals": len(elementary_cache),
            "row_pole": c, "column_pole": d, "row_denominator_bound": u,
            "column_denominator_bound": v, "displacement_denominator_valuation": displacement_valuation,
            "precision_certificate": certificate, "enumerated_matrix_nodes": 0,
            "scope": "Two external integer poles; joint product-pole kernels are not evaluated."}
