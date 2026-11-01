#!/usr/bin/env python3
"""Unit interval products with precision-bounded Newton computations."""

from __future__ import annotations

from functools import lru_cache


def _v2(value: int) -> int:
    assert value > 0
    return (value & -value).bit_length() - 1


def prefix_power_sums_mod(stop: int, degree: int, precision: int) -> tuple[int, ...]:
    """Sum j^r over 0 <= j < stop, modulo 2^precision, for 0 <= r <= degree."""
    if stop < 0 or degree < 0 or precision < 1:
        raise ValueError("Require nonnegative endpoint and degree, positive precision")
    modulus = 1 << precision
    binomials = [1]
    unit, exponent = 1, 0
    for k in range(1, degree + 2):
        if k > stop:
            binomials.append(0)
            continue
        numerator = stop - k + 1
        numerator_v, denominator_v = _v2(numerator), _v2(k)
        exponent += numerator_v - denominator_v
        assert exponent >= 0
        unit = unit * (numerator >> numerator_v) % modulus
        unit = unit * pow(k >> denominator_v, -1, modulus) % modulus
        binomials.append((unit << exponent) % modulus if exponent < precision else 0)
    factorials = [1]
    for k in range(1, degree + 1):
        factorials.append(factorials[-1] * k % modulus)
    stirling = [1]
    result = []
    for r in range(degree + 1):
        result.append(sum(stirling[j] * factorials[j] * binomials[j + 1]
                          for j in range(r + 1)) % modulus)
        following = [0] * (r + 2)
        for j in range(r + 2):
            following[j] = ((j * stirling[j] if j <= r else 0)
                            + (stirling[j - 1] if j else 0)) % modulus
        stirling = following
    return tuple(result)


@lru_cache(maxsize=128)
def interval_elementaries_mod(lower: int, upper: int, degree: int,
                              precision: int) -> tuple[int, ...]:
    """First elementary symmetric functions modulo 2^precision.

    Newton division loses at most v2(degree!) bits. Keeping those extra
    bits makes every division exact without constructing huge exact
    elementary symmetric functions.
    """
    if lower < 0 or upper < lower or degree < 0 or precision < 1:
        raise ValueError("Invalid interval, degree, or precision")
    extra = degree - degree.bit_count()
    working_precision = precision + extra
    high = prefix_power_sums_mod(upper, degree, working_precision)
    low = prefix_power_sums_mod(lower, degree, working_precision)
    working_modulus = 1 << working_precision
    powers = [(a - b) % working_modulus for a, b in zip(high, low)]
    elementary = [1]
    lost = 0
    for k in range(1, degree + 1):
        division_loss = _v2(k)
        lost += division_loss
        target_precision = working_precision - lost
        numerator_modulus = 1 << (target_precision + division_loss)
        numerator = sum((-1) ** (j - 1) * elementary[k - j] * powers[j]
                        for j in range(1, k + 1)) % numerator_modulus
        assert numerator % (1 << division_loss) == 0
        modulus = 1 << target_precision
        value = (numerator >> division_loss) * pow(k >> division_loss, -1, modulus) % modulus
        elementary.append(value)
    modulus = 1 << precision
    return tuple(value % modulus for value in elementary)


def unit_progression_jet(ring, interval, unit, step):
    """Same contract as interval_unit_jets.unit_progression_jet."""
    start, stop = interval
    assert 0 <= start <= stop and all(a % 2 == b % 2 == 0 for a, b in step)
    length = stop - start
    nonzero = [(index, coordinate) for index, pair in enumerate(step)
               for coordinate in pair if coordinate]
    if not nonzero or not length:
        return ring.power(unit, length)
    degree = min((ring.precision - 1) // min(_v2(value) for _, value in nonzero), length)
    formal_order = min(index for index, _ in nonzero)
    if formal_order:
        degree = min(degree, (ring.size - 1) // formal_order)
    increment = ring.multiply(step, ring.inverse(unit))
    elementary = interval_elementaries_mod(start, stop, degree, ring.precision)
    series, power = ring.constant(1), ring.constant(1)
    for k in range(1, degree + 1):
        power = ring.multiply(power, increment)
        series = ring.add(series, ring.scale(power, elementary[k]))
    return ring.multiply(ring.power(unit, length), series)
