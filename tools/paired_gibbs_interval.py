#!/usr/bin/env python3
"""Outward integer intervals for a fixed-scale paired Gibbs pressure.

Orders and physical energies are exact integers. Scale and trial ratio
are supplied as rational dyadic intervals. This certifies a displayed
pressure sign or root bracket, not a minimum over scales or heights.
"""
from __future__ import annotations


def log_positive(arithmetic, value):
    """Logarithm of an interval known to be at least one, with range reduction."""
    assert value[0] >= arithmetic.scale
    power = 0
    while value[0] >= (arithmetic.scale << (power+1)):
        power += 1
    reduced = arithmetic.divide_integer(value, 1 << power)
    logarithm = arithmetic.log_one_plus(arithmetic.subtract(reduced, arithmetic.one))
    return arithmetic.add(logarithm, arithmetic.multiply(arithmetic.rational(power),
                                                        arithmetic.log_one_plus(arithmetic.one)))


def partition(arithmetic, energy, scale):
    minimum = min(energy)
    total = arithmetic.zero
    for value in energy:
        total = arithmetic.add(total, arithmetic.exp_negative(
            arithmetic.multiply(arithmetic.rational(value-minimum), scale)))
    log_value = arithmetic.subtract(log_positive(arithmetic, total),
                                    arithmetic.multiply(arithmetic.rational(minimum), scale))
    return total, log_value


def pressure(arithmetic, profiles, energies_x, energies_y, scale, rho):
    rate = arithmetic.rational(2*rho.numerator, rho.denominator)
    result = arithmetic.zero
    for profile in profiles:
        result = arithmetic.add(result, partition(arithmetic, profile, scale)[1])
    physical = arithmetic.zero
    for energies in (energies_x, energies_y):
        physical = arithmetic.add(physical, partition(arithmetic, energies,
                                  arithmetic.divide(scale, rate))[1])
    return arithmetic.subtract(result, arithmetic.multiply(rate, physical))
