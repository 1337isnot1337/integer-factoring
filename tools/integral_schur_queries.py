#!/usr/bin/env python3
"""Streamed integral Schur factor queries, valid before every ring reduction.

These routines retain no partition-indexed cache.  They use a fixed
number of streamed partition/grid loops and polynomial-size integers.
The exact algorithms are research implementations, not practical large
integer factoring software.
"""

from collections import Counter
from itertools import product
from math import comb, factorial, prod

from check_prime_power_gram_involution import (
    character_finite_difference, power_sum_mode, stream_partitions,
)


def denominator(partition):
    return prod(part ** count * factorial(count)
                for part, count in Counter(partition).items())


def translated_power_sum(parameter, output):
    """Coefficient of p_output in exp(D)p_parameter, without storing it."""
    if len(parameter) != len(output) or sum(output) > sum(parameter):
        return 0
    variables = sorted(Counter(output).items())
    parameter_counts = Counter(parameter)
    accumulated = 0
    for point in product(*(range(count + 1) for _, count in variables)):
        value = 1
        for mode, count in parameter_counts.items():
            linear = sum((comb(mode, target) if target <= mode else 0) * coordinate
                         for (target, _), coordinate in zip(variables, point))
            value *= linear ** count
        sign_weight = prod((-1) ** (count - coordinate) * comb(count, coordinate)
                           for (_, count), coordinate in zip(variables, point))
        accumulated += sign_weight * value
    divisor = prod(factorial(count) for _, count in variables)
    assert accumulated % divisor == 0
    return accumulated // divisor


def schur_mode(parameter, output, input_):
    """Coefficient of s_output in s_parameter[|output|-|input|]s_input."""
    weight, input_weight = sum(parameter), sum(input_)
    common = factorial(weight) * factorial(input_weight)
    accumulated = 0
    for mu in stream_partitions(weight):
        first = character_finite_difference(parameter, mu)
        if not first:
            continue
        for nu in stream_partitions(input_weight):
            second = character_finite_difference(input_, nu)
            if not second:
                continue
            coefficient = first * second * (common // (denominator(mu) * denominator(nu)))
            for rho in stream_partitions(sum(output)):
                value = power_sum_mode(mu, rho, nu)
                if value:
                    accumulated += coefficient * value * character_finite_difference(output, rho)
    assert accumulated % common == 0
    return accumulated // common


def translated_angle_mode(parameter, index, output, input_):
    """One integral angle-mode entry of exp(D) omega(s_parameter).

The component of parameter weight w is multiplied by binom(h+w,index).
The outer sign (-1)^index is left to the crossing factorization.
"""
    weight, input_weight = sum(parameter), sum(input_)
    if index < 0 or index > 2 * weight:
        return 0
    common = factorial(weight) * factorial(input_weight)
    accumulated = 0
    for mu in stream_partitions(weight):
        first = (-1) ** (weight - len(mu)) * character_finite_difference(parameter, mu)
        if not first:
            continue
        for lowered_weight in range(max(len(mu), index - weight), weight + 1):
            angle = comb(weight + lowered_weight, index)
            if not angle:
                continue
            for lowered in stream_partitions(lowered_weight):
                if len(lowered) != len(mu):
                    continue
                translation = translated_power_sum(mu, lowered)
                if not translation:
                    continue
                for nu in stream_partitions(input_weight):
                    second = character_finite_difference(input_, nu)
                    if not second:
                        continue
                    coefficient = (first * angle * translation * second
                                   * (common // (denominator(mu) * denominator(nu))))
                    for rho in stream_partitions(sum(output)):
                        value = power_sum_mode(lowered, rho, nu)
                        if value:
                            accumulated += coefficient * value * character_finite_difference(output, rho)
    assert accumulated % common == 0
    return accumulated // common
