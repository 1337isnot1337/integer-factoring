#!/usr/bin/env python3
"""Integer-only threshold centers for the Nicolas descriptor-cover theorem.

This module implements the proved O(log(p+a)) workspace arithmetic step.
It does not instantiate the unknown numerical cutoff constants in Nicolas's
asymptotic theorem, and is not a practical polynomial divisor-max oracle.
"""
from __future__ import annotations


def reduced_log(numerator: int, denominator: int, precision: int) -> int:
    """A lower fixed-point approximation to log(numerator/denominator).

    The ratio must lie in [1,2].  The returned value is scaled by 2^precision.
    Every multiplication is truncated down immediately.  The atanh series
    has argument at most 1/3 and uses exactly `precision` terms.
    """
    if not (0 < denominator <= numerator <= 2 * denominator):
        raise ValueError("the reduced ratio must lie in [1,2]")
    if precision < 1:
        raise ValueError("precision must be positive")
    scale = 1 << precision
    z = ((numerator - denominator) * scale) // (numerator + denominator)
    z_squared = (z * z) // scale
    power = z
    total = 0
    for j in range(precision):
        total += power // (2 * j + 1)
        power = (power * z_squared) // scale
    return 2 * total


def integer_log(value: int, precision: int, log_two: int | None = None) -> int:
    """Lower fixed-point natural logarithm of a positive integer."""
    if value < 1:
        raise ValueError("value must be positive")
    exponent = value.bit_length() - 1
    if log_two is None:
        log_two = reduced_log(2, 1, precision)
    return exponent * log_two + reduced_log(value, 1 << exponent, precision)


def threshold_center(prime_bound: int, level: int) -> int:
    """Return an integer within 2 of p^(log(1+1/a)/log(2)).

    Primality of p is unnecessary for this arithmetic lemma.  The proof
    uses binary search with a fixed-point logarithmic sign test.  No real
    exponentiation, logarithm, or floating-point comparison is used.
    """
    p, a = prime_bound, level
    if p < 2 or a < 1:
        raise ValueError("p must be at least two and the level positive")
    length = p.bit_length() + a.bit_length() + 10
    precision = 32 * length
    scale = 1 << precision
    log_two = reduced_log(2, 1, precision)
    log_p = integer_log(p, precision, log_two)
    log_ratio = reduced_log(a + 1, a, precision)
    right = (log_p * log_ratio) // scale
    lo, hi = 1, p + 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        left = (integer_log(mid, precision, log_two) * log_two) // scale
        if left <= right:
            lo = mid
        else:
            hi = mid
    return lo


def is_prime(value: int) -> bool:
    """Elementary polynomial-time, logarithmic-space prime test."""
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def prime_count(value: int) -> int:
    return sum(is_prime(candidate) for candidate in range(2, value + 1))


def descriptor_centers(prime_bound: int, max_level: int) -> list[int]:
    """Convenience wrapper for finite audits; a list is not needed in the proof."""
    return [prime_count(threshold_center(prime_bound, level))
            for level in range(1, max_level + 1)]
