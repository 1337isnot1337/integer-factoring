#!/usr/bin/env python3
"""Recover the smallest prime from supplied suffix counts in one prime field.

This is a conditional wrapper. The normalized Bezout coefficient provider is
an explicit argument and is not implemented here. A short trial-division
budget makes counts in successively certified prime-only intervals exact.
"""

from __future__ import annotations

from math import isqrt

from prime_bezout_counting import is_prime, suffix_rectangle_contract


def prime_above_bits(bits: int) -> int:
    """Find one odd prime above bits by a polynomially bounded literal search."""
    if not isinstance(bits, int) or bits < 1:
        raise ValueError("Require a positive input bit length")
    candidate = max(3, bits + 1)
    if candidate % 2 == 0:
        candidate += 1
    bound = 16 * (bits + 1) ** 4
    while candidate <= bound:
        if is_prime(candidate):
            return candidate
        candidate += 2
    raise AssertionError("The elementary primorial existence bound was violated")


def one_prime_recovery_plan(n: int, prime: int | None = None) -> dict:
    if not isinstance(n, int) or n < 2:
        raise ValueError("Require an integer N >= 2")
    bits = n.bit_length()
    if prime is None:
        prime = prime_above_bits(bits)
    if not isinstance(prime, int) or prime < 3 or not is_prime(prime):
        raise ValueError("Use an odd prime characteristic")
    trial_bits = (bits + prime - 1) // prime
    initial_lower_bound = 1 << trial_bits
    trial_limit = min(isqrt(n), max(initial_lower_bound - 1, prime))
    return {
        "n": n, "input_bits": bits, "prime": prime,
        "trial_bits": trial_bits,
        "initial_lower_bound": initial_lower_bound,
        "trial_limit": trial_limit,
        "factor_multiplicity_bound": (bits - 1) // trial_bits,
        "coefficient_precision": 1,
        "provider_obligation": "Supply the normalized leading Bezout coefficient for every requested residue rectangle; no provider is implemented.",
    }


def smallest_prime_from_one_field(n: int, coefficient_oracle,
                                  prime: int | None = None) -> dict:
    """Conditional smallest-prime recovery with a tunable single field.

    Trial division is explicit, and may be exponential in input_bits/prime.
    The default prime is above input_bits, making this preprocessing polynomial.
    All later oracle queries use precisely the same prime and precision one.
    """
    plan = one_prime_recovery_plan(n, prime)
    prime = plan["prime"]
    trial_tests = 0
    for candidate in range(2, plan["trial_limit"] + 1):
        trial_tests += 1
        if n % candidate == 0:
            return {**plan, "smallest_prime_factor": candidate,
                    "trial_tests": trial_tests, "coefficient_queries": 0,
                    "suffix_queries": 0, "prime_intervals": [],
                    "scope": "Explicit preliminary trial division"}
    if plan["trial_limit"] == isqrt(n):
        return {**plan, "smallest_prime_factor": n,
                "trial_tests": trial_tests, "coefficient_queries": 0,
                "suffix_queries": 0, "prime_intervals": [],
                "scope": "Complete preliminary trial division"}
    assert n % prime and plan["factor_multiplicity_bound"] < prime

    coefficients = 0
    suffix_cache = {}
    intervals = []

    def suffix(threshold):
        nonlocal coefficients
        if threshold >= n:
            return 0
        if threshold not in suffix_cache:
            total = 0
            for residue in range(1, prime):
                contract = suffix_rectangle_contract(n, threshold, prime, residue)
                supplied = coefficient_oracle(contract)
                if not isinstance(supplied, int):
                    raise ValueError("Supply a real integer residue modulo the query prime")
                total += supplied
                coefficients += 1
            suffix_cache[threshold] = total % prime
        return suffix_cache[threshold]

    def result(factor):
        if factor < 2 or n % factor:
            raise ValueError("The supplied coefficients did not locate an integer factor")
        return {**plan, "smallest_prime_factor": factor,
                "trial_tests": trial_tests, "coefficient_queries": coefficients,
                "suffix_queries": len(suffix_cache), "prime_intervals": intervals,
                "scope": "Conditional on the supplied global leading coefficient provider"}

    lower = plan["initial_lower_bound"]
    while lower * lower <= n:
        upper = lower * lower
        base = suffix(lower)
        count = (base - suffix(upper)) % prime
        intervals.append({"lower": lower, "upper": upper,
                          "exact_prime_divisor_count": count})
        if count:
            left, right = lower, upper - 1
            while left < right:
                middle = (left + right) // 2
                if (base - suffix(middle + 1)) % prime:
                    right = middle
                else:
                    left = middle + 1
            return result(left)
        lower = upper
    return result(n)
