#!/usr/bin/env python3
"""Factoring reduction to supplied prime-adic leading inverse coefficients.

The coefficient provider is an explicit argument. No global inverse evaluator
is implemented. The wrapper performs only trial division by small primes,
contract construction, modular summation, and a monotone threshold search.
"""

from __future__ import annotations

from math import isqrt


def is_prime(value: int) -> bool:
    return value >= 2 and all(value % d for d in range(2, isqrt(value) + 1))


def quadratic_nonresidue(prime: int) -> int:
    if prime < 3 or not is_prime(prime):
        raise ValueError("Use an odd prime")
    return next(d for d in range(2, prime)
                if pow(d, (prime - 1) // 2, prime) == prime - 1)


def prime_rectangle_contract(n: int, prime: int, first_a: int,
                             first_b: int, size_exponent: int) -> dict:
    if (not isinstance(n, int) or n < 2 or not is_prime(prime) or prime == 2
            or n % prime == 0 or first_a < 1 or first_b < 1
            or size_exponent < 0 or first_a * first_b % prime != n % prime):
        raise ValueError("Require positive matched unit residues and a nonnegative size exponent")
    count = prime ** size_exponent
    derivative_order = prime * (count - 1) // (prime - 1) - size_exponent
    last_a = first_a + prime * (count - 1)
    last_b = first_b + prime * (count - 1)
    bound = max(abs(n - first_a * first_b), abs(n - last_a * last_b), 1)
    height, height_power = 0, 1
    while height_power <= bound:
        height += 1
        height_power *= prime
    height = max(height, size_exponent + 1)
    padding = height + 1
    return {"n": n, "prime": prime, "first_a": first_a, "first_b": first_b,
            "count": count, "size_exponent": size_exponent,
            "derivative_order": derivative_order, "height_bound": height,
            "padding_exponent": padding,
            "clearing_exponent": padding + 2 * derivative_order,
            "quadratic_nonresidue": quadratic_nonresidue(prime),
            "precision": 1, "requested_coefficients": 1,
            "requested_coefficient_degree": count - 1,
            "provider_obligation": "Supply [X^(m-1)] omega*ell^(K+2E)*V_t modulo ell, with omega^2=d and V_t*Q_t=1 modulo P. No provider is implemented."}


def suffix_rectangle_contract(n: int, threshold: int, prime: int,
                              row_residue: int) -> dict:
    if threshold < 2 or threshold > n or not 1 <= row_residue < prime:
        raise ValueError("Use a nonzero residue and a threshold between 2 and N")
    first_a = threshold + (row_residue - threshold) % prime
    column_residue = n * pow(row_residue, -1, prime) % prime
    first_b = 2 + (column_residue - 2) % prime
    size_exponent, count = 0, 1
    while count < n:
        size_exponent += 1
        count *= prime
    return prime_rectangle_contract(n, prime, first_a, first_b, size_exponent)


def small_prime_plan(n: int) -> dict:
    """Either find the smallest prime directly or supply a product exceeding N."""
    if n < 2:
        raise ValueError("Require N>=2")
    if n % 2 == 0:
        return {"smallest_prime_factor": 2, "primes": [], "product": 1}
    limit = 16 * n.bit_length() ** 2
    primes, product = [], 1
    for prime in range(3, limit + 1, 2):
        if not is_prime(prime):
            continue
        if n % prime == 0:
            return {"smallest_prime_factor": prime, "primes": primes,
                    "product": product, "prime_limit": limit}
        primes.append(prime)
        product *= prime
        if product > n:
            return {"primes": primes, "product": product, "prime_limit": limit}
    raise AssertionError("The proved elementary primorial bound was violated")


def smallest_prime_from_prime_coefficients(n: int, coefficient_oracle) -> dict:
    plan = small_prime_plan(n)
    if "smallest_prime_factor" in plan:
        return {**plan, "coefficient_queries": 0,
                "scope": "Trial division by the bounded prime list"}
    queries = 0

    def has_factor(threshold):
        nonlocal queries
        for prime in plan["primes"]:
            total = 0
            for residue in range(1, prime):
                contract = suffix_rectangle_contract(n, threshold, prime, residue)
                supplied = coefficient_oracle(contract)
                if not isinstance(supplied, int):
                    raise ValueError("Supply a real integer residue modulo the query prime")
                total += supplied
                queries += 1
            if total % prime:
                return True
        return False

    if not has_factor(2):
        return {**plan, "smallest_prime_factor": n, "coefficient_queries": queries,
                "scope": "Conditional on the supplied leading coefficient provider"}
    low, high = 2, n // 2 + 1
    while low + 1 < high:
        middle = (low + high) // 2
        if has_factor(middle):
            low = middle
        else:
            high = middle
    if n % low:
        raise ValueError("The supplied coefficients do not locate an integer divisor")
    return {**plan, "smallest_prime_factor": n // low,
            "largest_proper_divisor": low, "coefficient_queries": queries,
            "scope": "Conditional on the supplied leading coefficient provider"}
