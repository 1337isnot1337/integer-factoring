#!/usr/bin/env python3
"""Finite exact checks of the rough-input modular transfer.

The prefix oracle here is computed from an explicit divisor list. This
checks the preprocessing/recovery wrapper, not a faster factoring evaluator.
"""

from functools import lru_cache
from math import isqrt
import json


def ceil_div(first: int, second: int) -> int:
    return (first + second - 1) // second


def ceil_log2(value: int) -> int:
    return (value - 1).bit_length()


def least_factor(number: int) -> int:
    for candidate in range(2, isqrt(number) + 1):
        if number % candidate == 0:
            return candidate
    return number


def omega_and_tau(number: int) -> tuple[int, int]:
    omega, tau, candidate = 0, 1, 2
    while candidate * candidate <= number:
        exponent = 0
        while number % candidate == 0:
            number //= candidate
            exponent += 1
        if exponent:
            omega += exponent
            tau *= exponent + 1
        candidate += 1
    if number > 1:
        omega += 1
        tau *= 2
    return omega, tau


PRIMES = [3]


@lru_cache(maxsize=None)
def odd_primes(count: int) -> tuple[int, ...]:
    candidate = PRIMES[-1] + 2
    while len(PRIMES) < count:
        if all(candidate % prime for prime in range(3, isqrt(candidate) + 1, 2)):
            PRIMES.append(candidate)
        candidate += 2
    return tuple(PRIMES[:count])


def parity_support(exponents: list[int], modulus: int) -> set[int]:
    result = set()
    for exponent in exponents:
        residue = exponent % modulus
        if residue in result:
            result.remove(residue)
        else:
            result.add(residue)
    return result


def recover(number: int, threshold_bits: int) -> tuple[int, int, int]:
    bits = number.bit_length()
    for candidate in range(2, min(1 << threshold_bits, number) + 1):
        if number % candidate == 0:
            return candidate, 0, 0

    q = ceil_div(bits, threshold_bits)
    bound = 1 << q
    count = bits * (bound - 1) + 1
    moduli = odd_primes(count)
    assert moduli[-1] <= 4 * (count + 1) ** 2
    assert moduli[-1] <= 8 * count * ceil_log2(count + 2)
    divisors = [value for value in range(2, number + 1) if number % value == 0]
    assert len(divisors) <= bound
    prefix, query_count = 0, 0
    for remaining in range(bits - 1, -1, -1):
        zero_child = 2 * prefix
        permitted = [value for value in divisors if value >> remaining == zero_child]
        found = False
        for modulus in moduli:
            query_count += 1
            if parity_support(permitted, modulus):
                found = True
                break
        assert found == bool(permitted)
        prefix = zero_child if found else zero_child + 1
        assert any(value >> remaining == prefix for value in divisors)
    return prefix, query_count, moduli[-1]


def verify() -> dict:
    rough_cases = trivial_cases = 0
    for number in range(2, 8192):
        bits = number.bit_length()
        first = least_factor(number)
        omega, tau = omega_and_tau(number)
        for threshold_bits in range(1, bits + 1):
            if first <= 1 << threshold_bits:
                trivial_cases += 1
                continue
            q = ceil_div(bits, threshold_bits)
            assert threshold_bits * omega < bits
            assert tau <= 1 << omega <= 1 << q
            rough_cases += 1

    parameter_cases = half_cases = 0
    for bits in range(2, 100001):
        root = isqrt(bits)
        threshold = root if root * root == bits else root + 1
        assert ceil_div(bits, threshold) <= threshold
        assert threshold + ceil_div(bits, threshold) <= 2 * threshold
        parameter_cases += 1
        q = isqrt(bits.bit_length() - 1)
        threshold = ceil_div(bits, q)
        assert q >= 1 and ceil_div(bits, threshold) <= q
        bound = 1 << q
        count = bits * (bound - 1) + 1
        limit = 4 * (count + 1) ** 2
        assert bound <= bits
        assert count <= bits * bits
        assert limit <= 4 * bits ** 4
        half_cases += 1

    prime_supply_cases = 0
    for count in range(2, 100001):
        logarithm = ceil_log2(count + 2)
        limit = 8 * count * logarithm
        assert limit - 1 > (count + 2) * (2 * logarithm + 2)
        assert limit <= 1 << (2 * logarithm + 2)
        prime_supply_cases += 1
    actual_primes = odd_primes(6000)
    for count, prime in enumerate(actual_primes, 1):
        assert prime <= 8 * count * ceil_log2(count + 2)

    recovery_cases = prime_endpoints = cyclic_queries = max_prime = 0
    for number in range(2, 512):
        expected = least_factor(number)
        for threshold_bits in range(1, number.bit_length() + 1):
            found, queries, largest = recover(number, threshold_bits)
            assert found == expected
            recovery_cases += 1
            prime_endpoints += expected == number and queries > 0
            cyclic_queries += queries
            max_prime = max(max_prime, largest)

    # An individual cyclic reduction can vanish through cancellation.
    assert parity_support([3, 15], 3) == set()
    assert parity_support([3, 15], 5) != set()

    # The rough-input bound is not valid without the small-factor test.
    number, threshold = 210, 4
    q = ceil_div(number.bit_length(), threshold)
    omega, tau = omega_and_tau(number)
    assert tau > 1 << q
    assert least_factor(number) <= 1 << threshold

    return {
        "scope": "finite support, parameter, and recovery checks; arithmetic evaluator is separate",
        "rough_support_cases": rough_cases,
        "small_factor_cases": trivial_cases,
        "sqrt_overhead_parameter_cases": parameter_cases,
        "half_rate_polynomial_ring_parameter_cases": half_cases,
        "sharper_prime_supply_inequality_cases": prime_supply_cases,
        "sharper_prime_supply_actual_cases": len(actual_primes),
        "full_recovery_cases": recovery_cases,
        "prime_endpoint_cases_using_cyclic_queries": prime_endpoints,
        "cyclic_queries": cyclic_queries,
        "largest_prime_used": max_prime,
        "cyclic_cancellation_control": {
            "support": [3, 15],
            "vanishing_modulus": 3,
            "nonvanishing_modulus": 5,
        },
        "missing_pretrial_negative_control": {
            "number": number,
            "threshold_bits": threshold,
            "actual_divisor_count": tau,
            "invalid_unconditional_bound": 1 << q,
        },
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
