#!/usr/bin/env python3
"""Exact finite checks for sparse divisor nonzeroness in cyclic rings.

The reference divisor enumeration here is exhaustive, not a fast evaluator.
The all-input conditional transfer is proved in the companion research note.
"""

from itertools import combinations
import json
from math import comb, prod


def divisor_bound_bits(n: int) -> int:
    if n < 2:
        raise ValueError("use n >= 2")
    floor_log = n.bit_length() - 1
    ceil_log = (n - 1).bit_length()
    return min([n] + [(1 << cut) * ceil_log + (n + cut - 1) // cut
                      for cut in range(1, floor_log + 1)])


def first_odd_primes(count: int) -> list[int]:
    all_primes = [2]
    candidate = 3
    while len(all_primes) <= count:
        prime = True
        for divisor in all_primes:
            if divisor * divisor > candidate:
                break
            if candidate % divisor == 0:
                prime = False
                break
        if prime:
            all_primes.append(candidate)
        candidate += 2
    return all_primes[1:]


def cyclic_coefficients(exponents: tuple[int, ...], modulus: int,
                        characteristic: int) -> dict[int, int]:
    values: dict[int, int] = {}
    for exponent in exponents:
        residue = exponent % modulus
        updated = (values.get(residue, 0) + 1) % characteristic
        if updated:
            values[residue] = updated
        else:
            values.pop(residue, None)
    return values


def cyclic_product_two(first: int, second: int, dimension: int) -> int:
    """F2 coefficient vectors packed into integers; exact cyclic convolution."""
    mask = (1 << dimension) - 1
    if first & ~mask or second & ~mask:
        raise ValueError("coefficient vector exceeds ring dimension")
    value = 0
    for shift in range(dimension):
        if (second >> shift) & 1:
            rotated = ((first << shift) | (first >> (dimension - shift))) & mask
            value ^= rotated
    return value


def divisors(value: int) -> tuple[int, ...]:
    return tuple(candidate for candidate in range(1, value + 1)
                 if value % candidate == 0)


def verify() -> dict:
    maximum_n = 8
    maximum_k = 1 << divisor_bound_bits(maximum_n)
    maximum_m = maximum_n * (maximum_k - 1) + 1
    primes = first_odd_primes(maximum_m)

    supply_cases = 0
    for count, prime in enumerate(primes, 1):
        assert prime <= 4 * (count + 1) ** 2
        supply_cases += 1

    binomial_cases = 0
    for value in range(2, 130, 2):
        central = comb(value, value // 2)
        small_primes = [2] + [prime for prime in primes if prime <= value]
        assert central * (value + 1) >= 1 << value
        assert central <= value ** len(small_primes)
        remainder = central
        for prime in small_primes:
            exponent = 0
            while remainder % prime == 0:
                exponent += 1
                remainder //= prime
            assert prime ** exponent <= value
        assert remainder == 1
        binomial_cases += 1

    tau = [0] * 4096
    for divisor in range(1, len(tau)):
        for multiple in range(divisor, len(tau), divisor):
            tau[multiple] += 1
    divisor_cases = 0
    for value in range(2, len(tau)):
        n = value.bit_length()
        assert tau[value] <= 1 << divisor_bound_bits(n)
        divisor_cases += 1

    factored_cases = 0
    for width in range(1, 65):
        basis = ([2] + primes)[:width]
        for budget in range(1, 17):
            exponents = [max(1, budget // (index + 1)) for index in range(width)]
            value = prod(prime ** exponent for prime, exponent in zip(basis, exponents))
            divisor_count = prod(exponent + 1 for exponent in exponents)
            assert divisor_count <= 1 << divisor_bound_bits(value.bit_length())
            factored_cases += 1

    sparse_cases = bad_reductions = 0
    for n in range(2, 5):
        for length in range(1, 5):
            query_count = n * (length - 1) + 1
            family = primes[:query_count]
            for support in combinations(range(1 << n), length):
                anchor_product = prod(abs(value - support[0]) for value in support[1:])
                assert anchor_product >= 1
                if length > 1:
                    assert anchor_product < 1 << (n * (length - 1))
                for characteristic in (2, 3, 5):
                    bad = [prime for prime in family
                           if not cyclic_coefficients(support, prime, characteristic)]
                    assert len(bad) < len(family)
                    assert all(anchor_product % prime == 0 for prime in bad)
                    bad_reductions += len(bad)
                    sparse_cases += 1

    substitution_cases = multiplication_cases = 0
    for dimension in (3, 5, 7, 11, 13):
        for first in range(dimension):
            for second in range(dimension):
                assert cyclic_product_two(1 << first, 1 << second, dimension) == (
                    1 << ((first + second) % dimension))
                multiplication_cases += 1
        for n in range(1, 9):
            weights = [1 << pow(2, bit, dimension) for bit in range(n)]
            for value in range(1 << n):
                weighted = 1
                for bit, weight in enumerate(weights):
                    if (value >> bit) & 1:
                        weighted = cyclic_product_two(weighted, weight, dimension)
                assert weighted == 1 << (value % dimension)
                substitution_cases += 1

    prefix_queries = ring_queries = prefix_recovery_cases = prime_endpoints = 0

    def prefix_exists(value: int, n: int, free_bits: int, prefix: int,
                      characteristic: int, remove_one: bool = True) -> bool:
        nonlocal prefix_queries, ring_queries
        prefix_queries += 1
        support = tuple(divisor for divisor in divisors(value)
                        if divisor >> free_bits == prefix)
        proper = tuple(divisor for divisor in support if divisor > 1)
        count = n * ((1 << divisor_bound_bits(n)) - 1) + 1
        for prime in primes[:count]:
            ring_queries += 1
            coefficients = cyclic_coefficients(support, prime, characteristic)
            if remove_one and (1 >> free_bits) == prefix:
                residue = 1 % prime
                updated = (coefficients.get(residue, 0) - 1) % characteristic
                if updated:
                    coefficients[residue] = updated
                else:
                    coefficients.pop(residue, None)
            if remove_one:
                assert coefficients == cyclic_coefficients(proper, prime, characteristic)
            if coefficients:
                return True
        return False

    def recover(value: int, n: int, characteristic: int,
                remove_one: bool = True) -> int:
        prefix = 0
        for free_bits in range(n - 1, -1, -1):
            zero_child = 2 * prefix
            if prefix_exists(value, n, free_bits, zero_child, characteristic, remove_one):
                prefix = zero_child
            else:
                prefix = zero_child + 1
            eligible = tuple(divisor for divisor in divisors(value)
                             if divisor >> free_bits == prefix
                             and (divisor > 1 or not remove_one))
            assert eligible
        return prefix

    for characteristic, limit in ((2, 256), (3, 64), (5, 64)):
        for value in range(2, limit):
            result = recover(value, value.bit_length(), characteristic)
            assert result == min(divisor for divisor in divisors(value) if divisor > 1)
            assert all(result % candidate for candidate in range(2, result))
            prime_endpoints += result == value
            prefix_recovery_cases += 1

    full_factorization_cases = 0
    for value in range(2, 128):
        remainder, factors = value, []
        while remainder > 1:
            factor = recover(remainder, value.bit_length(), 2)
            assert remainder % factor == 0
            factors.append(factor)
            remainder //= factor
        assert prod(factors) == value
        assert factors == sorted(factors)
        full_factorization_cases += 1

    # One fixed ring and ordinary parity both miss the divisors {3,9} of 9.
    assert not cyclic_coefficients((3, 9), 3, 2)
    assert cyclic_coefficients((3, 9), 5, 2)
    assert len((3, 9)) % 2 == 0
    # Retaining divisor one makes greedy recovery return one instead of minFac(15).
    assert recover(15, 4, 2, remove_one=False) == 1
    # Two auxiliary witnesses per divisor would erase every F2 coefficient.
    assert not cyclic_coefficients((3, 3, 9, 9), 5, 2)

    return {
        "scope": "exact sparse-ring identities and finite prefix recovery with an exhaustive divisor oracle; no fast evaluator implementation or improved crossing family",
        "finite_divisor_bound_cases": divisor_cases,
        "constructed_prime_factorization_cases": factored_cases,
        "central_binomial_prime_supply_cases": binomial_cases,
        "prime_supply_cutoff_cases": supply_cases,
        "sparse_support_characteristic_cases": sparse_cases,
        "vanishing_reductions_checked_against_anchor_product": bad_reductions,
        "cyclic_monomial_multiplication_cases": multiplication_cases,
        "digit_weight_substitution_cases": substitution_cases,
        "smallest_prime_factor_recovery_cases": prefix_recovery_cases,
        "prime_input_endpoint_cases": prime_endpoints,
        "complete_prime_factorization_cases": full_factorization_cases,
        "prefix_queries": prefix_queries,
        "cyclic_ring_queries": ring_queries,
        "explicit_bound_examples": [
            {"input_bit_bound": n, "log2_K": divisor_bound_bits(n)}
            for n in (64, 256, 1024, 4096, 16384, 65536)
        ],
        "negative_controls": [
            "single modulus three misses X^3+X^9 over F2",
            "unweighted parity misses two permitted divisors",
            "omitting subtraction of X returns divisor one",
            "duplicate witnesses cancel all F2 coefficients",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
