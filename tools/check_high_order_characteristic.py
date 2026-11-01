#!/usr/bin/env python3
"""Exact finite checks for high-order small-characteristic selection.

Selection uses primality trial division and modular doubling. It never
factors the forbidden-prime certificate or a number 2^j-1.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from math import comb, gcd, isqrt, lcm
from pathlib import Path
from random import Random

from check_unitary_orbit_recovery import factor_reference, perfect_power, recover_nonzero


@lru_cache(maxsize=None)
def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def selection_bound(n: int, forbidden: int) -> int:
    assert n >= 1 and forbidden != 0
    h = 1 + n * (n + 1) // 2 + abs(forbidden).bit_length()
    return max(256, 4 * (h + 1))


@lru_cache(maxsize=None)
def select_prime(n: int, forbidden: int = 1) -> tuple[int, tuple[int, ...], int]:
    bound = selection_bound(n, forbidden)
    for prime in range(3, bound + 1, 2):
        if not is_prime(prime) or forbidden % prime == 0:
            continue
        powers = []
        current = 1
        for _ in range(n):
            current = (2 * current) % prime
            powers.append(current)
            if current == 1:
                break
        if len(powers) == n and 1 not in powers:
            assert prime > n + 1
            return prime, tuple(powers), bound
    raise AssertionError("The explicit prime-selection bound failed")


def exact_order(prime: int) -> int:
    assert is_prime(prime) and prime != 2
    current = 1
    for exponent in range(1, prime):
        current = (2 * current) % prime
        if current == 1:
            return exponent
    raise AssertionError("Fermat's theorem failed")


def check_primorial() -> dict[str, int]:
    current_lcm = 1
    primorial = 1
    binomial_cases = 0
    primorial_cases = 0
    for m in range(1, 1025):
        current_lcm = lcm(current_lcm, m)
        if is_prime(m):
            primorial *= m
        assert current_lcm >= 1 << (m - 1)
        assert current_lcm // primorial <= m ** isqrt(m)
        js = range(m) if m <= 96 else (0, 1, (m - 1) // 2, m - 2, m - 1)
        for j in js:
            assert current_lcm % (m * comb(m - 1, j)) == 0
            binomial_cases += 1
        if m >= 256:
            assert primorial * primorial >= 1 << (m - 2)
            primorial_cases += 1
    return {"lcm_sizes": 1024, "binomial_divisibility_cases": binomial_cases,
            "primorial_lower_bound_cases": primorial_cases}


def check_selection() -> dict[str, object]:
    rng = Random(20261030)
    cases = 0
    modular_powers = 0
    maximum_prime = 0
    maximum_bound = 0
    records = []
    for n in range(1, 257):
        forbidden_values = [1, -1, rng.getrandbits(128) | 1]
        if n <= 96:
            forbidden = 1
            for p in range(2, 2 * n + 1):
                if is_prime(p):
                    forbidden *= p
            forbidden_values.append(forbidden)
        for forbidden in forbidden_values:
            prime, powers, bound = select_prime(n, forbidden)
            assert prime <= bound and forbidden % prime != 0
            assert powers == tuple(pow(2, j, prime) for j in range(1, n + 1))
            assert exact_order(prime) > n
            cases += 1
            modular_powers += n
            maximum_prime = max(maximum_prime, prime)
            maximum_bound = max(maximum_bound, bound)
        if n in (1, 2, 8, 32, 64, 128, 256):
            prime, _, bound = select_prime(n)
            records.append({"n": n, "prime": prime, "order": exact_order(prime), "bound": bound})
    # Literal product certificate for small sizes, separate from selection.
    proof_products = 0
    for n in range(1, 33):
        forbidden = 2 * 3 * 5 * 7 * 11
        prime, _, bound = select_prime(n, forbidden)
        bad_product = 2 * abs(forbidden)
        for j in range(1, n + 1):
            bad_product *= (1 << j) - 1
        assert bad_product % prime != 0
        primorial = 1
        for p in range(2, bound + 1):
            if is_prime(p):
                primorial *= p
        assert primorial > bad_product
        proof_products += 1
    return {"selection_cases": cases, "certified_modular_powers": modular_powers,
            "maximum_selected_prime": maximum_prime, "maximum_search_bound": maximum_bound,
            "literal_product_cases": proof_products, "sample_selections": records}


@lru_cache(maxsize=None)
def canonical_roots(value: int) -> tuple[tuple[int], ...]:
    assert value > 1 and value % 2 == 1
    return tuple((z,) for z in range(2, value // 2 + 1) if (z * z - 1) % value == 0)


def factor_selected(value: int, prime: int) -> tuple[list[int], int]:
    output = []
    while value % 2 == 0:
        value //= 2
        output.append(2)
    pending = [] if value == 1 else [(value, 1)]
    queries = 0
    while pending:
        current, multiplicity = pending.pop()
        root = perfect_power(current)
        if root is not None:
            base, exponent = root
            pending.append((base, exponent * multiplicity))
            continue
        witnesses = canonical_roots(current)
        queries += 1
        if len(witnesses) % prime == 0:
            output.extend([current] * multiplicity)
            continue
        (z,), used = recover_nonzero(witnesses, current.bit_length(), 1, prime)
        queries += used
        divisor = gcd(current, z - 1)
        assert 1 < divisor < current and current % divisor == 0
        pending.extend(((divisor, multiplicity), (current // divisor, multiplicity)))
    return sorted(output), queries


def check_factoring(limit: int) -> dict[str, int]:
    maximum_queries = 0
    for value in range(2, limit + 1):
        prime, _, _ = select_prime(value.bit_length())
        result, queries = factor_selected(value, prime)
        assert result == factor_reference(value)
        maximum_queries = max(maximum_queries, queries)
    # The root count can exceed the selected prime; order, not cardinality,
    # is what prevents cancellation.
    large_value = 3 * 5 * 7 * 11 * 13 * 17
    prime, _, _ = select_prime(large_value.bit_length())
    witnesses = canonical_roots(large_value)
    assert len(witnesses) == 31 and len(witnesses) > prime
    assert len(witnesses) % prime != 0
    result, queries = factor_selected(large_value, prime)
    assert result == factor_reference(large_value)
    return {"consecutive_inputs": limit - 1, "maximum_factorization_queries": maximum_queries,
            "large_example_N": large_value, "large_example_prime": prime,
            "large_example_count": len(witnesses), "large_example_queries": queries}


def negative_controls() -> dict[str, object]:
    assert len(canonical_roots(105)) == 3
    assert len(canonical_roots(105)) % 3 == 0
    value = 255255
    assert value.bit_length() == 18 and 31 > value.bit_length()
    assert exact_order(31) == 5
    assert len(canonical_roots(value)) == 31
    assert len(canonical_roots(value)) % 31 == 0
    prime, _, _ = select_prime(value.bit_length(), 31)
    assert prime != 31 and len(canonical_roots(value)) % prime != 0
    return {
        "arbitrary_characteristic_cancels": {"N": 105, "prime": 3, "root_count": 3},
        "prime_larger_than_bitlength_is_insufficient": {"N": value, "n": 18, "prime": 31, "order": 5, "root_count": 31},
        "forbidden_prime_avoided": {"forbidden": 31, "selected": prime},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=4096)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"status": "PASS", "scope": "Exact finite checks, not a new factoring exponent",
              "primorial": check_primorial(), "selection": check_selection(),
              "factoring": check_factoring(args.limit), "negative_controls": negative_controls(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
