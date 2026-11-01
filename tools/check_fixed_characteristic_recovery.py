#!/usr/bin/env python3
"""Exact finite checks for recovery over any fixed prime characteristic.

Counting oracles enumerate small integer witnesses; they are not the
large-input clique evaluator. This checker uses only the standard library.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import comb, isqrt
from pathlib import Path
import random


def factor_reference(value: int) -> list[int]:
    result: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result.append(divisor)
            value //= divisor
        divisor += 1
    if value > 1:
        result.append(value)
    return result


@lru_cache(maxsize=None)
def divisors(value: int) -> tuple[int, ...]:
    low: list[int] = []
    high: list[int] = []
    for divisor in range(1, isqrt(value) + 1):
        if value % divisor == 0:
            low.append(divisor)
            if divisor * divisor != value:
                high.append(value // divisor)
    return tuple(low + high[::-1])


@lru_cache(maxsize=None)
def ordered_products(value: int, length: int) -> tuple[tuple[int, ...], ...]:
    if length == 1:
        return ((value,),)
    return tuple(
        (first, *tail)
        for first in divisors(value)
        for tail in ordered_products(value // first, length - 1)
    )


def rotations(values: tuple[int, ...]) -> list[tuple[int, ...]]:
    return [values[shift:] + values[:shift] for shift in range(len(values))]


def power_tuples(value: int, prime: int, exponent_bit: int) -> tuple[tuple[int, ...], ...]:
    power = prime**exponent_bit
    found: list[tuple[int, ...]] = []
    for base in divisors(value):
        if base == 1:
            continue
        raised = base**power
        if raised <= value and value % raised == 0:
            found.extend(ordered_products(base, prime - 1))
    return tuple(found)


def nontrivial_orbit_representatives(value: int, prime: int) -> tuple[tuple[int, ...], ...]:
    trivial = (1,) * (prime - 1) + (value,)
    return tuple(
        values
        for values in ordered_products(value, prime)
        if values != trivial and all(values < other for other in rotations(values)[1:])
    )


def pack(values: tuple[int, ...], width: int) -> int:
    encoded = 0
    for value in values:
        assert 0 <= value < 1 << width
        encoded = (encoded << width) | value
    return encoded


def unpack(encoded: int, width: int, length: int) -> tuple[int, ...]:
    mask = (1 << width) - 1
    return tuple((encoded >> (width * (length - index - 1))) & mask for index in range(length))


def recover_nonzero(
    witnesses: tuple[tuple[int, ...], ...], width: int, length: int, prime: int
) -> tuple[tuple[int, ...], int]:
    encoded = tuple(pack(witness, width) for witness in witnesses)
    assert len(set(encoded)) == len(encoded)
    parent = len(encoded) % prime
    assert parent != 0
    prefix = 0
    free = width * length
    calls = 0
    while free:
        free -= 1
        candidate = prefix * 2
        left = sum((value >> free) == candidate for value in encoded) % prime
        calls += 1
        prefix = candidate if left else candidate + 1
        child = sum((value >> free) == prefix for value in encoded) % prime
        assert child != 0
        parent = child
    result = unpack(prefix, width, length)
    assert result in witnesses
    return result, calls


def split_or_prime(value: int, prime: int) -> tuple[int | None, int, str]:
    width = value.bit_length()
    exponent_bit = 1
    power = prime
    calls = 0
    while power <= width - 1:
        witnesses = power_tuples(value, prime, exponent_bit)
        calls += 1
        if len(witnesses) % prime:
            component_width = (width + power - 1) // power
            found, used = recover_nonzero(witnesses, component_width, prime - 1, prime)
            divisor = next(component for component in found if component > 1)
            assert 1 < divisor < value and value % divisor == 0
            return divisor, calls + used, "power"
        exponent_bit += 1
        power *= prime
    witnesses = nontrivial_orbit_representatives(value, prime)
    calls += 1
    if len(witnesses) % prime:
        found, used = recover_nonzero(witnesses, width, prime, prime)
        divisor = next(component for component in found if component > 1)
        assert sum(component > 1 for component in found) >= 2
        assert 1 < divisor < value and value % divisor == 0
        return divisor, calls + used, "orbit"
    return None, calls, "prime"


def factor_by_counts(value: int, prime: int) -> tuple[list[int], int]:
    pending = [value]
    result: list[int] = []
    calls = 0
    while pending:
        current = pending.pop()
        divisor, used, _ = split_or_prime(current, prime)
        calls += used
        if divisor is None:
            result.append(current)
        else:
            pending.extend((divisor, current // divisor))
    return sorted(result), calls


def compare_bits(left: int, right: int, width: int, initial: int = 0) -> tuple[int, tuple[int, ...]]:
    """State -1=less, 0=equal, +1=greater; unique deterministic extension."""
    state = initial
    trace = [state]
    for position in range(width - 1, -1, -1):
        a = (left >> position) & 1
        b = (right >> position) & 1
        if state == 0:
            state = (a > b) - (a < b)
        trace.append(state)
    return state, tuple(trace)


def check_comparators() -> dict[str, int]:
    cases = 0
    for width in range(1, 8):
        for left, right in product(range(1 << width), repeat=2):
            state, trace = compare_bits(left, right, width)
            assert state == (left > right) - (left < right)
            assert len(trace) == width + 1
            cases += 1
    generator = random.Random(20261030)
    tuples = 0
    equal_coordinate_cases = 0
    for prime in (2, 3, 5, 7):
        for _ in range(300):
            values = tuple(generator.randrange(1, 32) for _ in range(prime))
            word = pack(values, 5)
            accepted = all(
                compare_bits(word, pack(other, 5), 5 * prime)[0] == -1
                for other in rotations(values)[1:]
            )
            assert accepted == all(values < other for other in rotations(values)[1:])
            if len(set(values)) < prime and len(set(values)) > 1:
                equal_coordinate_cases += 1
            tuples += 1
    # A repeated-coordinate representative has a unique strict least
    # rotation despite its equal first two components.
    assert min(rotations((2, 2, 3))) == (2, 2, 3)
    assert all((2, 2, 3) < other for other in rotations((2, 2, 3))[1:])
    return {"integer_comparisons": cases, "tuple_rotation_comparisons": tuples, "nonconstant_repeated_coordinate_tuples": equal_coordinate_cases}


def check_characteristic(prime: int, limit: int) -> dict[str, int]:
    power_cases = 0
    free_orbit_cases = 0
    representatives = 0
    max_queries = 0
    routes: Counter[str] = Counter()
    for value in range(2, limit + 1):
        width = value.bit_length()
        exponents = Counter(factor_reference(value))
        digit = 1
        power = prime
        all_zero = True
        while power <= width - 1:
            witnesses = power_tuples(value, prime, digit)
            count = 1
            for a in exponents.values():
                count *= comb(a // power + prime - 1, prime - 1)
            count -= 1
            assert count == len(witnesses)
            nonzero_digit = any((a // power) % prime for a in exponents.values())
            expected = prime - 1 if nonzero_digit else 0
            assert count % prime == expected
            all_zero &= count % prime == 0
            power_cases += 1
            digit += 1
            power *= prime
        p_free = all(a < prime for a in exponents.values())
        assert all_zero == p_free
        if p_free:
            ordered = ordered_products(value, prime)
            formula = 1
            for a in exponents.values():
                formula *= comb(a + prime - 1, prime - 1)
            assert len(ordered) == formula
            assert all(len(set(rotations(values))) == prime for values in ordered)
            canonical = {min(rotations(values)) for values in ordered}
            assert len(canonical) * prime == len(ordered)
            admitted = nontrivial_orbit_representatives(value, prime)
            assert len(admitted) == len(canonical) - 1 == formula // prime - 1
            if len(exponents) >= 2:
                assert len(admitted) % prime == prime - 1
            else:
                a = next(iter(exponents.values()))
                assert len(admitted) % prime == (pow(a, -1, prime) - 1) % prime
            assert (len(admitted) % prime != 0) == (factor_reference(value) != [value])
            free_orbit_cases += 1
            representatives += len(canonical)
        found, used = factor_by_counts(value, prime)
        assert found == factor_reference(value)
        max_queries = max(max_queries, used)
        _, _, route = split_or_prime(value, prime)
        routes[route] += 1
    return {
        "prime": prime,
        "integers": limit - 1,
        "power_digit_cases": power_cases,
        "p_free_orbit_cases": free_orbit_cases,
        "canonical_representatives": representatives,
        "maximum_factorization_queries": max_queries,
        "power_route_inputs": routes["power"],
        "orbit_route_inputs": routes["orbit"],
        "prime_declarations": routes["prime"],
    }


def check_binomial_digits() -> dict[str, int]:
    cases = 0
    for prime in (2, 3, 5, 7, 11, 13, 17, 19):
        for value in range(2001):
            residue = comb(value + prime - 1, prime - 1) % prime
            assert residue == int(value % prime == 0)
            cases += 1
        for exponent in range(1, prime):
            value = comb(exponent + prime - 1, prime - 1)
            assert value % prime == 0 and value % (prime * prime) != 0
            assert (value // prime) % prime == pow(exponent, -1, prime)
    return {"digit_identity_cases": cases, "primes": 8}


def check_budgets() -> dict[str, int]:
    cases = 0
    for prime in (2, 3, 5, 7, 11, 17):
        for width in list(range(2, 513)) + [10**k for k in range(3, 10)]:
            digit = 1
            power = prime
            while power <= width - 1:
                sizes = [(width + prime ** (digit - i) - 1) // prime ** (digit - i) for i in range(digit + 1)]
                requested = (prime - 2) * sizes[0] + (prime - 1) * sum(sizes[1:]) + width
                assert requested <= 4 * prime * width
                assert (prime - 1) * sizes[0] <= width + prime - 1
                assert (prime - 1) * prime * width <= prime * prime * width
                cases += 1
                digit += 1
                power *= prime
    return {"power_chain_cases": cases, "fixed_characteristics": 6}


def negative_controls() -> dict[str, object]:
    # The C_3 action is not free when the integer is a cube.
    cubes = ordered_products(8, 3)
    assert (2, 2, 2) in cubes and len(set(rotations((2, 2, 2)))) == 1
    assert len(cubes) == 10 and len(cubes) % 3 != 0
    # Composite rotation length has nonconstant tuples with shorter periods.
    periodic = (2, 3, 2, 3)
    assert len(set(rotations(periodic))) == 2 and 2 * 3 * 2 * 3 == 36
    assert 36 not in {value**4 for value in range(1, 4)}
    # Forgetting the initial comparison state falsely makes equal words less.
    assert compare_bits(13, 13, 4, initial=-1)[0] == -1
    assert compare_bits(13, 13, 4, initial=0)[0] == 0
    # A prime has its one trivial orbit; failing to subtract it gives a
    # false composite declaration over every characteristic checked here.
    counts = {}
    for prime in (2, 3, 5, 7):
        canonical = {min(rotations(values)) for values in ordered_products(11, prime)}
        assert len(canonical) == 1
        counts[str(prime)] = len(canonical)
    return {"nonfree_cube": {"N": 8, "p": 3, "ordered_count": 10}, "composite_rotation_length": {"tuple": periodic, "orbit_size": 2}, "bad_initial_comparison": {"left": 13, "right": 13, "initial_state": "less"}, "unremoved_prime_orbit_counts": counts}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=512)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert args.max_n >= 2
    report = {
        "status": "PASS",
        "scope": "Exact finite tuple, orbit, comparator and recovery checks; no asymptotic evaluator is executed.",
        "characteristics": [check_characteristic(prime, args.max_n) for prime in (2, 3, 5, 7)],
        "binomial_identities": check_binomial_digits(),
        "comparators": check_comparators(),
        "budgets": check_budgets(),
        "negative_controls": negative_controls(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
