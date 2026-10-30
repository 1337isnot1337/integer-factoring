#!/usr/bin/env python3
"""Check an elementary divisor-count bound with leading coefficient one."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path


def divisor_count_bits(n: int) -> int:
    """Return q such that tau(N) <= 2**q for every positive N < 2**n."""
    if n < 1:
        raise ValueError("n must be positive")
    floor_log = n.bit_length() - 1
    ceil_log = (n - 1).bit_length()
    bound = min(n - 1, (n + 1) // 2 + 1)
    for split in range(1, floor_log + 1):
        candidate = (1 << (split - 1)) * ceil_log + (n - 1) // split
        bound = min(bound, candidate)
    return bound


def old_divisor_count_bits(n: int) -> int:
    if n < 4:
        return n
    floor_log = n.bit_length() - 1
    return min(n, -(-(9 * n) // floor_log))


def divisor_sieve(limit: int) -> list[int]:
    """Independent direct definition: count all positive divisors."""
    counts = [0] * (limit + 1)
    for divisor in range(1, limit + 1):
        for multiple in range(divisor, limit + 1, divisor):
            counts[multiple] += 1
    return counts


def least_bucket_degree(n: int, q: int) -> int:
    target = n * ((1 << q) - 1)
    degree = 1
    order = 3
    while order <= target:
        degree += 1
        order *= 3
    assert order > target
    if degree > 1:
        assert order // 3 <= target
    return degree


def check_divisor_counts(limit: int) -> dict[str, int]:
    counts = divisor_sieve(limit)
    max_by_bits: dict[int, int] = {}
    exact_residue_checks = 0
    for number in range(1, limit + 1):
        n = number.bit_length()
        q = divisor_count_bits(n)
        assert counts[number] <= 1 << q
        # The strict +1 modulus margin is essential for exact counts.
        modulus = 1 << (q + 1)
        assert counts[number] % modulus == counts[number]
        exact_residue_checks += 1
        max_by_bits[n] = max(max_by_bits.get(n, 0), counts[number])

    split_checks = 0
    for n, max_count in max_by_bits.items():
        floor_log = n.bit_length() - 1
        ceil_log = (n - 1).bit_length()
        for split in range(1, floor_log + 1):
            candidate = (1 << (split - 1)) * ceil_log + (n - 1) // split
            assert max_count <= 1 << candidate
            split_checks += 1
    return {
        "independently_sieved_integers": limit,
        "all_split_candidates_on_sieved_bit_ranges": split_checks,
        "exact_count_modulus_checks": exact_residue_checks,
        "largest_sieved_divisor_count": max(counts),
    }


def check_parameters() -> dict[str, int]:
    invalid_inputs = 0
    for value in (0, -1, -100):
        try:
            divisor_count_bits(value)
        except ValueError:
            invalid_inputs += 1
        else:
            raise AssertionError("nonpositive input was accepted")
    assert divisor_count_bits(1) == 0
    assert divisor_count_bits(2) == 1
    assert divisor_count_bits(3) == 2

    old_comparisons = 0
    for n in range(1, 65537):
        q = divisor_count_bits(n)
        assert 0 <= q <= n - 1
        assert q <= old_divisor_count_bits(n)
        old_comparisons += 1

    coarse_inequalities = 0
    for floor_log in range(2, 8193):
        # Squared form avoids real square roots for odd floor_log.
        assert (floor_log * (floor_log + 1)) ** 2 <= 36 * (1 << floor_log)
        coarse_inequalities += 1

    asymptotic_inequalities = 0
    bit_boundaries = 0
    for exponent in range(1, 513):
        for n in ((1 << exponent) - 1, 1 << exponent, (1 << exponent) + 1):
            if n < 1:
                continue
            floor_log = n.bit_length() - 1
            ceil_log = (n - 1).bit_length()
            assert (1 << floor_log) <= n < (1 << (floor_log + 1))
            assert n <= (1 << ceil_log)
            if ceil_log:
                assert (1 << (ceil_log - 1)) < n
            q = divisor_count_bits(n)
            assert q <= old_divisor_count_bits(n)
            bit_boundaries += 1
            if floor_log >= 16:
                r = floor_log.bit_length()
                split = floor_log - 3 * r
                assert 1 <= split <= floor_log
                assert 1 << r >= floor_log + 1
                candidate = (1 << (split - 1)) * ceil_log + (n - 1) // split
                explicit_upper = Fraction(n, 2 * (floor_log + 1) ** 2) + Fraction(
                    n, floor_log - 3 * r
                )
                assert q <= candidate <= explicit_upper
                assert q >= (n - 1) // floor_log
                asymptotic_inequalities += 1
    return {
        "invalid_inputs_rejected": invalid_inputs,
        "comparisons_against_old_bound": old_comparisons,
        "exact_coarse_comparison_inequalities": coarse_inequalities,
        "bit_boundary_cases": bit_boundaries,
        "large_input_leading_coefficient_inequalities": asymptotic_inequalities,
    }


def check_recovery_parameters() -> dict[str, int]:
    cases = 0
    lengths = list(range(1, 2049)) + [1 << exponent for exponent in range(12, 17)]
    for n in lengths:
        q = divisor_count_bits(n)
        degree = least_bucket_degree(n, q)
        # q+ceil(log2 n)+1 is already sufficient even in base two,
        # hence also for base three.  This integer check includes q=0.
        assert degree <= q + (n - 1).bit_length() + 1
        old_degree = least_bucket_degree(n, old_divisor_count_bits(n))
        assert degree <= old_degree
        cases += 1

    rational_modulus_checks = 0
    for q in range(65):
        for denominator_product in (1, 2, 6, 30, 210, 1001):
            modulus = (1 << (q + 1)) * denominator_product + 1
            assert modulus > 1 << q
            assert gcd(modulus, denominator_product) == 1
            rational_modulus_checks += 1
    return {
        "strict_ternary_bucket_parameter_checks": cases,
        "rational_exact_count_modulus_checks": rational_modulus_checks,
    }


def comparison_table() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in (16, 32, 64, 128, 256, 512, 1024, 4096, 1_000_000, 1 << 40, 1 << 100):
        q = divisor_count_bits(n)
        previous = old_divisor_count_bits(n)
        rows.append(
            {
                "input_bits_n": n,
                "previous_count_bits": previous,
                "new_count_bits": q,
                "new_over_previous": str(Fraction(q, previous)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=200_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.limit < 2:
        parser.error("--limit must be at least two")
    results: dict[str, object] = {
        "status": "passed",
        "scope": "Uniform auxiliary divisor bound; no improvement of the unconditional factoring exponent.",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "Exact integer and rational arithmetic",
    }
    results.update(check_divisor_counts(args.limit))
    results.update(check_parameters())
    results.update(check_recovery_parameters())
    results["comparison_table"] = comparison_table()
    encoded = json.dumps(results, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
