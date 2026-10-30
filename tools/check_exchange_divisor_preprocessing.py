#!/usr/bin/env python3
"""Exact audit of next-prime exchange pruning and divisor maxima."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from math import prod
from pathlib import Path
from time import monotonic

from check_primorial_weighted_preprocessing import (
    divisor_sieve_maxima, independent_primes, ordinary_exponent_oracle,
)
from exchange_divisor_preprocessing import (
    exchange_descriptors, fast_exchange_maximum, maximum_exponent,
    streamed_exchange_descriptor,
)
from primorial_weighted_preprocessing import (
    descriptor_exponents, descriptor_values, weighted_descriptors,
)


def floor_root(value, degree):
    lo, hi = 0, value + 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** degree <= value:
            lo = mid
        else:
            hi = mid
    return lo


def factor_by_trial(value, primes):
    exponents = []
    for p in primes:
        if value == 1:
            break
        exponent = 0
        while value % p == 0:
            exponent += 1
            value //= p
        exponents.append(exponent)
    assert value == 1
    return exponents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    parser.add_argument("--sieve-bits", type=int, default=16)
    parser.add_argument("--streamed-bits", type=int, default=6)
    parser.add_argument("--large", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.sieve_bits <= 20:
        parser.error("sieve bits must lie between 1 and 20")
    if not 1 <= args.streamed_bits <= 9:
        parser.error("streamed bits must lie between 1 and 9")
    started = monotonic()
    primes = independent_primes(160)
    sieve_maxima = divisor_sieve_maxima(args.sieve_bits)
    table = []
    for n in range(1, args.sieve_bits + 1):
        restricted = fast_exchange_maximum(n)
        values = descriptor_values(restricted.descriptor)
        unrestricted = ordinary_exponent_oracle(n)
        assert values == sieve_maxima[n] == unrestricted[:2]
        table.append({"input_bits": n, "least_witness": values[0],
                      "maximum_divisors": values[1],
                      "restricted_patterns": restricted.feasible_candidates,
                      "unrestricted_patterns": unrestricted[2]})

    # Independent divisor sieve, retaining every strict record, not just powers of two.
    limit = 1 << args.sieve_bits
    counts = [0] * limit
    for divisor in range(1, limit):
        for multiple in range(divisor, limit, divisor):
            counts[multiple] += 1
    maximum = record_count = exponent_inequalities = 0
    records = []
    for value in range(1, limit):
        if counts[value] <= maximum:
            continue
        maximum = counts[value]
        exponents = factor_by_trial(value, primes)
        assert exponents == sorted(exponents, reverse=True)
        assert all(exponent > 0 for exponent in exponents)
        q = primes[len(exponents)]
        for p, exponent in zip(primes, exponents):
            assert p ** ((exponent + 1) // 2) < q
            exponent_inequalities += 1
        records.append({"integer": value, "divisors": maximum, "next_prime": q})
        record_count += 1

    replacement_checks = 0
    for n in range(2, 19):
        for descriptor in weighted_descriptors(n):
            exponents = descriptor_exponents(descriptor)
            if not exponents:
                continue
            value = prod(p ** e for p, e in zip(primes, exponents))
            count = prod(e + 1 for e in exponents)
            q = primes[len(exponents)]
            for i, (p, exponent) in enumerate(zip(primes, exponents)):
                r = (exponent + 1) // 2
                if p ** r < q:
                    continue
                replacement = value // (p ** r) * q
                new_exponents = exponents.copy()
                new_exponents[i] -= r
                new_count = 2 * prod(e + 1 for e in new_exponents)
                assert replacement < value
                assert new_count >= count
                assert new_count * (exponent + 1) == count * 2 * (exponent - r + 1)
                replacement_checks += 1

    cap_checks = profile_checks = 0
    for n in range(1, 21):
        seen = set()
        B = 4 * n * n
        L = (B - 1).bit_length()
        E = min(n - 1, 2 * L - 2)
        caps = [min(n - 1, floor_root(B - 1, (a + 1) // 2) - 1)
                for a in range(1, E + 1)]
        feasible_profiles = 0
        for descriptor in exchange_descriptors(n):
            assert descriptor not in seen
            seen.add(descriptor)
            heights = [j for j, multiplicity in descriptor for _ in range(multiplicity)]
            assert heights == sorted(heights, reverse=True)
            assert len(heights) <= E
            assert sum(heights) <= n - 1
            for a, height in enumerate(heights, 1):
                assert height <= caps[a-1]
                cap_checks += 1
            exponents = descriptor_exponents(descriptor)
            if exponents:
                q = primes[len(exponents)]
                assert len(heights) <= maximum_exponent(n, q)
                assert all(p ** ((e + 1) // 2) < q for p, e in zip(primes, exponents))
            value = prod(p ** e for p, e in zip(primes, exponents))
            feasible_profiles += value < 1 << n
            profile_checks += 1
        assert len(seen) <= prod(cap + 1 for cap in caps)
        assert feasible_profiles == fast_exchange_maximum(n).feasible_candidates

    streamed_table = []
    for n in range(1, args.streamed_bits + 1):
        streamed = streamed_exchange_descriptor(n)
        fast = fast_exchange_maximum(n)
        assert streamed.descriptor == fast.descriptor
        assert streamed.feasible_candidates == fast.feasible_candidates
        assert descriptor_values(streamed.descriptor) == sieve_maxima[n]
        streamed_table.append(asdict(streamed))

    comparisons = []
    for n in ([16, 32, 64, 128] if args.large else [16, 32, 64]):
        restricted = fast_exchange_maximum(n)
        values = descriptor_values(restricted.descriptor)
        unrestricted = ordinary_exponent_oracle(n)
        assert values == unrestricted[:2]
        comparisons.append({"input_bits": n, "least_witness": values[0],
                            "maximum_divisors": values[1],
                            "restricted_patterns": restricted.feasible_candidates,
                            "unrestricted_patterns": unrestricted[2],
                            "largest_visited_exponent": restricted.maximum_height_depth})

    invalid = 0
    for function in (fast_exchange_maximum, streamed_exchange_descriptor):
        for n in (0, -1, False, 1.5):
            try:
                function(n)
            except ValueError:
                invalid += 1
            else:
                raise AssertionError("invalid input accepted")

    files = [Path(__file__), Path(__file__).with_name("exchange_divisor_preprocessing.py"),
             Path(__file__).with_name("primorial_weighted_preprocessing.py")]
    result = {"status": "pass",
              "scope": "Exact next-prime exchange and implementation audit; TC0 arithmetic is a cited theorem, not run here; no factoring-score improvement.",
              "divisor_sieve_integers": limit - 1,
              "power_two_maxima_compared": args.sieve_bits,
              "strict_record_maxima_checked": record_count,
              "record_exponent_inequalities": exponent_inequalities,
              "literal_replacement_witnesses": replacement_checks,
              "height_profiles_checked": profile_checks,
              "height_cap_checks": cap_checks,
              "elementary_streamed_ranges": args.streamed_bits,
              "invalid_inputs_rejected": invalid,
              "comparison_table": comparisons,
              "small_maxima": table,
              "record_maxima": records,
              "streamed_results": streamed_table,
              "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
              "seconds": round(monotonic() - started, 3)}
    output = json.dumps(result, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(output)
    print(json.dumps({k: v for k, v in result.items() if k not in ("small_maxima", "record_maxima", "streamed_results")}, indent=2))


if __name__ == "__main__":
    main()
