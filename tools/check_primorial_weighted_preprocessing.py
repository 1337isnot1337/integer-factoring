#!/usr/bin/env python3
"""Independent checks for the primorial-weighted preprocessing theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import itertools
import json
from math import isqrt, prod
from pathlib import Path
from time import monotonic

from primorial_weighted_preprocessing import (
    descriptor_exponents, descriptor_product_bit, descriptor_values,
    fast_weighted_maximum, lower_weight, streamed_exact_descriptor,
    trial_prime, weighted_descriptors,
)


def independent_primes(count):
    bound = max(4, 4 * count * count)
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(bound) + 1):
        if sieve[p]:
            sieve[p*p:bound+1:p] = b"\x00" * ((bound-p*p)//p+1)
    return [p for p in range(bound + 1) if sieve[p]][:count]


def ordinary_exponent_oracle(n):
    """Separate prime-exponent prefix search; no differences or weights."""
    primes = independent_primes(n - 1)
    best_value = best_count = 1
    visited = 0

    def visit(index, previous, value, count):
        nonlocal best_value, best_count, visited
        visited += 1
        if count > best_count or (count == best_count and value < best_value):
            best_count, best_value = count, value
        if index == len(primes):
            return
        for exponent in range(1, previous + 1):
            value *= primes[index]
            if value >= 1 << n:
                break
            visit(index + 1, exponent, value, count * (exponent + 1))

    visit(0, n - 1, 1, 1)
    return best_value, best_count, visited


def weighted_dp_count(budget):
    coefficients = [1] + [0] * budget
    for j in range(1, budget + 1):
        weight = max(j, j * ((j + 1).bit_length() - 1) // 4)
        for degree in range(weight, budget + 1):
            coefficients[degree] += coefficients[degree - weight]
    return sum(coefficients)


def divisor_sieve_maxima(bits):
    limit = 1 << bits
    counts = [0] * limit
    for divisor in range(1, limit):
        for multiple in range(divisor, limit, divisor):
            counts[multiple] += 1
    maximum = witness = 0
    answer = {}
    next_bits = 1
    for value in range(1, limit):
        if counts[value] > maximum:
            maximum, witness = counts[value], value
        if value == (1 << next_bits) - 1:
            answer[next_bits] = (witness, maximum)
            next_bits += 1
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    parser.add_argument("--sieve-bits", type=int, default=16)
    parser.add_argument("--streamed-bits", type=int, default=6)
    args = parser.parse_args()
    if not 1 <= args.sieve_bits <= 20:
        parser.error("sieve bits must lie between 1 and 20")
    if not 1 <= args.streamed_bits <= 9:
        parser.error("streamed bits must lie between 1 and 9")
    started = monotonic()

    primes = independent_primes(512)
    primorial = 1
    for j, p in enumerate(primes, 1):
        primorial *= p
        weight = lower_weight(j)
        assert weight >= j
        assert (1 << weight) <= primorial
        # Exact form of 16*w_j >= j*log2(j+1).
        assert (1 << (16 * weight)) >= (j + 1) ** j
        assert p <= 4 * j * j
        if j <= 64:
            assert trial_prime(j) == p

    descriptor_checks = 0
    support_checks = 0
    count_table = []
    for n in range(1, 25):
        seen = set()
        feasible = 0
        max_support = 0
        for descriptor in weighted_descriptors(n):
            assert descriptor not in seen
            seen.add(descriptor)
            assert all(j > k for (j, _), (k, _) in zip(descriptor, descriptor[1:]))
            assert sum(m * lower_weight(j) for j, m in descriptor) <= n - 1
            exponents = descriptor_exponents(descriptor)
            recovered = tuple((j, exponents[j-1] - (exponents[j] if j < len(exponents) else 0))
                              for j in range(len(exponents), 0, -1)
                              if exponents[j-1] > (exponents[j] if j < len(exponents) else 0))
            assert recovered == descriptor
            assert sum(exponents) == sum(j * m for j, m in descriptor) <= n - 1
            value = prod(p ** e for p, e in zip(primes, exponents))
            column_product = 1
            for j, m in descriptor:
                column_product *= prod(primes[:j]) ** m
            assert value == column_product
            feasible += value < 1 << n
            max_support = max(max_support, len(descriptor))
            # Distinct indices must cost at least the first r lower weights.
            assert sum(lower_weight(j) for j, _ in descriptor) >= sum(lower_weight(j) for j in range(1, len(descriptor)+1))
            descriptor_checks += 1
            support_checks += 1
        assert len(seen) == weighted_dp_count(n - 1)
        ordinary = ordinary_exponent_oracle(n)
        assert feasible == ordinary[2]
        if n in (1, 4, 8, 12, 16, 20, 24):
            count_table.append({"input_bits": n, "surrogate_patterns": len(seen),
                                "feasible_patterns": feasible, "maximum_support": max_support})

    sieved = divisor_sieve_maxima(args.sieve_bits)
    maxima = []
    for n in range(1, args.sieve_bits + 1):
        result = fast_weighted_maximum(n)
        values = descriptor_values(result.descriptor)
        ordinary = ordinary_exponent_oracle(n)
        assert values == sieved[n] == ordinary[:2]
        assert result.feasible_patterns_checked == ordinary[2]
        maxima.append({**asdict(result), "least_witness": values[0], "maximum_divisors": values[1]})

    streamed = []
    for n in range(1, args.streamed_bits + 1):
        result = streamed_exact_descriptor(n)
        expected = fast_weighted_maximum(n)
        assert result == expected
        # Consume generated bits here only for verification, outside the space claim.
        value = sum(descriptor_product_bit(result.descriptor, n, bit) << bit for bit in range(n))
        count = sum(descriptor_product_bit(result.descriptor, n, bit, True) << bit for bit in range(n))
        assert (value, count) == descriptor_values(result.descriptor)
        streamed.append({**asdict(result), "least_witness": value, "maximum_divisors": count})

    bit_checks = 0
    for n in range(1, 7):
        for descriptor in weighted_descriptors(n):
            if sum(j * m for j, m in descriptor) > 3:
                continue
            values = descriptor_values(descriptor)
            for kind, value in enumerate(values):
                for bit in range(value.bit_length() + 2):
                    observed = descriptor_product_bit(descriptor, n, bit, bool(kind))
                    assert observed == ((value >> bit) & 1)
                    bit_checks += 1

    lower_families = []
    for n in (32, 64, 128, 512, 1024, 4096):
        k = isqrt((n - 1) // (4 * n.bit_length()))
        primorials = [prod(primes[:j]) for j in range(1, k + 1)]
        maximal = prod(primorials)
        assert maximal < 1 << n
        family = {prod(P for P, bit in zip(primorials, bits) if bit)
                  for bits in itertools.product((0, 1), repeat=k)}
        assert len(family) == 1 << k
        assert all(value < 1 << n for value in family)
        lower_families.append({"input_bits": n, "support_limit": k,
                               "distinct_feasible_patterns": len(family),
                               "largest_product_bits": maximal.bit_length()})

    prefix_counts = []
    for n in (16, 32, 64):
        value, count, visited = ordinary_exponent_oracle(n)
        assert visited <= weighted_dp_count(n - 1)
        prefix_counts.append({"input_bits": n, "canonical_patterns": visited,
                              "surrogate_patterns": weighted_dp_count(n - 1),
                              "least_witness": value, "maximum_divisors": count})

    invalid = 0
    for function in (fast_weighted_maximum, streamed_exact_descriptor):
        for n in (0, -1, False, 1.5):
            try:
                function(n)
            except ValueError:
                invalid += 1
            else:
                raise AssertionError("invalid input accepted")

    files = (Path(__file__), Path(__file__).with_name("primorial_weighted_preprocessing.py"))
    result = {
        "status": "pass",
        "scope": "Exact finite enumeration, arithmetic, and lower-weight checks; no factoring-score improvement.",
        "primorial_weight_cases": 512,
        "trial_prime_comparisons": 64,
        "weighted_descriptors_checked": descriptor_checks,
        "weighted_count_ranges": 24,
        "support_budget_checks": support_checks,
        "divisor_sieve_integers": (1 << args.sieve_bits) - 1,
        "divisor_sieve_ranges": args.sieve_bits,
        "streamed_complete_ranges": args.streamed_bits,
        "descriptor_product_bits_checked": bit_checks,
        "invalid_inputs_rejected": invalid,
        "pattern_count_table": count_table,
        "ordinary_prefix_count_table": prefix_counts,
        "maxima": maxima,
        "streamed_maxima": streamed,
        "lower_bound_families": lower_families,
        "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
        "seconds": round(monotonic() - started, 3),
    }
    output = json.dumps(result, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(output)
    print(json.dumps({k: v for k, v in result.items() if k not in ("maxima", "streamed_maxima")}, indent=2))


if __name__ == "__main__":
    main()
