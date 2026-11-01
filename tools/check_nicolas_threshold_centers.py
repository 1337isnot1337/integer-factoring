#!/usr/bin/env python3
"""Exact rational-interval controls for the polynomial descriptor-cover proof.

The finite highly-composite checks do not certify numerical global Nicolas
constants. The unconditional complexity theorem uses fixed existential
constants, as stated separately in the research note.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from nicolas_threshold_centers import (
    threshold_center, prime_count, is_prime, integer_log, reduced_log,
)


@lru_cache(maxsize=None)
def reduced_log_interval(numerator: int, denominator: int, terms: int = 32):
    """Independent, exact rational interval; no fixed-point routine is used."""
    assert denominator <= numerator <= 2 * denominator
    z = Q(numerator - denominator, numerator + denominator)
    square = z * z
    power = z
    total = Q(0)
    for j in range(terms):
        total += 2 * power / (2 * j + 1)
        power *= square
    # The positive omitted denominators are all at least 2*terms+1.
    tail = 2 * power / ((2 * terms + 1) * (1 - square))
    return total, total + tail


@lru_cache(maxsize=None)
def integer_log_interval(value: int):
    shift = value.bit_length() - 1
    lo, hi = reduced_log_interval(value, 1 << shift)
    l2, h2 = reduced_log_interval(2, 1)
    return lo + shift * l2, hi + shift * h2


def difference_interval(p: int, a: int, candidate: int):
    lp, hp = integer_log_interval(p)
    la, ha = reduced_log_interval(a + 1, a)
    lj, hj = integer_log_interval(candidate)
    l2, h2 = reduced_log_interval(2, 1)
    return lj * l2 - hp * ha, hj * h2 - lp * la


def audit_centers():
    cases = []
    for p in (2, 3, 5, 7, 11, 17, 31, 61, 127, 257, 1021, 4093, 65521):
        for a in (1, 2, 3, 4, 5, 7, 11, 16, 23, 31, 47, 64):
            center = threshold_center(p, a)
            # Monotonicity of the logarithmic test gives |center-theta|<2.
            if center - 2 >= 1:
                assert difference_interval(p, a, center - 2)[1] < 0
            if center + 2 <= p:
                assert difference_interval(p, a, center + 2)[0] > 0
            cases.append((p, a, center))
    return cases


def audit_log_errors():
    count = 0
    for value in range(1, 129):
        for precision in (16, 32, 64):
            # Exact interval precision is much finer than the tested fixed point.
            lower, upper = integer_log_interval(value)
            estimate = Q(integer_log(value, precision), 1 << precision)
            assert estimate <= lower
            # Bound from the proof: per reduced log <= (8F+4)*2^-F.
            limit = Q((value.bit_length() + 1) * (8 * precision + 4),
                      1 << precision)
            assert upper - estimate <= limit
            count += 1
    return count


def trial_exponents(value):
    result = []
    p = 2
    while value > 1:
        if is_prime(p):
            exponent = 0
            while value % p == 0:
                value //= p
                exponent += 1
            result.append((p, exponent))
        p += 1
    return result


def audit_finite_records(bits=16):
    limit = 1 << bits
    divisors = [0] * limit
    for divisor in range(1, limit):
        for multiple in range(divisor, limit, divisor):
            divisors[multiple] += 1
    best = 0
    records = []
    level_checks = 0
    maximum_offset = 0
    for value in range(1, limit):
        if divisors[value] <= best:
            continue
        best = divisors[value]
        factors = trial_exponents(value)
        if not factors:
            records.append({"integer": value, "divisors": best, "offsets": []})
            continue
        assert all(e > 0 for _, e in factors)
        assert [e for _, e in factors] == sorted((e for _, e in factors), reverse=True)
        p = factors[-1][0]
        offsets = []
        for level in range(3, factors[0][1] + 2):
            height = sum(e >= level for _, e in factors)
            center = threshold_center(p, level)
            offset = height - prime_count(center)
            offsets.append(offset)
            maximum_offset = max(maximum_offset, abs(offset))
            assert abs(offset) <= 3
            level_checks += 1
        records.append({"integer": value, "divisors": best,
                        "largest_prime": p, "offsets_above_level_2": offsets})
    return records, level_checks, maximum_offset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    centers = audit_centers()
    logs = audit_log_errors()
    records, levels, offset = audit_finite_records()
    result = {
        "scope": "Exact arithmetic controls and finite record coverage; no numerical global Nicolas cutoff is certified.",
        "arithmetic": {"center_cases": len(centers), "logarithm_error_checks": logs,
                       "center_tuples": centers},
        "finite_records": {"integers_sieved": (1 << 16) - 1,
                           "strict_records": len(records), "high_level_checks": levels,
                           "illustrative_cutoff": 2, "maximum_absolute_offset": offset,
                           "records": records},
        "source_sha256": hashlib.sha256(Path(__file__).with_name(
            "nicolas_threshold_centers.py").read_bytes()).hexdigest(),
    }
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
