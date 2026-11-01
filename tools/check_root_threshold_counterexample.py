#!/usr/bin/env python3
"""Exact counterexample to detecting composites by whole-root-interval parity."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from math import isqrt, prod
from pathlib import Path


PRIMES = (61, 71, 380117, 175531457527)
EXPONENTS = (1, 2, 2, 2)
INPUT = 1368960876177490538878515064474553061781
HISTOGRAM = {2: 10, 3: 6, 4: 2, 5: 2, 7: 2, 10: 2, 21: 2}


def prime_by_trial_division(value: int) -> tuple[bool, int]:
    if value < 2:
        return False, 0
    if value % 2 == 0:
        return value == 2, 1
    if value % 3 == 0:
        return value == 3, 2
    tests = 2
    for divisor in range(5, isqrt(value) + 1, 6):
        tests += 1
        if value % divisor == 0:
            return False, tests
        if divisor + 2 <= isqrt(value):
            tests += 1
            if value % (divisor + 2) == 0:
                return False, tests
    return True, tests


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    primality = []
    for prime in PRIMES:
        valid, checks = prime_by_trial_division(prime)
        assert valid
        primality.append({"prime": prime, "trial_divisions": checks, "tested_through": isqrt(prime)})
    assert prod(prime**exponent for prime, exponent in zip(PRIMES, EXPONENTS)) == INPUT
    divisors = sorted(
        prod(prime**exponent for prime, exponent in zip(PRIMES, exponents))
        for exponents in product(*(range(exponent + 1) for exponent in EXPONENTS))
    )
    assert len(divisors) == len(set(divisors)) == 54
    assert all(INPUT % divisor == 0 for divisor in divisors)
    bins: Counter[int] = Counter()
    records = []
    for divisor in divisors[1:]:
        exponent = 1
        while divisor ** (exponent + 1) <= INPUT:
            exponent += 1
        assert divisor**exponent <= INPUT < divisor ** (exponent + 1)
        if exponent >= 2:
            bins[exponent] += 1
            records.append({"divisor": divisor, "maximum_power": exponent})
    assert dict(sorted(bins.items())) == HISTOGRAM
    assert all(count % 2 == 0 for count in bins.values())
    root_counts = {}
    for exponent in range(2, 23):
        count = sum(divisor**exponent <= INPUT for divisor in divisors[1:])
        assert count == sum(multiplicity for power, multiplicity in bins.items() if power >= exponent)
        assert count % 2 == 0
        root_counts[str(exponent)] = count
    assert PRIMES[0] ** 22 > INPUT
    # The rigorous power-divisor query succeeds on this input: its bases
    # are the divisors of 71*380117*175531457527, with one excluded.
    square_bases = [divisor for divisor in divisors[1:] if INPUT % (divisor * divisor) == 0]
    assert len(square_bases) == 7
    report = {
        "status": "PASS",
        "scope": "A finite counterexample for whole-root-interval parity; arbitrary prefix algorithms are not excluded.",
        "N": INPUT,
        "input_bits": INPUT.bit_length(),
        "factorization": [[prime, exponent] for prime, exponent in zip(PRIMES, EXPONENTS)],
        "primality_certification": primality,
        "divisor_count": len(divisors),
        "maximum_power_histogram": dict(sorted(bins.items())),
        "root_interval_counts": root_counts,
        "divisor_bin_certificate": records,
        "square_divisor_bases_excluding_one": len(square_bases),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
