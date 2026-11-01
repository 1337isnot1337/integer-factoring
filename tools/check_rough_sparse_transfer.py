#!/usr/bin/env python3
"""Exact controls for rough-input sparse divisor-polynomial transfer.

The finite evaluator enumerates divisors and is only a reference oracle.
It is not an implementation of the hypothesized fast arithmetic circuits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from math import isqrt, prod
from pathlib import Path


def ceil_sqrt(n):
    r = isqrt(n)
    return r + (r * r < n)


def trial_prefix(value, exponent):
    limit = 1 << exponent
    for divisor in range(2, min(limit, value) + 1):
        if value % divisor == 0:
            return divisor
    return None


def trial_factors(value):
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            value //= divisor
            factors.append(divisor)
        divisor += 1
    if value > 1:
        factors.append(value)
    return factors


@lru_cache(maxsize=None)
def odd_primes(count):
    result = []
    candidate = 3
    while len(result) < count:
        if all(candidate % divisor for divisor in range(2, isqrt(candidate) + 1)):
            result.append(candidate)
        candidate += 2
    return tuple(result)


@lru_cache(maxsize=None)
def reference_divisors(value):
    return tuple(divisor for divisor in range(1, value + 1) if value % divisor == 0)


def cyclic_nonzero(exponents, modulus, characteristic):
    coefficients = Counter(exponent % modulus for exponent in exponents)
    return any(value % characteristic for value in coefficients.values())


def rough_min_factor(value, exponent, characteristic, statistics):
    n = value.bit_length()
    direct = trial_prefix(value, exponent)
    statistics['trial_calls'] += 1
    if direct is not None:
        statistics['trial_successes'] += 1
        return direct
    factors = trial_factors(value)
    assert all(p > 1 << exponent for p in factors)
    height = (n - 1) // exponent
    support_bound = 1 << height
    assert exponent * len(factors) <= n - 1
    assert len(reference_divisors(value)) <= support_bound
    statistics['rough_calls'] += 1
    statistics['rough_support_checks'] += 1
    family_size = n * (support_bound - 1) + 1
    primes = odd_primes(family_size)
    assert primes[-1] <= 4 * (family_size + 1) ** 2
    eligible = tuple(d for d in reference_divisors(value) if d > 1)
    k, prefix = n, 0
    while k:
        lower_prefix = 2 * prefix
        subset = tuple(d for d in eligible if d >> (k - 1) == lower_prefix)
        detected = False
        for modulus in primes:
            statistics['cyclic_reductions'] += 1
            if cyclic_nonzero(subset, modulus, characteristic):
                detected = True
                break
        assert detected == bool(subset)
        prefix = lower_prefix if detected else lower_prefix + 1
        k -= 1
        assert any(d >> k == prefix for d in eligible)
        statistics['prefix_queries'] += 1
    assert prefix == min(eligible)
    return prefix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    statistics = Counter()
    samples = []
    for characteristic in (2, 3, 5):
        for value in range(2, 2049):
            n = value.bit_length()
            exponent = ceil_sqrt(n)
            result = rough_min_factor(value, exponent, characteristic, statistics)
            expected = trial_factors(value)[0]
            assert result == expected
            statistics['least_factor_cases'] += 1
            if value in (97, 221, 899, 1763, 2021, 2047):
                samples.append({'integer': value, 'characteristic': characteristic,
                                'sieve_exponent': exponent, 'least_factor': result,
                                'support_bound_if_rough': 1 << ((n - 1) // exponent)})
    # Every integer threshold is allowed, including thresholds exceeding sqrt(N).
    for value in range(2, 257):
        for exponent in range(1, value.bit_length() + 1):
            result = rough_min_factor(value, exponent, 2, statistics)
            assert result == trial_factors(value)[0]
            statistics['tradeoff_cases'] += 1
    # Complete prime factorization with the same original threshold.
    for value in range(2, 513):
        exponent = ceil_sqrt(value.bit_length())
        remaining = value
        factors = []
        while remaining > 1:
            prime = rough_min_factor(remaining, exponent, 2, statistics)
            factors.append(prime)
            remaining //= prime
        assert factors == trial_factors(value)
        assert prod(factors) == value
        statistics['complete_factorizations'] += 1
    # The support bound is used only after the small-factor branch failed.
    assert len(reference_divisors(36)) > 1 << (((36).bit_length() - 1) // 3)
    assert trial_prefix(36, 3) == 2
    # A prime above the sieve threshold must survive and be returned by prefixes.
    assert trial_prefix(97, ceil_sqrt((97).bit_length())) is None
    # Integer-only schedule inequalities for a broad parameter range.
    for n in range(2, 100001):
        exponent = ceil_sqrt(n)
        height = (n - 1) // exponent
        assert exponent * exponent >= n and (exponent - 1) ** 2 < n
        assert height <= exponent
        assert exponent + height <= 2 * exponent
        statistics['schedule_checks'] += 1
    result = {'scope': 'Finite exact controls for a conditional transfer; divisor enumeration is a reference oracle, not the assumed fast circuit evaluator.',
              'counts': dict(statistics), 'samples': samples,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
