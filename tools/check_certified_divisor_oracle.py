#!/usr/bin/env python3
"""Independent certification for arbitrary-bound divisor queries and records."""
from __future__ import annotations

import argparse
from bisect import bisect_left
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from math import prod
from pathlib import Path

from certified_divisor_oracle import maximum_below, is_highly_composite
from adaptive_divisor_preprocessing import adaptive_maximum
from check_adaptive_divisor_preprocessing import (
    integer_log_interval, rational_log_interval, marginal_interval,
    independent_primes,
)

def check_by_exponent_enumeration(result):
    """Rebuild the terminal cover, then enumerate in independent prime order."""
    if result.exclusive_bound == 2:
        return 1
    stage = result.stages[-1]
    radius = Q(stage['radius_numerator'], stage['radius_denominator'])
    expected = stage['height_windows']
    upper_support = expected[0][1]
    primes = independent_primes(upper_support + 2)
    windows = []
    for level in range(1, len(expected) + 1):
        mandatory = possible = 0
        for p in primes:
            lo, hi = marginal_interval(p, level, result.epsilon_numerator,
                                       result.epsilon_denominator)
            if lo > radius:
                break
            possible += 1
            if hi < -radius:
                mandatory += 1
        windows.append((mandatory, min(result.input_bits - 1, possible)))
    assert windows == [tuple(pair) for pair in expected]
    limit = result.exclusive_bound
    best_value = best_divisors = 1
    leaves = 0
    benefit_scale = 1 << 64
    superior_exponents = []
    remaining = result.superior_integer
    for p in primes:
        exponent = 0
        while remaining % p == 0:
            remaining //= p
            exponent += 1
        superior_exponents.append(exponent)
    assert remaining == 1
    cost_cache = {}

    def lower_cost(index, exponent):
        key = index, exponent
        if key in cost_cache:
            return cost_cache[key]
        optimal = superior_exponents[index]
        lopt, hopt = integer_log_interval(optimal + 1)
        le, he = integer_log_interval(exponent + 1)
        lp, hp = integer_log_interval(primes[index])
        difference = exponent - optimal
        lo = lopt - he + Q(result.epsilon_numerator, result.epsilon_denominator) * difference * (lp if difference >= 0 else hp)
        cost = max(0, (lo.numerator * benefit_scale) // lo.denominator)
        cost_cache[key] = cost
        return cost

    for support in range(windows[0][0], windows[0][1] + 1):
        if support == 0:
            if all(lo == 0 for lo, _ in windows):
                leaves += 1
            continue
        next_prime = primes[support]
        ranges = []
        for index in range(1, support + 1):
            lower = max([1] + [a for a, (lo, hi) in enumerate(windows, 1) if lo >= index])
            upper = min([len(windows)] + [a - 1 for a, (lo, hi) in enumerate(windows, 1) if hi < index])
            exchange_cap = 0
            while primes[index - 1] ** (exchange_cap + 1) < next_prime:
                exchange_cap += 1
            upper = min(upper, 2 * exchange_cap)
            ranges.append((lower, upper))
        if any(lo > hi for lo, hi in ranges):
            continue
        # Every unvisited support prime must occur at least once.
        suffix = [1] * (support + 1)
        for index in range(support - 1, -1, -1):
            suffix[index] = suffix[index + 1] * primes[index]
        selected = []

        fixed_ones_start = support
        while fixed_ones_start and ranges[fixed_ones_start - 1] == (1, 1):
            fixed_ones_start -= 1
        one_cost_suffix = [0] * (support + 1)
        for index in range(support - 1, -1, -1):
            one_cost_suffix[index] = one_cost_suffix[index + 1] + lower_cost(index, 1)

        def visit(index, cap, value, divisors, partial_benefit):
            nonlocal best_value, best_divisors, leaves
            if index == support:
                heights = [sum(e >= a for e in selected) for a in range(1, len(windows) + 1)]
                assert all(lo <= h <= hi for h, (lo, hi) in zip(heights, windows))
                leaves += 1
                if divisors > best_divisors or (divisors == best_divisors and value < best_value):
                    best_value, best_divisors = value, divisors
                return
            if index >= fixed_ones_start:
                assert cap >= 1
                product = value * suffix[index]
                total_benefit = partial_benefit + one_cost_suffix[index]
                if product >= limit or total_benefit * radius.denominator > radius.numerator * benefit_scale:
                    return
                count = divisors * (1 << (support - index))
                leaves += 1
                if count > best_divisors or (count == best_divisors and product < best_value):
                    best_value, best_divisors = product, count
                return
            lower, upper = ranges[index]
            for exponent in range(lower, min(upper, cap) + 1):
                product = value * primes[index] ** exponent
                if product * suffix[index + 1] >= limit:
                    break
                next_benefit = partial_benefit + lower_cost(index, exponent)
                if next_benefit * radius.denominator > radius.numerator * benefit_scale:
                    continue
                selected.append(exponent)
                visit(index + 1, exponent, product, divisors * (exponent + 1), next_benefit)
                selected.pop()

        visit(0, len(windows), 1, 1, 0)
    assert (best_value, best_divisors) == (result.least_witness, result.maximum_divisors)
    return leaves


def check_certificate(result):
    if result.exclusive_bound == 2:
        assert result.least_witness == result.maximum_divisors == 1
        return 0
    primes = independent_primes(result.stages[-1]['height_windows'][0][1] + 20)
    reference = []
    remaining = result.superior_integer
    for p in primes:
        exponent = 0
        while remaining % p == 0:
            remaining //= p
            exponent += 1
        if exponent:
            reference.append(exponent)
        else:
            assert remaining == 1
            break
    assert remaining == 1
    comparisons = 0
    for p, exponent in zip(primes, reference):
        assert marginal_interval(p, exponent, result.epsilon_numerator,
                                 result.epsilon_denominator)[1] <= 0
        assert marginal_interval(p, exponent + 1, result.epsilon_numerator,
                                 result.epsilon_denominator)[0] > 0
        comparisons += 2
    assert marginal_interval(primes[len(reference)], 1, result.epsilon_numerator,
                             result.epsilon_denominator)[0] > 0
    assert result.superior_integer < result.exclusive_bound
    lower_value = 1
    for p in primes:
        exponent = 0
        for level in range(1, 2 * result.input_bits + 2):
            lo, hi = marginal_interval(p, level, result.epsilon_numerator - 1,
                                       result.epsilon_denominator)
            if hi <= 0:
                exponent = level
            elif lo > 0:
                break
            else:
                if ((level + 1) ** result.epsilon_denominator >=
                        level ** result.epsilon_denominator * p ** (result.epsilon_numerator - 1)):
                    exponent = level
                else:
                    break
        lower_value *= p ** exponent
        if lower_value >= result.exclusive_bound:
            break
        assert exponent > 0
    assert lower_value >= result.exclusive_bound
    eps = Q(result.epsilon_numerator, result.epsilon_denominator)
    nlo, nhi = integer_log_interval(result.superior_integer)
    dlo, dhi = integer_log_interval(result.maximum_divisors)
    tlo, thi = integer_log_interval(result.superior_divisors)
    xlo, xhi = integer_log_interval(result.exclusive_bound)
    upper = eps * (xhi - nlo) + thi - dlo
    final = result.stages[-1]
    radius = Q(final['radius_numerator'], final['radius_denominator'])
    assert upper <= radius
    assert result.least_witness < result.exclusive_bound
    assert result.least_witness == prod(p ** e for p, e in zip(primes, result.exponents))
    assert result.maximum_divisors == prod(e + 1 for e in result.exponents)
    assert result.divisor_parameter_bits == (result.maximum_divisors - 1).bit_length()
    return comparisons + 1


def sieve_controls():
    limit = 1 << 16
    divisors = [0] * limit
    for divisor in range(1, limit):
        for multiple in range(divisor, limit, divisor):
            divisors[multiple] += 1
    records = []
    best = 0
    for value in range(1, limit):
        if divisors[value] > best:
            records.append(value)
            best = divisors[value]
    bounds = sorted(set(range(2, 258)) |
                    {x + offset for x in records for offset in (-1, 0, 1, 2)
                     if 2 <= x + offset <= limit} |
                    {limit * i // 64 for i in range(1, 65)})
    for bound in bounds:
        actual = maximum_below(bound)
        expected = records[bisect_left(records, bound) - 1]
        assert (actual.least_witness, actual.maximum_divisors) == (expected, divisors[expected])
    assert not is_highly_composite(0)
    for value in range(1, 257):
        assert is_highly_composite(value) == (value in records)
    # Tie-sensitive boundaries: 12 is a strict record, 18 ties its divisor
    # count, and the first record 1 is accepted under the stated convention.
    assert is_highly_composite(1) and is_highly_composite(12)
    assert not is_highly_composite(18)
    return {'integers_sieved': limit - 1, 'strict_records': len(records),
            'arbitrary_bound_queries': len(bounds), 'membership_queries': 260}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--large', action='store_true')
    parser.add_argument('--progress', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls = sieve_controls()
    cases = [16, 32, 64, 128, 256, 512]
    if args.large:
        cases += [1024, 2048, 4096, 8192]
    table = []
    marginal_checks = 0
    for n in cases:
        result = maximum_below(1 << n)
        if n <= 2048:
            previous = adaptive_maximum(n)
            assert (result.least_witness, result.maximum_divisors) == (
                previous.least_witness, previous.maximum_divisors)
            earlier_profiles = previous.visited_profiles
        else:
            earlier_profiles = None
        marginal_checks += check_certificate(result)
        independent_profiles = check_by_exponent_enumeration(result)
        table.append({'input_bits': n, 'least_witness': result.least_witness,
                      'maximum_divisors': result.maximum_divisors,
                      'divisor_parameter_bits': result.divisor_parameter_bits,
                      'pruned_feasible_profiles': result.visited_profiles,
                      'pruned_search_nodes': result.total_search_nodes,
                      'unpruned_feasible_profiles': earlier_profiles,
                      'independent_terminal_profiles': independent_profiles,
                      'benefit_stages': len(result.stages),
                      'terminal_radius': [result.stages[-1]['radius_numerator'],
                                          result.stages[-1]['radius_denominator']]})
        if args.progress:
            print('checked', n, 'independent profiles', independent_profiles,
                  file=__import__('sys').stderr, flush=True)
    for bound in (99991, 100000, 100001, 10**12 + 39, (1 << 128) - 17):
        result = maximum_below(bound)
        marginal_checks += check_certificate(result)
        check_by_exponent_enumeration(result)
    report = {'scope': 'Exact arbitrary-bound controls, strict-record membership tests, and independent benefit certificates; no timing measurement is used as a proof.',
              'sieve_controls': controls, 'large_case_count': len(cases),
              'superior_marginal_checks': marginal_checks,
              'additional_arbitrary_bound_certificates': 5,
              'comparisons': table,
              'oracle_sha256': hashlib.sha256(Path(__file__).with_name(
                  'certified_divisor_oracle.py').read_bytes()).hexdigest()}
    text = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
