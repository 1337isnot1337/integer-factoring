#!/usr/bin/env python3
"""Independent exact controls for the adaptive-benefit divisor oracle."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from adaptive_divisor_preprocessing import BenefitSearch, adaptive_maximum
from exchange_divisor_preprocessing import fast_exchange_maximum
from primorial_weighted_preprocessing import descriptor_values
from check_primorial_weighted_preprocessing import (
    divisor_sieve_maxima, ordinary_exponent_oracle,
)


@lru_cache(maxsize=None)
def rational_log_interval(numerator, denominator, terms=40):
    assert denominator <= numerator <= 2 * denominator
    z = Q(numerator - denominator, numerator + denominator)
    z2 = z * z
    power = z
    total = Q(0)
    for j in range(terms):
        total += 2 * power / (2 * j + 1)
        power *= z2
    tail = 2 * power / ((2 * terms + 1) * (1 - z2))
    return total, total + tail


@lru_cache(maxsize=None)
def integer_log_interval(value):
    exponent = value.bit_length() - 1
    shift = max(0, exponent - 96)
    top = value >> shift
    lo = rational_log_interval(top, 1 << (exponent - shift))[0]
    hi = rational_log_interval(top + int(shift > 0),
                               1 << (exponent - shift))[1]
    l2, h2 = rational_log_interval(2, 1)
    return exponent * l2 + lo, exponent * h2 + hi


def marginal_interval(p, level, numerator, denominator):
    lp, hp = integer_log_interval(p)
    lr, hr = rational_log_interval(level + 1, level)
    return Q(numerator, denominator) * lp - hr, Q(numerator, denominator) * hp - lr


def independent_primes(count):
    result = []
    value = 2
    while len(result) < count:
        if all(value % d for d in range(2, __import__('math').isqrt(value) + 1)):
            result.append(value)
        value += 1
    return result


def check_superior(result):
    if result.input_bits == 1:
        return 0
    primes = independent_primes(len(result.exponents) + 20)
    value = result.superior_integer
    superior_exponents = []
    for p in primes:
        exponent = 0
        while value % p == 0:
            value //= p
            exponent += 1
        if exponent:
            superior_exponents.append(exponent)
        elif value == 1:
            break
    assert value == 1
    checks = 0
    numerator, denominator = result.epsilon_numerator, result.epsilon_denominator
    for p, exponent in zip(primes, superior_exponents):
        assert marginal_interval(p, exponent, numerator, denominator)[1] <= 0
        assert marginal_interval(p, exponent + 1, numerator, denominator)[0] > 0
        checks += 2
    assert marginal_interval(primes[len(superior_exponents)], 1,
                             numerator, denominator)[0] > 0
    checks += 1
    # Independently reconstruct enough of the immediately lower mesh point
    # to establish that its superior product reaches the strict input bound.
    lower_value = 1
    for p in primes:
        exponent = 0
        for level in range(1, 2 * result.input_bits + 2):
            interval = marginal_interval(p, level, numerator - 1, denominator)
            if interval[1] <= 0:
                exponent = level
            elif interval[0] > 0:
                break
            else:
                # An exact power comparison settles an unusually close case.
                if (level + 1) ** denominator >= level ** denominator * p ** (numerator - 1):
                    exponent = level
                else:
                    break
        lower_value *= p ** exponent
        if lower_value >= 1 << result.input_bits:
            break
        assert exponent > 0
    assert lower_value >= 1 << result.input_bits
    return checks


def check_terminal_certificate(result):
    if result.input_bits == 1:
        return
    eps = Q(result.epsilon_numerator, result.epsilon_denominator)
    n_lower, n_upper = integer_log_interval(result.superior_integer)
    d_lower, d_upper = integer_log_interval(result.maximum_divisors)
    t_lower, t_upper = integer_log_interval(result.superior_divisors)
    l2, h2 = rational_log_interval(2, 1)
    upper = eps * result.input_bits * h2 - d_lower + t_upper - eps * n_lower
    stage = result.stages[-1]
    radius = Q(stage['radius_numerator'], stage['radius_denominator'])
    assert upper <= radius
    assert stage['certified']


def check_by_exponent_enumeration(result):
    """Rebuild the terminal cover, then enumerate in independent prime order."""
    if result.input_bits == 1:
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
    limit = 1 << result.input_bits
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

        def visit(index, cap, value, divisors, partial_benefit):
            nonlocal best_value, best_divisors, leaves
            if index == support:
                heights = [sum(e >= a for e in selected) for a in range(1, len(windows) + 1)]
                assert all(lo <= h <= hi for h, (lo, hi) in zip(heights, windows))
                leaves += 1
                if divisors > best_divisors or (divisors == best_divisors and value < best_value):
                    best_value, best_divisors = value, divisors
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


def check_increment_comparisons():
    engine = BenefitSearch(8)
    comparisons = 0
    for numerator in range(1, engine.denominator + 1, 7):
        for p in (2, 3, 5, 11, 31, 67):
            for exponent in (1, 2, 3, 5, 9, 17):
                got = engine.profitable_increment(p, exponent, numerator)
                expected = ((exponent + 1) ** engine.denominator >=
                            exponent ** engine.denominator * p ** numerator)
                assert got == expected
                comparisons += 1
    assert engine.profitable_increment(2, 1, engine.denominator)
    assert not engine.profitable_increment(3, 1, engine.denominator)
    return comparisons + 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--large', action='store_true')
    parser.add_argument('--progress', action='store_true')
    args = parser.parse_args()
    sieve = divisor_sieve_maxima(16)
    table = []
    cases = list(range(1, 33)) + [48, 64, 96, 128, 192, 256]
    if args.large:
        cases += [512, 1024, 2048]
    checks = 0
    for n in cases:
        result = adaptive_maximum(n)
        value, count = result.least_witness, result.maximum_divisors
        exchange = None
        if n <= 256:
            exchange = fast_exchange_maximum(n)
            assert (value, count) == descriptor_values(exchange.descriptor)
        if n <= 16:
            assert (value, count) == sieve[n]
        elif n <= 24:
            assert (value, count) == ordinary_exponent_oracle(n)[:2]
        assert value < 1 << n
        primes = independent_primes(len(result.exponents))
        assert value == __import__('math').prod(p ** e for p, e in zip(primes, result.exponents))
        assert count == __import__('math').prod(e + 1 for e in result.exponents)
        assert result.divisor_parameter_bits == (count - 1).bit_length()
        checks += check_superior(result)
        check_terminal_certificate(result)
        independent_leaves = check_by_exponent_enumeration(result)
        if args.progress:
            print('checked', n, 'independent profiles', independent_leaves, file=__import__('sys').stderr, flush=True)
        table.append({'input_bits': n, 'least_witness': value, 'maximum_divisors': count,
                      'adaptive_feasible_profiles': result.visited_profiles,
                      'adaptive_search_nodes': result.total_search_nodes,
                      'independent_terminal_profiles': independent_leaves,
                      'exchange_feasible_profiles': exchange.feasible_candidates if exchange else None,
                      'benefit_stages': len(result.stages),
                      'terminal_radius': ([result.stages[-1]['radius_numerator'], result.stages[-1]['radius_denominator']]
                                         if result.stages else None)})
    result = {
        'scope': 'Exact oracle comparisons and independent rational-interval certificates; no timing measurement certifies the asymptotic bound.',
        'input_cases': len(cases), 'superior_marginal_checks': checks,
        'exact_power_comparison_checks': check_increment_comparisons(),
        'divisor_sieve_integers': (1 << 16) - 1, 'comparisons': table,
        'oracle_sha256': hashlib.sha256(Path(__file__).with_name('adaptive_divisor_preprocessing.py').read_bytes()).hexdigest(),
    }
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
