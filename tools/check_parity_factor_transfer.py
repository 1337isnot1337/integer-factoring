#!/usr/bin/env python3
"""Exact finite controls for the scalar-parity factoring transfer.

Witness queries use explicit integer enumeration as a reference oracle.
This is not an implementation of a fast clique-parity evaluator.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from heapq import heapify, heappop, heappush
import json
from math import isqrt, prod
from pathlib import Path


def integer_root(value, exponent):
    assert value >= 0 and exponent >= 1
    lo, hi = 0, 1 << ((value.bit_length() + exponent - 1) // exponent)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** exponent <= value:
            lo = mid
        else:
            hi = mid - 1
    return lo


def least_prime_sieve(limit):
    least = list(range(limit + 1))
    if limit >= 1:
        least[1] = 1
    for p in range(2, isqrt(limit) + 1):
        if least[p] == p:
            for multiple in range(p * p, limit + 1, p):
                if least[multiple] == multiple:
                    least[multiple] = p
    return least


def factor_from_sieve(value, least):
    result = []
    while value > 1:
        result.append(least[value])
        value //= least[value]
    return result


@lru_cache(maxsize=None)
def witness_support(value, kind, level=0):
    if kind == 'lower':
        return tuple(a for a in range(2, isqrt(value - 1) + 1)
                     if value % a == 0)
    assert kind == 'power' and level >= 1
    power = 1 << level
    return tuple(a for a in range(2, integer_root(value, power) + 1)
                 if value % (a ** power) == 0)


def recover_odd_witness(support, width, statistics):
    assert len(support) % 2
    prefix = 0
    for remaining in range(width - 1, -1, -1):
        zero = 2 * prefix
        parity = sum((a >> remaining) == zero for a in support) % 2
        statistics['prefix_queries'] += 1
        prefix = zero if parity else zero + 1
        assert sum((a >> remaining) == prefix for a in support) % 2
    assert prefix in support
    return prefix


def find_proper_factor(value, padded_width, statistics):
    statistics['factor_calls'] += 1
    top_level = (padded_width - 1).bit_length() - 1
    for level in range(1, top_level + 1):
        support = witness_support(value, 'power', level)
        statistics['power_root_queries'] += 1
        if len(support) % 2:
            factor = recover_odd_witness(support, padded_width, statistics)
            assert 1 < factor < value and value % factor == 0
            statistics['power_successes'] += 1
            return factor
    support = witness_support(value, 'lower')
    statistics['lower_root_queries'] += 1
    if len(support) % 2:
        factor = recover_odd_witness(support, padded_width, statistics)
        assert 1 < factor < value and value % factor == 0
        statistics['lower_successes'] += 1
        return factor
    return None


def parity_factorization(value, statistics):
    width = value.bit_length()
    pending = [value]
    result = []
    calls = 0
    while pending:
        current = pending.pop()
        factor = find_proper_factor(current, width, statistics)
        calls += 1
        if factor is None:
            result.append(current)
        else:
            pending.extend((factor, current // factor))
    result.sort()
    assert calls <= 2 * (width - 1) - 1
    return result


def binary_intervals(lower, upper, width):
    """Disjoint (remaining bits, upper prefix) intervals covering [lower,upper]."""
    if lower > upper:
        return []
    result = []
    current = lower
    while current <= upper:
        alignment = ((current & -current).bit_length() - 1) if current else width
        room = (upper - current + 1).bit_length() - 1
        remaining = min(alignment, room)
        result.append((remaining, current >> remaining))
        current += 1 << remaining
    return result


def parameters(width):
    root = 1
    while 4 * root ** 4 < width:
        root *= 2
    return root, 2 * root, 2 * root ** 2, 8 * root ** 2 + 1, 4 * root ** 4


def multiplication_mass(width):
    root, gates, digit, message, capacity = parameters(width)
    # Per-vertex upper bounds: incident bits plus local-setting entropy.
    classes = [
        (2 * gates, 3 * gates * message + capacity // gates),
        (gates, 5 * gates * message),
        (gates, 4 * gates * message),
        (gates * gates, 5 * message + digit),
        (4 * capacity, 4),
    ]
    total = sum(count * mass for count, mass in classes)
    maximum = max(mass for _, mass in classes)
    assert total == 20 * gates ** 2 * message + 20 * capacity
    assert capacity <= 16 * width
    assert total <= 3200 * width
    assert maximum ** 4 <= 512 ** 4 * width ** 3
    return total, maximum, classes


def power_widths(width, level):
    power = 1 << level
    start = (width + power - 1) // power
    squares = [start << j for j in range(level)]
    return start, squares + [max(width, power * start)]


def carry_product(left, right, width):
    _, _, digit, message, capacity = parameters(width)
    count = digit
    base = 1 << digit
    assert left < (1 << width) and right < (1 << width)
    x = [(left >> (digit * j)) & (base - 1) for j in range(count)]
    y = [(right >> (digit * j)) & (base - 1) for j in range(count)]
    coefficients = [0] * (2 * count)
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            coefficients[i + j] += xi * yj
    carry = 0
    result = 0
    for j, coefficient in enumerate(coefficients):
        assert coefficient < (1 << (message - 1)) + 1
        carry, output_digit = divmod(carry + coefficient, base)
        assert carry < 1 << message
        result += output_digit << (digit * j)
    assert carry == 0 and result == left * right
    assert result < 1 << (2 * capacity)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    statistics = Counter()
    least = least_prime_sieve(65535)

    # Independent exponent formulas establish which witness family must be odd.
    for value in range(2, 65536):
        factors = factor_from_sieve(value, least)
        exponents = tuple(Counter(factors).values())
        maximum = max(exponents)
        if maximum >= 2:
            level = maximum.bit_length() - 1
            power = 1 << level
            count = prod(exponent // power + 1 for exponent in exponents) - 1
            assert count % 2 == 1
            assert level <= (value.bit_length() - 1).bit_length() - 1
            statistics['repeated_prime_completeness'] += 1
        elif len(factors) >= 2:
            count = (1 << (len(factors) - 1)) - 1
            assert count % 2 == 1
            statistics['squarefree_completeness'] += 1
        else:
            assert least[value] == value
            statistics['prime_completeness'] += 1
        statistics['completeness_cases'] += 1

    # Exact arithmetic-predicate supports, with padded recursive calls.
    for value in range(2, 16385):
        found = parity_factorization(value, statistics)
        assert found == factor_from_sieve(value, least)
        assert prod(found) == value
        statistics['complete_factorizations'] += 1
        width = value.bit_length()
        for level in range(1, (width - 1).bit_length()):
            exponents = Counter(factor_from_sieve(value, least)).values()
            predicted = prod(e // (1 << level) + 1 for e in exponents) - 1
            assert len(witness_support(value, 'power', level)) == predicted
            statistics['power_count_controls'] += 1

    # Parity descent may return a larger composite witness, which is valid.
    example_support = witness_support(36, 'power', 1)
    recovered = recover_odd_witness(example_support, 6, statistics)
    assert example_support == (2, 3, 6) and recovered == 6
    # Retaining the full intermediate power prevents a false witness for prime17.
    assert pow(3, 4, 1 << 5) == 17 and 3 ** 4 != 17
    assert 3 not in witness_support(17, 'power', 2)

    for width in range(2, 15):
        limit = 1 << width
        for lower, upper in ((2, limit - 1), (3, limit - 3),
                             (limit // 3, 2 * limit // 3), (2, 1)):
            intervals = binary_intervals(lower, upper, width)
            assert len(intervals) <= 2 * width
            reconstructed = []
            for remaining, prefix in intervals:
                reconstructed.extend(range(prefix << remaining,
                                           (prefix + 1) << remaining))
            assert reconstructed == list(range(lower, upper + 1))
            statistics['interval_controls'] += 1

    for width in range(2, 10001):
        for level in range(1, (width - 1).bit_length()):
            start, widths = power_widths(width, level)
            power = 1 << level
            assert start * power < 2 * width
            assert max(widths) < 2 * width
            assert sum(widths) < 4 * width
            total_mass = 0
            maximum_mass = 0
            for current in widths:
                total, maximum, _ = multiplication_mass(current)
                total_mass += total
                maximum_mass = max(maximum_mass, maximum)
            assert total_mass <= 12800 * width
            assert maximum_mass ** 4 <= 1024 ** 4 * width ** 3
            statistics['power_width_controls'] += 1

    # Check the deterministic least-loaded-bucket mass bound on actual classes.
    for width in (5, 17, 33, 65, 129):
        for level in range(1, (width - 1).bit_length()):
            _, widths = power_widths(width, level)
            masses = []
            for current in widths:
                _, _, classes = multiplication_mass(current)
                masses.extend(mass for count, mass in classes for _ in range(count))
            for buckets in (1, 3, 17, 257):
                heap = [(0, index) for index in range(buckets)]
                heapify(heap)
                for mass in masses:
                    load, index = heappop(heap)
                    heappush(heap, (load + mass, index))
                maximum_load = max(load for load, _ in heap)
                assert maximum_load <= (sum(masses) + buckets - 1) // buckets + max(masses)
                statistics['greedy_cluster_controls'] += 1

    # The large-cutoff estimates use only exact integer arithmetic.
    for groups in range(2, 65):
        for width in (groups ** 40, (groups + 1) ** 40 - 1, groups ** 80):
            clusters = 3 * groups ** 8
            fourth = integer_root(width ** 3, 4)
            fourth += fourth ** 4 < width ** 3
            total = 12800 * width
            fit = (total + clusters - 1) // clusters + 1024 * fourth
            label = max(1, (groups * total + clusters * fit) // groups ** 3)
            entropy = (2 * clusters // groups) * fit
            assert groups ** 2 * label <= 32000 * width
            assert groups * entropy <= 38000 * width
            assert clusters * groups * groups <= 3 * width
            statistics['constructor_schedule_controls'] += 1

    for width in range(1, 9):
        bound = min(1 << width, 16)
        for left in range(bound):
            for right in range(bound):
                carry_product(left, right, width)
                statistics['free_output_carry_controls'] += 1
    for width in (17, 32, 64, 127, 256):
        bound = (1 << width) - 1
        for left, right in ((bound, bound), (bound, 1),
                            (bound // 3, bound // 5), (0, bound)):
            carry_product(left, right, width)
            statistics['free_output_carry_controls'] += 1

    result = {
        'scope': 'Finite exact controls for scalar-parity transfer; witness enumeration is a reference oracle, not a fast evaluator.',
        'counts': dict(statistics),
        'nonleast_witness_control': {'integer': 36, 'power': 2,
                                    'support': list(example_support), 'recovered': recovered},
        'overflow_negative_control': {'integer': 17, 'root': 3, 'power': 4,
                                      'full_power': 81, 'truncated_power': 17},
        'source_revision': 'c91cc3f4947caca8977e81a5d14244f73b5ae2b1',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
