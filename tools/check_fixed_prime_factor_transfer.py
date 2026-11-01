#!/usr/bin/env python3
"""Exact controls for explicit fixed-characteristic scalar factoring witnesses.

Tuple enumeration is an independent finite reference oracle. It does not
implement a faster scalar clique evaluator.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
from math import comb, isqrt, prod
from pathlib import Path
import sys

# An integrated checkout uses the adjacent, already published parity checker.
# The fallback permits independent staging before the parent integrates this file.
try:
    import check_parity_factor_transfer as parity_controls
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'algorithm-parity' / 'tools'))
    import check_parity_factor_transfer as parity_controls

integer_root = parity_controls.integer_root
least_prime_sieve = parity_controls.least_prime_sieve
factor_from_sieve = parity_controls.factor_from_sieve
multiplication_mass = parity_controls.multiplication_mass


@lru_cache(maxsize=None)
def divisors(value):
    low, high = [], []
    for candidate in range(1, isqrt(value) + 1):
        if value % candidate == 0:
            low.append(candidate)
            if candidate * candidate != value:
                high.append(value // candidate)
    return tuple(low + high[::-1])


@lru_cache(maxsize=None)
def ordered_factorizations(value, length):
    if length == 1:
        return ((value,),)
    return tuple((first,) + suffix for first in divisors(value)
                 for suffix in ordered_factorizations(value // first, length - 1))


def strict_cyclic_minimum(values):
    return all(values < values[shift:] + values[:shift]
               for shift in range(1, len(values)))


@lru_cache(maxsize=None)
def cyclic_support(value, characteristic):
    trivial = (1,) * (characteristic - 1) + (value,)
    return tuple(values for values in ordered_factorizations(value, characteristic)
                 if strict_cyclic_minimum(values) and values != trivial)


@lru_cache(maxsize=None)
def power_support(value, characteristic, level):
    power = characteristic ** level
    result = []
    for product_value in range(2, integer_root(value, power) + 1):
        if value % (product_value ** power) == 0:
            result.extend(ordered_factorizations(product_value, characteristic - 1))
    return tuple(result)


def power_levels(width, characteristic):
    level, power = 1, characteristic
    while power <= width - 1:
        yield level
        level += 1
        power *= characteristic


def encode_tuple(values, width):
    result = 0
    for value in values:
        result = (result << width) | value
    return result


def recover_nonzero_tuple(support, width, characteristic, statistics):
    assert len(support) % characteristic
    length = len(support[0])
    words = tuple(encode_tuple(values, width) for values in support)
    prefix = 0
    for remaining in range(length * width - 1, -1, -1):
        zero = 2 * prefix
        residue = sum(word >> remaining == zero for word in words) % characteristic
        statistics['prefix_queries'] += 1
        prefix = zero if residue else zero + 1
        assert sum(word >> remaining == prefix for word in words) % characteristic
    mask = (1 << width) - 1
    result = tuple((prefix >> (width * index)) & mask
                   for index in range(length - 1, -1, -1))
    assert result in support
    return result


def find_factor(value, width, characteristic, statistics):
    statistics['factor_calls'] += 1
    for level in power_levels(width, characteristic):
        support = power_support(value, characteristic, level)
        statistics['power_queries'] += 1
        if len(support) % characteristic:
            values = recover_nonzero_tuple(support, width, characteristic, statistics)
            result = next(a for a in values if a > 1)
            assert 1 < result < value and value % result == 0
            statistics['power_successes'] += 1
            return result
    support = cyclic_support(value, characteristic)
    statistics['cyclic_queries'] += 1
    if len(support) % characteristic:
        values = recover_nonzero_tuple(support, width, characteristic, statistics)
        result = next(a for a in values if 1 < a < value)
        assert value % result == 0
        statistics['cyclic_successes'] += 1
        return result
    return None


def scalar_factorization(value, characteristic, statistics):
    width = value.bit_length()
    pending, result = [value], []
    calls = 0
    while pending:
        current = pending.pop()
        factor = find_factor(current, width, characteristic, statistics)
        calls += 1
        if factor is None:
            result.append(current)
        else:
            pending.extend((factor, current // factor))
    assert calls <= 2 * width - 3
    return sorted(result)


def composed_widths(width, characteristic, level):
    power = characteristic ** level
    start = (width + power - 1) // power
    aggregate = (characteristic - 1) * start
    widths = [i * start for i in range(1, characteristic - 1)]
    copy_nodes = 0
    for step in range(level):
        current = aggregate * characteristic ** step
        widths.extend(i * current for i in range(1, characteristic))
        copy_nodes += characteristic * current
    widths.append(max(width, power * aggregate))
    return start, widths, copy_nodes


def comparison_trace(left, right, bits):
    state = 0  # equal, less, greater
    trace = []
    for position in range(bits - 1, -1, -1):
        a, b = (left >> position) & 1, (right >> position) & 1
        if state == 0 and a != b:
            state = 1 if a < b else 2
        trace.append(state)
    assert (state == 1) == (left < right)
    assert (state == 0) == (left == right)
    return trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    statistics = Counter()
    least = least_prime_sieve(65535)
    samples = []

    for characteristic in (2, 3, 5, 7):
        # Formula proof independently covers every input in the interval.
        for value in range(2, 65536):
            exponents = tuple(Counter(factor_from_sieve(value, least)).values())
            maximum = max(exponents)
            if maximum >= characteristic:
                level, power = 0, 1
                while power * characteristic <= maximum:
                    level += 1
                    power *= characteristic
                assert level in tuple(power_levels(value.bit_length(), characteristic))
                count = prod(comb(e // power + characteristic - 1, characteristic - 1)
                             for e in exponents) - 1
                assert count % characteristic == characteristic - 1
                statistics['power_completeness'] += 1
            else:
                ordered = prod(comb(e + characteristic - 1, characteristic - 1)
                               for e in exponents)
                assert ordered % characteristic == 0
                count = ordered // characteristic - 1
                if len(exponents) >= 2:
                    assert count % characteristic == characteristic - 1
                else:
                    expected = (pow(exponents[0], -1, characteristic) - 1) % characteristic
                    assert count % characteristic == expected
                assert bool(count % characteristic) == (least[value] != value)
                statistics['cyclic_completeness'] += 1
            statistics['completeness_cases'] += 1

        for value in range(2, 1025):
            result = scalar_factorization(value, characteristic, statistics)
            assert result == factor_from_sieve(value, least)
            assert prod(result) == value
            statistics['complete_factorizations'] += 1

        for value in range(2, 513):
            exponents = tuple(Counter(factor_from_sieve(value, least)).values())
            ordered_count = prod(comb(e + characteristic - 1, characteristic - 1)
                                 for e in exponents)
            fixed = int(all(e % characteristic == 0 for e in exponents))
            assert (ordered_count - fixed) % characteristic == 0
            predicted = (ordered_count - fixed) // characteristic - 1
            assert len(cyclic_support(value, characteristic)) == predicted
            statistics['cyclic_count_controls'] += 1
            for level in power_levels(value.bit_length(), characteristic):
                power = characteristic ** level
                predicted = prod(comb(e // power + characteristic - 1, characteristic - 1)
                                 for e in exponents) - 1
                assert len(power_support(value, characteristic, level)) == predicted
                statistics['power_count_controls'] += 1

        for value in (4, 8, 12, 16, 27, 36, 64, 81, 128, 221):
            support = cyclic_support(value, characteristic)
            samples.append({'integer': value, 'characteristic': characteristic,
                            'cyclic_witnesses': len(support),
                            'cyclic_residue': len(support) % characteristic,
                            'factors': scalar_factorization(value, characteristic, Counter())})

        for width in range(2, 2001):
            for level in power_levels(width, characteristic):
                start, widths, copies = composed_widths(width, characteristic, level)
                assert sum(widths) <= 2 * characteristic ** 2 * width
                assert max(widths) <= 2 * characteristic * width
                assert copies < 2 * characteristic * width
                total = 4 * copies
                maximum = 4
                for current in widths:
                    mass, high, _ = multiplication_mass(current)
                    total += mass
                    maximum = max(maximum, high)
                assert total <= 6500 * characteristic ** 2 * width
                assert maximum ** 4 <= (1024 * characteristic) ** 4 * width ** 3
                statistics['power_mass_controls'] += 1
            cycle_widths = [i * width for i in range(1, characteristic)]
            total = sum(multiplication_mass(w)[0] for w in cycle_widths)
            total += 18 * characteristic * (characteristic - 1) * width
            maximum = max(10, *(multiplication_mass(w)[1] for w in cycle_widths))
            assert total <= 6500 * characteristic ** 2 * width
            assert maximum ** 4 <= (1024 * characteristic) ** 4 * width ** 3
            statistics['cyclic_mass_controls'] += 1

        for groups in range(2, 33):
            for width in (groups ** 40, (groups + 1) ** 40 - 1):
                clusters = 3 * groups ** 8
                root = integer_root(width ** 3, 4)
                root += root ** 4 < width ** 3
                total = 6500 * characteristic ** 2 * width
                fit = (total + clusters - 1) // clusters + 1024 * characteristic * root
                label = max(1, (groups * total + clusters * fit) // groups ** 3)
                entropy = (2 * clusters // groups) * fit
                assert groups ** 2 * label <= 20000 * characteristic ** 2 * width
                assert groups * entropy <= 20000 * characteristic ** 2 * width
                statistics['constructor_schedule_controls'] += 1

    # Comparison chains are deterministic, including equality and leading zeros.
    for bits in range(1, 7):
        for left in range(1 << bits):
            for right in range(1 << bits):
                comparison_trace(left, right, bits)
                statistics['comparison_trace_controls'] += 1
    assert not strict_cyclic_minimum((3, 3, 3))
    assert strict_cyclic_minimum((1, 1, 27))
    assert (1, 1, 27) not in cyclic_support(27, 3)
    # Prime characteristic is required for the full-length orbit argument.
    assert not strict_cyclic_minimum((1, 2, 1, 2))
    assert (1, 2, 1, 2)[2:] + (1, 2, 1, 2)[:2] == (1, 2, 1, 2)

    result = {
        'scope': 'Finite exact controls for linear-scale scalar transfer in each fixed prime characteristic; tuple enumeration is a reference oracle.',
        'characteristics': [2, 3, 5, 7],
        'counts': dict(statistics),
        'samples': samples,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'parity_checker_sha256': hashlib.sha256(Path(parity_controls.__file__).read_bytes()).hexdigest(),
    }
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
