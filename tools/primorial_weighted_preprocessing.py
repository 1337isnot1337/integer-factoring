#!/usr/bin/env python3
"""Exact divisor-count preprocessing using prime-exponent differences.

The ordinary-integer mode is practical.  The streamed descriptor mode
implements the smaller-space specification without storing candidate products,
divisor counts, prime tables, or memoized product bits.  Python interpreter
overhead is not a finite-tape measurement; the companion note proves the bound.

The balanced streamed multiplication primitive follows the repository's
exact_divisor_maximum.py at dd28d375ea65863356c67b051f5d4f217562475d.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from math import isqrt, prod
from time import monotonic


Descriptor = tuple[tuple[int, int], ...]


def validate_bits(n: int) -> None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("input bit bound must be a positive integer")


def lower_weight(j: int) -> int:
    if not isinstance(j, int) or isinstance(j, bool) or j < 1:
        raise ValueError("support index must be a positive integer")
    return max(j, j * ((j + 1).bit_length() - 1) // 4)


def weighted_descriptors(n: int):
    """Every descriptor with sum m_j w_j <= n-1, in a unique DFS order."""
    validate_bits(n)
    prefix: list[tuple[int, int]] = []

    def visit(maximum_index: int, remaining: int):
        yield tuple(prefix)
        for j in range(min(maximum_index, remaining), 0, -1):
            weight = lower_weight(j)
            for multiplicity in range(1, remaining // weight + 1):
                prefix.append((j, multiplicity))
                yield from visit(j - 1, remaining - multiplicity * weight)
                prefix.pop()

    yield from visit(n - 1, n - 1)


def exponent_at(descriptor: Descriptor, i: int) -> int:
    """Exponent of the i-th prime, with prime positions numbered from one."""
    if i < 1:
        raise ValueError("prime position must be positive")
    return sum(m for j, m in descriptor if j >= i)


def descriptor_exponents(descriptor: Descriptor) -> list[int]:
    return [exponent_at(descriptor, i) for i in range(1, max((j for j, _ in descriptor), default=0) + 1)]


def trial_prime(i: int) -> int:
    """Regenerate one prime without a retained list or sieve."""
    if i < 1:
        raise ValueError("prime position must be positive")
    seen = 0
    for candidate in range(2, 4 * i * i + 1):
        if all(candidate % divisor for divisor in range(2, isqrt(candidate) + 1)):
            seen += 1
            if seen == i:
                return candidate
    raise AssertionError("elementary prime bound failed")


def streamed_product_bit(factor_count, factor_at, leaf_bits, bit):
    """Recompute one product bit with bounded carries and no memoization."""
    if factor_count < 0 or leaf_bits < 1 or bit < 0:
        raise ValueError("invalid streamed-product parameters")

    def get(lo, hi, k):
        count = hi - lo
        if count == 0:
            return int(k == 0)
        if k >= count * leaf_bits:
            return 0
        if count == 1:
            value = factor_at(lo)
            if value < 0 or value.bit_length() > leaf_bits:
                raise ValueError("factor exceeds declared bit bound")
            return (value >> k) & 1
        mid = (lo + hi) // 2
        left_bits = (mid - lo) * leaf_bits
        right_bits = (hi - mid) * leaf_bits
        carry = 0
        digit = 0
        for column in range(k + 1):
            total = carry
            for a in range(max(0, column - right_bits + 1), min(column, left_bits - 1) + 1):
                if get(lo, mid, a):
                    total += get(mid, hi, column - a)
            digit = total & 1
            carry = total >> 1
        return digit

    return get(0, factor_count, bit)


def descriptor_product_bit(descriptor: Descriptor, n: int, bit: int, divisor_count: bool = False) -> int:
    """One bit of the integer or its divisor count, from its descriptor."""
    if divisor_count:
        count = max((j for j, _ in descriptor), default=0)
        return streamed_product_bit(count, lambda a: exponent_at(descriptor, a + 1) + 1,
                                    max(1, n.bit_length()), bit)
    count = sum(j * m for j, m in descriptor)
    leaf_bits = max(1, (4 * n * n).bit_length())

    def factor_at(a):
        previous = 0
        for i in range(1, max((j for j, _ in descriptor), default=0) + 1):
            previous += exponent_at(descriptor, i)
            if a < previous:
                return trial_prime(i)
        raise AssertionError("factor index outside compressed descriptor")

    return streamed_product_bit(count, factor_at, leaf_bits, bit)


def compare_descriptors(left: Descriptor, right: Descriptor, n: int, divisor_count: bool) -> int:
    bound = n if divisor_count else max(1, (n - 1) * (4 * n * n).bit_length())
    for bit in range(bound - 1, -1, -1):
        a = descriptor_product_bit(left, n, bit, divisor_count)
        b = descriptor_product_bit(right, n, bit, divisor_count)
        if a != b:
            return 1 if a > b else -1
    return 0


@dataclass
class DescriptorMaximum:
    input_bits: int
    descriptor: Descriptor
    surrogate_patterns_checked: int
    feasible_patterns_checked: int
    maximum_support: int


def streamed_exact_descriptor(n: int) -> DescriptorMaximum:
    """Complete smaller-space search; return a descriptor, never a product."""
    validate_bits(n)
    best: Descriptor = ()
    visited = feasible = support = 0
    bound = max(1, (n - 1) * (4 * n * n).bit_length())
    for descriptor in weighted_descriptors(n):
        visited += 1
        support = max(support, len(descriptor))
        if any(descriptor_product_bit(descriptor, n, bit) for bit in range(n, bound)):
            continue
        feasible += 1
        comparison = compare_descriptors(descriptor, best, n, True)
        if comparison > 0 or (comparison == 0 and compare_descriptors(descriptor, best, n, False) < 0):
            best = descriptor
    return DescriptorMaximum(n, best, visited, feasible, support)


def descriptor_values(descriptor: Descriptor) -> tuple[int, int]:
    """Practical helper: materializes full integers, outside the space claim."""
    exponents = descriptor_exponents(descriptor)
    return (prod(trial_prime(i + 1) ** e for i, e in enumerate(exponents)),
            prod(e + 1 for e in exponents))


def fast_weighted_maximum(n: int) -> DescriptorMaximum:
    """Same enumeration with ordinary-integer feasibility and comparisons."""
    validate_bits(n)
    best: Descriptor = ()
    best_value = best_count = 1
    visited = feasible = support = 0
    primes = [trial_prime(i) for i in range(1, n)]
    for descriptor in weighted_descriptors(n):
        visited += 1
        support = max(support, len(descriptor))
        exponents = descriptor_exponents(descriptor)
        value = prod(p ** e for p, e in zip(primes, exponents))
        if value >= 1 << n:
            continue
        feasible += 1
        count = prod(e + 1 for e in exponents)
        if count > best_count or (count == best_count and value < best_value):
            best = descriptor
            best_value, best_count = value, count
    return DescriptorMaximum(n, best, visited, feasible, support)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bits", type=int, nargs="+")
    parser.add_argument("--streamed", action="store_true", help="use descriptor-only arithmetic; intentionally slow")
    args = parser.parse_args()
    started = monotonic()
    function = streamed_exact_descriptor if args.streamed else fast_weighted_maximum
    results = [asdict(function(n)) for n in args.bits]
    print(json.dumps({"mode": "streamed-descriptor" if args.streamed else "ordinary-integer",
                      "results": results, "seconds": round(monotonic() - started, 3)}, indent=2))


if __name__ == "__main__":
    main()
