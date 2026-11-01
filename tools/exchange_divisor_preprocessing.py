#!/usr/bin/env python3
"""Exact maximum divisor counts with next-prime exchange restrictions.

The practical oracle uses ordinary integer arithmetic.  The --streamed mode
uses the elementary descriptor-only multiplication routine, and has the weaker
2**O(log(n)**2) arithmetic bound.  The sharper small-space theorem uses the
classical uniform iterated-multiplication result, not an implementation hidden
inside this program.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from math import prod
from time import monotonic

from primorial_weighted_preprocessing import (
    compare_descriptors, descriptor_product_bit, trial_prime, validate_bits,
)


def maximum_exponent(n: int, next_prime: int) -> int:
    return min(n - 1, 2 * ((next_prime - 1).bit_length()) - 2)


def height_cap(k: int, a: int, next_prime: int) -> int:
    """Number of the first k primes with p**ceil(a/2) < next_prime."""
    power = (a + 1) // 2
    count = 0
    for i in range(1, k + 1):
        if trial_prime(i) ** power >= next_prime:
            break
        count += 1
    return count


def heights_to_descriptor(heights):
    """Transpose Ferrers columns into nonzero exponent differences."""
    answer = []
    for height in heights:
        if answer and answer[-1][0] == height:
            answer[-1] = (height, answer[-1][1] + 1)
        else:
            answer.append((height, 1))
    return tuple(answer)


def exchange_descriptors(n):
    """All restricted profiles; exact size feasibility is checked separately.

    This recursive implementation has O(log(n)**2) working bits in a
    finite-tape translation.  The theorem's packed box enumerator avoids
    per-frame copies of the remaining budget.
    """
    validate_bits(n)
    yield ()
    for k in range(1, n):
        q = trial_prime(k + 1)
        end = maximum_exponent(n, q)
        heights = [k]

        def visit(a, previous, remaining):
            yield heights_to_descriptor(heights)
            if a > end:
                return
            cap = min(previous, remaining, height_cap(k, a, q))
            for height in range(1, cap + 1):
                heights.append(height)
                yield from visit(a + 1, height, remaining - height)
                heights.pop()

        yield from visit(2, k, n - 1 - k)


@dataclass
class ExchangeMaximum:
    input_bits: int
    descriptor: tuple[tuple[int, int], ...]
    candidates_checked: int
    feasible_candidates: int
    maximum_height_depth: int


def streamed_exchange_descriptor(n: int) -> ExchangeMaximum:
    validate_bits(n)
    best = ()
    visited = feasible = depth = 0
    bound = max(1, (n - 1) * (4 * n * n).bit_length())
    for descriptor in exchange_descriptors(n):
        visited += 1
        depth = max(depth, sum(m for _, m in descriptor))
        if any(descriptor_product_bit(descriptor, n, bit) for bit in range(n, bound)):
            continue
        feasible += 1
        comparison = compare_descriptors(descriptor, best, n, True)
        if comparison > 0 or (comparison == 0 and compare_descriptors(descriptor, best, n, False) < 0):
            best = descriptor
    return ExchangeMaximum(n, best, visited, feasible, depth)


def fast_exchange_maximum(n: int) -> ExchangeMaximum:
    """Practical DFS with full integer products; polynomial workspace."""
    validate_bits(n)
    primes = [trial_prime(i) for i in range(1, n + 1)]
    primorials = [1]
    for p in primes:
        primorials.append(primorials[-1] * p)
    limit = 1 << n
    best_value = best_count = 1
    best = ()
    visited = 1
    max_depth = 0
    for k in range(1, n):
        if primorials[k] >= limit:
            break
        q = primes[k]
        end = maximum_exponent(n, q)
        caps = [0] * (end + 1)
        for a in range(1, end + 1):
            power = (a + 1) // 2
            caps[a] = sum(p ** power < q for p in primes[:k])
        heights = [k]

        def visit(a, previous, remaining, value, count):
            nonlocal best, best_value, best_count, visited, max_depth
            visited += 1
            max_depth = max(max_depth, len(heights))
            if count > best_count or (count == best_count and value < best_value):
                best = heights_to_descriptor(heights)
                best_value, best_count = value, count
            if a > end:
                return
            for height in range(1, min(previous, remaining, caps[a]) + 1):
                candidate = value * primorials[height]
                if candidate >= limit:
                    break
                divisor = a ** height
                assert count % divisor == 0
                next_count = (count // divisor) * ((a + 1) ** height)
                heights.append(height)
                visit(a + 1, height, remaining - height, candidate, next_count)
                heights.pop()

        visit(2, k, n - 1 - k, primorials[k], 1 << k)
    return ExchangeMaximum(n, best, visited, visited, max_depth)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bits", type=int, nargs="+")
    parser.add_argument("--streamed", action="store_true")
    args = parser.parse_args()
    function = streamed_exchange_descriptor if args.streamed else fast_exchange_maximum
    started = monotonic()
    results = [asdict(function(n)) for n in args.bits]
    print(json.dumps({"mode": "elementary-streamed" if args.streamed else "ordinary-integer",
                      "results": results, "seconds": round(monotonic() - started, 3)}, indent=2))


if __name__ == "__main__":
    main()
