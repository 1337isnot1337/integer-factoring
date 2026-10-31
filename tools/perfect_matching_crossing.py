"""Polynomial crossings on the alphabet of perfect matchings.

The exact entry formula comes from a fixed-degree matching sector of
the neutral source kernel.  No favorable pressure bound is assumed.
"""

from functools import lru_cache
from itertools import product
from math import comb, factorial
from collections import Counter

from neutral_reflected_crossing import binomial_polynomial


@lru_cache(None)
def perfect_matchings(vertices):
    """Sorted tuples of edges, with vertices supplied as a tuple."""
    vertices = tuple(vertices)
    if not vertices:
        return ((),)
    if len(vertices) % 2 or len(set(vertices)) != len(vertices):
        raise ValueError("use an even number of distinct vertices")
    first, tail = vertices[0], vertices[1:]
    return tuple(((first, second),) + rest
                 for second in tail
                 for rest in perfect_matchings(tuple(v for v in tail if v != second)))


def polynomial_entry(x, y, u, v):
    """Binary polynomial, encoded by its coefficient bits."""
    x, y, u, v = map(frozenset, (x, y, u, v))
    if x - y != v - u or y - x != u - v:
        return 0
    common = x & y
    retained = len(common & u & v)
    return binomial_polynomial(len(common) - retained) << retained


def dense_seed(pairs):
    labels = perfect_matchings(tuple(range(1, 2 * pairs + 1)))
    q = len(labels)
    sets = tuple(frozenset(label) for label in labels)
    values = [[0, 0] for _ in range(q ** 4)]
    for x, y, u, v in product(range(q), repeat=4):
        value = polynomial_entry(sets[x], sets[y], sets[u], sets[v])
        assert value.bit_count() % 2 == int(x == v and y == u)
        values[x + q * y + q * q * u + q ** 3 * v] = [value, 0]
    return {"constructor": "perfect_matching_sector", "pairs": pairs,
            "dimensions": [q] * 4, "labels": labels, "values": values,
            "denominator_degree": 0, "polynomial_degree": pairs,
            "strict_endpoint": True}


def common_block(pairs):
    labels = perfect_matchings(tuple(range(2 * pairs)))
    sets = tuple(frozenset(label) for label in labels)
    return [[binomial_polynomial(pairs - len(a & b)) << len(a & b)
             for b in sets] for a in sets]


def reflected_parts(matching, period):
    return tuple(sorted(period - 1 - (1 << i) - (1 << j)
                        for i, j in matching))


def matching_count(pairs):
    result = 1
    for k in range(1, pairs + 1):
        result *= 2 * k - 1
    return result


def disjoint_pair_count(pairs):
    """Ordered pairs of perfect matchings with no common edge."""
    return matching_count(pairs) * sum(
        (-1) ** k * comb(pairs, k) * matching_count(pairs - k)
        for k in range(pairs + 1))


def row_block_multiplicities(pairs):
    return {r: comb(2 * pairs, 2 * r) * disjoint_pair_count(pairs - r)
            for r in range(pairs + 1)}


def column_zero_counts(pairs):
    """Complete zero Smith multiplicities, by a labeled-cycle recurrence."""
    values = [Counter({0: 1})]
    for n in range(1, pairs + 1):
        current = Counter({order + 1: (2 * n - 1) * count
                           for order, count in values[n - 1].items()})
        for k in range(2, n + 1):
            coefficient = factorial(2 * n - 1) // (2 * factorial(2 * n - 2 * k))
            cycle_order = 2 * (k & -k)
            for order, count in values[n - k].items():
                current[order] += coefficient * count
                current[order + cycle_order] += coefficient * count
        assert sum(current.values()) == matching_count(n) ** 2
        values.append(current)
    return dict(sorted(values[-1].items()))


def binary_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def toeplitz_profile(matrix, maximum_order):
    """Exact profile from ranks modulo t^k, with a stated order bound.

    For a common block, G(w)G(w^-1)=I proves maximum_order=r. The
    final full slope is checked, rather than treating missing pivots
    as zeros. This implementation uses binary integers only.
    """
    n = len(matrix)
    coefficient_rows = [
        [sum(((value >> power) & 1) << j for j, value in enumerate(row))
         for row in matrix]
        for power in range(max((value.bit_length() for row in matrix for value in row), default=0))]
    ranks = [0]
    counts = Counter()
    old_slope = 0
    for depth in range(1, maximum_order + 2):
        rows = []
        for time in range(depth):
            for i in range(n):
                row = 0
                for power in range(min(time + 1, len(coefficient_rows))):
                    row ^= coefficient_rows[power][i] << ((time - power) * n)
                rows.append(row)
        ranks.append(binary_rank(rows))
        slope = ranks[-1] - ranks[-2]
        counts[depth - 1] = slope - old_slope
        old_slope = slope
    assert old_slope == n, "the asserted finite order bound was insufficient"
    assert all(count >= 0 for count in counts.values())
    return {order: count for order, count in counts.items() if count}, ranks


def row_zero_counts(pairs):
    counts = Counter()
    for r, multiplicity in row_block_multiplicities(pairs).items():
        if multiplicity:
            block, _ = toeplitz_profile(common_block(r), r)
            for order, count in block.items():
                counts[order] += multiplicity * count
    assert sum(counts.values()) == matching_count(pairs) ** 2
    return dict(sorted(counts.items()))
