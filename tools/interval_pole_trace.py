#!/usr/bin/env python3
"""Exact dyadic reference for a global interval pole trace; no faster bound."""

from dataclasses import dataclass


def parity_interval(interval, parity):
    lower, upper = interval
    return -((parity - lower) // 2), -((parity - upper) // 2)


@dataclass(frozen=True)
class BilinearRectangle:
    """Zeros of A*u*v+B*u+C*v+D on two nonnegative integer intervals."""

    a: int
    b: int
    c: int
    d: int
    left: tuple
    right: tuple

    def evaluate(self, u, v):
        return self.a * u * v + self.b * u + self.c * v + self.d

    def children(self):
        """The two parity classes on which the bilinear expression is even."""
        assert self.a % 2 == 0 and self.b % 2 and self.c % 2
        for e in (0, 1):
            f = (self.b * e + self.d) % 2
            yield BilinearRectangle(
                2 * self.a, self.b + self.a * f, self.c + self.a * e,
                (self.a * e * f + self.b * e + self.c * f + self.d) // 2,
                parity_interval(self.left, e), parity_interval(self.right, f))


def pole_count(number, left, right):
    """Count factor pairs by exact global recursion, with its tree size exposed.

    This enumerative reference costs O(min(len(left),len(right))*poly(n));
    it does not implement a fast normalized resolvent oracle.
    """
    assert number >= 1 and number % 2
    assert 1 <= left[0] <= left[1] and 1 <= right[0] <= right[1]
    root = BilinearRectangle(2, 1, 1, (1 - number) // 2,
                             parity_interval(left, 1), parity_interval(right, 1))
    stats = {'nodes': 0, 'max_depth': 0, 'max_coefficient_bits': 0,
             'enclosure_prunes': 0, 'corner_solves': 0, 'singleton_solves': 0}

    def visit(state, depth):
        stats['nodes'] += 1
        stats['max_depth'] = max(stats['max_depth'], depth)
        stats['max_coefficient_bits'] = max(
            stats['max_coefficient_bits'],
            *(abs(x).bit_length() for x in (state.a, state.b, state.c, state.d)))
        l, u = state.left
        v, w = state.right
        if l == u or v == w:
            return 0
        # All original factors are positive, so both partial slopes are positive.
        if state.evaluate(l, v) > 0 or state.evaluate(u - 1, w - 1) < 0:
            stats['enclosure_prunes'] += 1
            return 0
        if state.evaluate(l, v) == 0 or state.evaluate(u - 1, w - 1) == 0:
            stats['corner_solves'] += 1
            return 1
        if u - l == 1:
            stats['singleton_solves'] += 1
            slope, constant = state.a * l + state.c, state.b * l + state.d
            return int((-constant) % slope == 0 and v <= -constant // slope < w)
        if w - v == 1:
            stats['singleton_solves'] += 1
            slope, constant = state.a * v + state.b, state.c * v + state.d
            return int((-constant) % slope == 0 and l <= -constant // slope < u)
        return sum(visit(child, depth + 1) for child in state.children())

    return visit(root, 0), stats


def pole_filter_bits(number, left, right, precision):
    """A safe O(n+p)-bit displacement scale for every pair in the rectangle."""
    assert precision >= 1 and 1 <= left[0] <= left[1] and 1 <= right[0] <= right[1]
    if left[0] == left[1] or right[0] == right[1]:
        return precision
    bound = max(abs(number - left[0] * right[0]),
                abs(number - (left[1] - 1) * (right[1] - 1)))
    return precision + bound.bit_length()


def factorial_vandermonde_content(length):
    """v_2(prod_{0<=i<j<length}(j-i)), in polynomial binary-input cost."""
    assert length >= 0
    result, stride = 0, 2
    while stride < length:
        quotient, remainder = divmod(length, stride)
        result += stride * quotient * (quotient - 1) // 2 + quotient * remainder
        stride *= 2
    return result


def interval_vandermonde_content(intervals):
    """Exact content for a fixed number of ordered disjoint integer intervals."""
    result = sum(factorial_vandermonde_content(u - l) for l, u in intervals)
    for i, (lower, upper) in enumerate(intervals):
        assert lower <= upper
        for start, stop in intervals[i + 1:]:
            assert upper <= start
            result += (factorial_vandermonde_content(stop - lower)
                       + factorial_vandermonde_content(start - upper)
                       - factorial_vandermonde_content(stop - upper)
                       - factorial_vandermonde_content(start - lower))
    return result


def square_pole_domains(number, left, right):
    """Pad the shorter factor domain above N; at most two intervals per side."""
    assert number >= 1
    domains = []
    for lower, upper in (left, right):
        assert 1 <= lower <= upper
        domains.append([(min(lower, number + 1), min(upper, number + 1))])
    size = max(domain[0][1] - domain[0][0] for domain in domains)
    for domain in domains:
        missing = size - (domain[0][1] - domain[0][0])
        if missing:
            domain.append((number + 1, number + 1 + missing))
    return tuple(tuple(domain) for domain in domains)
