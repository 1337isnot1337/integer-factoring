#!/usr/bin/env python3
"""Characteristic-two matching crossings with monomial self-dual labels."""

from functools import lru_cache

from fermionic_piece_polynomial import polynomial_product


def matching_labels(vertices):
    """Enumerate perfect matchings as sorted tuples of two-bit edge masks."""
    if vertices < 0 or vertices % 2:
        raise ValueError("The number of physical vertices must be even.")

    def generate(remaining):
        if not remaining:
            yield ()
            return
        first = remaining & -remaining
        rest = remaining ^ first
        candidates = rest
        while candidates:
            second = candidates & -candidates
            candidates ^= second
            for tail in generate(rest ^ second):
                yield tuple(sorted((first | second,) + tail))

    return tuple(generate((1 << vertices) - 1))


def validate_labels(labels, vertices):
    full = (1 << vertices) - 1
    for label in labels:
        used = 0
        for edge in label:
            if edge <= 0 or edge.bit_count() != 2 or edge & used or edge & ~full:
                raise ValueError("Each label must be a perfect matching.")
            used |= edge
        if used != full:
            raise ValueError("Each label must cover every physical vertex.")


def lifted_polynomial(operator_a, operator_b, output, input_, vertices,
                      *, kernel="intersection"):
    """Return the packed six-block Pfaffian.

    The L radix is 2*vertices+1 and every physical-variable radix is5.
    Each argument is a perfect matching on the same physical vertices.
    kernel='identity' is an exploratory comparison, not the copy fold.
    """
    validate_labels((operator_a, operator_b, output, input_), vertices)
    if kernel not in ("intersection", "disjointness", "identity"):
        raise ValueError("Unknown annihilation kernel.")
    radix = 2 * vertices + 1
    places = [radix * 5**i for i in range(vertices)]

    @lru_cache(None)
    def monomial(mask):
        return 1 << sum(places[i] for i in range(vertices) if mask >> i & 1)

    @lru_cache(None)
    def binomial(mask):
        answer = 1
        for i in range(vertices):
            if mask >> i & 1:
                answer ^= answer << places[i]
        return answer

    def annihilation(first, second):
        if kernel == "intersection":
            return bool(first & second)
        if kernel == "disjointness":
            return not first & second
        return first == second

    particles = (
        [("U", edge) for edge in output] + [("B", edge) for edge in operator_b]
        + [("A", edge) for edge in operator_a] + [("Y", edge) for edge in input_]
    )
    count = len(particles)
    entries = [[0] * count for _ in range(count)]
    for i, (left, e) in enumerate(particles):
        for j in range(i + 1, count):
            right, f = particles[j]
            pair = left, right
            value = 0
            if pair == ("U", "B") and e == f:
                value = polynomial_product(3, binomial(f))
            elif pair == ("U", "A") and e == f:
                value = polynomial_product(6, polynomial_product(monomial(e), binomial(f)))
            elif pair == ("U", "Y") and e == f:
                value = monomial(f) << 1
            elif pair == ("B", "A") and annihilation(e, f):
                value = monomial(f) << 1
            elif pair == ("B", "Y") and annihilation(e, f):
                value = polynomial_product(6, polynomial_product(monomial(f), binomial(e)))
            elif pair == ("A", "Y") and annihilation(e, f):
                value = polynomial_product(3, binomial(e))
            entries[i][j] = value

    @lru_cache(None)
    def pfaffian(remaining):
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        value = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            other = bit.bit_length() - 1
            if entries[first][other]:
                value ^= polynomial_product(entries[first][other], pfaffian(rest ^ bit))
        return value

    return pfaffian((1 << count) - 1)
