#!/usr/bin/env python3
"""Homogeneous matching-label Wick polynomials over F_2.

The all-even-particle unit dual construction is proved in
research/matching-gram-unipotence.md. These raw monomial polynomials have
the same two ranks as the dual-normalized pieces, by an invertible change
of the x label. No asymptotic rank bound is asserted.
"""

from functools import lru_cache

from fermionic_piece_polynomial import polynomial_product
from matching_code_lift import edge_polynomial


def homogeneous_polynomial(vertices, x, v, u, y, collapsed=False):
    """Return the independent-bit polynomial, or its collapsed t curve.

    The collapse sets all Z_i=Z and L=Z^2, then t=Z^8. It is implemented
    directly with F=1+s,G=s and extracts the even s powers. Every matching
    diagram has even F and G degrees, so this extraction is exact.
    """
    groups = ([('U', e) for e in u] + [('B', e) for e in v]
              + [('A', e) for e in x] + [('Y', e) for e in y])
    edges = {}
    for i, (first, e) in enumerate(groups):
        for j in range(i + 1, len(groups)):
            second, f = groups[j]
            pair = first, second
            value = 0
            if collapsed:
                if pair == ('U', 'B') and e == f:
                    value = 3  # F
                elif pair == ('U', 'A') and e == f:
                    value = 6  # FG
                elif pair == ('U', 'Y') and e == f:
                    value = 2  # G
                elif pair == ('B', 'A') and not (e & f):
                    value = 2
                elif pair == ('B', 'Y') and not (e & f):
                    value = 6
                elif pair == ('A', 'Y') and not (e & f):
                    value = 3
            else:
                if pair == ('U', 'B') and e == f:
                    value = edge_polynomial(vertices, 0, 0, e)
                    value ^= value << 1
                elif pair == ('U', 'A') and e == f:
                    value = edge_polynomial(vertices, 1, e, e)
                    value ^= value << 1
                elif pair == ('U', 'Y') and e == f:
                    value = edge_polynomial(vertices, 1, e)
                elif pair == ('B', 'A') and not (e & f):
                    value = edge_polynomial(vertices, 1, f)
                elif pair == ('B', 'Y') and not (e & f):
                    value = edge_polynomial(vertices, 1, f, e)
                    value ^= value << 1
                elif pair == ('A', 'Y') and not (e & f):
                    value = edge_polynomial(vertices, 0, 0, e)
                    value ^= value << 1
            if value:
                edges[i, j] = value

    @lru_cache(None)
    def pfaffian(remaining):
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        answer = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            edge = edges.get((first, bit.bit_length() - 1), 0)
            if edge:
                answer ^= polynomial_product(edge, pfaffian(rest ^ bit))
        return answer

    answer = pfaffian((1 << len(groups)) - 1)
    if not collapsed:
        return answer
    assert all(not ((answer >> exponent) & 1)
               for exponent in range(1, answer.bit_length(), 2))
    return sum(((answer >> (2 * exponent)) & 1) << exponent
               for exponent in range((answer.bit_length() + 1) // 2))
