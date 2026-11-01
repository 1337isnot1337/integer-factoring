#!/usr/bin/env python3
"""Exact complete-piece ranks for reflected complementary-pair labels.

The reduction is proved in research/paired-mode-johnson.md.  This module
uses subset-intersection blocks, not exponentially larger four-label
tensor tables.  It makes no improved asymptotic rate claim.
"""

from functools import lru_cache
from itertools import combinations
from math import comb

from check_modular_crossing import rank_two
from fermionic_piece_polynomial import binomial_polynomial, polynomial_product


@lru_cache(None)
def subsets(size, weight):
    return tuple(sum(1 << i for i in selected)
                 for selected in combinations(range(size), weight))


@lru_cache(None)
def trinomial_power(exponent):
    """Bitset for (1+Z+Z^2)^exponent in F_2[Z]."""
    if not exponent:
        return 1
    bit = exponent & -exponent
    previous = trinomial_power(exponent ^ bit)
    return previous ^ (previous << bit) ^ (previous << (2 * bit))


@lru_cache(None)
def entry_numerator(inside_distance, outside_distance, outside_weight):
    return polynomial_product(binomial_polynomial(inside_distance + outside_distance),
                              trinomial_power(outside_weight - outside_distance))


def complete_entry(q, x, y, u, v):
    """Return numerator bitset and Laurent offset for one full tensor entry."""
    assert all(label.bit_count() == q for label in (x, y, u, v))
    xor = x ^ y
    if xor != v ^ u:
        return 0, 0
    inside_weight = xor.bit_count() // 2
    outside_weight = q - inside_weight
    inside_distance = inside_weight - (x & v & xor).bit_count()
    outside_distance = outside_weight - (x & v & ~xor).bit_count()
    return entry_numerator(inside_distance, outside_distance, outside_weight), -outside_weight


def coefficient_blocks(pair_count, q, inside_weight):
    """Coefficient matrices for one fixed XOR support of size 2*inside_weight."""
    j = inside_weight
    k = q - j
    basis = tuple((a, b) for a in subsets(2 * j, j)
                  for b in subsets(pair_count - 2 * j, k))
    dimension = len(basis)
    result = {shift: [0] * dimension for shift in range(-k, q + 1)}
    for row, (a, b) in enumerate(basis):
        for column, (u, v) in enumerate(basis):
            polynomial = entry_numerator(j - (a & u).bit_count(),
                                         k - (b & v).bit_count(), k)
            while polynomial:
                bit = polynomial & -polynomial
                polynomial ^= bit
                shift = bit.bit_length() - 1 - k
                result[shift][row] |= 1 << column
    return result


def rank_profile(pair_count, q):
    if not 0 < q < pair_count:
        raise ValueError("use a nontrivial constant-weight label set")
    row_ranks = {shift: 0 for shift in range(-q, q + 1)}
    blocks = []
    for j in range(min(q, pair_count - q) + 1):
        matrices = coefficient_blocks(pair_count, q, j)
        ranks = {shift: rank_two(rows) for shift, rows in matrices.items()}
        multiplicity = comb(pair_count, 2 * j)
        for shift, rank in ranks.items():
            row_ranks[shift] += multiplicity * rank
        blocks.append({"inside_weight": j,
                       "dimension": len(next(iter(matrices.values()))),
                       "multiplicity": multiplicity,
                       "ranks": ranks})
    profile = [[shift, row_ranks[shift], row_ranks[-shift]] for shift in row_ranks]
    return {"pair_count": pair_count, "occupied_pairs": q,
            "label_dimension": comb(pair_count, q), "profile": profile,
            "maximum_rank_product": max(r * c for _, r, c in profile),
            "blocks": blocks}


def paired_labels(pair_count, q, mode_box):
    if mode_box < 2 * pair_count or mode_box & (mode_box - 1):
        raise ValueError("choose a power-of-two box with at least pair_count complementary pairs")
    pairs = tuple((2 * i + 1, 2 * mode_box - (2 * i + 1)) for i in range(pair_count))
    return tuple(tuple(sorted(part for i in range(pair_count) if mask >> i & 1 for part in pairs[i]))
                 for mask in subsets(pair_count, q))


def expand_power(polynomial, step, offset):
    result = {}
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        result[step * (offset + bit.bit_length() - 1)] = 1
    return result
