#!/usr/bin/env python3
"""Exact folded Johnson-Gram ranks via elementary Frobenius block splitting."""

from collections import defaultdict
from functools import lru_cache
from math import comb

from check_modular_crossing import rank_two


@lru_cache(None)
def frobenius_weights(highest):
    """B_m is a tensor product of 2D factors with these Frobenius weights."""
    if highest < 0:
        raise ValueError("highest weight must be nonnegative")
    if not highest:
        return ()
    if highest & 1:
        return (1,) + tuple(2 * w for w in frobenius_weights((highest - 1) // 2))
    return (1, 1) + tuple(2 * w for w in frobenius_weights((highest - 2) // 2))


@lru_cache(None)
def tensor_with_natural(highest):
    """(block highest weight, determinant twist, multiplicity) for V*B_m."""
    if not highest:
        return ((1, 0, 1),)
    if highest & 1:
        return ((highest + 1, 0, 1),)
    result = defaultdict(int)
    result[highest - 1, 1] += 2
    for block, twist, multiplicity in tensor_with_natural((highest - 2) // 2):
        result[2 * block + 1, 2 * twist] += multiplicity
    return tuple((block, twist, multiplicity)
                 for (block, twist), multiplicity in sorted(result.items(), reverse=True))


@lru_cache(None)
def tensor_power_blocks(power):
    """Constructive decomposition of V^tensor power, valid over every F2 algebra."""
    if power < 0:
        raise ValueError("tensor power must be nonnegative")
    result = {(0, 0): 1}
    for degree in range(power):
        following = defaultdict(int)
        for (block, twist), multiplicity in result.items():
            for next_block, next_twist, next_multiplicity in tensor_with_natural(block):
                following[next_block, twist + next_twist] += multiplicity * next_multiplicity
        result = following
        assert all(block + 2 * twist == degree + 1 for block, twist in result)
    return tuple((block, twist, multiplicity)
                 for (block, twist), multiplicity in sorted(result.items(), reverse=True))


@lru_cache(None)
def central_block(highest):
    """Weighted subsets indexing the zero-weight subspace of B_m."""
    if highest & 1:
        return (), ()
    weights = frobenius_weights(highest)
    sums = [0] * (1 << len(weights))
    for mask in range(1, len(sums)):
        bit = mask & -mask
        sums[mask] = sums[mask ^ bit] + weights[bit.bit_length() - 1]
    masks = tuple(mask for mask, value in enumerate(sums) if 2 * value == highest)
    intersections = tuple(tuple(sums[a & b] for b in masks) for a in masks)
    return masks, intersections


def binomial_parity(top, bottom):
    if bottom < 0:
        return 0
    if top < 0:
        top = bottom - top - 1
    return int(bottom & ~top == 0)


@lru_cache(None)
def folded_gram_rank(doubles, half_singles, winding):
    """Rank of binom(Q+|A intersect B|-1,Q-i) on p-subsets of a 2p-set."""
    i, p, q = doubles, half_singles, winding
    if min(i, p, q) < 0 or q < i:
        return 0
    answer = 0
    for highest, twist, multiplicity in tensor_power_blocks(2 * p):
        order = q - i - twist
        if order < 0:
            continue
        _, intersections = central_block(highest)
        rows = [sum(binomial_parity(q - 1 + overlap, order) << k
                    for k, overlap in enumerate(row)) for row in intersections]
        answer += multiplicity * rank_two(rows)
    return answer


@lru_cache(None)
def johnson_binomial_rank(half_universe, top_offset, order, intersection_sign=1):
    """Rank of binom(top_offset +/- |A intersect B|,order) on p-subsets."""
    p = half_universe
    if p < 0 or intersection_sign not in (-1, 1):
        raise ValueError("nonnegative half-universe and sign +/-1 are required")
    if order < 0:
        return 0
    answer = 0
    for highest, twist, multiplicity in tensor_power_blocks(2 * p):
        degree = order - twist
        if degree < 0:
            continue
        _, intersections = central_block(highest)
        offset = top_offset - (twist if intersection_sign == -1 else 0)
        rows = [sum(binomial_parity(offset + intersection_sign * overlap, degree) << k
                    for k, overlap in enumerate(row)) for row in intersections]
        answer += multiplicity * rank_two(rows)
    return answer


@lru_cache(None)
def complete_component_rank(length, doubles, half_singles, signed_shift):
    """All source angles combined inside one folded-component Gram matrix."""
    i, p, v = doubles, half_singles, signed_shift
    k = i + p
    if min(i, p) < 0 or k > length:
        return 0
    order = length + v - k - i
    if order < 0 or order > length - i:
        return 0
    return johnson_binomial_rank(p, length - i, order, -1)


@lru_cache(maxsize=4)
def residue_occupancy_counts(modulus, length):
    """Exact counts by (double count, half-single count, residue energy)."""
    if modulus < 2 or modulus & (modulus - 1):
        raise ValueError("modulus must be a power of two")
    if not 0 <= length <= modulus // 2:
        raise ValueError("length must lie between zero and the number of low modes")
    counts = {(0, 0, 0): 1}
    for residue in range(1, modulus, 2):
        following = defaultdict(int)
        for (doubles, singles, energy), count in counts.items():
            following[doubles, singles, energy] += count
            if 2 * doubles + singles + 2 <= 2 * length:
                following[doubles + 1, singles, energy + 2 * residue] += count
            if 2 * doubles + singles + 1 <= 2 * length:
                following[doubles, singles + 1, energy + residue] += count
        counts = following
    return tuple((i, singles // 2, energy, count)
                 for (i, singles, energy), count in counts.items() if singles % 2 == 0)


def coaccessible_gram_dimensions(modulus, length, delta=0, maximum_index=None):
    """Gram dimensions at K=delta+M*t after imposing W=delta+M*k."""
    maximum_index = 2 * length if maximum_index is None else maximum_index
    result = [0] * (maximum_index + 1)
    for i, p, energy, count in residue_occupancy_counts(modulus, length):
        k = i + p
        if energy != delta + modulus * k:
            continue
        for index in range(k, maximum_index + 1):
            result[index] += count * folded_gram_rank(i, p, index - k)
    return tuple(result)


def unweighted_gram_dimensions(modes, length, maximum_index=None):
    """Uniform upper bound obtained by dropping the residue-energy constraint."""
    if not 0 <= length <= modes:
        raise ValueError("length must lie between zero and the number of modes")
    maximum_index = 2 * length if maximum_index is None else maximum_index
    result = [0] * (maximum_index + 1)
    for i in range(length + 1):
        for p in range(min(length - i, (modes - i) // 2) + 1):
            k = i + p
            count = comb(modes, i) * comb(modes - i, 2 * p)
            for index in range(k, maximum_index + 1):
                result[index] += count * folded_gram_rank(i, p, index - k)
    return tuple(result)


def complete_gram_channels(modulus, length, parameter_weight):
    """Direct complete-piece channels, after combining all source angles."""
    delta = 2 * parameter_weight - length * modulus
    if delta > 0:
        raise ValueError("the parameter shell must lie at or below central energy")
    shifts = range(-(2 * parameter_weight // modulus), length + 1)
    rows = {v: 0 for v in shifts}
    columns = {v: 0 for v in shifts}
    for i, p, energy, count in residue_occupancy_counts(modulus, length):
        k = i + p
        row_active = energy == delta + modulus * k
        column_active = energy == modulus * k
        if not row_active and not column_active:
            continue
        for v in shifts:
            if row_active:
                rows[v] += count * complete_component_rank(length, i, p, v)
            if column_active:
                columns[v] += count * complete_component_rank(length, i, p, -v)
    return tuple((v, rows[v], columns[v]) for v in shifts)


def unweighted_complete_gram_channels(modes, length):
    """Uniform direct bounds, dropping residue energy but retaining all angles."""
    if not 0 <= length <= modes:
        raise ValueError("length must lie between zero and the number of modes")
    rows = [0] * (2 * length + 1)
    for i in range(length + 1):
        for p in range(min(length - i, (modes - i) // 2) + 1):
            count = comb(modes, i) * comb(modes - i, 2 * p)
            for v in range(-length, length + 1):
                rows[v + length] += count * complete_component_rank(length, i, p, v)
    return tuple((v, rows[v + length], rows[length - v])
                 for v in range(-length, length + 1))
