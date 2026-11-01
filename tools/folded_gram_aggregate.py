#!/usr/bin/env python3
"""Fast exact aggregation of the unweighted complete folded-Gram channels.

The local residue module is 1 + V + det(V).  In exterior degree 2k, its
torus character is sum_r C(N,r) C(N,2k-r) X^r.  The elementary Frobenius
blocks from folded_gram_ranks.py therefore aggregate all double/single
residue patterns before binary elimination.
"""

from functools import lru_cache
from math import comb

from check_modular_crossing import rank_two
from folded_gram_ranks import binomial_parity, central_block, frobenius_weights


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def decompose_character(character):
    degree = len(character) - 1
    if degree & 1:
        raise ValueError("an even exterior degree is required")
    k = degree // 2
    remaining = list(character)
    result = []
    for highest in range(2 * k, -1, -2):
        twist = k - highest // 2
        multiplicity = remaining[twist]
        if not multiplicity:
            continue
        assert multiplicity > 0
        result.append((highest, twist, multiplicity))
        for exponent, value in enumerate(block_character(highest)):
            remaining[twist + exponent] -= multiplicity * value
        assert min(remaining) >= 0
    assert not any(remaining)
    return tuple(result)


@lru_cache(None)
def block_character(highest):
    result = [1]
    for weight in frobenius_weights(highest):
        following = result + [0] * weight
        for degree, value in enumerate(result):
            following[degree + weight] += value
        result = following
    return tuple(result)


@lru_cache(maxsize=128)
def exterior_blocks(modes, half_degree):
    """(highest weight, determinant twist, multiplicity) in exterior grade2k."""
    k = half_degree
    if not 0 <= k <= modes:
        raise ValueError("half-degree must lie between zero and the number of modes")
    character = [choose(modes, r) * choose(modes, 2 * k - r) for r in range(2 * k + 1)]
    result = decompose_character(character)
    assert sum(multiplicity * len(central_block(highest)[0])
               for highest, _, multiplicity in result) == comb(modes, k) ** 2
    return tuple(result)


@lru_cache(maxsize=2)
def gaussian_distributions(modes):
    """Coefficients of all Gaussian binomials [N choose r]_q."""
    polynomials = [[1]]
    cap = modes // 2
    for n in range(1, modes + 1):
        if len(polynomials) <= min(n, cap):
            polynomials.append([])
        for r in range(min(n, cap), 0, -1):
            degree = r * (n - r)
            current = polynomials[r] + [0] * (degree + 1 - len(polynomials[r]))
            offset = n - r
            for exponent, value in enumerate(polynomials[r - 1]):
                current[exponent + offset] += value
            polynomials[r] = current
    low = [tuple(poly) for poly in polynomials]
    return tuple(low[min(r, modes - r)] for r in range(modes + 1))


def subset_pair_energy_count(modes, first_size, second_size, energy):
    """[X^energy] E_first(X) E_second(X) for low odd-mode subsets."""
    if not 0 <= first_size <= modes or not 0 <= second_size <= modes:
        return 0
    remaining = energy - first_size ** 2 - second_size ** 2
    if remaining < 0 or remaining & 1:
        return 0
    degree = remaining // 2
    distributions = gaussian_distributions(modes)
    first, second = distributions[first_size], distributions[second_size]
    lower = max(0, degree - len(second) + 1)
    upper = min(degree, len(first) - 1)
    return sum(first[j] * second[degree - j] for j in range(lower, upper + 1))


@lru_cache(maxsize=1024)
def weighted_exterior_blocks(modes, half_degree, energy):
    k = half_degree
    lower_half = [subset_pair_energy_count(modes, r, 2 * k - r, energy)
                  for r in range(k + 1)]
    character = lower_half + lower_half[-2::-1]
    result = decompose_character(character)
    assert sum(multiplicity * len(central_block(highest)[0])
               for highest, _, multiplicity in result) == subset_pair_energy_count(
                   modes, k, k, energy)
    return result


def shell_dimension(modulus, length, weight):
    modes = modulus // 2
    remaining = weight - length ** 2
    if not 0 <= length <= modes or remaining < 0 or remaining & 1:
        return 0
    polynomial = gaussian_distributions(modes)[length]
    degree = remaining // 2
    return polynomial[degree] if degree < len(polynomial) else 0


@lru_cache(None)
def merged_block_rank(length, half_degree, highest, signed_shift):
    k, m, v = half_degree, highest, signed_shift
    if m & 1 or not 0 <= m <= 2 * k <= 2 * length:
        raise ValueError("an even highest weight at most 2k<=2L is required")
    twist = k - m // 2
    order = length + v - k - twist
    offset = length - twist
    if order < 0 or order > offset:
        return 0
    _, intersections = central_block(m)
    rows = [sum(binomial_parity(offset - overlap, order) << j
                for j, overlap in enumerate(row)) for row in intersections]
    return rank_two(rows)


def unweighted_complete_gram_channels(modes, length):
    if not 0 <= length <= modes:
        raise ValueError("length must lie between zero and the number of modes")
    rows = [0] * (2 * length + 1)
    for k in range(length + 1):
        for highest, _, multiplicity in exterior_blocks(modes, k):
            for shift in range(-length, length + 1):
                rows[shift + length] += multiplicity * merged_block_rank(
                    length, k, highest, shift)
    return tuple((shift, rows[shift + length], rows[length - shift])
                 for shift in range(-length, length + 1))


def complete_gram_channels(modulus, length, parameter_weight):
    """Exact residue-energy-aware complete channels using aggregated blocks."""
    if modulus < 2 or modulus & (modulus - 1):
        raise ValueError("modulus must be a power of two")
    modes = modulus // 2
    if not 0 <= length <= modes:
        raise ValueError("length must lie between zero and the number of modes")
    delta = 2 * parameter_weight - length * modulus
    if delta > 0:
        raise ValueError("the parameter shell must lie at or below central energy")
    shifts = range(-(2 * parameter_weight // modulus), length + 1)
    rows = {v: 0 for v in shifts}
    columns = {v: 0 for v in shifts}
    for k in range(length + 1):
        for highest, _, multiplicity in weighted_exterior_blocks(modes, k, k * modulus + delta):
            for v in shifts:
                rows[v] += multiplicity * merged_block_rank(length, k, highest, v)
        for highest, _, multiplicity in weighted_exterior_blocks(modes, k, k * modulus):
            for v in shifts:
                columns[v] += multiplicity * merged_block_rank(length, k, highest, -v)
    return tuple((v, rows[v], columns[v]) for v in shifts)
