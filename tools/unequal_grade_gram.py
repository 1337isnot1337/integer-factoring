#!/usr/bin/env python3
"""Complete Gram channels for unequal operator and state exterior grades."""

from functools import lru_cache
from math import isqrt, log2

from check_modular_crossing import rank_two
from folded_gram_aggregate import (
    block_character, shell_dimension,
    subset_pair_energy_count,
)
from folded_gram_ranks import binomial_parity, frobenius_weights


@lru_cache(None)
def weight_block(highest, degree):
    """A possibly noncentral torus-weight block of the explicit B_m module."""
    if degree < 0 or degree > highest:
        return ()
    weights = frobenius_weights(highest)
    sums = [0] * (1 << len(weights))
    for mask in range(1, len(sums)):
        bit = mask & -mask
        sums[mask] = sums[mask ^ bit] + weights[bit.bit_length() - 1]
    masks = [mask for mask, value in enumerate(sums) if value == degree]
    return tuple(tuple(sums[a & b] for b in masks) for a in masks)


def decompose_character(character):
    total = len(character) - 1
    remaining = list(character)
    result = []
    for highest in range(total, -1, -2):
        twist = (total - highest) // 2
        multiplicity = remaining[twist]
        if not multiplicity:
            continue
        assert multiplicity > 0
        result.append((highest, twist, multiplicity))
        for degree, value in enumerate(block_character(highest)):
            remaining[twist + degree] -= multiplicity * value
        assert min(remaining) >= 0
    assert not any(remaining)
    return tuple(result)


@lru_cache(maxsize=4096)
def weighted_exterior_blocks(modes, total_degree, energy):
    lower_half = [subset_pair_energy_count(modes, r, total_degree - r, energy)
                  for r in range(total_degree // 2 + 1)]
    mirror = lower_half[-2::-1] if total_degree % 2 == 0 else lower_half[::-1]
    return decompose_character(lower_half + mirror)


@lru_cache(None)
def component_block_rank(operator_length, state_length, creators, highest, signed_shift):
    length, states, c, m, v = operator_length, state_length, creators, highest, signed_shift
    delta = states - length
    total = 2 * c + delta
    if total < 0 or m > total or (total - m) & 1:
        return 0
    twist = (total - m) // 2
    degree = c - twist
    if degree < 0 or degree > m:
        return 0
    order = length + v - c - twist
    offset = length - twist
    if order < 0 or order > offset:
        return 0
    intersections = weight_block(m, degree)
    rows = [sum(binomial_parity(offset - overlap, order) << j
                for j, overlap in enumerate(row)) for row in intersections]
    return rank_two(rows)


def complete_gram_channels(modulus, operator_length, state_length, operator_energy, state_energy):
    """Both complete cut bounds; return shifts in units of the modulus."""
    if modulus < 2 or modulus & (modulus - 1):
        raise ValueError("the modulus must be a power of two")
    modes = modulus // 2
    length, states = operator_length, state_length
    if not 0 <= length <= modes or not 0 <= states <= modes:
        raise ValueError("both grades must lie between zero and the number of modes")
    j, energy = operator_energy, state_energy
    h = length * modulus - j
    if j > h:
        raise ValueError("the operator energy must lie at or below its central shell")
    delta = states - length
    minimum_shift = max(-length, -((energy + j) // modulus))
    maximum_shift = min(length, (energy + h) // modulus)
    shifts = range(minimum_shift, maximum_shift + 1)
    rows = {v: 0 for v in shifts}
    columns = {v: 0 for v in shifts}
    for creators in range(max(0, -delta), length + 1):
        total_degree = 2 * creators + delta
        row_energy = energy + j - length * modulus + creators * modulus
        col_energy = energy - j + creators * modulus
        for highest, _, multiplicity in weighted_exterior_blocks(modes, total_degree, row_energy):
            for v in shifts:
                rows[v] += multiplicity * component_block_rank(
                    length, states, creators, highest, v)
        for highest, _, multiplicity in weighted_exterior_blocks(modes, total_degree, col_energy):
            for v in shifts:
                columns[v] += multiplicity * component_block_rank(
                    length, states, creators, highest, -v)
    return tuple((v, rows[v], columns[v]) for v in shifts)


def score_case(modulus, operator_length, state_length, operator_energy, state_energy):
    operators = shell_dimension(modulus, operator_length, operator_energy)
    states = shell_dimension(modulus, state_length, state_energy)
    dimension = min(operators, states)
    if dimension < 2:
        raise ValueError("both domains must have at least two labels")
    profile = complete_gram_channels(
        modulus, operator_length, state_length, operator_energy, state_energy)
    cap = dimension ** 2
    total = 0
    for _, r, c in profile:
        product = min(cap, r) * min(cap, c)
        root = isqrt(product)
        total += root + int(root * root < product)
    return {
        "mode_cutoff": modulus, "operator_length": operator_length,
        "state_length": state_length, "operator_energy": operator_energy,
        "state_energy": state_energy, "operator_labels": operators,
        "state_labels": states, "selected_labels": dimension,
        "integer_rank_sum_upper_bound": total,
        "rho_of_upper_bound": log2(total) / (2 * log2(dimension)),
        "complete_uncapped_channels": [list(piece) for piece in profile],
    }
