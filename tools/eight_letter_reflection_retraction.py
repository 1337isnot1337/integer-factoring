#!/usr/bin/env python3
"""An explicit nonorthogonal eight-letter restriction of the M=8 reflection."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_block_smith_crossings import reverse_polynomial, smith_profile
from check_cancellation_module_automata import pmul
from reflection_dual_rail_gaussian import gauged_kernel


AMBIENT_HEIGHTS = tuple(sum((7 - 2 * bit) for bit in range(4)
                          if word >> bit & 1) for word in range(16))
CODE = (0, 1, 2, 4, 9, 8, 10, 12)
HEIGHTS = tuple(AMBIENT_HEIGHTS[word] for word in CODE)


def halve_exponents(value):
    """Remove the common ramification, rejecting any odd exponent."""
    if any(value >> exponent & 1 for exponent in range(1, value.bit_length(), 2)):
        raise ArithmeticError("The reflection kernel has an odd exponent.")
    return sum(1 << (exponent // 2)
               for exponent in range(0, value.bit_length(), 2)
               if value >> exponent & 1)


def shift_polynomial(value, shift):
    if shift >= 0:
        return value << shift
    if value & ((1 << -shift) - 1):
        raise ArithmeticError("The restored tensor has a Laurent pole.")
    return value >> -shift


@lru_cache(None)
def ambient_normalized_tensor():
    """Complete normalized zero tensor, with all 65,536 entries checked."""
    zero = tuple(tuple(halve_exponents(value) for value in row)
                 for row in gauged_kernel(8, "zero"))
    infinity = tuple(tuple(halve_exponents(value) for value in row)
                     for row in gauged_kernel(8, "infinity"))

    def permute(index):
        leg, mode = divmod(index, 4)
        return (0, 2, 1, 3)[leg] * 4 + mode

    assert all(infinity[i][j] == zero[permute(i)][permute(j)]
               for i, j in product(range(16), repeat=2))

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        i = first.bit_length() - 1
        rest, remaining, result = mask ^ first, mask ^ first, 0
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length() - 1
            if zero[i][j]:
                result ^= pmul(zero[i][j], pfaffian(rest ^ chosen))
        return result

    values = tuple(pfaffian(mask) for mask in range(1 << 16))
    for mask, value in enumerate(values):
        x, y, u, v = ((mask >> (4 * leg)) & 15 for leg in range(4))
        assert value.bit_count() % 2 == int(x == v and y == u)
    return values


def displayed_maps(nonorthogonal=True):
    embedding = tuple(1 << word for word in CODE)
    retraction = list(embedding)
    if nonorthogonal:
        retraction[4] ^= 1 << 6
    return embedding, tuple(retraction), HEIGHTS


def homogeneous_restriction(embedding, retraction, heights, keep_matrices=False):
    """Construct both cuts independently for homogeneous binary LP=I maps.

    Integer masks encode embedding columns and retraction rows. The return
    profiles include neither of the cancelling common shifts +16 and -16.
    """
    embedding, retraction, heights = map(tuple, (embedding, retraction, heights))
    q = len(embedding)
    if not q or len(retraction) != q or len(heights) != q:
        raise ValueError("Embedding, retraction, and heights must have equal lengths.")
    for maps in (embedding, retraction):
        for mask, height in zip(maps, heights):
            if not 0 < mask < 1 << 16:
                raise ValueError("Every map must be a nonzero ambient binary vector.")
            if any(AMBIENT_HEIGHTS[word] != height
                   for word in range(16) if mask >> word & 1):
                raise ValueError("The physical maps must preserve height.")
    if not all((retraction[i] & embedding[j]).bit_count() % 2 == int(i == j)
               for i, j in product(range(q), repeat=2)):
        raise ValueError("The retraction identity LP=I failed.")
    support = lambda mask: tuple(word for word in range(16) if mask >> word & 1)
    ps, ls = tuple(map(support, embedding)), tuple(map(support, retraction))
    ambient = ambient_normalized_tensor()
    zero = [[0] * (q * q) for _ in range(q * q)]
    infinity = [[0] * (q * q) for _ in range(q * q)]
    raw_entries, normalized_entries = [], []
    for x, y, u, v in product(range(q), repeat=4):
        at_zero, at_infinity = 0, 0
        for a, b, c, d in product(ls[x], ls[y], ps[u], ps[v]):
            at_zero ^= ambient[a | b << 4 | c << 8 | d << 12]
            # The independent natural infinity kernel swaps physical y and u.
            at_infinity ^= ambient[a | c << 4 | b << 8 | d << 12]
        zero[q * x + y][q * u + v] = at_zero
        infinity[q * x + u][q * y + v] = at_infinity
        raw = shift_polynomial(at_zero, 16 - heights[x] - heights[y])
        assert raw.bit_count() % 2 == int(x == v and y == u)
        raw_entries.append(raw)
        normalized_entries.append((u, v, at_infinity))
    degree = max(value.bit_length() - 1 for value in raw_entries)
    for raw, (u, v, normalized) in zip(raw_entries, normalized_entries):
        restored = shift_polynomial(reverse_polynomial(raw, degree),
                                    heights[u] + heights[v] + 16 - degree)
        assert restored == normalized
    profiles, precisions = [], []
    for matrix in (zero, infinity):
        matrix_degree = max(value.bit_length() - 1 for row in matrix for value in row)
        orders, precision = smith_profile(matrix, matrix_degree)
        assert len(orders) == q * q
        profiles.append(tuple(orders))
        precisions.append(precision)
    result = {
        "embedding": embedding, "retraction": retraction,
        "heights": heights, "type_energies": tuple(4 * height for height in heights),
        "actual_polynomial_degree": degree,
        "zero_orders": profiles[0], "infinity_orders": profiles[1],
        "smith_precisions": precisions,
        "literal_endpoint_entries": q**4,
        "literal_infinity_restoration_entries": q**4,
    }
    if keep_matrices:
        result["matrices"] = (zero, infinity)
    return result
