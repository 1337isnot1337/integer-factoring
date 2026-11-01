#!/usr/bin/env python3
"""Universal exact local tensor for complementary fixed-length shells.

For four L-part low-mode parameter/state labels of common weight j,
with partner label b=M-v, this returns t**L times the complete crossing
profile, where t=z**M. Arithmetic is exact in F_2[t].
"""

from functools import lru_cache
from itertools import product


def multiply(first, second):
    result = 0
    while second:
        bit = second & -second
        result ^= first << (bit.bit_length() - 1)
        second ^= bit
    return result


@lru_cache(None)
def local_tensor(x, y, u, v):
    """A two-mode complementary pair, with occupation masks in 0..3."""
    vertices = ([(0, bit) for bit in range(2) if u >> bit & 1]
                + [(1, bit) for bit in range(2) if x >> bit & 1]
                + [(2, bit) for bit in range(2) if v >> bit & 1]
                + [(3, bit) for bit in range(2) if y >> bit & 1])
    size = len(vertices)
    if size & 1:
        return 0
    matrix = [[0] * size for _ in range(size)]
    for i, (first, a) in enumerate(vertices):
        for k in range(i + 1, size):
            second, b = vertices[k]
            if (first, second) == (0, 1) and a == b:
                matrix[i][k] = 0b110
            elif (first, second) == (1, 3) and a != b:
                matrix[i][k] = 0b11
            elif (first, second) == (1, 2) and a == b:
                matrix[i][k] = 0b10
            elif (first, second) == (0, 2) and a != b:
                matrix[i][k] = 1
            elif (first, second) == (2, 3) and a == b:
                matrix[i][k] = 1
            elif (first, second) == (0, 3) and a == b:
                matrix[i][k] = 1

    @lru_cache(None)
    def pfaffian(mask):
        if mask == 0:
            return 1
        first_bit = mask & -mask
        first = first_bit.bit_length() - 1
        remaining = mask ^ first_bit
        result = 0
        candidates = remaining
        while candidates:
            second_bit = candidates & -candidates
            second = second_bit.bit_length() - 1
            value = matrix[first][second]
            if value:
                result ^= multiply(value, pfaffian(remaining ^ second_bit))
            candidates ^= second_bit
        return result

    return pfaffian((1 << size) - 1)


LOCAL_TENSOR = tuple(local_tensor(x, y, u, v) for x, y, u, v in product(range(4), repeat=4))


def pair_occupations(parts, modulus):
    if modulus < 4 or modulus & (modulus - 1):
        raise ValueError("modulus must be a power of two at least four")
    result = [0] * (modulus // 4)
    for part in parts:
        if not 0 < part < modulus or part % 2 == 0:
            raise ValueError("all modes must be odd and below the modulus")
        representative = min(part, modulus - part)
        index = (representative - 1) // 2
        bit = int(part > modulus // 2)
        if result[index] >> bit & 1:
            raise ValueError("mode sets must have distinct parts")
        result[index] |= 1 << bit
    return tuple(result)


def folded_profile(x, y, u, v):
    """The arguments are complementary-pair occupation tuples."""
    if not len(x) == len(y) == len(u) == len(v):
        raise ValueError("all occupation tuples must have the same length")
    result = 1
    for a, b, c, d in zip(x, y, u, v):
        value = LOCAL_TENSOR[((a * 4 + b) * 4 + c) * 4 + d]
        if not value:
            return 0
        result = multiply(result, value)
    return result


def unfold_profile(value, modulus, length, parameter_weight):
    """Convert t**L F to z**(2j) F, matching the native Pfaffian."""
    shift = 2 * parameter_weight - length * modulus
    result = 0
    while value:
        bit = value & -value
        degree = (bit.bit_length() - 1) * modulus + shift
        if degree < 0:
            raise AssertionError("a forbidden negative profile degree survived")
        result |= 1 << degree
        value ^= bit
    return result
