"""Polynomial-size characteristic-two blocks for symmetric Boolean operators.

The module formula uses SL2 tiltings, not a characteristic-zero decomposition.
See research/sl2-coefficient-ranks.md for sources and the proof.
"""

from functools import lru_cache
from math import comb

from check_modular_crossing import rank_two
from check_modular_recombination import insert, represent
from modular_pfaffian import one_plus_z_power


@lru_cache(maxsize=None)
def tilting_weights(n):
    v = n + 1
    leading = v.bit_length() - 1
    return tuple(
        [1 << i for i in range(leading)]
        + [1 << i for i in range(leading) if v >> i & 1]
    )


@lru_cache(maxsize=None)
def tilting_character(n):
    result = {0: 1}
    for weight in tilting_weights(n):
        following = {}
        for exponent, count in result.items():
            for shifted in (exponent - weight, exponent + weight):
                following[shifted] = following.get(shifted, 0) + count
        result = following
    assert max(result) == n and result[n] == 1
    return result


def decompose_character(character):
    """The caller supplies the character of a known tilting module."""
    remaining, answer = dict(character), []
    for n in range(max(remaining), -1, -1):
        multiplicity = remaining.get(n, 0)
        if multiplicity:
            assert multiplicity > 0
            answer.append((n, multiplicity))
            for exponent, count in tilting_character(n).items():
                remaining[exponent] = remaining.get(exponent, 0) - multiplicity * count
                assert remaining[exponent] >= 0
    assert not any(remaining.values())
    return tuple(answer)


@lru_cache(maxsize=None)
def tensor_tiltings(m):
    character = {m - 2 * j: comb(m, j) for j in range(m + 1)}
    result = decompose_character(character)
    assert sum(mu * (1 << len(tilting_weights(n))) for n, mu in result) == 1 << m
    return result


@lru_cache(maxsize=None)
def reciprocal_character(t):
    result = {0: 1}
    for _ in range(t):
        following = {}
        for exponent, count in result.items():
            for shift, coefficient in ((-2, 3), (-1, 1), (0, 8), (1, 1), (2, 3)):
                following[exponent + shift] = (
                    following.get(exponent + shift, 0) + count * coefficient
                )
        result = following
    return result


@lru_cache(maxsize=None)
def reciprocal_tiltings(t):
    result = decompose_character(reciprocal_character(t))
    assert sum(mu * (1 << len(tilting_weights(n))) for n, mu in result) == 16**t
    return result


@lru_cache(maxsize=None)
def subset_weights(n):
    weights = tilting_weights(n)
    result = [0] * (1 << len(weights))
    for subset in range(1, len(result)):
        bit = subset & -subset
        result[subset] = result[subset ^ bit] + weights[bit.bit_length() - 1]
    return tuple(result)


@lru_cache(maxsize=4096)
def tilting_rank(n, code):
    size = 1 << len(tilting_weights(n))
    if code & 1:
        return size
    if not code:
        return 0
    weights, rows = subset_weights(n), [0] * size
    for row in range(size):
        free = (size - 1) ^ row
        subset = free
        while True:
            if code >> weights[subset] & 1:
                rows[row] |= 1 << (row | subset)
            if not subset:
                break
            subset = (subset - 1) & free
    return rank_two(rows)


def rank_on_tiltings(decomposition, code):
    return sum(
        mu * tilting_rank(n, code & ((1 << (n + 1)) - 1)) for n, mu in decomposition
    )


def boolean_rank(m, code):
    """Rank of sum code[k]*e_k(N_1,...,N_m) on the Boolean algebra."""
    assert m >= 0 and 0 <= code < 1 << (m + 1)
    return rank_on_tiltings(tensor_tiltings(m), code)


def reciprocal_rank(t, coefficient_mask):
    assert t >= 1 and 0 <= coefficient_mask < 1 << (2 * t + 1)
    code = sum(
        ((coefficient_mask & (one_plus_z_power(k) << (2 * t - k))).bit_count() % 2) << k
        for k in range(2 * t + 1)
    )
    return rank_on_tiltings(reciprocal_tiltings(t), code)


def reverse_bits(value, length):
    return sum(((value >> j) & 1) << (length - 1 - j) for j in range(length))


def reciprocal_optimum(t):
    dimension, basis, selected = 2 * t + 1, {}, []
    row_ranks = [reciprocal_rank(t, mask) for mask in range(1 << dimension)]
    target = (1 << dimension) - 1
    for cost, mask in sorted(
        (row_ranks[m] * row_ranks[reverse_bits(m, dimension)], m)
        for m in range(1, 1 << dimension)
    ):
        if insert(basis, mask, 1 << len(selected)):
            selected.append(mask)
        hit = represent(basis, target)
        if hit is not None:
            witness = [
                {
                    "coefficient_mask": hex(m),
                    "row_rank": row_ranks[m],
                    "column_rank": row_ranks[reverse_bits(m, dimension)],
                }
                for j, m in enumerate(selected)
                if hit >> j & 1
            ]
            return {
                "word_length": t,
                "binary_width": 2 * t,
                "exact_full_span_cutoff": cost,
                "witness": witness,
            }
    raise AssertionError("crossing missing from coefficient span")
