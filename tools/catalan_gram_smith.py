#!/usr/bin/env python3
"""Exact truncated-series Smith profiles for full Catalan cup Gram matrices."""

from __future__ import annotations

from functools import lru_cache

from catalan_cup_crossing import cup_basis


def loop_count(first, second):
    """Count components in the union of two perfect matchings."""
    seen = set()
    loops = 0
    for start in range(len(first)):
        if start in seen:
            continue
        loops += 1
        seen.add(start)
        pending = [start]
        while pending:
            vertex = pending.pop()
            for following in (first[vertex], second[vertex]):
                if following not in seen:
                    seen.add(following)
                    pending.append(following)
    return loops


def rank_binary_rows(rows):
    """Rank packed row vectors over F2; each Python int is one row."""
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


@lru_cache(maxsize=8)
def gram_coefficients(length):
    """Packed coefficient matrices of M=(F^T F)/(q^r delta)."""
    matchings, _ = cup_basis(length)
    dimension = len(matchings)
    max_degree = length // 2 - 1
    coefficients = [[0] * dimension for _ in range(max_degree + 1)]
    for row, first in enumerate(matchings):
        for column, second in enumerate(matchings):
            degree = loop_count(first, second) - 1
            coefficients[degree][row] |= 1 << column
    return tuple(tuple(rows) for rows in coefficients)


def truncated_kernel_dimension(length, precision, coefficients=None):
    """Kernel dimension of M on F2[delta]/(delta^precision)."""
    assert length >= 2 and length % 2 == 0 and precision >= 1
    if coefficients is None:
        coefficients = gram_coefficients(length)
    dimension = len(coefficients[0])
    max_degree = len(coefficients) - 1
    rows = []
    for output_degree in range(precision):
        for row in range(dimension):
            packed = 0
            for input_degree in range(output_degree + 1):
                matrix_degree = output_degree - input_degree
                if matrix_degree <= max_degree:
                    packed |= coefficients[matrix_degree][row] << (
                        input_degree * dimension
                    )
            rows.append(packed)
    return precision * dimension - rank_binary_rows(rows)


def smith_valuation_profile(length, precision=8):
    """Return the finite, exact delta-adic Smith profile of the Gram matrix.

    For Smith exponents s_i of M, the truncated kernel has dimension
    sum_i min(k, s_i). Its first difference counts exponents >= k.
    The Gram matrix delta*M shifts every exponent up by one.
    """
    assert precision >= 2
    dimension = len(cup_basis(length)[0])
    coefficients = gram_coefficients(length)
    kernels = [0]
    for k in range(1, precision + 1):
        kernels.append(truncated_kernel_dimension(length, k, coefficients))
    at_least = [dimension]
    at_least.extend(kernels[k] - kernels[k - 1] for k in range(1, precision + 1))

    counts = {}
    for gram_exponent in range(1, precision + 1):
        count = at_least[gram_exponent - 1] - at_least[gram_exponent]
        if count:
            counts[gram_exponent] = count
    return {
        "length": length,
        "dimension": dimension,
        "precision": precision,
        "kernel_dimensions": kernels[1:],
        "gram_exponent_counts": counts,
        "unresolved_at_least": at_least[precision],
    }


def dyck_path_valuation_distribution(semilength):
    """Count return paths by the delta-valuation of their spin-path norm."""
    assert semilength >= 0
    current = {0: {0: 1}}
    for step in range(2 * semilength):
        following = {}
        for height, histogram in current.items():
            up_count = (step + height) // 2
            down_count = (step - height) // 2
            if up_count < semilength and height < semilength:
                target = following.setdefault(height + 1, {})
                for valuation, count in histogram.items():
                    target[valuation] = target.get(valuation, 0) + count
            if down_count < semilength and height > 0:
                even_index = height if height % 2 == 0 else height + 1
                local_valuation = (even_index & -even_index) - 1
                target = following.setdefault(height - 1, {})
                for valuation, count in histogram.items():
                    total = valuation + local_valuation
                    target[total] = target.get(total, 0) + count
        current = following
    return current.get(0, {})


def maximum_path_norm_delta_valuation(semilength):
    """Closed-form maximum over the dyadic record-height profiles."""
    assert semilength >= 1
    maximum = 0
    exponent = 1
    while (1 << exponent) - 1 <= semilength:
        height = (1 << exponent) - 1
        profile_sum = (exponent - 1) * (1 << exponent) + 1
        maximum = max(maximum, profile_sum + (semilength - height) * height)
        exponent += 1
    return maximum


def max_path_norm_delta_valuation_dp(semilength):
    """Independent max-plus walk recurrence for the path-norm valuation."""
    assert semilength >= 0
    current = {0: 0}
    for step in range(2 * semilength):
        following = {}
        for height, valuation in current.items():
            up_count = (step + height) // 2
            down_count = (step - height) // 2
            if up_count < semilength and height < semilength:
                following[height + 1] = max(
                    following.get(height + 1, -1), valuation
                )
            if down_count < semilength and height > 0:
                even_index = height if height % 2 == 0 else height + 1
                local_valuation = (even_index & -even_index) - 1
                following[height - 1] = max(
                    following.get(height - 1, -1), valuation + local_valuation
                )
        current = following
    return current.get(0, -1)
