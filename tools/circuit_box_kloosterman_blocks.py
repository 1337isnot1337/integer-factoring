"""Exact interval-weighted dyadic Kloosterman frequency blocks.

The output is an integer. No complex approximation, cyclotomic field, or
enumeration of the frequencies inside one block is used.
"""

from __future__ import annotations

from math import gcd


def _floor_sum(count: int, modulus: int, slope: int, offset: int) -> int:
    """Sum floor((slope*i+offset)/modulus), including signed coefficients."""
    if count < 0 or modulus <= 0:
        raise ValueError("invalid floor-sum domain")
    answer = 0
    while True:
        quotient, slope = divmod(slope, modulus)
        answer += quotient * count * (count - 1) // 2
        quotient, offset = divmod(offset, modulus)
        answer += quotient * count
        upper = slope * count + offset
        if upper < modulus:
            return answer
        count, offset = divmod(upper, modulus)
        modulus, slope = slope, modulus


def _cyclic_pieces(start: int, length: int, modulus: int) -> list[tuple[int, int]]:
    if not length:
        return []
    start %= modulus
    end = start + length
    if end <= modulus:
        return [(start, end)]
    return [(start, modulus), (0, end - modulus)]


def affine_interval_pairs(
    first: tuple[int, int],
    second: tuple[int, int],
    slope: int,
    target: int,
    modulus: int,
) -> int:
    """Count (a,b) in two signed intervals with a+slope*b=target mod modulus."""
    left_a, right_a = first
    left_b, right_b = second
    if right_a < left_a or right_b < left_b or modulus <= 0:
        raise ValueError("invalid affine interval domain")
    if gcd(slope, modulus) != 1:
        raise ValueError("the affine slope must be a unit")
    qa, ra = divmod(right_a - left_a, modulus)
    qb, rb = divmod(right_b - left_b, modulus)
    answer = qa * qb * modulus + qa * rb + qb * ra
    for start_a, stop_a in _cyclic_pieces(left_a, ra, modulus):
        for start_b, stop_b in _cyclic_pieces(left_b, rb, modulus):
            count = stop_b - start_b
            coefficient = -slope
            constant = target - slope * start_b
            # Number of residues in [start_a,stop_a). The common floor
            # terms in the two prefix counts cancel.
            answer += _floor_sum(
                count, modulus, coefficient, constant + modulus - start_a
            ) - _floor_sum(
                count, modulus, coefficient, constant + modulus - stop_a
            )
    return answer


def odd_square_roots(value: int, exponent: int) -> list[int]:
    """All odd square roots modulo 2**exponent, without a residue scan."""
    if exponent < 0:
        raise ValueError("negative root precision")
    if exponent == 0:
        return [0]
    if not value & 1:
        return []
    roots = [1]
    modulus = 2
    for _ in range(1, exponent):
        next_modulus = 2 * modulus
        roots = [
            candidate
            for root in roots
            for candidate in (root, root + modulus)
            if (candidate * candidate - value) % next_modulus == 0
        ]
        modulus = next_modulus
        if not roots:
            break
    return roots


def stationary_nodes(number: int, slope: int, exponent: int) -> list[int]:
    """At most eight stationary representatives for K(u,u*slope*number)."""
    if exponent < 1 or not (number & slope & 1):
        raise ValueError("a positive precision and two odd parameters are required")
    low = exponent // 2
    high = exponent - low
    step = 1 << low
    return [
        candidate
        for root in odd_square_roots(number * slope, low)
        for lift in range(1 << (high - low))
        for candidate in [root + step * lift]
        if candidate & 1
    ]


def _condition_interval(
    interval: tuple[int, int], residue: int, modulus: int
) -> tuple[int, int]:
    left, right = interval
    return (
        (left - residue + modulus - 1) // modulus,
        (right - residue + modulus - 1) // modulus,
    )


def weighted_frequency_block(
    number: int,
    first: tuple[int, int],
    second: tuple[int, int],
    exponent: int,
    slope_residue: int,
    slope_precision: int | None = None,
) -> dict:
    """Sum all odd u and c=c0 mod 2**h of Ahat(u) Bhat(cu) K(u,cu*N).

    Here q=2**exponent, Ahat(u)=sum_{a in first} exp(-2*pi*i*u*a/q),
    and K(u,v)=sum_{x odd mod q} exp(2*pi*i*(u*x+v/x)/q).
    The block contains (q/2)*2**(exponent-h) frequency pairs. Its exact
    integer value is evaluated by at most sixteen affine pair counts.
    """
    if exponent < 1:
        raise ValueError("positive conductor exponent required")
    low = exponent // 2
    minimum = max(1, low)
    if slope_precision is None:
        slope_precision = minimum
    if not minimum <= slope_precision <= exponent:
        raise ValueError("slope precision must be between max(1,floor(m/2)) and m")
    if not (number & slope_residue & 1):
        raise ValueError("number and slope residue must be odd")
    if first[1] < first[0] or second[1] < second[0]:
        raise ValueError("interval endpoints are reversed")
    modulus = 1 << exponent
    slope_modulus = 1 << slope_precision
    slope = slope_residue % slope_modulus
    fibre_exponent = exponent - slope_precision
    fibre = 1 << fibre_exponent
    roots = stationary_nodes(number, slope, exponent)
    answer = 0
    active = 0
    for root in roots:
        inverse_image = number * pow(root, -1, modulus) % modulus
        row_residue = root % fibre
        column_residue = inverse_image % fibre
        rows = _condition_interval(first, row_residue, fibre)
        columns = _condition_interval(second, column_residue, fibre)
        if rows[0] == rows[1] or columns[0] == columns[1]:
            continue
        active += 1
        target = (
            root + slope * inverse_image - row_residue - slope * column_residue
        ) // fibre
        full = affine_interval_pairs(rows, columns, slope, target, slope_modulus)
        half = affine_interval_pairs(rows, columns, slope, target, slope_modulus // 2)
        answer += modulus * full - (modulus // 2) * half
    answer <<= low + fibre_exponent
    return {
        "value": answer,
        "conductor_exponent": exponent,
        "slope_precision": slope_precision,
        "slope_residue": slope,
        "stationary_nodes": roots,
        "active_nodes": active,
        "affine_pair_queries": 2 * active,
        "slope_count": fibre,
        "frequency_pair_count": (modulus // 2) * fibre,
    }


def count_inverse_rectangle(
    number: int,
    first: tuple[int, int],
    second: tuple[int, int],
    exponent: int,
) -> dict:
    """Complete reference count by shared frequency blocks.

    This intentionally retains the charged O(2**(m/2)*poly(m)) outer
    block cost. It is not a new all-input factoring bound.
    """
    if exponent < 1 or not number & 1:
        raise ValueError("positive precision and odd target required")
    if first[1] < first[0] or second[1] < second[0]:
        raise ValueError("interval endpoints are reversed")
    modulus = 1 << exponent
    rows = _condition_interval(first, 1, 2)
    columns = _condition_interval(second, 1, 2)
    numerator = 2 * modulus * (rows[1] - rows[0]) * (columns[1] - columns[0])
    blocks = nodes = queries = 0
    for conductor in range(2, exponent + 1):
        precision = max(1, conductor // 2)
        for slope in range(1, 1 << precision, 2):
            block = weighted_frequency_block(
                number, first, second, conductor, slope, precision
            )
            numerator += (1 << (exponent - conductor)) * block["value"]
            blocks += 1
            nodes += len(block["stationary_nodes"])
            queries += block["affine_pair_queries"]
    count, remainder = divmod(numerator, modulus * modulus)
    if remainder:
        raise AssertionError("Fourier normalization did not produce an integer")
    return {
        "count": count,
        "conductor_exponent": exponent,
        "block_count": blocks,
        "stationary_nodes_total": nodes,
        "affine_pair_queries": queries,
    }
