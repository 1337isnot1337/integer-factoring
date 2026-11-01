#!/usr/bin/env python3
"""Balanced pole-kernel padding and polynomial-cost normalized inverse entries.

Domains are lists of (first positive odd node, number of nodes) progressions
with spacing two. The numerical domain sizes are never enumerated here.
"""

from __future__ import annotations

from math import gcd

from interval_unit_jets import QuadraticJets, unit_progression_jet


def valuation_two(value: int) -> int:
    if not value:
        raise ValueError("Zero has no finite valuation")
    value = abs(value)
    return (value & -value).bit_length() - 1


def domain_size(domain: tuple[tuple[int, int], ...]) -> int:
    return sum(length for _, length in domain)


def balanced_padding(n: int, left: tuple[int, int], right: tuple[int, int]) -> tuple:
    """Pad the shorter nonempty odd progression without losing residue balance."""
    if n < 1 or n % 2 == 0:
        raise ValueError("N must be positive and odd")
    for first, count in (left, right):
        if first < 1 or first % 2 == 0 or count < 1 or first + 2 * (count - 1) > n:
            raise ValueError("Inputs must be nonempty odd progressions contained in [1,N]")
    count = max(left[1], right[1])
    modulus = 1 << ((2 * count - 1).bit_length())

    def pad(interval: tuple[int, int]) -> tuple[tuple[int, int], ...]:
        first, length = interval
        if length == count:
            return (interval,)
        tail = first + 2 * length
        multiple = max(0, (n - tail) // modulus + 1)
        start = tail + multiple * modulus
        assert start > n
        return (interval, (start, count - length))

    return pad(left), pad(right), modulus


def residue_count(domain: tuple[tuple[int, int], ...], residue: int, level: int) -> int:
    modulus = 1 << level
    residue %= modulus
    if residue % 2 == 0:
        return 0
    index_modulus = modulus // 2
    total = 0
    for first, length in domain:
        root = ((residue - first) // 2) % index_modulus
        if root < length:
            total += 1 + (length - 1 - root) // index_modulus
    return total


def inverse_entry_valuation(n: int, left: tuple, right: tuple, a: int, b: int, k: int) -> int:
    """Cofactor-energy valuation, using only O(k) residue capacity queries."""
    result = 0
    modulus = 2
    for level in range(1, k + 1):
        row_residue = a % modulus
        column_residue = n * pow(b, -1, modulus) % modulus
        if row_residue != column_residue:
            def imbalance(residue: int) -> int:
                return (residue_count(left, residue, level)
                        - residue_count(right, n * pow(residue, -1, modulus), level))
            result += 1 - imbalance(row_residue) + imbalance(column_residue)
        modulus *= 2
    return result


def affine_product(ring: QuadraticJets, first: int, step: int, length: int,
                   shift_bits: int | None, *, omit_zero: bool = False) -> tuple[int, tuple]:
    """Return valuation and unit of prod(first+step*i+2^shift_bits*w).

    If shift_bits is None the factors are rational integers. Their zero,
    when present, must be explicitly omitted for a derivative product.
    """
    if step == 0 or length < 0:
        raise ValueError("Require a nonzero progression step and nonnegative length")
    if length == 0:
        return 0, ring.constant(1)
    pieces = [(0, length)]
    exponent = 0
    product = ring.constant(1)
    if (-first) % step == 0 and 0 <= -first // step < length:
        root = -first // step
        pieces = [(0, root), (root + 1, length)]
        if not omit_zero:
            if shift_bits is None:
                raise ValueError("An integer zero factor must be omitted")
            exponent += shift_bits
            product = ring.constant(0, 1)
    common = valuation_two(gcd(first, step))
    base, increment = first >> common, step >> common
    shifted = None if shift_bits is None else shift_bits - common
    for lower, upper in pieces:
        if lower == upper:
            continue
        exponent += common * (upper - lower)
        if increment % 2 == 0:
            assert base % 2
            if shifted is not None:
                assert shifted >= 1
            unit = ring.constant(base, 0 if shifted is None else 1 << shifted)
            product = ring.multiply(product, unit_progression_jet(
                ring, (lower, upper), unit, ring.constant(increment)))
            continue
        bound = max(abs(base + increment * lower), abs(base + increment * (upper - 1))).bit_length()
        covered = 0
        for v in range(bound):
            stride = 1 << (v + 1)
            root = ((1 << v) - base) * pow(increment, -1, stride) % stride
            start = -((root - lower) // stride)
            stop = -((root - upper) // stride)
            if start == stop:
                continue
            odd_base = (base + increment * root) >> v
            if shifted is not None:
                assert shifted > v
            unit = ring.constant(odd_base, 0 if shifted is None else 1 << (shifted - v))
            factor = unit_progression_jet(ring, (start, stop), unit, ring.constant(2 * increment))
            product = ring.multiply(product, factor)
            exponent += v * (stop - start)
            covered += stop - start
        assert covered == upper - lower
    assert product[0][0] % 2 or product[0][1] % 2
    return exponent, product


def normalized_inverse_entry(n: int, left: tuple, right: tuple,
                             a: int, b: int, k: int, precision: int) -> dict:
    """Value and unit of (M_k^-1)[b,a] without enumerating either domain."""
    if n < 1 or n % 2 == 0 or precision < 1:
        raise ValueError("Require positive odd N and positive unit precision")
    for domain, point in ((left, a), (right, b)):
        previous = 0
        member = False
        for first, length in domain:
            if first <= previous or first % 2 == 0 or length < 1:
                raise ValueError("Require ordered disjoint positive odd progressions")
            previous = first + 2 * (length - 1)
            member = member or (first <= point <= previous and (point - first) % 2 == 0)
        if not member:
            raise ValueError("The requested row and column must lie in their domains")
    count = domain_size(left)
    if count != domain_size(right) or count < 1:
        raise ValueError("Domains must have equal positive size")
    maximum_row = left[-1][0] + 2 * (left[-1][1] - 1)
    maximum_column = right[-1][0] + 2 * (right[-1][1] - 1)
    bound = max(abs(n - left[0][0] * right[0][0]), abs(n - maximum_row * maximum_column))
    if k < max(1, bound.bit_length()):
        raise ValueError("Displacement must separate every nonzero product difference")
    ring = QuadraticJets(precision, 0)
    valuation = 0
    unit = ring.constant(1)

    def include(exponent_unit: tuple, sign: int) -> None:
        nonlocal valuation, unit
        exponent, factor = exponent_unit
        valuation += sign * exponent
        unit = ring.multiply(unit, factor if sign == 1 else ring.inverse(factor))

    for first, length in left:
        include(affine_product(ring, n - b * first, -2 * b, length, k), 1)
        include(affine_product(ring, a - first, -2, length, None, omit_zero=True), -1)
    for first, length in right:
        include(affine_product(ring, n - a * first, -2 * a, length, k), 1)
        include(affine_product(ring, b - first, -2, length, None, omit_zero=True), -1)
    include(affine_product(ring, n - a * b, 1, 1, k), -1)
    unit = ring.multiply(unit, ring.inverse(ring.power(ring.constant(n, 1 << k), count - 1)))
    valuation -= k
    unit = ring.multiply(unit, ring.inverse(ring.constant(0, 1)))
    assert valuation == inverse_entry_valuation(n, left, right, a, b, k)
    return {"valuation": valuation, "unit": list(unit[0]), "precision": precision}
