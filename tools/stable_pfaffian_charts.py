#!/usr/bin/env python3
"""Complete bivariate coordinate-chart bounds for the stable Pfaffian."""

from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb, isqrt, log2

from check_modular_crossing import rank_two
from pfaffian_particle_hole import one_particle_matrix
from stable_pfaffian_bounds import rectangular_geometry, stable_multiplicities


def encode_linear(form, stride):
    return (form & 1) | (((form >> 1) & 1) << 1) | (((form >> 2) & 1) << stride)


def linear_power_product(factors, stride):
    value = 1
    for form, exponent in factors:
        while exponent and value:
            bit = exponent & -exponent
            result = value if form & 1 else 0
            if form & 2:
                result ^= value << bit
            if form & 4:
                result ^= value << (stride * bit)
            value = result
            exponent ^= bit
    return value


@lru_cache(None)
def chart_block_profile(first, second, difference, offset, p, q, r,
                        ancillary_row, ancillary_column, nonmonomial_twist=0):
    geometry = rectangular_geometry(first, second, difference, offset,
                                     ancillary_row, ancillary_column)
    if geometry is None:
        return (), 0
    ur, uc, vr, _, first_intersections, second_intersections = geometry
    row_dimension = len(first_intersections) * len(second_intersections)
    col_dimension = len(first_intersections[0]) * len(second_intersections[0])
    stride = first + second + nonmonomial_twist + 5
    values = {}
    for w in {value for row in first_intersections for value in row}:
        for z in {value for row in second_intersections for value in row}:
            exponents = (first - ur - uc + w, uc - w, vr - z, z)
            assert min(exponents) >= 0
            values[w, z] = linear_power_product(
                ((p ^ r, nonmonomial_twist), *zip((r, p, q, p ^ q ^ r), exponents)), stride)
    rows = {}
    for a, c in product(range(len(first_intersections)), range(len(second_intersections))):
        row_index = a * len(second_intersections) + c
        for b, e in product(range(len(first_intersections[0])), range(len(second_intersections[0]))):
            value = values[first_intersections[a][b], second_intersections[c][e]]
            col_index = b * len(second_intersections[0]) + e
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                key = degree % stride, degree // stride
                if key not in rows:
                    rows[key] = [0] * row_dimension
                rows[key][row_index] ^= 1 << col_index
                value ^= bit
    return tuple(sorted((a, b, rank_two(matrix)) for (a, b), matrix in rows.items())), max(row_dimension, col_dimension)


@lru_cache(None)
def chart_row_bounds(pairs, length, p, q, r, all_ancillas=True):
    if p.bit_count() & 1 or q.bit_count() & 1 or not r.bit_count() & 1:
        raise ValueError("the chart must evaluate to crossing at x=y=1")
    g = one_particle_matrix(*(encode_linear(form, 8) for form in (p, q, r)))
    ancillas = [(a, b, (value.bit_length() - 1) % 8, (value.bit_length() - 1) // 8)
                for a, row in enumerate(g) for b, value in enumerate(row)
                if value and value & (value - 1) == 0]
    if not all_ancillas:
        ancillas = [(1, 0, 0, 0)]
    delta = p ^ r
    monomial = delta in (1, 2, 4)
    candidates, largest = [], 0
    for ar, ac, ax, ay in ancillas:
        rows = defaultdict(int)
        for first, second, twist, multiplicity in stable_multiplicities(pairs):
            shift_x = twist if delta == 2 else 0
            shift_y = twist if delta == 4 else 0
            nonmonomial = 0 if monomial else twist
            profile, dimension = chart_block_profile(first, second, 0, 2 * length - 2 * pairs,
                                                      p, q, r, ar, ac, nonmonomial)
            largest = max(largest, dimension)
            for x, y, rank in profile:
                x, y = x + shift_x - ax, y + shift_y - ay
                if x >= 0 and y >= 0 and x + y <= length:
                    rows[x, y] += multiplicity * rank
        candidates.append(rows)
    return tuple((x, y, min(rows.get((x, y), 0) for rows in candidates))
                 for x in range(length + 1) for y in range(length + 1 - x)), tuple(ancillas), largest


def score_chart(pairs, length, p, q, r, label_dimension=None, all_ancillas=True):
    dimension = comb(2 * pairs, length) if label_dimension is None else label_dimension
    first, row_ancillas, rlargest = chart_row_bounds(pairs, length, p, q, r, all_ancillas)
    second, col_ancillas, clargest = chart_row_bounds(pairs, length, q, p, r, all_ancillas)
    total, profile = 0, []
    for (a, b, row), (ca, cb, col) in zip(first, second):
        assert a == ca and b == cb
        row, col = min(dimension ** 2, row), min(dimension ** 2, col)
        product_ = row * col
        root = isqrt(product_)
        total += root + int(root * root < product_)
        if row or col:
            profile.append([a, b, row, col])
    return {"complement_pairs": pairs, "length": length, "labels": dimension,
            "chart_p_q_r": [p, q, r], "profile": profile,
            "row_ancillary_monomials": [list(term) for term in row_ancillas],
            "column_ancillary_monomials": [list(term) for term in col_ancillas],
            "largest_explicit_block": max(rlargest, clargest),
            "integer_rank_sum_upper_bound": total,
            "rho_of_upper_bound": log2(total) / (2 * log2(dimension))}
