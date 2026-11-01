#!/usr/bin/env python3
"""Complete coefficient bounds for endpoint-preserving stable Pfaffian pencils."""

from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb, isqrt, log2

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from pfaffian_particle_hole import one_particle_matrix
from stable_pfaffian_bounds import (
    power_polynomial, rectangular_geometry, stable_multiplicities,
)


def monomial_degree(value):
    return value.bit_length() - 1 if value and value & (value - 1) == 0 else None


@lru_cache(None)
def pencil_block_profile(first, second, difference, offset, p, q, r,
                         ancillary_row, ancillary_column, nonmonomial_twist=0):
    geometry = rectangular_geometry(first, second, difference, offset,
                                     ancillary_row, ancillary_column)
    if geometry is None:
        return (), 0
    ur, uc, vr, _, first_intersections, second_intersections = geometry
    row_dimension = len(first_intersections) * len(second_intersections)
    col_dimension = len(first_intersections[0]) * len(second_intersections[0])
    factor = power_polynomial(p ^ r, nonmonomial_twist)
    values = {}
    for w in {value for row in first_intersections for value in row}:
        for z in {value for row in second_intersections for value in row}:
            exponents = (first - ur - uc + w, uc - w, vr - z, z)
            assert min(exponents) >= 0
            value = factor
            for polynomial, exponent in zip((r, p, q, p ^ q ^ r), exponents):
                value = multiply(value, power_polynomial(polynomial, exponent))
            values[w, z] = value
    rows = {}
    for a, c in product(range(len(first_intersections)), range(len(second_intersections))):
        row_index = a * len(second_intersections) + c
        for b, e in product(range(len(first_intersections[0])), range(len(second_intersections[0]))):
            value = values[first_intersections[a][b], second_intersections[c][e]]
            col_index = b * len(second_intersections[0]) + e
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                if degree not in rows:
                    rows[degree] = [0] * row_dimension
                rows[degree][row_index] ^= 1 << col_index
                value ^= bit
    return tuple(sorted((degree, rank_two(matrix)) for degree, matrix in rows.items())), max(row_dimension, col_dimension)


@lru_cache(None)
def pencil_row_bounds(pairs, length, p, q, r, all_ancillas=True):
    if p.bit_count() & 1 or q.bit_count() & 1 or not r.bit_count() & 1:
        raise ValueError("the pencil must evaluate to (p,q,r)=(0,0,1) at t=1")
    g = one_particle_matrix(p, q, r)
    ancillas = [(a, b, monomial_degree(value)) for a, row in enumerate(g)
                for b, value in enumerate(row) if monomial_degree(value) is not None]
    if not all_ancillas:
        ancillas = [(1, 0, 0)]
    delta_degree = monomial_degree(p ^ r)
    maximum_degree = length * max(p.bit_length() - 1, q.bit_length() - 1, r.bit_length() - 1)
    candidates, largest = [], 0
    for ar, ac, ancillary_degree in ancillas:
        rows = defaultdict(int)
        for first, second, twist, multiplicity in stable_multiplicities(pairs):
            nonmonomial = twist if delta_degree is None else 0
            shift = 0 if delta_degree is None else delta_degree * twist
            profile, dimension = pencil_block_profile(first, second, 0, 2 * length - 2 * pairs,
                                                       p, q, r, ar, ac, nonmonomial)
            largest = max(largest, dimension)
            for degree, rank in profile:
                original_degree = degree + shift - ancillary_degree
                if 0 <= original_degree <= maximum_degree:
                    rows[original_degree] += multiplicity * rank
        candidates.append(rows)
    return tuple((degree, min(rows.get(degree, 0) for rows in candidates))
                 for degree in range(maximum_degree + 1)), tuple(ancillas), largest


def score_pencil(pairs, length, p, q, r, label_dimension=None, all_ancillas=True):
    dimension = comb(2 * pairs, length) if label_dimension is None else label_dimension
    first, first_ancillas, largest_first = pencil_row_bounds(pairs, length, p, q, r, all_ancillas)
    second, second_ancillas, largest_second = pencil_row_bounds(pairs, length, q, p, r, all_ancillas)
    assert [degree for degree, _ in first] == [degree for degree, _ in second]
    cap, total, profile = dimension ** 2, 0, []
    for (degree, row), (_, col) in zip(first, second):
        row, col = min(cap, row), min(cap, col)
        product_ = row * col
        root = isqrt(product_)
        total += root + int(root * root < product_)
        if row or col:
            profile.append([degree, row, col])
    return {
        "complement_pairs": pairs, "length": length, "labels": dimension,
        "pencil_p_q_r": [p, q, r], "profile": profile,
        "row_ancillary_monomials": [list(term) for term in first_ancillas],
        "column_ancillary_monomials": [list(term) for term in second_ancillas],
        "largest_explicit_block": max(largest_first, largest_second),
        "integer_rank_sum_upper_bound": total,
        "rho_of_upper_bound": log2(total) / (2 * log2(dimension)),
    }
