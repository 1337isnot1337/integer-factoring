#!/usr/bin/env python3
"""Exact full-grade Pfaffian cut ranks from polynomially many local types."""

from functools import lru_cache
from itertools import product
from math import comb

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from pfaffian_particle_hole import local_exterior_table, one_particle_matrix


G_GRADES = ((1, 0), (0, 1), (2, 1), (1, 2))
H_GRADES = ((0, 0), (2, 0), (1, 1), (1, 1), (0, 2), (2, 2))
H_STATES = ((0, 0), (3, 0), (1, 2), (2, 1), (0, 3), (3, 3))


@lru_cache(None)
def type_basis(odd_sites, middle_sites, operator_grade, state_grade):
    grades = (G_GRADES,) * odd_sites + (H_GRADES,) * middle_sites
    @lru_cache(None)
    def recurse(position, first, second):
        remaining = len(grades) - position
        if first < 0 or second < 0 or first > 2 * remaining or second > 2 * remaining:
            return ()
        if not remaining:
            return ((),) if first == second == 0 else ()
        return tuple((index,) + suffix for index, (a, b) in enumerate(grades[position])
                     for suffix in recurse(position + 1, first - a, second - b))
    return recurse(0, operator_grade, state_grade)


def type_row_ranks(scalars, odd_sites, middle_sites, operator_grade, state_grade,
                   p=6, q=3, r=2):
    basis = type_basis(odd_sites, middle_sites, operator_grade - scalars, state_grade - scalars)
    if not basis:
        return {}, 0
    factor = 1
    for _ in range(scalars):
        factor = multiply(factor, p ^ r)
    if not factor:
        return {}, len(basis)
    g = one_particle_matrix(p, q, r)
    table = local_exterior_table(p, q, r)
    h = tuple(tuple(table[((x * 4 + y) * 4 + u) * 4 + v]
                    for v, u in H_STATES) for x, y in H_STATES)
    matrices = (g,) * odd_sites + (h,) * middle_sites
    cache, rows = {}, {}
    for row, first in enumerate(basis):
        for col, second in enumerate(basis):
            entries = tuple(sorted(matrix[a][b] for matrix, a, b in zip(matrices, first, second)))
            if entries and entries[0] == 0:
                continue
            value = cache.get(entries)
            if value is None:
                value = factor
                for entry in entries:
                    value = multiply(value, entry)
                cache[entries] = value
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                if degree not in rows:
                    rows[degree] = [0] * len(basis)
                rows[degree][row] ^= 1 << col
                value ^= bit
    return {degree: rank_two(matrix) for degree, matrix in rows.items()}, len(basis)


def full_grade_row_profile(pairs, operator_grade, state_grade=None, p=6, q=3, r=2):
    state_grade = operator_grade if state_grade is None else state_grade
    ranks, types = {}, []
    dimension_sum = 0
    for scalars in range(pairs + 1):
        for odd_sites in range(pairs - scalars + 1):
            middle_sites = pairs - scalars - odd_sites
            multiplicity = comb(pairs, scalars) * comb(pairs - scalars, odd_sites) * 2 ** (scalars + odd_sites)
            local, dimension = type_row_ranks(scalars, odd_sites, middle_sites,
                                             operator_grade, state_grade, p, q, r)
            dimension_sum += multiplicity * dimension
            if dimension:
                types.append({"scalar_sites": scalars, "odd_sites": odd_sites,
                              "middle_sites": middle_sites, "multiplicity": multiplicity,
                              "matrix_dimension": dimension,
                              "coefficient_ranks": sorted(local.items())})
            for degree, rank in local.items():
                ranks[degree] = ranks.get(degree, 0) + multiplicity * rank
    expected = comb(2 * pairs, operator_grade) * comb(2 * pairs, state_grade)
    assert dimension_sum == expected
    return {"complement_pairs": pairs, "operator_grade": operator_grade,
            "state_grade": state_grade, "matrix_dimension": expected,
            "row_ranks": sorted(ranks.items()), "types": types,
            "largest_type_matrix": max((row["matrix_dimension"] for row in types), default=0),
            "type_entries_evaluated": sum(row["matrix_dimension"] ** 2 for row in types)}


def full_grade_source_profile(pairs, operator_grade, state_grade=None):
    report = full_grade_row_profile(pairs, operator_grade, state_grade)
    rows = dict(report["row_ranks"])
    # On full state grades, mode complementation preserves the domain.
    # Swapping p,q in the homogeneous degree-L source reverses degree 2L.
    report["profile"] = [[degree, rows.get(degree, 0), rows.get(2 * operator_grade - degree, 0)]
                         for degree in range(2 * operator_grade + 1)
                         if rows.get(degree, 0) or rows.get(2 * operator_grade - degree, 0)]
    return report
