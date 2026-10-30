#!/usr/bin/env python3
"""Stable one-ancilla GL2xGL2 bounds for actual full-grade Pfaffian cuts."""

from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import comb, isqrt, log2

from check_homogeneous_central_symbols import binomial_two
from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from folded_gram_ranks import frobenius_weights, tensor_power_blocks


HG_INTERTWINER_COLUMNS = (
    4352, 8704, 17408, 34816, 4132, 66056, 1049728, 2392064,
    16, 74240, 4384, 131072, 64, 296960, 17536, 524288,
    1, 1065984, 2, 2131968, 4356, 4194304, 8712, 8388608,
)


@lru_cache(None)
def weighted_masks(highest, degree):
    if degree < 0 or degree > highest:
        return (), ()
    weights = frobenius_weights(highest)
    sums = [0] * (1 << len(weights))
    for mask in range(1, len(sums)):
        bit = mask & -mask
        sums[mask] = sums[mask ^ bit] + weights[bit.bit_length() - 1]
    return tuple(mask for mask, value in enumerate(sums) if value == degree), tuple(sums)


@lru_cache(None)
def power_polynomial(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


@lru_cache(None)
def stable_multiplicities(pairs):
    """G tensor (g+h)^tensor(2P), decomposed with constant binary bases."""
    counts = defaultdict(int)
    for chosen in range(2 * pairs + 1):
        outer = comb(2 * pairs, chosen)
        for first, first_twist, first_count in tensor_power_blocks(chosen + 1):
            for second, second_twist, second_count in tensor_power_blocks(2 * pairs - chosen + 1):
                twist = first_twist + second_twist
                assert first + second + 2 * twist == 2 * pairs + 2
                counts[first, second, twist] += outer * first_count * second_count
    return tuple((m, n, d, count) for (m, n, d), count in sorted(counts.items()))


@lru_cache(None)
def rectangular_block_profile(first, second, grade_difference, grade_offset):
    """Coefficient ranks before the common determinant twist t^(2d).

    Ancillary G[1,0]=1 has row grade (0,1) and column grade (1,0).
    Thus the first factor has adjacent row/column weight slices, while
    the second has one common slice.
    """
    if (grade_difference + 1 + first) & 1 or (grade_offset - 1 + second) & 1:
        return (), 0, 0
    row_weight = (grade_difference + 1 + first) // 2
    col_weight = row_weight - 1
    second_weight = (grade_offset - 1 + second) // 2
    row_first, sums_first = weighted_masks(first, row_weight)
    col_first, _ = weighted_masks(first, col_weight)
    second_masks, sums_second = weighted_masks(second, second_weight)
    row_dimension = len(row_first) * len(second_masks)
    col_dimension = len(col_first) * len(second_masks)
    if not row_dimension or not col_dimension:
        return (), row_dimension, col_dimension
    first_intersections = [[sums_first[a & b] for b in col_first] for a in row_first]
    second_intersections = [[sums_second[a & b] for b in second_masks] for a in second_masks]
    values = {}
    for w in {x for row in first_intersections for x in row}:
        for z in {x for row in second_intersections for x in row}:
            exponent = col_weight + second_weight - w - z
            assert exponent >= 0
            values[w, z] = multiply(binomial_two(exponent), power_polynomial(7, z)) << (first - row_weight)
    rows = {}
    for a, c in product(range(len(row_first)), range(len(second_masks))):
        row_index = a * len(second_masks) + c
        for b, e in product(range(len(col_first)), range(len(second_masks))):
            value = values[first_intersections[a][b], second_intersections[c][e]]
            col_index = b * len(second_masks) + e
            while value:
                bit = value & -value
                degree = bit.bit_length() - 1
                if degree not in rows:
                    rows[degree] = [0] * row_dimension
                rows[degree][row_index] ^= 1 << col_index
                value ^= bit
    return tuple(sorted((degree, rank_two(matrix)) for degree, matrix in rows.items())), row_dimension, col_dimension


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def stable_row_bounds(pairs, operator_grade, state_grade=None):
    states = operator_grade if state_grade is None else state_grade
    difference, offset = states - operator_grade, operator_grade + states - 2 * pairs
    rows = defaultdict(int)
    row_dimension = col_dimension = 0
    largest = 0
    active = 0
    for first, second, twist, multiplicity in stable_multiplicities(pairs):
        profile, rdim, cdim = rectangular_block_profile(first, second, difference, offset)
        row_dimension += multiplicity * rdim
        col_dimension += multiplicity * cdim
        largest = max(largest, rdim, cdim)
        if profile:
            active += 1
        for degree, rank in profile:
            rows[degree + 2 * twist] += multiplicity * rank
    # These are the exact dimensions of the augmented row/column weight spaces.
    odd_grades = ((1, 0), (0, 1), (2, 1), (1, 2))
    expected_row = sum(choose(2 * pairs, operator_grade - a) * choose(2 * pairs, states + 1 - b)
                       for a, b in odd_grades)
    expected_col = sum(choose(2 * pairs, operator_grade + 1 - a) * choose(2 * pairs, states - b)
                       for a, b in odd_grades)
    assert row_dimension == expected_row and col_dimension == expected_col
    assert all(0 <= degree <= 2 * operator_grade for degree in rows)
    return {
        "complement_pairs": pairs, "operator_grade": operator_grade, "state_grade": states,
        "row_bounds": sorted(rows.items()), "largest_explicit_block": largest,
        "active_GL2_product_blocks": active,
        "augmented_row_dimension": row_dimension, "augmented_column_dimension": col_dimension,
    }


def score_case(pairs, length, label_dimension=None):
    report = stable_row_bounds(pairs, length)
    dimension = comb(2 * pairs, length) if label_dimension is None else label_dimension
    cap = dimension ** 2
    rows = dict(report["row_bounds"])
    profile = []
    total = 0
    for degree in range(2 * length + 1):
        first, second = min(cap, rows.get(degree, 0)), min(cap, rows.get(2 * length - degree, 0))
        product_ = first * second
        root = isqrt(product_)
        total += root + int(root * root < product_)
        profile.append([degree, first, second])
    report.update({"labels": dimension, "profile": profile,
                   "integer_rank_sum_upper_bound": total,
                   "rho_of_upper_bound": log2(total) / (2 * log2(dimension))})
    return report


@lru_cache(None)
def rectangular_geometry(first, second, grade_difference, grade_offset, ancillary_row=1, ancillary_column=0):
    """Weight intersections for a fixed constant ancillary G entry."""
    grades = ((1, 0), (0, 1), (2, 1), (1, 2))
    ar, br = grades[ancillary_row]
    ac, bc = grades[ancillary_column]
    numerators = (grade_difference + br - ar + first,
                  grade_difference + bc - ac + first,
                  grade_offset + ar + br - 2 + second,
                  grade_offset + ac + bc - 2 + second)
    if any(value & 1 for value in numerators):
        return None
    ur, uc, vr, vc = (value // 2 for value in numerators)
    first_rows, fsums = weighted_masks(first, ur)
    first_cols, _ = weighted_masks(first, uc)
    second_rows, ssums = weighted_masks(second, vr)
    second_cols, _ = weighted_masks(second, vc)
    if not first_rows or not first_cols or not second_rows or not second_cols:
        return None
    return (ur, uc, vr, vc,
            tuple(tuple(fsums[a & b] for b in first_cols) for a in first_rows),
            tuple(tuple(ssums[a & b] for b in second_cols) for a in second_rows))


@lru_cache(None)
def centered_block_profile(first, second, grade_difference, grade_offset,
                           ancillary_row=1, ancillary_column=0):
    """Exact bivariate jet ranks at p=1+x, q=1+y, r=1, before x^d."""
    geometry = rectangular_geometry(first, second, grade_difference, grade_offset,
                                     ancillary_row, ancillary_column)
    if geometry is None:
        return (), 0
    _, uc, vr, _, first_intersections, second_intersections = geometry
    row_dimension = len(first_intersections) * len(second_intersections)
    col_dimension = len(first_intersections[0]) * len(second_intersections[0])
    stride = first + second + 3
    values = {}
    for w in {x for row in first_intersections for x in row}:
        for z in {x for row in second_intersections for x in row}:
            exponents = (uc - w, vr - z, z)
            assert min(exponents) >= 0
            value = 1
            # Frobenius factors multiply each binomial/trinomial in O(log P)
            # shifts; exponents stay below the injective monomial stride.
            for which, exponent in enumerate(exponents):
                while exponent:
                    bit = exponent & -exponent
                    if which == 0:
                        value ^= value << bit
                    elif which == 1:
                        value ^= value << (stride * bit)
                    else:
                        value ^= (value << bit) ^ (value << (stride * bit))
                    exponent ^= bit
            values[w, z] = value
    rows = {}
    for a, c in product(range(len(first_intersections)), range(len(second_intersections))):
        row_index = a * len(second_intersections) + c
        for b, e in product(range(len(first_intersections[0])), range(len(second_intersections[0]))):
            value = values[first_intersections[a][b], second_intersections[c][e]]
            col_index = b * len(second_intersections[0]) + e
            while value:
                bit = value & -value
                exponent = bit.bit_length() - 1
                key = exponent % stride, exponent // stride
                if key not in rows:
                    rows[key] = [0] * row_dimension
                rows[key][row_index] ^= 1 << col_index
                value ^= bit
    return tuple(sorted((a, b, rank_two(matrix)) for (a, b), matrix in rows.items())), max(row_dimension, col_dimension)


def centered_row_bounds(pairs, length, state_grade=None, ancillary_row=1, ancillary_column=0):
    states = length if state_grade is None else state_grade
    rows = defaultdict(int)
    largest = 0
    for first, second, twist, multiplicity in stable_multiplicities(pairs):
        profile, dimension = centered_block_profile(first, second, states - length,
                                                    length + states - 2 * pairs,
                                                    ancillary_row, ancillary_column)
        largest = max(largest, dimension)
        for a, b, rank in profile:
            rows[a + twist, b] += multiplicity * rank
    ancillary_parameter_grade = (1, 0, 2, 1)[ancillary_row]
    assert all(a + b <= length + ancillary_parameter_grade for a, b in rows)
    return rows, largest


def centered_score_case(pairs, length, label_dimension=None, all_ancillas=True):
    dimension = comb(2 * pairs, length) if label_dimension is None else label_dimension
    # All six entries are identically one after specializing r=1.
    ancillas = ((1, 0), (1, 1), (1, 2), (1, 3), (0, 0), (0, 2)) if all_ancillas else ((1, 0),)
    candidates = [centered_row_bounds(pairs, length, ancillary_row=a, ancillary_column=b)
                  for a, b in ancillas]
    bounds = {(a, b): min(min(dimension ** 2, rows.get((a, b), 0)) for rows, _ in candidates)
              for a in range(length + 1) for b in range(length + 1 - a)}
    total, profile = 0, []
    for a, b in sorted(bounds):
        first, second = bounds[a, b], bounds[b, a]
        product_ = first * second
        root = isqrt(product_)
        total += root + int(root * root < product_)
        if first or second:
            profile.append([a, b, first, second])
    return {"complement_pairs": pairs, "length": length, "labels": dimension,
            "ancillary_entries": [list(pair) for pair in ancillas],
            "largest_explicit_block": max(size for _, size in candidates),
            "complete_centered_jet_bounds": profile,
            "integer_rank_sum_upper_bound": total,
            "rho_of_upper_bound": log2(total) / (2 * log2(dimension))}
