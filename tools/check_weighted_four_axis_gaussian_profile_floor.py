#!/usr/bin/env python3
"""Exact finite controls for the weighted four-axis Gaussian profile floor.

The universal statement is proved in the accompanying note. These controls
check explicit shared-minor indexing, polynomial degree and valuation bounds,
ordered-prefix comparisons, and finite correlated-type entropy inequalities.
All polynomial calculations keep the four coordinate variables independent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import random

from check_four_axis_gaussian_profile_floor import (
    MASK, ONE, ZERO, compositions, decode, encode, multiply,
)


def determinant(matrix):
    """Subset dynamic program for a determinant in characteristic two."""
    size = len(matrix)
    values = {0: ONE}
    for row in range(size):
        following = {}
        for mask, polynomial in values.items():
            for column in range(size):
                if mask >> column & 1 or not matrix[row][column]:
                    continue
                term = multiply(polynomial, matrix[row][column])
                next_mask = mask | (1 << column)
                following[next_mask] = following.get(next_mask, ZERO) ^ term
        values = {mask: value for mask, value in following.items() if value}
    return values.get((1 << size) - 1, ZERO)


def polynomial_kernel(heights, seed):
    rng = random.Random(seed)
    modes = len(heights[0])
    size = 16 * modes
    kernel = [[ZERO for _ in range(size)] for _ in range(size)]
    pools = {}
    for left, right in combinations(range(size), 2):
        v, i = divmod(left, modes)
        w, j = divmod(right, modes)
        endpoint = (v ^ w == 15) and i == j
        cap = heights[v][i] + heights[w][j]
        if cap < 0:
            assert not endpoint
            continue
        if not endpoint and (cap == 0 or rng.randrange(7)):
            continue
        if cap not in pools:
            pools[cap] = [encode(exponents) for degree in range(cap + 1)
                          for exponents in compositions(degree)]
        pool = pools[cap]
        count = (3 if len(pool) >= 3 and rng.randrange(4) == 0 else 1) if endpoint else 2
        polynomial = frozenset(rng.sample(pool, count))
        kernel[left][right] = kernel[right][left] = polynomial
    return kernel


def shared_minor_controls():
    settings = [
        ("one-mode-positive-heights", [[1 + (v & 1)] for v in range(16)], 107),
        ("one-mode-negative-heights", [[-1 if v & 1 else 2] for v in range(16)], 109),
        ("two-uniform-modes", [[1, 2] for _ in range(16)], 113),
        ("two-modes-zero-average", [[-1 if v & 1 else 1, 1 + ((v >> 1) & 1)]
                                    for v in range(16)], 127),
    ]
    determinant_count = transposed_count = negative_cases = nonhomogeneous_cases = 0
    records = []
    for name, heights, seed in settings:
        modes = len(heights[0])
        size = 16 * modes
        kernel = polynomial_kernel(heights, seed)
        average = [Fraction(sum(row[j] for row in heights), 16) for j in range(modes)]
        assert all(value >= 0 for value in average)
        negative = any(value < 0 for row in heights for value in row)
        negative_cases += negative
        nonhomogeneous = False
        for left, right in product(range(size), repeat=2):
            v, i = divmod(left, modes)
            w, j = divmod(right, modes)
            polynomial = kernel[left][right]
            assert polynomial == kernel[right][left]
            assert len(polynomial) % 2 == int((v ^ w == 15) and i == j)
            assert not polynomial or left != right
            degrees = {sum(decode(monomial)) for monomial in polynomial}
            assert all(degree <= heights[v][i] + heights[w][j] for degree in degrees)
            nonhomogeneous |= len(degrees) > 1
        assert nonhomogeneous
        nonhomogeneous_cases += nonhomogeneous
        crossblocks = {}
        for axis in range(4):
            rows = [index for index in range(size) if not (index // modes >> axis & 1)]
            columns = [index for index in range(size) if index // modes >> axis & 1]
            crossblocks[axis] = (rows, columns,
                                 [[kernel[left][right] for right in columns] for left in rows])
        subset_records = []
        for count in range(1, modes + 1):
            for chosen_modes in combinations(range(modes), count):
                budget = sum(heights[v][j] for v in range(16) for j in chosen_modes)
                upper = [None] * 4
                pair_records = []
                for a, b in combinations(range(4), 2):
                    orientations = []
                    for sigma in (0, 1):
                        rows = [v * modes + j for v in range(16) for j in chosen_modes
                                if not (v >> a & 1) and (v >> b & 1) == sigma]
                        columns = [(index // modes ^ 15) * modes + index % modes for index in rows]
                        minor = [[kernel[left][right] for right in columns] for left in rows]
                        endpoint = [[len(value) % 2 for value in row] for row in minor]
                        assert endpoint == [[int(i == j) for j in range(4 * count)]
                                            for i in range(4 * count)]
                        for axis in (a, b):
                            ambient_rows, ambient_columns, matrix = crossblocks[axis]
                            if axis == b and sigma:
                                extracted = [[matrix[ambient_rows.index(right)][ambient_columns.index(left)]
                                              for left in rows] for right in columns]
                                assert extracted == [list(row) for row in zip(*minor)]
                                transposed_count += 1
                            else:
                                extracted = [[matrix[ambient_rows.index(left)][ambient_columns.index(right)]
                                              for right in columns] for left in rows]
                                assert extracted == minor
                        degree_cap = sum(heights[index // modes][index % modes]
                                         for index in rows + columns)
                        assert 0 <= degree_cap <= MASK  # Prevent packed-exponent overflow.
                        polynomial = determinant(minor)
                        assert polynomial and len(polynomial) % 2 == 1
                        exponents = [decode(monomial) for monomial in polynomial]
                        actual_degree = max(map(sum, exponents))
                        assert actual_degree <= degree_cap
                        va = min(exponent[a] for exponent in exponents)
                        vb = min(exponent[b] for exponent in exponents)
                        assert va + vb <= actual_degree <= degree_cap
                        for axis, valuation in ((a, va), (b, vb)):
                            upper[axis] = valuation if upper[axis] is None else min(upper[axis], valuation)
                        orientations.append({"orientation": sigma, "rows": rows,
                                             "antipodal_columns": columns,
                                             "degree_cap": degree_cap,
                                             "actual_degree": actual_degree,
                                             "determinant_monomials": len(polynomial),
                                             "valuations": [va, vb]})
                        determinant_count += 1
                    assert sum(case["degree_cap"] for case in orientations) == budget
                    pair_records.append({"axes": [a, b], "orientations": orientations})
                assert all(value is not None for value in upper)
                assert sum(upper) <= budget
                subset_records.append({"modes": list(chosen_modes),
                                       "shared_minor_upper_bounds_for_Delta_4k": upper,
                                       "sum_upper_bounds": sum(upper), "budget": budget,
                                       "pairs": pair_records})
        records.append({"name": name, "mode_count": modes, "heights": heights,
                        "averaged_heights": [[value.numerator, value.denominator] for value in average],
                        "negative_individual_heights": negative,
                        "nonhomogeneous_entries": nonhomogeneous,
                        "mode_subsets": subset_records})
    assert negative_cases >= 2
    return {"kernels": len(settings), "negative_height_kernels": negative_cases,
            "nonhomogeneous_kernels": nonhomogeneous_cases,
            "exact_shared_determinants": determinant_count,
            "transposed_minor_extractions": transposed_count, "cases": records}


def prefix_controls():
    rng = random.Random(131)
    cases = products = 0
    settings = [(0, 1), (1, 3), (1, 2, 4), (0, 2, 5)]
    for heights in settings:
        reference = [height for height in heights for _ in range(4)]
        size = len(reference)
        candidates = [reference, [0] * size]
        for _ in range(600):
            total = rng.randrange(sum(reference) + 1)
            sequence = [0] * size
            for _ in range(total):
                sequence[rng.randrange(size)] += 1
            candidates.append(sorted(sequence))
        for sequence in candidates:
            if any(sum(sequence[:4 * k]) > sum(reference[:4 * k])
                   for k in range(1, len(heights) + 1)):
                continue
            assert all(sum(sequence[:j]) <= sum(reference[:j]) for j in range(1, size + 1))
            for t in (Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)):
                left = right = Fraction(1)
                for u, v in zip(sequence, reference):
                    left *= 1 + t**u
                    right *= 1 + t**v
                assert left >= right
                products += 1
            cases += 1
    assert cases > 1000
    return {"ordered_sequence_cases": cases, "exact_partition_product_comparisons": products,
            "reference_height_lists": settings}


def entropy_controls():
    cases = equalities = 0
    alphabet_cases = 0
    alphabets = [(0, 1, 2, 3), (0, 1, 2), (0, 3), (1, 2)]
    for heights in ((0, 2), (1, 3), (2, 2)):
        for alphabet in alphabets:
            for t in (Fraction(1, 2), Fraction(2, 3)):
                weights = [t**sum(heights[j] for j in range(2) if letter >> j & 1)
                           for letter in alphabet]
                partition = sum(weights)
                full_partition = (1 + t**heights[0]) * (1 + t**heights[1])
                assert partition <= full_partition
                for denominator in range(2, 9):
                    for counts in compositions(denominator, len(alphabet)):
                        # Exponentiating M*(log Z-H(pi)-E(pi)*log t) avoids logs.
                        left = partition**denominator
                        right = Fraction(denominator**denominator)
                        for count, weight in zip(counts, weights):
                            left *= count**count  # Python's 0**0 is the required factor 1.
                            right *= weight**count
                        assert left >= right
                        expected_equality = all(Fraction(count, denominator) == weight / partition
                                                for count, weight in zip(counts, weights))
                        assert (left == right) == expected_equality
                        equalities += expected_equality
                        cases += 1
                alphabet_cases += 1
    return {"alphabet_height_scale_cases": alphabet_cases,
            "exact_correlated_type_comparisons": cases, "equality_cases": equalities,
            "identity": "Z^M product(k_b^k_b) >= M^M product(w_b^k_b)"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    dependency = Path(__file__).with_name("check_four_axis_gaussian_profile_floor.py")
    result = {"status": "passed",
              "scope": "Exact finite controls; the universal weighted profile floor is proved in the note.",
              "shared_minors": shared_minor_controls(),
              "ordered_prefixes": prefix_controls(), "entropy": entropy_controls(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "polynomial_helpers_sha256": hashlib.sha256(dependency.read_bytes()).hexdigest()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"],
                      "kernels": result["shared_minors"]["kernels"],
                      "shared_determinants": result["shared_minors"]["exact_shared_determinants"],
                      "transposed_minors": result["shared_minors"]["transposed_minor_extractions"],
                      "ordered_prefix_cases": result["ordered_prefixes"]["ordered_sequence_cases"],
                      "partition_product_comparisons": result["ordered_prefixes"]["exact_partition_product_comparisons"],
                      "correlated_type_comparisons": result["entropy"]["exact_correlated_type_comparisons"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
