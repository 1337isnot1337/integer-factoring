#!/usr/bin/env python3
"""Exact dyadic affine equivalence of one-species typed Schur tensors.

All arithmetic is in F2[w]. No floating pressure calculation is used.
This verifies finite instances of the accompanying all-size identity;
it does not claim an improved integer-factoring exponent.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank
from check_schur_downset_crossing import binomial_polynomial, sparse_kernel
from check_schur_gaussian_reduction import polynomial_determinant, polynomial_multiply
from weighted_smith_oracle import WeightedSmithOracle


ROOT = Path(__file__).resolve().parent.parent


def substitute_power(value, scale):
    result = 0
    while value:
        bit = value & -value
        result ^= 1 << (scale * (bit.bit_length() - 1))
        value ^= bit
    return result


def polynomial_power(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_multiply(result, value)
        value = polynomial_multiply(value, value)
        exponent >>= 1
    return result


def one_species_kernel(indices, ambient):
    """Construct (3) directly, independently of sparse_kernel extraction."""
    indices = tuple(indices)
    n = len(indices)
    matrix = [[0] * (2 * n) for _ in range(2 * n)]
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            if a & ~b == 0:
                q = binomial_polynomial(a, 2)
                matrix[i][j] = q ^ (q << 1)
            if b & ~a == 0:
                matrix[i][n + j] = 1 << (2 * ambient)
            if i == j:
                matrix[n + i][j] = 1
            h = 1 << (2 * (ambient - 1 - (a | b)))
            matrix[n + i][n + j] = h ^ (h << 1)
    return matrix


def selected_minor(matrix, rows, columns):
    if len(rows) != len(columns):
        return 0
    return polynomial_determinant([[matrix[i][j] for j in columns] for i in rows])


def check_blocks(indices, ambient, scale):
    enlarged = tuple(scale * a + scale - 1 for a in indices)
    n = len(indices)
    small = one_species_kernel(indices, ambient)
    large = one_species_kernel(enlarged, scale * ambient)
    changed = [[substitute_power(v, scale) for v in row] for row in small]
    lam = binomial_polynomial(scale - 1)
    full = sparse_kernel(enlarged, scale * ambient)
    for i in range(2 * n):
        for j in range(2 * n):
            assert large[i][j] == full[(2 + i // n) * n + i % n][j]
            left = polynomial_multiply(lam if j >= n else 1, large[i][j])
            right = polynomial_multiply(lam if i < n else 1, changed[i][j])
            assert left == right
    return small, large, changed, lam


def check_all_minors(indices, ambient, scale):
    small, large, changed, lam = check_blocks(indices, ambient, scale)
    n = len(indices)
    powers = [polynomial_power(lam, r) for r in range(n + 1)]
    count = 0
    unequal_charge = 0
    nontrivial_factor = 0
    for size in range(2 * n + 1):
        for rows in combinations(range(2 * n), size):
            rv = sum(i < n for i in rows)
            for columns in combinations(range(2 * n), size):
                rx = sum(j >= n for j in columns)
                actual = selected_minor(large, rows, columns)
                expected = selected_minor(changed, rows, columns)
                assert polynomial_multiply(powers[rx], actual) == polynomial_multiply(powers[rv], expected)
                count += 1
                unequal_charge += int(rv != rx)
                nontrivial_factor += int(actual != expected)
    assert nontrivial_factor > 0
    return {"indices": indices, "ambient": ambient, "scale": scale,
            "all_minors_checked": count, "unequal_charge_minors": unequal_charge,
            "minors_retaining_nontrivial_boundary_factor": nontrivial_factor}


def local_entries(matrix, n):
    result = {}
    for x, y, u, v in product(range(1 << n), repeat=4):
        rows = [i for i in range(n) if v >> i & 1]
        rows += [n + i for i in range(n) if y >> i & 1]
        columns = [i for i in range(n) if u >> i & 1]
        columns += [n + i for i in range(n) if x >> i & 1]
        result[x, y, u, v] = selected_minor(matrix, rows, columns)
    return result


def type_words(length, counts):
    labels = []
    choices = [tuple(combinations(range(length), k)) for k in counts]
    for occupied in product(*choices):
        word = [0] * length
        for a, positions in enumerate(occupied):
            for j in positions:
                word[j] |= 1 << a
        labels.append(tuple(word))
    return labels


def coefficient_ranks(entries, dimension):
    maximum = max(entries.values(), default=0).bit_length()
    result = {}
    for degree in range(maximum):
        row_cut = [0] * (dimension * dimension)
        column_cut = [0] * (dimension * dimension)
        for (x, y, u, v), value in entries.items():
            if value >> degree & 1:
                row_cut[x * dimension + y] |= 1 << (u * dimension + v)
                column_cut[x * dimension + u] |= 1 << (y * dimension + v)
        if any(row_cut):
            result[degree] = (binary_rank(row_cut), binary_rank(column_cut))
    return result


def check_type(indices, ambient, scale, length, counts):
    small, large, _, _ = check_blocks(indices, ambient, scale)
    n = len(indices)
    small_local = local_entries(small, n)
    large_local = local_entries(large, n)
    labels = type_words(length, counts)
    entries_small, entries_large = {}, {}
    for xi, yi, ui, vi in product(range(len(labels)), repeat=4):
        x, y, u, v = (labels[i] for i in (xi, yi, ui, vi))
        p = q = 1
        for states in zip(x, y, u, v):
            p = polynomial_multiply(p, small_local[states])
            q = polynomial_multiply(q, large_local[states])
        assert q == substitute_power(p, scale)
        assert p.bit_count() % 2 == int(xi == vi and yi == ui)
        entries_small[xi, yi, ui, vi] = p
        entries_large[xi, yi, ui, vi] = q
    small_ranks = coefficient_ranks(entries_small, len(labels))
    large_ranks = coefficient_ranks(entries_large, len(labels))
    assert large_ranks == {scale * j: ranks for j, ranks in small_ranks.items()}
    varying = sum(len({mask.bit_count() for mask in word}) > 1 for word in labels)
    assert varying > 0
    return {"indices": indices, "ambient": ambient, "scale": scale,
            "copies": length, "counts": counts, "labels": len(labels),
            "typed_entries_checked": len(entries_small),
            "words_with_varying_local_occupation_count": varying,
            "smaller_nonzero_coefficient_ranks": small_ranks,
            "larger_nonzero_coefficient_degrees": sorted(large_ranks)}


def canonical_reduction(indices, ambient):
    e = min(((a + 1) & -(a + 1)).bit_length() - 1 for a in indices)
    scale = 1 << e
    assert ambient % scale == 0
    reduced = tuple((a + 1) // scale - 1 for a in indices)
    assert any(a % 2 == 0 for a in reduced)
    return reduced, ambient // scale, scale


def check_gauge_transport(indices, ambient, scale):
    small, large, _, _ = check_blocks(indices, ambient, scale)
    n = len(indices)
    base_zero = ([0] * n + [1 + 2 * a for a in indices]
                 + [0] * n + [-2 * ambient + 1 + 2 * a for a in indices])
    base_infinity = ([2 * a - n + 1 for a in indices] + [-n] * n
                     + [2 * a - n + 1 for a in indices] + [2 * ambient - n] * n)
    checks = 0
    for perturbation in (0, 1, 3):
        zero = [value + perturbation * (j % 3) for j, value in enumerate(base_zero)]
        infinity = [value + perturbation * ((j % 5) - 2)
                    for j, value in enumerate(base_infinity)]
        new_zero = [scale * value for value in zero]
        correction = [scale - 1] * n + [0] * (2 * n) + [1 - scale] * n
        new_infinity = [scale * value + shift for value, shift in zip(infinity, correction)]
        for matrix, heights in ((small, zero), (large, new_zero)):
            for block in (0, 1):
                for i in range(block * n, (block + 1) * n):
                    for j in range(block * n, (block + 1) * n):
                        value = matrix[i][j]
                        if value:
                            order = (value & -value).bit_length() - 1
                            assert order + heights[i] + heights[2 * n + j] >= 0
        for at_infinity, rows, columns in (
                (False, tuple(range(n)), tuple(range(n, 2 * n))),
                (False, tuple(range(n, 2 * n)), tuple(range(n))),
                (True, tuple(range(2 * n)), tuple(range(2 * n)))):
            old_heights, new_heights = ((infinity, new_infinity) if at_infinity
                                        else (zero, new_zero))
            old_block = [[small[i][j] for j in columns] for i in rows]
            new_block = [[large[i][j] for j in columns] for i in rows]
            old_orders, _ = WeightedSmithOracle(old_block, infinity=at_infinity).profile(
                [old_heights[i] for i in rows], [old_heights[2 * n + j] for j in columns])
            new_orders, _ = WeightedSmithOracle(new_block, infinity=at_infinity).profile(
                [new_heights[i] for i in rows], [new_heights[2 * n + j] for j in columns])
            assert new_orders == tuple(scale * order for order in old_orders)
            checks += 1
        for mode in range(n):
            old_physical = sum(zero[group * n + mode] + infinity[group * n + mode]
                               for group in range(4))
            new_physical = sum(new_zero[group * n + mode] + new_infinity[group * n + mode]
                               for group in range(4))
            assert new_physical == scale * old_physical
    return {"indices": indices, "ambient": ambient, "scale": scale,
            "exact_weighted_profile_transports": checks,
            "independent_height_perturbations": 3}


def check_full_kernel_defect(indices, ambient, scale):
    n = len(indices)
    enlarged = tuple(scale * a + scale - 1 for a in indices)
    small = sparse_kernel(indices, ambient)
    large = sparse_kernel(enlarged, scale * ambient)
    # R and C have lambda exponents 0,1,1,0 and 0,-1,-1,0.
    row_exponents = (0, 1, 1, 0)
    column_exponents = (0, -1, -1, 0)
    lam = binomial_polynomial(scale - 1)
    mismatches_without_defect = 0
    checked = 0
    for a in range(4):
        for b in range(4):
            exponent = row_exponents[a] + column_exponents[b]
            defect = a == 1 and b in (0, 3)
            # zeta = w^(2(L-1)) / lambda^2.
            if defect:
                exponent -= 2
            for i in range(n):
                for j in range(n):
                    target = substitute_power(small[a * n + i][b * n + j], scale)
                    if defect:
                        target <<= 2 * (scale - 1)
                    actual = large[a * n + i][b * n + j]
                    left = polynomial_multiply(polynomial_power(lam, max(0, -exponent)), actual)
                    right = polynomial_multiply(polynomial_power(lam, max(0, exponent)), target)
                    assert left == right
                    checked += 1
                    if defect:
                        bare = substitute_power(small[a * n + i][b * n + j], scale)
                        raw_exponent = row_exponents[a] + column_exponents[b]
                        lhs = polynomial_multiply(polynomial_power(lam, max(0, -raw_exponent)), actual)
                        rhs = polynomial_multiply(polynomial_power(lam, max(0, raw_exponent)), bare)
                        mismatches_without_defect += int(lhs != rhs)
    assert mismatches_without_defect > 0
    return {"indices": indices, "ambient": ambient, "scale": scale,
            "full_kernel_entries_checked": checked,
            "entries_failing_if_scalar_defect_is_omitted": mismatches_without_defect}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    block_count = 0
    for indices in ((0,), (0, 1), (0, 2), (1, 3), (0, 1, 2), (0, 2, 5)):
        ambient = 1 << max(indices).bit_length()
        for scale in (2, 4, 8):
            check_blocks(indices, ambient, scale)
            block_count += 1
    minors = [check_all_minors(*case) for case in
              (((0,), 1, 2), ((0, 2), 4, 4), ((1, 3), 4, 2), ((0, 1, 2), 4, 2))]
    types = [check_type(*case) for case in
             (((0, 1), 2, 4, 2, (1, 1)), ((0, 2), 4, 2, 3, (1, 2)))]
    reductions = []
    for indices, ambient in (((0, 2, 5), 8), ((1, 3, 5), 8), ((3, 11, 23), 32),
                             ((2 ** 128 - 1,), 2 ** 128)):
        reduced, reduced_ambient, scale = canonical_reduction(indices, ambient)
        assert tuple(scale * a + scale - 1 for a in reduced) == indices
        reductions.append({"indices": indices, "ambient": ambient, "scale": scale,
                           "reduced_indices": reduced, "reduced_ambient": reduced_ambient})
    report = {"status": "exact finite replay of an all-size polynomial tensor identity; no improved exponent",
              "polynomial_block_cases": block_count, "minor_checks": minors, "typed_checks": types,
              "canonical_reductions": reductions,
              "gauge_transport": [check_gauge_transport(*case) for case in
                                   (((0,), 1, 4), ((0, 2), 4, 4), ((0, 1, 2), 4, 2))],
              "two_species_defect": check_full_kernel_defect((0, 2), 4, 4),
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
