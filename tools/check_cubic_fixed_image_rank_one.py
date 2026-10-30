#!/usr/bin/env python3
"""Universal fixed-image rank-one barrier, with an explicit F8 control witness."""

import argparse
import hashlib
import json
from pathlib import Path

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import INVERSE, TABLE, matrices
from explore_cubic_mixed_gf8 import cases, mixed_basis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(7)
    baseline = matrices(family, [[0] * 7 for _ in range(7)])
    pairs = [7 * x + y for x in range(7) for y in range(x + 1, 7) if y - x not in (1, 6)]
    oriented = pairs + [7 * (i % 7) + i // 7 for i in pairs]
    assert len(pairs) == 14 and len(set(oriented)) == 28
    baseline_entries = 0
    for cut, indices in zip(baseline, (pairs, oriented)):
        for k, matrix in enumerate(cut):
            for a, row in enumerate(indices):
                for b, column in enumerate(indices):
                    assert ((matrix[row] >> (3 * column)) & 7) == int(k == 0 and a == b)
                    baseline_entries += 1

    mixed_entries = symmetry_entries = 0
    for k in range(5):
        for y in range(7):
            for v in range(7):
                parent = family.parent[k][y][v]
                for j in range(7):
                    for l in range(7):
                        assert ((parent[j] >> (7 + l)) & 1) == 0
                        assert ((parent[7 + j] >> l) & 1) == 0
                        mixed_entries += 2
                        assert ((parent[7 + j] >> (7 + l)) & 1) == (
                            (family.parent[k][v][y][7 + l] >> (7 + j)) & 1)
                        symmetry_entries += 1

    # For arbitrary u and symmetric K, the column update is a*b^T, where
    # a_(x,u)=u_x*u_u. Check b^T*a=0 as a polynomial in formal u variables,
    # separately for every basis matrix of the seven-by-seven symmetric space.
    symmetric_basis_checks = 0
    for i in range(7):
        for j in range(i, 7):
            monomials = set()
            for index in oriented:
                v, y = divmod(index, 7)
                active = ((y == i and v == j) or (i != j and y == j and v == i))
                if active:
                    monomial = tuple(sorted((v, y)))
                    if monomial in monomials:
                        monomials.remove(monomial)
                    else:
                        monomials.add(monomial)
            assert not monomials
            symmetric_basis_checks += 1
    assert symmetric_basis_checks == 28
    groups = {}
    for index, pair in enumerate(pairs):
        groups.setdefault(pair % 7, []).append(index)
    assert sorted(map(len, groups.values())) == [1, 2, 3, 4, 4]

    # A fully numerical control of the universal construction for the dense
    # rank-one pair in the complete mixed-space archive.
    name, first, second = cases()[3]
    assert name == "dense_rank_one_1_and_z"
    rows, columns, selected = mixed_basis(family, first, second)
    image = [1, 2, 4, 3, 6, 7, 5]
    covector = [1, 1, 2, 2, 4, 4, 3]
    assert first == [[TABLE[a][b] for b in covector] for a in image]
    assert len(rows) == 8
    eliminators, kept = [], []
    for group in groups.values():
        pivot = group[0]
        for index in group[1:]:
            scalar = TABLE[image[pairs[index] // 7]][INVERSE[image[pairs[pivot] // 7]]]
            eliminators.append({"index": index, "pivot": pivot, "pivot_coefficient": scalar})
            kept.append(index)
    assert len(kept) == 9
    row_checks = 0
    for generator, matrix in enumerate(rows):
        for a, row in enumerate(eliminators):
            for b, column in enumerate(kept):
                value = ((matrix[pairs[row["index"]]] >> (3 * pairs[column])) & 7) ^ TABLE[
                    row["pivot_coefficient"]][(matrix[pairs[row["pivot"]]] >> (3 * pairs[column])) & 7]
                assert value == int(generator == 0 and a == b)
                row_checks += 1
    a_vector = [TABLE[image[index // 7]][image[index % 7]] for index in oriented]
    pivot = next(i for i, value in enumerate(a_vector) if value)
    b_vectors, column_checks = [], 0
    for generator, matrix in enumerate(columns):
        b_vector = [TABLE[INVERSE[a_vector[pivot]]][
            ((matrix[oriented[pivot]] >> (3 * index)) & 7) ^ int(generator == 0 and j == pivot)]
            for j, index in enumerate(oriented)]
        for i, row in enumerate(oriented):
            for j, column in enumerate(oriented):
                assert ((matrix[row] >> (3 * column)) & 7) == (
                    int(generator == 0 and i == j) ^ TABLE[a_vector[i]][b_vector[j]])
                column_checks += 1
        inner = 0
        for a, b in zip(a_vector, b_vector):
            inner ^= TABLE[a][b]
        assert inner == 0
        b_vectors.append(b_vector)

    assert 9 * 28 == 252
    assert 252 ** 6 > 7 ** 17 and 17 ** 2 > 2 * 12 ** 2
    result = {
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "explore_cubic_mixed_gf8.py",
                                       "explore_cubic_complement_gf8.py", "explore_cubic_complement_outputs.py")},
        "labels": 7, "baseline_restriction_coefficient_entries_checked": baseline_entries,
        "parent_low_high_zero_entries_checked": mixed_entries,
        "parent_high_high_swap_symmetry_entries_checked": symmetry_entries,
        "formal_symmetric_matrix_basis_checks": symmetric_basis_checks,
        "baseline_row_indices": pairs, "column_indices": oriented,
        "row_image_groups": groups, "row_image_dimension_upper_bound": 5,
        "row_scalar_compression_dimension_lower_bound": 9,
        "column_common_triangular_dimension": 28,
        "tensor_power_rank_lower_bases": [9, 28],
        "tensor_power_rank_cost_base_squared": 252,
        "recorded_dense_case": {
            "selected_raw_coefficient_indices": selected, "image_vector": image,
            "covector": covector, "row_eliminators": eliminators,
            "row_right_inverse_coordinate_indices": kept,
            "row_compression_coefficient_entries_checked": row_checks,
            "column_common_image_vector": a_vector, "column_image_pivot": pivot,
            "column_update_covectors": b_vectors,
            "column_factorization_coefficient_entries_checked": column_checks,
            "all_column_update_covectors_annihilate_common_image": True,
        },
        "scope": "Every characteristic-two field, every fixed image vector u, arbitrary complementary maps C=u*v^T with varying v, their full combined coefficient space, and every tensor power; also the larger class baseline plus u_x*u_u*K(y,v) for arbitrary symmetric K.",
        "fixed_image_required_within_each_local_coefficient_space": True,
        "full_tensor_span_excluded": True,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "recorded_dense_case"}), flush=True)


if __name__ == "__main__":
    main()
