#!/usr/bin/env python3
"""All-field, all-power closure for mixtures of coordinate-image outputs."""

import argparse
import hashlib
import json
from pathlib import Path

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import matrices
from explore_cubic_mixed_gf8 import cases, mixed_basis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(7)
    baseline = matrices(family, [[0] * 7 for _ in range(7)])
    pairs = [(x, y) for x in range(7) for y in range(x + 1, 7) if y - x not in (1, 6)]
    row_pairs = sorted(pairs, reverse=True)
    column_pairs = pairs + [(y, x) for x, y in pairs]
    row_indices = [7 * x + y for x, y in row_pairs]
    column_indices = [7 * x + y for x, y in column_pairs]
    baseline_checks = 0
    for cut, indices in zip(baseline, (row_indices, column_indices)):
        for k, matrix in enumerate(cut):
            for a, row in enumerate(indices):
                for b, column in enumerate(indices):
                    assert ((matrix[row] >> (3 * column)) & 7) == int(k == 0 and a == b)
                    baseline_checks += 1

    # Basis of the enlarged update space: delta_(x,i)*delta_(w,i)*K_i(y,v),
    # with all 49 matrix units K_i permitted, without a symmetry assumption.
    row_checks = column_checks = 0
    update_positions = set()
    for image in range(7):
        for first in range(7):
            for second in range(7):
                for a, (x, y) in enumerate(row_pairs):
                    for b, (v, w) in enumerate(row_pairs):
                        value = int(x == w == image and y == first and v == second)
                        if value:
                            assert a < b and x > v
                            update_positions.add((a, b))
                        row_checks += 1
                for x, w in column_pairs:
                    for v, y in column_pairs:
                        assert not (x == w == image and y == first and v == second)
                        column_checks += 1
    assert len(update_positions) == 10

    name, first, second = cases()[4]
    assert name == "two_distinct_changed_output_rows"
    rows, columns, selected = mixed_basis(family, first, second)
    assert len(rows) == 7
    control_checks = 0
    for cut_name, cut, indices in (("row", rows, row_indices), ("column", columns, column_indices)):
        for k, matrix in enumerate(cut):
            for a, row in enumerate(indices):
                for b, column in enumerate(indices):
                    value = (matrix[row] >> (3 * column)) & 7
                    if a == b:
                        assert value == int(k == 0)
                    elif cut_name == "column" or a > b:
                        assert value == 0
                    control_checks += 1
    assert 14 * 28 == 392 > 7 ** 3 and 3 ** 2 > 2 * 2 ** 2
    result = {
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "explore_cubic_mixed_gf8.py",
                                       "explore_cubic_complement_gf8.py", "explore_cubic_complement_outputs.py")},
        "labels": 7, "baseline_restriction_coefficients_checked": baseline_checks,
        "universal_update_matrix_units": 343,
        "universal_row_restriction_coefficients_checked": row_checks,
        "universal_column_restriction_coefficients_checked": column_checks,
        "row_triangular_order": row_indices, "column_scalar_compression_indices": column_indices,
        "universal_strict_upper_update_positions": sorted(update_positions),
        "finite_case_selected_raw_coefficient_indices": selected,
        "finite_case_restriction_coefficients_checked": control_checks,
        "tensor_power_rank_lower_bases": [14, 28], "tensor_power_rank_cost_base_squared": 392,
        "scope": "Every characteristic-two field; arbitrary mixtures of all coefficient directions from complementary maps C=e_i*v^T, with both i and v varying; every tensor power. More generally, baseline plus delta_(x,w)*K_x(y,v) with arbitrary K_x.",
        "does_not_cover_simultaneous_nonzero_rows_in_one_complementary_map": True,
        "full_tensor_span_excluded": True,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
