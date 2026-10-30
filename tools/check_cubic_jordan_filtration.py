#!/usr/bin/env python3
"""All-power Jordan restrictions and an exact failed full-leg degeneration."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import INVERSE, matrices, scale
from explore_cubic_mixed_gf8 import cases, mixed_basis
from check_cubic_mixed_gf8 import Coordinates, flatten, reconstruct


def rational_rank(vectors):
    pivots = {}
    for source in vectors:
        row = list(map(Fraction, source))
        for pivot, old in sorted(pivots.items()):
            scalar = row[pivot]
            if scalar:
                row = [a - scalar * b for a, b in zip(row, old)]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is not None:
            scalar = row[pivot]
            pivots[pivot] = [value / scalar for value in row]
    return len(pivots)


def leg_weight(entry):
    x, y, v, u = entry
    result = [0] * 14
    result[x] += 1
    result[7 + y] += 1
    result[v] -= 1
    result[7 + u] -= 1
    return result


def projected_weight(entry):
    first = leg_weight(entry)
    return [first[i] + first[7 + i] for i in range(7)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(7)
    name, first, second = cases()[2]
    assert name == "jordan_shift_1_and_z"
    rows, columns, selected = mixed_basis(family, first, second)
    assert len(rows) == 8 and selected == [0, 1, 2, 3, 4, 6, 7, 8]
    coordinate_map = Coordinates(rows)
    zero = [[0] * 7 for _ in range(7)]
    values_zero = matrices(family, zero)
    values_one = matrices(family, first)
    values_z = matrices(family, second)
    polynomial_coordinates = []
    binary_entries = 0
    for k in range(5):
        formal = []
        for cut in range(2):
            constant = flatten(values_zero[cut][k])
            difference_one = flatten(values_one[cut][k]) ^ constant
            difference_z = flatten(values_z[cut][k]) ^ constant
            quadratic = scale(difference_z ^ scale(difference_one, 2, 2401), INVERSE[6], 2401)
            linear = difference_one ^ quadratic
            assert linear == 0
            formal.append((constant, linear, quadratic))
            for value in formal[-1]:
                assert all(((value >> (3 * entry)) & 7) <= 1 for entry in range(2401))
                binary_entries += 2401
        coordinates = []
        for degree in range(3):
            packed = formal[0][degree]
            matrix = [(packed >> (147 * i)) & ((1 << 147) - 1) for i in range(49)]
            code = coordinate_map.express(matrix)
            assert code is not None
            assert reconstruct(columns, code) == formal[1][degree]
            coordinates.append(code)
        polynomial_coordinates.append(coordinates)

    # Unordered nonedges of the logical seven-cycle, with both orientations.
    row_indices = [7 * x + y for x in range(7) for y in range(x + 1, 7)
                   if y - x not in (1, 6)]
    generic = row_indices + [7 * (i % 7) + i // 7 for i in row_indices]
    column_indices = sorted(generic, key=lambda i: (i // 7 + i % 7, i // 7, i % 7))
    assert len(row_indices) == 14 and len(set(column_indices)) == 28
    row_entries = column_entries = 0
    column_edges = []
    for a, row in enumerate(row_indices):
        for b, column in enumerate(row_indices):
            for generator, matrix in enumerate(rows):
                value = (matrix[row] >> (3 * column)) & 7
                assert value == int(generator == 0 and a == b)
                row_entries += 1
    for a, row in enumerate(column_indices):
        for b, column in enumerate(column_indices):
            values = [(matrix[row] >> (3 * column)) & 7 for matrix in columns]
            column_entries += len(values)
            if a == b:
                assert values == [1] + [0] * 7
            elif any(values):
                assert a < b
                old_weight = row // 7 + row % 7
                new_weight = column // 7 + column % 7
                assert new_weight == old_weight + 2
                column_edges.append({"row": row, "column": column, "coefficients": values})
    assert len(column_edges) == 20

    support = {}
    for x in range(7):
        for y in range(7):
            for v in range(7):
                for u in range(7):
                    values = [(matrix[7 * x + y] >> (3 * (7 * v + u))) & 7 for matrix in rows]
                    if any(values):
                        support[(x, y, v, u)] = values
    assert len(support) == 121

    # All full-leg weights have the form alpha_x+beta_y-alpha_v-beta_u.
    # Opposite inequalities first force alpha_i-beta_i to be constant.
    first_pairs, first_equations = [], []
    for j in range(1, 7):
        a, b = (0, j, j, 0), (j, 0, 0, j)
        assert a in support and b in support
        wa, wb = leg_weight(a), leg_weight(b)
        assert wa == [-value for value in wb]
        first_equations.append(wa)
        first_pairs.append([a, b])
    assert rational_rank(first_equations) == 6

    # After alpha=beta+constant, five opposite pairs force all adjacent sums
    # beta_j+beta_(j+1) to agree. The remaining beta vectors alternate by parity.
    second_pairs, second_equations = [], []
    for j in range(1, 6):
        a, b = (0, 1, j + 1, j), (j, j + 1, 1, 0)
        assert a in support and b in support
        wa, wb = projected_weight(a), projected_weight(b)
        assert wa == [-value for value in wb]
        second_equations.append(wa)
        second_pairs.append([a, b])
    assert rational_rank(second_equations) == 5
    remaining_gauges = [([1] * 7, [0] * 7), ([1] * 7, [1] * 7),
                        ([i % 2 for i in range(7)], [i % 2 for i in range(7)])]
    for entry in support:
        weight = leg_weight(entry)
        for alpha, beta in remaining_gauges:
            assert sum(a * b for a, b in zip(weight, alpha + beta)) == 0
    assert rational_rank([leg_weight(entry) for entry in support]) == 11
    assert 14 * 28 == 392 > 7 ** 3 and 3 ** 2 > 2 * 2 ** 2
    result = {
        "source_sha256": {source: hashlib.sha256(Path(__file__).with_name(source).read_bytes()).hexdigest()
                          for source in (Path(__file__).name, "explore_cubic_mixed_gf8.py",
                                         "explore_cubic_complement_gf8.py", "explore_cubic_complement_outputs.py",
                                         "check_cubic_mixed_gf8.py")},
        "labels": 7, "dimension": 8, "selected_raw_coefficient_indices": selected,
        "formal_scalar_family_descends_to_F2": True,
        "quadratic_scalar_family_coordinates": polynomial_coordinates,
        "binary_formal_coefficient_entries_checked": binary_entries,
        "separating_functional_in_selected_basis": [1] + [0] * 7,
        "row_scalar_compression_indices": row_indices,
        "row_scalar_compression_dimension": 14,
        "row_compression_coefficient_entries_checked": row_entries,
        "column_triangular_order": column_indices,
        "column_triangular_dimension": 28,
        "column_triangular_coefficient_entries_checked": column_entries,
        "column_strict_upper_edges": column_edges,
        "tensor_power_rank_lower_bases": [14, 28],
        "tensor_power_rank_cost_base_squared": 392,
        "full_support_entries": len(support),
        "full_leg_weight_constraint_rank": 11,
        "first_opposite_weight_witness_pairs": first_pairs,
        "second_opposite_weight_witness_pairs_after_alpha_minus_beta_constant": second_pairs,
        "feasible_full_leg_weight_gauges": remaining_gauges,
        "every_feasible_nonnegative_full_leg_weight_vanishes_on_support": True,
        "full_tensor_span_excluded": True,
        "scope": "Every characteristic-two field and every tensor power of the full formal coefficient space for C=lambda*J, with J the seven-label strict shift.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "column_strict_upper_edges"}), flush=True)


if __name__ == "__main__":
    main()
