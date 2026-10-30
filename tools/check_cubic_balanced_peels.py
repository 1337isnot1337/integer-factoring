#!/usr/bin/env python3
"""Exact 177-dimensional normal form and structured balanced-tensor peel checks."""

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_fields import rank4
from check_cubic_varying_image_escape import rank_pair, xor_matrices


def dyad(first, second):
    """Column flattening is vec(first)*vec(second)^T for symmetric matrices."""
    rows = []
    for x in range(7):
        a = (first >> (7 * x)) & 127
        for y in range(7):
            rows.append(sum(a << (7 * v) for v in range(7) if (second >> (7 * y + v)) & 1))
    columns = [second if (first >> i) & 1 else 0 for i in range(49)]
    return rows, columns


def flatten(matrix):
    return sum(row << (49 * i) for i, row in enumerate(matrix))


def build_coordinates(basis):
    pivots = {}
    for i, matrix in enumerate(basis):
        value, coordinate = flatten(matrix), 1 << i
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value, coordinate
                break
            old, tag = pivots[pivot]
            value ^= old
            coordinate ^= tag
        assert value
    return pivots


def coordinates(matrix, pivots):
    value, coordinate = flatten(matrix), 0
    while value:
        pivot = value.bit_length() - 1
        assert pivot in pivots
        old, tag = pivots[pivot]
        value ^= old
        coordinate ^= tag
    return coordinate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(Path("evidence/cubic-varying-image-escape.json").read_text())
    family = ComplementOutputs(7)
    baseline = family.matrices([0] * 7)
    diagonal_units = [1 << (8 * i) for i in range(7)]
    identity_vector = sum(diagonal_units)
    edges = [(i, j, (1 << (7 * i + j)) | (1 << (7 * j + i))) for i, j in combinations(range(7), 2)]
    normal, names = [], []
    for i in range(7):
        for j in range(7):
            normal.append(dyad(diagonal_units[i], diagonal_units[j]))
            names.append(["D_D", i, j])
    for i, j, edge in edges:
        for vertex in (i, j):
            normal.append(dyad(edge, diagonal_units[vertex]))
            names.append(["A_D", i, j, vertex])
            normal.append(dyad(diagonal_units[vertex], edge))
            names.append(["D_A", vertex, i, j])
        normal.append(dyad(edge, edge))
        names.append(["A_A", i, j])
    assert len(normal) == 154
    for i, j, edge in edges:
        normal.append(dyad(edge, identity_vector))
        names.append(["A_I", i, j])
    endpoint = ([1 << i for i in range(49)], [1 << i for i in range(49)])
    normal.append(endpoint)
    names.append(["E"])
    normal.append((baseline[0][1], baseline[1][1]))
    names.append(["baseline_A"])
    assert len(normal) == 177
    pivots = build_coordinates([tensor[0] for tensor in normal])
    normal_columns = [flatten(tensor[1]) for tensor in normal]
    raw_checked, physical_coordinate_pivots = 0, {}
    for record in source["rank_one_output_maps"]:
        rows, columns = family.matrices(record["complement_rows"])
        for row, column in zip(rows, columns):
            coordinate = coordinates(row, pivots)
            reduced = coordinate
            while reduced:
                pivot = reduced.bit_length() - 1
                if pivot not in physical_coordinate_pivots:
                    physical_coordinate_pivots[pivot] = reduced
                    break
                reduced ^= physical_coordinate_pivots[pivot]
            reconstructed = 0
            for i in range(177):
                if (coordinate >> i) & 1:
                    reconstructed ^= normal_columns[i]
            assert reconstructed == flatten(column)
            raw_checked += 1
    assert raw_checked == 770
    assert len(physical_coordinate_pivots) == 177

    swap = [1 << (7 * (i % 7) + i // 7) for i in range(49)]
    rank_one = [identity_vector if i // 7 == i % 7 else 0 for i in range(49)]
    h, n = (rank_one, swap), (swap, rank_one)
    balanced = tuple(xor_matrices(e, first, second) for e, first, second in zip(endpoint, h, n))
    assert rank_pair(balanced) == [22, 22]

    # Verify the operator identity for all basis directions of an alternating
    # A. Linearity then proves it for arbitrary alternating A over any field.
    operator_columns = 0
    for alternating in [0] + [edge for _, _, edge in edges]:
        peel = dyad(identity_vector ^ alternating, identity_vector)
        residual = [xor_matrices(a, b) for a, b in zip(balanced, peel)]
        for v in range(7):
            for w in range(7):
                expected = 1 << (7 * v + w)
                for x in range(7):
                    if (alternating >> (7 * x + w)) & 1:
                        expected ^= 1 << (7 * x + v)
                if v == w:
                    expected ^= identity_vector
                actual = sum(((residual[0][row] >> (7 * v + w)) & 1) << row for row in range(49))
                assert actual == expected
                operator_columns += 1
    assert operator_columns == 1078
    assert 6 * 7 // 2 + 1 == 22 and 49 - 22 == 27
    assert 7 * 27 * 21 == 63 ** 2 and (10 ** 2) * 7 > 22 ** 2

    # Full rank controls for every edge at t=0,1 and at a nonbinary F4 value.
    edge_controls = []
    for i, j, edge in edges:
        first = dyad(edge, edge)
        direction = dyad(edge, identity_vector)
        records = []
        for scalar, expected_peel, expected_residual in ((0, [4, 1], [25, 21]),
                                                        (1, [12, 1], [34, 21]),
                                                        (2, [14, 1], [35, 21])):
            peel_pair, residual_pair = [], []
            for cut in range(2):
                peel_rows, residual_rows = [], []
                for row in range(49):
                    peel_value = residual_value = 0
                    for column in range(49):
                        value = ((first[cut][row] >> column) & 1) ^ (scalar if (direction[cut][row] >> column) & 1 else 0)
                        peel_value |= value << (2 * column)
                        residual_value |= (value ^ ((balanced[cut][row] >> column) & 1)) << (2 * column)
                    peel_rows.append(peel_value)
                    residual_rows.append(residual_value)
                peel_pair.append(peel_rows)
                residual_pair.append(residual_rows)
            peel_ranks = [rank4(matrix, 49) for matrix in peel_pair]
            residual_ranks = [rank4(matrix, 49) for matrix in residual_pair]
            assert peel_ranks == expected_peel and residual_ranks == expected_residual
            assert peel_ranks == [rank4(matrix, 49, low_pivot=True) for matrix in peel_pair]
            assert residual_ranks == [rank4(matrix, 49, low_pivot=True) for matrix in residual_pair]
            records.append({"parameter": scalar, "peel_ranks": peel_ranks, "residual_ranks": residual_ranks})
        edge_controls.append({"edge": [i, j], "controls": records})
    # The symbolic 2x2 kernel equation has coefficients d+tb,c+td,b+ta,a+tc.
    for basis in range(4):
        p = [[int(basis == 2 * i + j) for j in range(2)] for i in range(2)]
        a, b, c, d = p[0][0], p[0][1], p[1][0], p[1][1]
        for scalar in (0, 1):
            actual = [[p[i][j] ^ p[j][i] ^ p[1 - j][1 - i] ^ (scalar & p[j][1 - i])
                       for j in range(2)] for i in range(2)]
            assert actual == [[d ^ (scalar & b), c ^ (scalar & d)],
                              [b ^ (scalar & a), a ^ (scalar & c)]]
    result = {
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "check_cubic_varying_image_escape.py",
                                       "explore_cubic_complement_outputs.py", "explore_cubic_complement_fields.py")},
        "normal_form_dimension": 177, "normal_form_basis": names,
        "physical_raw_directions_reconstructed_in_both_cuts": raw_checked,
        "normal_form": "Q154 + Sym(out) tensor I(input) + span(E,baseline_A)",
        "alternating_peel_operator_basis_columns_checked": operator_columns,
        "alternating_peel_residual_row_rank_lower_bound": 27,
        "alternating_peel_residual_column_rank_lower_bound": 21,
        "alternating_peel_cost_lower_bound": "10*sqrt(7), strictly greater than22",
        "single_edge_parameter_classes": {"zero": {"peel": [4, 1], "residual": [25, 21]},
                                          "one": {"peel": [12, 1], "residual": [34, 21]},
                                          "other": {"peel": [14, 1], "residual": [35, 21]}},
        "single_edge_rank_controls": edge_controls,
        "all_field_claim_scope": "Only the two stated peel families; the normal form holds over all characteristic-two extensions.",
        "general_cost_36_optimality_claim": False,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("normal_form_basis", "single_edge_rank_controls")}), flush=True)


if __name__ == "__main__":
    main()
