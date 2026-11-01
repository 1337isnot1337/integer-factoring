#!/usr/bin/env python3
"""Symbolic all-field, all-power rank bound for the mixed cyclic cubic space."""

import argparse
import hashlib
import json
from pathlib import Path

from explore_cubic_complement_outputs import ComplementOutputs
from explore_cubic_complement_gf8 import INVERSE, TABLE, matrices, scale
from explore_cubic_mixed_gf8 import cases, mixed_basis
from check_cubic_mixed_gf8 import Coordinates, flatten, reconstruct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    family = ComplementOutputs(7)
    _, first, second = cases()[0]
    rows, columns, selected = mixed_basis(family, first, second)
    assert len(rows) == 7 and selected == [0, 1, 2, 3, 6, 7, 8]
    coordinates = Coordinates(rows)
    zero = [[0] * 7 for _ in range(7)]
    zero_rows, zero_columns = matrices(family, zero)
    zero_coordinates = [coordinates.express(matrix) for matrix in zero_rows]
    assert all(value is not None for value in zero_coordinates)
    for column, coordinate in zip(zero_columns, zero_coordinates):
        assert reconstruct(columns, coordinate) == flatten(column)

    # Each coefficient matrix for C=lambda*S is quadratic in lambda. Recover
    # its three polynomial coefficients from lambda=0,1,z, and check that they
    # descend to F2 while remaining in the seven-dimensional selected space.
    value_one = matrices(family, first)
    value_z = matrices(family, second)
    polynomial_coordinates = []
    binary_entries = 0
    for k in range(5):
        coefficient_pair = []
        for cut, values_zero, values_one, values_z in zip(
                (rows, columns), (zero_rows, zero_columns), value_one, value_z):
            constant = flatten(values_zero[k])
            first_difference = flatten(values_one[k]) ^ constant
            second_difference = flatten(values_z[k]) ^ constant
            quadratic = scale(second_difference ^ scale(first_difference, 2, 2401), INVERSE[6], 2401)
            linear = first_difference ^ quadratic
            assert linear == 0
            coefficient_pair.append((constant, linear, quadratic))
            for value in coefficient_pair[-1]:
                assert all(((value >> (3 * entry)) & 7) <= 1 for entry in range(2401))
                binary_entries += 2401
        entry_coordinates = []
        for degree, row_flat in enumerate(coefficient_pair[0]):
            row_matrix = [(row_flat >> (147 * i)) & ((1 << 147) - 1) for i in range(49)]
            coordinate = coordinates.express(row_matrix)
            assert coordinate is not None
            assert reconstruct(columns, coordinate) == coefficient_pair[1][degree]
            entry_coordinates.append(coordinate)
        polynomial_coordinates.append(entry_coordinates)

    powers = [1]
    for _ in range(6):
        powers.append(TABLE[powers[-1]][2])
    assert len(set(powers)) == 7 and TABLE[powers[-1]][2] == 1
    for a in range(7):
        total = 0
        for b in range(7):
            total ^= powers[(a * b) % 7]
        assert total == int(a == 0)
    orbit_pairs = ((0,), (1, 6), (2, 5), (3, 4))
    block_of = {value: index for index, block in enumerate(orbit_pairs) for value in block}
    e = [1, 0, 0, 0, 0, 0, 0]
    t = [0, 0, 0, 0, 5, 5, 0]
    v = [0, 0, 0, 1, 0, 0, 4]
    equivariance_entries = symbolic_entries = triangular_entries = 0
    blocks = []
    for cut_name, cut in zip(("row", "column"), (rows, columns)):
        for matrix in cut:
            for a in range(49):
                x, y = divmod(a, 7)
                moved_a = 7 * ((x + 1) % 7) + (y + 1) % 7
                for b in range(49):
                    u, w = divmod(b, 7)
                    moved_b = 7 * ((u + 1) % 7) + (w + 1) % 7
                    assert ((matrix[a] >> (3 * b)) & 7) == ((matrix[moved_a] >> (3 * moved_b)) & 7)
                    equivariance_entries += 1
        for q in range(7):
            block = [[[0] * 7 for _ in range(7)] for _ in range(7)]
            for k in range(7):
                for l in range(7):
                    for i, matrix in enumerate(cut):
                        value = 0
                        for shift in range(7):
                            value ^= TABLE[(matrix[k] >> (3 * (7 * shift + (shift + l) % 7))) & 7][
                                powers[(-q * shift) % 7]]
                        block[k][l][i] = value
                        symbolic_entries += 1
                    if block_of[k] != block_of[l]:
                        assert block[k][l] == [0] * 7
            for k in (2, 3, 4, 5):
                assert block[k][k] == e
                if cut_name == "row":
                    expected = [TABLE[value][powers[(-q * k) % 7]] for value in t]
                else:
                    expected = [TABLE[value][powers[(-q * (k + 1)) % 7]] for value in v]
                assert block[k][7 - k] == expected
            if cut_name == "row":
                assert block[6][1] == [TABLE[value][powers[q]] for value in t]
                for k in (2, 3):
                    for l in (2, 3):
                        assert block[k][l] == (e if k == l else [0] * 7)
            else:
                c = powers[-q % 7]
                for k in (2, 3):
                    # P has columns (rho^-k, 1) and (0, 1).
                    base = [[powers[-q * k % 7], 0], [1, 1]]
                    inverse = [[powers[q * k % 7], 0], [powers[q * k % 7], 1]]
                    indices = (k, 7 - k)
                    for i in range(7):
                        transformed = [[0, 0], [0, 0]]
                        for a in range(2):
                            for b in range(2):
                                for u in range(2):
                                    for w in range(2):
                                        transformed[a][b] ^= TABLE[inverse[a][u]][
                                            TABLE[block[indices[u]][indices[w]][i]][base[w][b]]]
                                triangular_entries += 1
                        diagonal = e[i] ^ TABLE[c][v[i]]
                        assert transformed == [[diagonal, TABLE[c][v[i]]], [0, diagonal]]
            blocks.append({"cut": cut_name, "charge": q, "entry_coefficient_vectors": block})

    # Exact comparison needs no decimal or irrational evaluation.
    assert 21 * 26 == 546 > 7 ** 3 and 3 ** 2 > (2 ** 2) * 2
    assert 14 * 4 * 6 == 336
    assert 336 ** 35 > 7 ** 99 and 99 ** 2 > 8 * 35 ** 2
    result = {
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, "explore_cubic_mixed_gf8.py",
                                       "explore_cubic_complement_gf8.py", "explore_cubic_complement_outputs.py",
                                       "check_cubic_mixed_gf8.py")},
        "labels": 7, "dimension": 7, "selected_raw_coefficient_indices": selected,
        "separating_functional": e, "row_generic_parameter": t, "column_generic_parameter": v,
        "zero_section_coordinates": zero_coordinates,
        "quadratic_scalar_family_coordinates": polynomial_coordinates,
        "formal_scalar_family_descends_to_F2": True,
        "binary_formal_coefficient_entries_checked": binary_entries,
        "cyclic_equivariance_entries_checked": equivariance_entries,
        "symbolic_block_coefficients_checked": symbolic_entries,
        "simultaneous_triangularization_coefficient_entries_checked": triangular_entries,
        "fourier_blocks": blocks,
        "rank_lower_bounds_when_separating_functional_nonzero": [21, 26],
        "rank_product_lower_bound": 546,
        "row_scalar_compression_dimension": 14,
        "column_simultaneously_upper_triangular_restriction_dimension": 28,
        "column_diagonal_functional": "e + rho^-1 v, four copies per seventh root rho",
        "tensor_power_grid_size_per_coordinate": 7,
        "tensor_power_polynomial_degree_per_coordinate": 1,
        "tensor_power_minimum_grid_support_base": 6,
        "tensor_power_column_rank_lower_base": 24,
        "tensor_power_rank_cost_base_squared": 336,
        "scope": "Every field of characteristic two; every tensor power of the full formal coefficient space for all scalar cyclic complementary maps C=lambda*S.",
        "strict_seed_excluded": True, "full_tensor_span_excluded": True,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "fourier_blocks"}), flush=True)


if __name__ == "__main__":
    main()
