#!/usr/bin/env python3
"""Literal checks of the graded Ward rank-descent formula."""

from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import all_length_labels, evaluate_local
from explore_symplectic_boundary_codes import boundary_code, half_orientation_code


@lru_cache(None)
def labels(pairs, grade):
    if not 0 <= grade <= 2 * pairs:
        return ()
    return all_length_labels(pairs, grade)


@lru_cache(None)
def local_centered(pairs):
    stride = 2 * pairs + 1
    return evaluate_local((3, 1 | (1 << stride), 1, 1, 1, 1))


@lru_cache(None)
def pattern_value(pairs, pattern):
    local = local_centered(pairs)
    value = 1
    for index in pattern:
        value = multiply(value, local[index])
        if not value:
            break
    return value


def entry(pairs, x, y, u, v):
    pattern = tuple(sorted(((a * 4 + b) * 4 + c) * 4 + d
                           for a, b, c, d in zip(x, y, u, v)))
    return pattern_value(pairs, pattern)


@lru_cache(None)
def rectangle_profile(pairs, ra, ry, cv, cu):
    first, second, third, fourth = (labels(pairs, grade) for grade in (ra, ry, cv, cu))
    if not first or not second or not third or not fourth:
        return {}
    row_count, col_count = len(first) * len(second), len(third) * len(fourth)
    stride = 2 * pairs + 1
    matrices = {}
    for x, y in product(range(len(first)), range(len(second))):
        row = x * len(second) + y
        for v, u in product(range(len(third)), range(len(fourth))):
            value = entry(pairs, first[x], second[y], fourth[u], third[v])
            column = v * len(fourth) + u
            while value:
                bit = value & -value
                exponent = bit.bit_length() - 1
                key = exponent % stride, exponent // stride
                if key not in matrices:
                    matrices[key] = [0] * row_count
                matrices[key][row] ^= 1 << column
                value ^= bit
    assert all(a + b <= ra for a, b in matrices)
    return {key: rank_two(matrix) for key, matrix in matrices.items()}


def supports(value):
    result = []
    while value:
        bit = value & -value
        result.append(bit.bit_length() - 1)
        value ^= bit
    return result


def actual_code_profile(code):
    pairs, length, d = code["pairs"], code["length"], code["labels"]
    ambient, projections = code["ambient_occupations"], code["output_indices"]
    injection = [supports(vector) for vector in code["input_vectors"]]
    matrices, column_matrices = {}, {}
    stride = 2 * pairs + 1
    for x, y, u, v in product(range(d), repeat=4):
        value = 0
        for yy, vv in product(injection[y], injection[v]):
            value ^= entry(pairs, ambient[projections[x]], ambient[yy],
                           ambient[projections[u]], ambient[vv])
        while value:
            bit = value & -value
            exponent = bit.bit_length() - 1
            key = exponent % stride, exponent // stride
            if key not in matrices:
                matrices[key] = [0] * (d * d)
                column_matrices[key] = [0] * (d * d)
            matrices[key][x * d + y] ^= 1 << (v * d + u)
            column_matrices[key][x * d + u] ^= 1 << (v * d + y)
            value ^= bit
    return {key: (rank_two(matrix), rank_two(column_matrices[key]))
            for key, matrix in matrices.items()}


def partial_form_code():
    pairs, length = 3, 2
    ambient = labels(pairs, length)
    index = {word: i for i, word in enumerate(ambient)}
    return {"pairs": pairs, "length": length, "labels": 1,
            "ambient_occupations": ambient,
            "input_vectors": [(1 << index[(3, 0, 0)]) ^ (1 << index[(0, 3, 0)])],
            "output_indices": [index[(0, 3, 0)]]}


def audit_case(code, forms, name):
    pairs, length = code["pairs"], code["length"]
    actual = actual_code_profile(code)
    report = []
    for (a, b), (row, column) in sorted(actual.items()):
        if b < forms:
            assert row == column == 0
            continue
        row_bound = column_bound = 0
        for j in range(min(forms, b - forms) + 1):
            clifford_count = forms - j
            ra = length - forms - j
            channel_bound = (1 << clifford_count) * comb(pairs, clifford_count)
            for h in range(clifford_count + 1):
                ry = length - 3 * forms + j + 2 * h
                key = a, b - forms - j
                row_bound += rectangle_profile(pairs, ra, ry, length, length).get(key, 0)
                column_bound += channel_bound * rectangle_profile(
                    pairs, ra, length, length, ry).get((key[1], key[0]), 0)
        assert row <= row_bound and column <= column_bound, (
            name, a, b, row, column, row_bound, column_bound)
        report.append([a, b, row, column, row_bound, column_bound])
    assert not any(b < forms and (row or column) for (a, b), (row, column) in actual.items())
    return {"name": name, "pairs": pairs, "grade": length, "forms": forms,
            "logical_labels": code["labels"],
            "jets_actual_row_actual_column_descent_row_descent_column": report}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    first = boundary_code(2, 2, True)
    second = boundary_code(3, 3, True)
    cases = ((first, 1, "global_P2_L2"),
             (second, 1, "global_P3_L3"),
             (half_orientation_code(second), 1, "global_P3_L3_half"),
             (partial_form_code(), 1, "partial_two_of_three_pairs"))
    result = {
        "status": "Literal complete coefficient ranks satisfy both graded Ward descent bounds.",
        "scope": "Finite controls for the all-size ordinary factorization; no score improvement asserted.",
        "cases": [audit_case(code, forms, name) for code, forms, name in cases],
        "literal_rectangle_profiles_constructed": rectangle_profile.cache_info().currsize,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()

