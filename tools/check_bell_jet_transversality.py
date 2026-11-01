#!/usr/bin/env python3
"""Exact tangent-space and seventh/eighth Bell-jet controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path


def multiply(first, second):
    raw = 0
    for bit in range(2):
        if second >> bit & 1:
            raw ^= first << bit
    if raw & 4:
        raw ^= 7
    return raw


def inverse(value):
    assert value
    return next(candidate for candidate in (1, 2, 3)
                if multiply(value, candidate) == 1)


def matrix_inverse(matrix):
    q = len(matrix)
    rows = [list(row) + [int(i == j) for j in range(q)]
            for i, row in enumerate(matrix)]
    for column in range(q):
        pivot = next(i for i in range(column, q) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scalar = inverse(rows[column][column])
        rows[column] = [multiply(scalar, value) for value in rows[column]]
        for i in range(q):
            if i != column and rows[i][column]:
                scalar = rows[i][column]
                rows[i] = [first ^ multiply(scalar, second)
                           for first, second in zip(rows[i], rows[column])]
    return [row[q:] for row in rows]


def packed_scale(row, scalar):
    first, second = row
    return ((0, 0), (first, second), (second, first ^ second),
            (first ^ second, first))[scalar]


def packed_rank(rows):
    pivots = {}
    for first, second in rows:
        while first | second:
            pivot = (first | second).bit_length() - 1
            coefficient = ((first >> pivot) & 1) | (((second >> pivot) & 1) << 1)
            if pivot not in pivots:
                pivots[pivot] = packed_scale((first, second), inverse(coefficient))
                break
            other = packed_scale(pivots[pivot], coefficient)
            first ^= other[0]
            second ^= other[1]
    return len(pivots)


def encode(coefficients):
    first = second = 0
    for position, value in coefficients.items():
        first ^= (value & 1) << position
        second ^= ((value >> 1) & 1) << position
    return first, second


def tangent_equations(matrix, cut):
    q = len(matrix)
    second = matrix_inverse(matrix)
    alpha = [matrix[i % q][i // q] for i in range(q * q)]
    beta = [second[i % q][i // q] for i in range(q * q)]
    row_pivot = next(i for i, value in enumerate(alpha) if value)
    column_pivot = next(i for i, value in enumerate(beta) if value)
    alpha = [multiply(value, inverse(alpha[row_pivot])) for value in alpha]
    beta = [multiply(value, inverse(beta[column_pivot])) for value in beta]

    def position(i, j):
        if cut == 0:
            return i + q * q * j
        x, u, y, v = i % q, i // q, j % q, j // q
        return x + q * y + q * q * u + q**3 * v

    equations = []
    for i, j in product(range(q * q), repeat=2):
        if i == row_pivot or j == column_pivot:
            continue
        coefficients = {}
        for row, column, value in (
            (i, j, 1), (row_pivot, j, alpha[i]),
            (i, column_pivot, beta[j]),
            (row_pivot, column_pivot, multiply(alpha[i], beta[j])),
        ):
            p = position(row, column)
            coefficients[p] = coefficients.get(p, 0) ^ value
        equations.append(encode(coefficients))
    return equations


def identity(q):
    return [[int(i == j) for j in range(q)] for i in range(q)]


def upper(q, shift):
    result = identity(q)
    for i in range(q - 1):
        result[i][i + 1] = 1 + (i + shift) % 3
    return result


def power(value, exponent):
    result = 1
    for _ in range(exponent):
        result = multiply(result, value)
    return result


def hasse(pair, point, order):
    result = 0
    for degree in range(max(pair[0].bit_length(), pair[1].bit_length())):
        if degree >= order and comb(degree, order) % 2:
            coefficient = ((pair[0] >> degree) & 1) | (((pair[1] >> degree) & 1) << 1)
            result ^= multiply(coefficient, power(point, degree - order))
    return result


def flatten(values, cut):
    result = [[0] * 16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        i, j = ((x + 4 * y, u + 4 * v) if cut == 0
                else (x + 4 * u, y + 4 * v))
        result[i][j] = values[x + 4 * y + 16 * u + 64 * v]
    return result


def ordinary_rank(matrix):
    return packed_rank([encode({j: value for j, value in enumerate(row)})
                        for row in matrix])


def projected_rank(leading, next_matrix):
    assert ordinary_rank(leading) == 1
    rp, cp = next((i, j) for i, row in enumerate(leading)
                  for j, value in enumerate(row) if value)
    pivot_inverse = inverse(leading[rp][cp])
    projected = []
    for i in range(16):
        if i == rp:
            continue
        row = []
        ai = multiply(leading[i][cp], pivot_inverse)
        for j in range(16):
            if j == cp:
                continue
            bj = multiply(leading[rp][j], pivot_inverse)
            row.append(next_matrix[i][j]
                       ^ multiply(ai, next_matrix[rp][j])
                       ^ multiply(bj, next_matrix[i][cp])
                       ^ multiply(multiply(ai, bj), next_matrix[rp][cp]))
        projected.append(row)
    return ordinary_rank(projected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    alternating = [[int(j == (i ^ 1)) for j in range(4)] for i in range(4)]
    cases = [
        ("two-letter boundary", identity(2), identity(2), 3),
        ("three-letter identity", identity(3), identity(3), 0),
        ("three-letter nonsymmetric", upper(3, 0), upper(3, 1), 0),
        ("four-letter identity", identity(4), identity(4), 0),
        ("four-letter alternating", alternating, identity(4), 0),
        ("four-letter nonsymmetric", upper(4, 0), upper(4, 1), 0),
        ("five-letter nonsymmetric", upper(5, 1), upper(5, 2), 0),
    ]
    dimensions = []
    for label, row, column, expected in cases:
        equations = tangent_equations(row, 0) + tangent_equations(column, 1)
        rank = packed_rank(equations)
        dimension = len(row)**4 - rank
        assert dimension == expected
        dimensions.append({"case": label, "alphabet": len(row),
                           "equations": len(equations), "rank": rank,
                           "intersection_dimension": dimension})

    from bell_jet_crt_crossings import construct, scalar_data
    from search_entrywise_reflection_jets import bell_values
    scalar = scalar_data()
    scalar_controls = []
    for point in (2, 3):
        leading = hasse(scalar["kernel"], point, 7)
        normalized = {key: multiply(hasse(value, point, 7), inverse(leading))
                      for key, value in scalar.items()}
        assert normalized == {"endpoint": 1, "row": 2, "column": 3, "kernel": 1}
        scalar_controls.append({"point": point, "normalized_seventh_coefficients": normalized})

    jet_entries = 0
    jet_controls = []
    for label, a, b, _ in cases[3:6]:
        row, column = bell_values(a, 0), bell_values(b, 1)
        carrier = []
        for index in range(256):
            x, y, u, v = ((index // 4**leg) % 4 for leg in range(4))
            carrier.append(int(x == v and y == u)
                           ^ multiply(2, row[index]) ^ multiply(3, column[index]))
        seed = construct(a, b, carrier)
        ranks = []
        for point, cut in ((2, 0), (3, 1)):
            for order in range(1, 8):
                assert all(hasse(value, point, order) == 0 for value in seed["values"])
                jet_entries += 256
            leading = flatten([hasse(value, point, 0) for value in seed["values"]], cut)
            eighth = flatten([hasse(value, point, 8) for value in seed["values"]], cut)
            ranks.append(projected_rank(leading, eighth))
            jet_entries += 512
        if label == "four-letter identity":
            assert ranks == [6, 6]
        if label == "four-letter alternating":
            assert ranks == [4, 4]
        jet_controls.append({"case": label, "degree": seed["degree"],
                             "eighth_schur_ranks": ranks})
    names = (Path(__file__).name, "bell_jet_crt_crossings.py",
             "entrywise_reflection_constraints.py", "search_entrywise_reflection_jets.py")
    report = {
        "status": "passed",
        "scope": "exact tangent intersections and normalized seventh jets; complete eighth-order controls at three relative Bell forms",
        "tangent_cases": dimensions,
        "scalar_controls": scalar_controls,
        "complete_hasse_entries": jet_entries,
        "jet_controls": jet_controls,
        "source_sha256": {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "source_sha256"}))


if __name__ == "__main__":
    main()
