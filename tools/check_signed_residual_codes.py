#!/usr/bin/env python3
"""Signed exterior residual, all-field rank spaces, and an F27 phase repair.

This is a finite endpoint certificate. It does not assert a complete
odd-characteristic polynomial factoring construction or a pressure gain.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def matrix_product(left, right):
    return [[sum(a * b for a, b in zip(row, column)) for column in zip(*right)]
            for row in left]


def matrix_sum(left, right):
    return [[a + b for a, b in zip(row, other)] for row, other in zip(left, right)]


def exterior_matrix(matrix):
    n = len(matrix)

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        result, remaining, position = 0, columns, 0
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length() - 1
            result += (-1) ** position * matrix[i][j] * minor(rows ^ bit, columns ^ cb)
            position += 1
        return result

    return [[minor(i, j) for j in range(1 << n)] for i in range(1 << n)]


def build_residual():
    matrix = [[int(i == j) for j in range(4)] for i in range(4)]
    matrix[0][3] = matrix[2][1] = 1
    inverse = [[int(i == j) for j in range(4)] for i in range(4)]
    inverse[0][3] = inverse[2][1] = -1
    gaussian, gaussian_inverse = exterior_matrix(matrix), exterior_matrix(inverse)
    creators, annihilators = [], []
    for a in range(4):
        creation, annihilation = [[0] * 16 for _ in range(16)], [[0] * 16 for _ in range(16)]
        for mask in range(16):
            sign = (-1) ** ((mask & ((1 << a) - 1)).bit_count())
            if mask >> a & 1:
                annihilation[mask ^ (1 << a)][mask] = sign
            else:
                creation[mask | (1 << a)][mask] = sign
        creators.append(creation)
        annihilators.append(annihilation)
    left = matrix_product(creators[1], creators[3])
    right = matrix_product(annihilators[2], annihilators[0])
    perturbation = matrix_product(matrix_product(left, gaussian), right)
    relative = matrix_product(perturbation, gaussian_inverse)
    assert matrix_product(relative, relative) == [[-a for a in row] for row in relative]
    identity = [[int(i == j) for j in range(16)] for i in range(16)]
    assert matrix_product(gaussian, gaussian_inverse) == identity
    return matrix_sum(gaussian, perturbation), gaussian, relative


def flatten(residual, cut):
    answer = [[0] * 16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        row, column = ((4 * x + y, 4 * u + v) if cut == 0 else
                       (4 * x + u, 4 * y + v) if cut == 1 else
                       (4 * x + v, 4 * y + u))
        answer[row][column] = residual[v + 4 * y][x + 4 * u]
    return answer


def prime_rank(matrix, prime):
    work = [[a % prime for a in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [a * inverse % prime for a in work[rank]]
        for i in range(rank + 1, len(work)):
            factor = work[i][column]
            if factor:
                work[i] = [(a - factor * b) % prime for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


def integer_determinant(matrix):
    work = [row[:] for row in matrix]
    previous, sign, n = 1, 1, len(work)
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if work[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            work[k], work[pivot] = work[pivot], work[k]
            sign = -sign
        value = work[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = value * work[i][j] - work[i][k] * work[k][j]
                assert numerator % previous == 0
                work[i][j] = numerator // previous
            work[i][k] = 0
        previous = value
    return sign * work[-1][-1]


def check_crossing_cut_spaces(residual):
    matrix = flatten(residual, 2)
    columns = [list(column) for column in zip(*matrix)]
    for column in columns:
        a = [column[4 * i:4 * i + 4] for i in range(4)]
        assert a[0][2] == a[0][3] == a[1][0] == a[3][0] == 0
        assert a[0][0] == a[1][1] == a[2][2]
        assert a[0][1] == -a[2][3] and a[2][0] == a[3][1]
    for row in matrix:
        a = [row[4 * i:4 * i + 4] for i in range(4)]
        assert a[0][1] == a[0][3] == a[2][0] == a[3][0] == 0
        assert a[0][0] == a[1][1] == a[2][2]
        assert a[0][2] == -a[1][3] and a[1][0] == a[3][2]
    selected_columns = (0, 2, 4, 6, 9, 11, 13, 15)
    selected_rows = (0, 1, 6, 7, 8, 9, 14, 15)
    column_coordinates = (0, 1, 8, 9, 6, 7, 14, 15)
    row_coordinates = (0, 4, 2, 6, 9, 13, 11, 15)
    column_det = integer_determinant([[columns[j][i] for j in selected_columns]
                                      for i in column_coordinates])
    row_det = integer_determinant([[matrix[j][i] for j in selected_rows]
                                   for i in row_coordinates])
    assert abs(column_det) == abs(row_det) == 1
    return {"selected_columns": selected_columns, "selected_rows": selected_rows,
            "column_coordinate_determinant": column_det,
            "row_coordinate_determinant": row_det}


def check_symbolic_phase(residual):
    choices = (((1, 0),), ((0, 0), (3, 1)))
    records = []
    for x, y, u, v in product(range(2), repeat=4):
        coefficients = {}
        for (xx, a), (yy, b), (uu, c), (vv, d) in product(
                choices[x], choices[y], choices[u], choices[v]):
            key = a + c, b + d
            coefficients[key] = coefficients.get(key, 0) + residual[vv + 4 * yy][xx + 4 * uu]
        coefficients = {k: a for k, a in coefficients.items() if a}
        expected = {}
        if x == v and y == u:
            expected[(0, 0)] = 1
            if x != u:
                expected[(1, 1)] = 1
            elif x == 1:
                expected[(1, 1)] = 2
        assert coefficients == expected
        records.append({"labels": [x, y, u, v],
                        "coefficients": [[*key, value] for key, value in sorted(coefficients.items())]})
    return records


def field_add(a, b):
    result, place = 0, 1
    for _ in range(3):
        result += ((a % 3 + b % 3) % 3) * place
        a //= 3
        b //= 3
        place *= 3
    return result


def field_multiply(a, b):
    aa, bb = [], []
    for _ in range(3):
        aa.append(a % 3)
        bb.append(b % 3)
        a //= 3
        b //= 3
    coefficients = [0] * 5
    for i, j in product(range(3), repeat=2):
        coefficients[i + j] += aa[i] * bb[j]
    # alpha^3 = 1-alpha^2, from alpha^3+alpha^2-1=0.
    for k in (4, 3):
        value = coefficients[k] % 3
        coefficients[k - 3] += value
        coefficients[k - 1] -= value
    return sum((coefficients[i] % 3) * 3 ** i for i in range(3))


def field_power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = field_multiply(result, a)
        a = field_multiply(a, a)
        exponent >>= 1
    return result


def check_ternary_phase_repair():
    assert all((x ** 3 + x ** 2 - 1) % 3 for x in range(3))
    alpha = 3
    assert field_add(field_add(field_power(alpha, 3), field_power(alpha, 2)), 2) == 0
    assert all(field_power(a, 26) == 1 for a in range(1, 27))
    conjugates = [field_power(alpha, 3 ** i) for i in range(3)]
    phases = [[1, field_add(1, a), field_add(1, a), field_add(1, field_multiply(2, a))]
              for a in conjugates]
    norms = []
    for j in range(4):
        value = 1
        for phase in phases:
            value = field_multiply(value, phase[j])
        assert value == 1
        norms.append(value)
    # Repeating the same logical core pair in the three copies gives identity.
    for output, input_ in product(range(4), repeat=2):
        value = 1
        for phase in phases:
            value = field_multiply(value, phase[input_] if output == input_ else 0)
        assert value == int(output == input_)
    return {"modulus": "alpha^3+alpha^2-1 over F3", "encoding": "base3 coefficient digits",
            "alpha": alpha, "conjugates": conjugates, "phase_diagonals": phases,
            "phase_norms": norms, "grouped_endpoint_entries": 16}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    residual, gaussian, relative = build_residual()
    ranks = {}
    for prime in (2, 3, 5, 7):
        cut_ranks = [prime_rank(flatten(residual, cut), prime) for cut in range(3)]
        assert cut_ranks[1:] == [12, 8]
        assert prime_rank(relative, prime) == 4
        ranks[str(prime)] = cut_ranks
    result = {"status": "signed residual endpoint identities only; no odd-characteristic factoring claim",
              "crossing_cut_spaces": check_crossing_cut_spaces(residual),
              "cut_ranks": ranks, "symbolic_phase_entries": check_symbolic_phase(residual),
              "F27_phase_repair": check_ternary_phase_repair(),
              "integer_residual": residual,
              "source_sha256": {"tools/check_signed_residual_codes.py":
                                  sha256(Path(__file__).read_bytes()).hexdigest()}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
