#!/usr/bin/env python3
"""Exact block Smith profiles for ordered polynomial crossings over F_2."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import check_cancellation_module_automata as local_tools


def series_inverse(value, precision):
    assert value & 1
    inverse, available = 1, 1
    while available < precision:
        available = min(2 * available, precision)
        inverse = local_tools.pmul(inverse, local_tools.pmul(value, inverse)) & ((1 << available) - 1)
    assert local_tools.pmul(value, inverse) & ((1 << precision) - 1) == 1
    return inverse


def smith_at_precision(matrix, precision):
    mask = (1 << precision) - 1
    work = [[value & mask for value in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    orders = []
    for k in range(columns):
        best = None
        for i in range(k, rows):
            for j in range(k, columns):
                if work[i][j]:
                    order = local_tools.porder(work[i][j])
                    if best is None or order < best[0]:
                        best = order, i, j
                        if not order:
                            break
            if best is not None and not best[0]:
                break
        if best is None:
            return None
        order, pi, pj = best
        orders.append(order)
        work[k], work[pi] = work[pi], work[k]
        if pj != k:
            for row in work:
                row[k], row[pj] = row[pj], row[k]
        inverse = series_inverse(work[k][k] >> order, precision)
        for j in range(k, columns):
            work[k][j] = local_tools.pmul(work[k][j], inverse) & mask
        assert work[k][k] == 1 << order
        support = [j for j in range(k + 1, columns) if work[k][j]]
        for i in range(k + 1, rows):
            factor = work[i][k] >> order
            if factor:
                for j in support:
                    work[i][j] ^= local_tools.pmul(factor, work[k][j]) & mask
            work[i][k] = 0
    assert orders == sorted(orders)
    return tuple(orders)


def smith_profile(matrix, degree):
    """All visible minimum pivots certify exact orders; never infer infinity."""
    maximum = len(matrix[0]) * degree + 1
    precision = min(maximum, max(2, 2 * degree + 2))
    while True:
        orders = smith_at_precision(matrix, precision)
        if orders is not None:
            return orders, precision
        assert precision < maximum, "generic rank deficiency or invalid supplied degree bound"
        precision = min(2 * precision, maximum)


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            column = row.bit_length() - 1
            if column not in pivots:
                pivots[column] = row
                break
            row ^= pivots[column]
    return len(pivots)


def toeplitz_rank(matrix, precision):
    rows, columns = len(matrix), len(matrix[0])
    result = []
    for output_degree in range(precision):
        for i in range(rows):
            row = 0
            for input_degree in range(output_degree + 1):
                coefficient = output_degree - input_degree
                for j in range(columns):
                    if matrix[i][j] >> coefficient & 1:
                        row |= 1 << (input_degree * columns + j)
            result.append(row)
    return binary_rank(result)


def polynomial_matrix_product(left, right):
    dimension = len(left)
    answer = [[0] * dimension for _ in range(dimension)]
    for i, j, k in product(range(dimension), repeat=3):
        answer[i][k] ^= local_tools.pmul(left[i][j], right[j][k])
    return answer


def block_crossing(local, q, length):
    dimension = len(next(iter(local.values())))
    alphabet = q ** length
    answer = {}
    for labels in product(range(alphabet), repeat=4):
        matrix = [[int(i == j) for j in range(dimension)] for i in range(dimension)]
        current = labels
        for _ in range(length):
            digits = tuple(value % q for value in current)
            current = tuple(value // q for value in current)
            matrix = polynomial_matrix_product(matrix, local[digits])
        for a, b in product(range(dimension), repeat=2):
            assert matrix[a][b].bit_count() % 2 == int(labels[0] == labels[3]
                                                      and labels[1] == labels[2] and a == b)
        answer[labels] = matrix
    return answer


def flatten(local, q, cut, side=0):
    dimension = len(next(iter(local.values())))
    n = q * q
    matrix = [[0] * n for _ in range(n * dimension * dimension)]
    for (x, y, u, v), values in local.items():
        left, right = (q * x + y, q * u + v) if cut == 0 else (q * x + u, q * y + v)
        physical, other = (left, right) if side == 0 else (right, left)
        for a, b in product(range(dimension), repeat=2):
            matrix[(physical * dimension + a) * dimension + b][other] = values[a][b]
    return matrix


def reverse_polynomial(value, degree):
    assert value.bit_length() <= degree + 1
    return sum(1 << (degree - k) for k in range(degree + 1) if value >> k & 1)


def polynomial_profiles(local, q, degree):
    profiles = {}
    counts = Counter()
    for cut in range(2):
        matrix = flatten(local, q, cut)
        zero, zero_precision = smith_profile(matrix, degree)
        reversed_matrix = [[reverse_polynomial(value, degree) for value in row] for row in matrix]
        infinity, infinity_precision = smith_profile(reversed_matrix, degree)
        for polynomial, orders in ((matrix, zero), (reversed_matrix, infinity)):
            for precision in range(1, max(orders) + 2):
                actual = toeplitz_rank(polynomial, precision)
                assert actual == sum(max(0, precision - order) for order in orders)
                counts["independent_toeplitz_ranks"] += 1
        profiles[cut] = {"zero": zero, "infinity_orders": infinity,
                         "degrees": tuple(degree - value for value in infinity),
                         "zero_precision": zero_precision, "infinity_precision": infinity_precision}
    return profiles, counts


def word_sums(profile, s):
    return (sum(s ** value for value in profile["degrees"]),
            sum(s ** (-value) for value in profile["zero"]))


def check_family(name, local, degree, maximum_block=3):
    records, profiles = [], {}
    counts = Counter()
    points = tuple(map(Fraction, ("1", "5/4", "3/2", "2", "3", "4", "5")))
    for length in range(1, maximum_block + 1):
        q = 2 ** length
        block = block_crossing(local, 2, length)
        profiles[length], checked = polynomial_profiles(block, q, degree * length)
        counts.update(checked)
        counts["local_matrix_endpoint_entries"] += len(block) * len(next(iter(local.values()))) ** 2
        best = None
        for s, first_cut in product(points, range(2)):
            U, _ = word_sums(profiles[length][first_cut], s)
            _, V = word_sums(profiles[length][1 - first_cut], s)
            candidate = U * V
            if best is None or candidate < best[0]:
                best = candidate, s, first_cut
        records.append({"block_length": length, "alphabet": q, "B": str(best[0]),
                        "s": str(best[1]), "infinity_cut": best[2],
                        "profiles": {str(cut): {key: list(value) if isinstance(value, tuple) else value
                                               for key, value in profile.items()}
                                     for cut, profile in profiles[length].items()}})
    for a in range(1, maximum_block):
        for b in range(1, maximum_block - a + 1):
            for cut in range(2):
                left, right, whole = profiles[a][cut], profiles[b][cut], profiles[a + b][cut]
                zero_sum = sorted(x + y for x in left["zero"] for y in right["zero"])
                degree_sum = sorted((x + y for x in left["degrees"] for y in right["degrees"]), reverse=True)
                assert all(x >= y for x, y in zip(whole["zero"], zero_sum))
                assert all(x <= y for x, y in zip(whole["degrees"], degree_sum))
                counts["coordinatewise_smith_comparisons"] += 2 * len(zero_sum)
                for s in points:
                    Ua, Va = word_sums(left, s)
                    Ub, Vb = word_sums(right, s)
                    Uab, Vab = word_sums(whole, s)
                    assert Uab <= Ua * Ub
                    assert Vab <= Va * Vb
                    counts["exact_submultiplicativity_checks"] += 2
    if name == "auxiliary elementary":
        assert all(Fraction(record["B"]) == Fraction(256, 27) ** record["block_length"] for record in records)
    if name == "ordered pairing" and maximum_block >= 3:
        assert records[0]["B"] == "16" and records[1]["B"] == "256"
        assert Fraction(records[2]["B"]) == Fraction(1241625872629, 320000000)
        assert Fraction(records[2]["B"]) < 16 ** 3
    return {"name": name, "degree": degree, "records": records, **dict(counts)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-block", type=int, default=3)
    args = parser.parse_args()
    cases = [
        check_family("auxiliary elementary", *local_tools.elementary(2), args.maximum_block),
        check_family("ordered pairing", *local_tools.pairing(), args.maximum_block),
    ]
    # A unit term beyond a supplied jet must never be declared exact zero.
    assert smith_at_precision([[1, 0], [0, 1 << 7]], 7) is None
    assert smith_profile([[1, 0], [0, 1 << 7]], 7)[0] == (0, 7)
    report = {"status": "pass", "field": 2, "maximum_block": args.maximum_block,
              "cases": cases, "truncated_zero_negative_control": True,
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "local_constructor_sha256": sha256(Path(local_tools.__file__).read_bytes()).hexdigest(),
              "scope": "Exact fixed-block Smith improvement and submultiplicative profile checks; no improved factoring exponent."}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
