#!/usr/bin/env python3
"""Exact natural-profile and unrestricted-gauge controls for weighted Pascal flags."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from cube_symmetric_exterior_blocks import series_smith_orders
from tensor_pascal_box_crossings import (
    exterior_tensor, matrix_product, polynomial_power, reverse_polynomial,
    weighted_infinity_matrix,
)
from weighted_pascal_flags import (
    change_flag, clifford_data, inverse_kernel, natural_heights,
    transpose, weighted_kernel,
)


def polynomial_shift(value, shift):
    if shift >= 0:
        return value << shift
    assert not value & ((1 << -shift) - 1)
    return value >> -shift


def exact_weighted_profile(matrix, shifts, row_heights, column_heights, infinity=False):
    n, entries, minimum, maximum = len(matrix), [], 0, 0
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if not value:
                continue
            shift = shifts[i][j]
            if infinity:
                shift = -shift - (value.bit_length() - 1)
                value = reverse_polynomial(value)
            shift += row_heights[i] + column_heights[j]
            valuation = (value & -value).bit_length() - 1
            minimum = min(minimum, shift + valuation)
            maximum = max(maximum, shift + value.bit_length() - 1)
            entries.append((i, j, value, shift))
    shifted = [[0] * n for _ in range(n)]
    for i, j, value, shift in entries:
        shifted[i][j] = polynomial_shift(value, shift - minimum)
    bound = n * (maximum - minimum)
    precision = min(bound + 1, max(64, 2 * (maximum - minimum) + 1))
    while True:
        try:
            return [value + minimum for value in series_smith_orders(shifted, precision)]
        except ArithmeticError:
            if precision >= bound + 1:
                raise
            precision = min(2 * precision, bound + 1)


def subset_sums(values):
    result = [0]
    for value in values:
        result += [item + value for item in result]
    return sorted(result)


def check_kernel(weights, outer):
    seed = weighted_kernel(weights, outer)
    n, matrix, pascal = seed["n"], seed["matrix"], seed["pascal"]
    first = [row[:n] for row in matrix[:n]]
    second = [row[n:] for row in matrix[n:]]
    target = [[polynomial_power(3, 2 * outer) if entry else 0 for entry in row] for row in pascal]
    assert matrix_product(first, second) == target
    assert matrix_product(matrix, inverse_kernel(seed)) == [[int(i == j) for j in range(2 * n)]
                                                           for i in range(2 * n)]
    for i, row in enumerate(second):
        for j, value in enumerate(row):
            order = (value & -value).bit_length() - 1
            predicted = 2 * sum(weight for k, weight in enumerate(seed["low_bits"])
                                if not (i | j) >> k & 1)
            assert order == predicted
            gauged = order + seed["zero_energies"][i] + seed["zero_energies"][j] - 2 * outer
            assert gauged == 2 * sum(weight for k, weight in enumerate(seed["low_bits"])
                                     if i & j & (1 << k))
    shifted, shift = weighted_infinity_matrix(seed)
    expected = sorted(seed["energies"] + [-value for value in seed["energies"]])
    assert [value + shift for value in series_smith_orders(shifted, max(expected) - shift + 1)] == expected
    heights = natural_heights(seed)
    cut_checks = []
    for place in (0, 1):
        values, shifts, rows, columns = clifford_data(seed, place)
        row = [heights[place][i] for i in rows]
        column = [heights[place][i] for i in columns]
        orders = exact_weighted_profile(values, shifts, row + [-x for x in row],
                                        column + [-x for x in column], bool(place))
        energies = seed["energies"] if place else seed["zero_energies"]
        assert orders == sorted((energies + [-x for x in energies]) * 2)
        cut_checks.append(dict(Counter(orders)))
    return {"weights": weights, "outer": outer, "modes": n,
            "zero_energies": seed["zero_energies"], "infinity_energies": seed["energies"],
            "clifford_profiles": cut_checks}


def original_fock_profile(seed, tensor, heights, place):
    n, q = seed["n"], len(tensor["alphabet"])
    physical = [[sum(heights[leg * n + bit] for bit in range(n) if symbol >> bit & 1)
                 for symbol in range(q)] for leg in range(4)]
    matrix = [[0] * (q * q) for _ in range(q * q)]
    shifts = [[0] * (q * q) for _ in range(q * q)]
    for (b, x, v, y), value in tensor["values"].items():
        row, column = (b + q * v, x + q * y) if not place else (b + q * x, v + q * y)
        matrix[row][column] = value
    if not place:
        row_heights = [physical[0][b] + physical[2][v] for v in range(q) for b in range(q)]
        column_heights = [physical[1][x] + physical[3][y] for y in range(q) for x in range(q)]
    else:
        row_heights = [physical[0][b] + physical[1][x] for x in range(q) for b in range(q)]
        column_heights = [physical[2][v] + physical[3][y] for y in range(q) for v in range(q)]
    return exact_weighted_profile(matrix, shifts, row_heights, column_heights, bool(place))


def check_arbitrary_gauges(weights, outer, perturbation, flag=False):
    seed = weighted_kernel(weights, outer)
    n = seed["n"]
    if flag:
        constant = [[int(i == j) for j in range(n)] for i in range(n)]
        constant[0][n - 1] = 1
        seed = change_flag(seed, constant)
    tensor = exterior_tensor(seed)
    heights = natural_heights(seed)
    result = []
    for place in (0, 1):
        changed = [value + perturbation * (((5 * index + 3 * place) % 7) - 3)
                   for index, value in enumerate(heights[place])]
        values, shifts, rows, columns = clifford_data(seed, place)
        row, column = [changed[i] for i in rows], [changed[i] for i in columns]
        orders = exact_weighted_profile(values, shifts, row + [-x for x in row],
                                        column + [-x for x in column], bool(place))
        assert orders == [-value for value in orders[::-1]]
        positive = orders[2 * n:]
        twice_mean = (2 * outer * n if place == 0 else 0) + sum(changed)
        twice_minimum = twice_mean - sum(positive)
        assert twice_minimum % 2 == 0
        minimum = twice_minimum // 2
        expected = [minimum + value for value in subset_sums(positive)]
        actual = original_fock_profile(seed, tensor, changed, place)
        assert expected == actual
        result.append({"place": place, "orders_compared": len(actual),
                       "minimum_order": min(actual), "maximum_order": max(actual)})
    return {"weights": weights, "outer": outer, "perturbation": perturbation,
            "flag_changed": flag, "cuts": result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    kernels = [check_kernel(weights, outer) for weights, outer in
               (([1], 2), ([1], 4), ([3], 4), ([1, 1], 4), ([1, 2], 4),
                ([1, 3], 8), ([2, 3], 8), ([1, 2, 4], 8), ([1, 2, 3], 8),
                ([1, 3, 5], 16), ([1, 2, 4, 7], 16), ([1, 2, 4, 8], 16))]
    gauges = [check_arbitrary_gauges(weights, outer, perturbation, flag)
              for weights, outer, perturbation, flag in
              (([1], 2, 1, False), ([3], 4, 2, False), ([3], 4, 1, True),
               ([1, 1], 4, 1, False), ([1, 3], 8, 2, False), ([1, 3], 8, 1, True))]
    sources = [Path(__file__), Path(__file__).with_name("weighted_pascal_flags.py"),
               Path(__file__).with_name("tensor_pascal_box_crossings.py"),
               Path(__file__).with_name("cube_symmetric_exterior_blocks.py")]
    report = {"status": "pass", "kernel_cases": kernels, "unrestricted_gauge_cases": gauges,
              "factoring_score": "1/2",
              "source_sha256": {path.name: sha256(path.read_bytes()).hexdigest() for path in sources}}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(json.dumps({"status": "pass", "kernel_cases": len(kernels),
                      "unrestricted_gauge_cases": len(gauges),
                      "original_fock_orders_compared": sum(cut["orders_compared"]
                          for case in gauges for cut in case["cuts"])}, indent=2))


if __name__ == "__main__":
    main()
