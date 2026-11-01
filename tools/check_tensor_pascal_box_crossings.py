#!/usr/bin/env python3
"""Independent finite controls for the tensor-product Pascal box theorem."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from cube_symmetric_exterior_blocks import series_smith_orders
from tensor_pascal_box_crossings import (
    box_kernel, exterior_tensor, matrix_product, physical_profiles,
    polynomial_power, small_box_code, three_quarter_certificate,
    weighted_infinity_matrix,
)


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            position = row.bit_length() - 1
            if position in pivots:
                row ^= pivots[position]
            else:
                pivots[position] = row
                break
    return len(pivots)


def toeplitz_rank(matrix, length):
    n = len(matrix)
    rows = []
    for output in range(length):
        for i in range(n):
            row = 0
            for source in range(output + 1):
                degree = output - source
                for j in range(n):
                    if matrix[i][j] >> degree & 1:
                        row ^= 1 << (source * n + j)
            rows.append(row)
    return binary_rank(rows)


def subset_sums(values):
    sums = [0]
    for value in values:
        sums += [item + value for item in sums]
    return sorted(sums)


def one_particle_case(coordinates, side):
    seed = box_kernel(coordinates, side)
    n, outer, matrix = seed["n"], seed["outer"], seed["matrix"]
    left = [row[:n] for row in matrix[:n]]
    right = [row[n:] for row in matrix[n:]]
    expected = [[polynomial_power(3, 2 * outer) if value else 0 for value in row]
                for row in seed["pascal"]]
    assert matrix_product(left, right) == expected
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            endpoint = value.bit_count() & 1
            target = seed["pascal"][i][j - n] if i < n <= j else int(i == j + n)
            assert endpoint == target
    weighted, shift = weighted_infinity_matrix(seed)
    expected_orders = sorted([value - shift for energy in seed["energies"]
                              for value in (energy - 2 * outer, energy)])
    orders = series_smith_orders(weighted, max(expected_orders) + 1)
    assert orders == expected_orders
    controls = 0
    if n <= 4:
        for length in sorted({1} | {order + 1 for order in orders}):
            assert toeplitz_rank(weighted, length) == sum(max(0, length - order) for order in orders)
            controls += 1
    return {"coordinates": coordinates, "side": side, "modes": n,
            "infinity_orders": dict(Counter(order + shift for order in orders)),
            "toeplitz_controls": controls}


def physical_case(coordinates, side):
    seed = box_kernel(coordinates, side)
    tensor = exterior_tensor(seed)
    report = physical_profiles(seed, tensor)
    expected_zero = subset_sums(seed["energies"] * 2)
    expected_infinity = subset_sums(seed["energies"] + [-value for value in seed["energies"]])
    assert report["profiles"] == [expected_zero, expected_infinity]
    n = seed["n"]
    energies = [sum(seed["energies"][j] for j in range(n) if symbol >> j & 1)
                for symbol in tensor["alphabet"]]
    endpoint_controls = 0
    for (b, x, v, y), value in tensor["values"].items():
        if energies[x] == energies[v]:
            assert value.bit_count() & 1 == int(b == y and x == v)
            endpoint_controls += 1
    return {"coordinates": coordinates, "side": side, "physical_letters": 1 << n,
            "complete_profile_orders": sum(map(len, report["profiles"])),
            "equal_energy_endpoint_controls": endpoint_controls}


def stochastic_controls():
    counts = 0
    for side in (2, 4, 8, 16, 32, 64):
        one = [Fraction(2 * j + 1 - side, side) for j in range(side)]
        two = [(a + b) / 2 for a, b in product(one, repeat=2)]
        for r in range(1, side // 2 + 1):
            threshold = Fraction(2 * r - 1, side)
            single = Fraction(sum(abs(value) <= threshold for value in one), side)
            averaged = Fraction(sum(abs(value) <= threshold for value in two), side * side)
            formula = Fraction((side - 2 * r) * (2 * r - 1), side * side)
            assert averaged - single == formula and formula >= 0
            counts += 1
    return counts


def coherent_controls():
    seed = box_kernel(2, 2)
    tensor = exterior_tensor(seed)
    choices = [{index: (index + 1) % 6} for index in range(5)]
    choices += [{index: choice for index in range(5)} for choice in range(6)]
    reports = []
    for code in choices:
        reduced = small_box_code(tensor, code)
        report = physical_profiles(seed, reduced)
        certificate = three_quarter_certificate(report["profiles"], report["total_heights"])
        reports.append({"choices": code, "physical_letters": len(reduced["alphabet"]),
                        "complete_profile_orders": sum(map(len, report["profiles"])),
                        "polynomial_degree": certificate["degree"],
                        "minimum_translated_coefficient": min(certificate["translated_coefficients"])})
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [one_particle_case(k, side) for k, side in
             ((1, 1), (1, 2), (1, 4), (1, 8), (2, 1), (2, 2),
              (2, 4), (2, 8), (4, 1), (4, 2), (8, 1))]
    physical = [physical_case(k, side) for k, side in ((1, 1), (1, 2), (1, 4), (2, 2))]
    coherent = coherent_controls()
    stochastic = stochastic_controls()
    sources = [Path(__file__), Path(__file__).with_name("tensor_pascal_box_crossings.py"),
               Path(__file__).with_name("cube_symmetric_exterior_blocks.py")]
    report = {"status": "pass", "one_particle_cases": cases, "physical_cases": physical,
              "coherent_cases": coherent, "stochastic_cdf_controls": stochastic,
              "factoring_score": "1/2",
              "source_sha256": {path.name: sha256(path.read_bytes()).hexdigest() for path in sources}}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(json.dumps({"status": "pass", "one_particle_cases": len(cases),
                      "physical_cases": len(physical), "coherent_cases": len(coherent),
                      "stochastic_cdf_controls": stochastic}, indent=2))


if __name__ == "__main__":
    main()
