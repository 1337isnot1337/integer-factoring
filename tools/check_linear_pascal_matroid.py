#!/usr/bin/env python3
"""Exact projected-basis, zero-spectrum, and original-cut controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_weighted_pascal_flags import (
    exact_weighted_profile, original_fock_profile, subset_sums,
)
from linear_pascal_matroid import (
    binary_rank, linear_images, linear_pascal_kernel, projected_basis_weight,
)
from tensor_pascal_box_crossings import exterior_tensor, matrix_product, polynomial_power
from weighted_pascal_flags import clifford_data, inverse_kernel, natural_heights


def polynomial_order(value):
    return (value & -value).bit_length() - 1 if value else None


def check_rank_weights(rows, translations=True):
    dimension, length = len(rows), 1 << len(rows)
    weights = [projected_basis_weight(rows, mask) for mask in range(length)]
    for mask, value in enumerate(weights):
        candidates = [chosen for chosen in range(length)
                      if chosen.bit_count() == mask.bit_count()
                      and binary_rank([rows[j] & mask for j in range(dimension)
                                       if chosen >> j & 1]) == mask.bit_count()]
        assert value == min(candidates)
        for other in range(length):
            assert value + weights[other] <= weights[mask | other] + weights[mask & other]
            if mask & other == mask:
                assert value <= weights[other]
    images = linear_images(rows)
    count = 0
    for translation in range(length) if translations else (0,):
        for q in range(length):
            literal = 0
            for i in range(length):
                if i & q == q:
                    literal ^= polynomial_power(3, images[i] ^ translation)
            assert polynomial_order(literal) == weights[length - 1 - q]
            count += 1
    zero = [2 * length - 1 - 2 * weights[length - 1 - i] for i in range(length)]
    assert all(zero[i] + zero[length - 1 - i] >= 2 * length for i in range(length))
    return count


def complete_small_controls():
    result = []
    for dimension in (1, 2, 3):
        count, subcube_checks = 0, 0
        for rows in product(range(1 << dimension), repeat=dimension):
            if binary_rank(rows) != dimension:
                continue
            subcube_checks += check_rank_weights(rows)
            count += 1
        expected = 1
        for j in range(dimension):
            expected *= (1 << dimension) - (1 << j)
        assert count == expected
        result.append({"binary_dimension": dimension, "all_invertible_maps": count,
                       "literal_affine_subcube_orders_checked": subcube_checks})
    return result


def check_kernel(rows, full_fock=False):
    seed = linear_pascal_kernel(rows)
    n, matrix, pascal = seed["n"], seed["matrix"], seed["pascal"]
    first = [row[:n] for row in matrix[:n]]
    second = [row[n:] for row in matrix[n:]]
    target = [[polynomial_power(3, 2 * n) if value else 0 for value in row] for row in pascal]
    assert matrix_product(first, second) == target
    assert matrix_product(matrix, inverse_kernel(seed)) == [
        [int(i == j) for j in range(2 * n)] for i in range(2 * n)]
    for i, row in enumerate(second):
        for j, value in enumerate(row):
            predicted = 2 * seed["projected_basis_weights"][n - 1 - (i | j)]
            assert polynomial_order(value) == predicted
            assert predicted + seed["zero_energies"][i] + seed["zero_energies"][j] - 2 * n >= 0
    heights = natural_heights(seed)
    values, shifts, row_indices, column_indices = clifford_data(seed, 0)
    row = [heights[0][i] for i in row_indices]
    column = [heights[0][i] for i in column_indices]
    actual = exact_weighted_profile(values, shifts, row + [-x for x in row],
                                    column + [-x for x in column])
    expected = sorted((seed["zero_energies"] + [-x for x in seed["zero_energies"]]) * 2)
    assert actual == expected
    infinity = exact_weighted_profile(matrix, [[0] * (2 * n) for _ in matrix],
                                     heights[1][2 * n:], heights[1][:2 * n], True)
    record = {"binary_rows": list(rows), "length": n, "images": seed["images"],
              "projected_basis_weights": seed["projected_basis_weights"],
              "zero_energies": seed["zero_energies"], "zero_clifford_orders": actual,
              "candidate_infinity_exponents": seed["exponents"],
              "actual_infinity_orders_at_reference_gauge": infinity,
              "original_fock_orders_checked": 0}
    if full_fock:
        tensor = exterior_tensor(seed)
        actual_zero = original_fock_profile(seed, tensor, heights[0], 0)
        actual_infinity = original_fock_profile(seed, tensor, heights[1], 1)
        assert actual_zero == subset_sums(seed["zero_energies"] * 2)
        assert actual_infinity == subset_sums(infinity)
        record["original_fock_orders_checked"] = len(actual_zero) + len(actual_infinity)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    complete = complete_small_controls()
    rng = random.Random(202611010543)
    larger = []
    for dimension in (4, 4, 5, 5, 6):
        while True:
            rows = [rng.randrange(1 << dimension) for _ in range(dimension)]
            if binary_rank(rows) == dimension:
                break
        larger.append({"binary_rows": rows,
                       "literal_subcube_orders_checked": check_rank_weights(rows, False)})
    samples = [([1], True)]
    samples += [(list(rows), True) for rows in product(range(4), repeat=2)
                if binary_rank(rows) == 2]
    samples += [([3, 6, 4], False), ([7, 3, 5], False), ([4, 2, 1], False)]
    samples += [(case["binary_rows"], False) for case in larger[:3]]
    kernels = [check_kernel(rows, full_fock) for rows, full_fock in samples]
    example = next(case for case in kernels if case["binary_rows"] == [3, 2])
    assert example["zero_energies"] == [1, 5, 5, 7]
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), Path(__file__).with_name("linear_pascal_matroid.py"),
               Path(__file__).with_name("check_weighted_pascal_flags.py"),
               Path(__file__).with_name("weighted_pascal_flags.py"),
               Path(__file__).with_name("tensor_pascal_box_crossings.py"),
               Path(__file__).with_name("cube_symmetric_exterior_blocks.py"),
               root / "research/linear-pascal-matroid-energies.md"]
    report = {"status": "pass", "complete_small_controls": complete,
              "larger_rank_controls": larger, "kernel_cases": kernels,
              "factoring_score": "1/2",
              "source_sha256": {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest()
                                for path in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "all_invertible_small_maps": sum(
        case["all_invertible_maps"] for case in complete),
        "literal_subcube_orders_checked": sum(case["literal_affine_subcube_orders_checked"]
            for case in complete) + sum(case["literal_subcube_orders_checked"] for case in larger),
        "kernel_cases": len(kernels), "original_fock_orders_checked": sum(
            case["original_fock_orders_checked"] for case in kernels)}, indent=2))


if __name__ == "__main__":
    main()
