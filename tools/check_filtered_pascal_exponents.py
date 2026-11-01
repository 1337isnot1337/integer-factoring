#!/usr/bin/env python3
"""Exact classification and actual-cut controls for filtered Pascal exponents."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import random

from check_weighted_pascal_flags import (
    exact_weighted_profile, original_fock_profile, subset_sums,
)
from filtered_pascal_exponents import (
    decode_triangular_map, filtered_pascal_kernel, triangular_map, upset_orders,
)
from tensor_pascal_box_crossings import exterior_tensor, matrix_product, polynomial_power
from weighted_pascal_flags import clifford_data, inverse_kernel, natural_heights


def preserves_flag(images):
    length = len(images)
    return all(order is None or order >= length - 1 - q
               for q, order in enumerate(upset_orders(images)))


def is_classified(images):
    try:
        decode_triangular_map(images)
        return True
    except ValueError:
        return False


def complete_small_classification():
    result = []
    for length in (2, 4):
        accepted, total = 0, 0
        for images in product(range(length), repeat=length):
            actual, predicted = preserves_flag(images), is_classified(images)
            assert actual == predicted, images
            accepted += actual
            total += 1
        dimension = length.bit_length() - 1
        assert accepted == 1 << (length + dimension - 1)
        result.append({"length": length, "all_maps": total, "accepted": accepted})
    accepted = 0
    for images in permutations(range(8)):
        actual = preserves_flag(images)
        assert actual == is_classified(images), images
        accepted += actual
        if actual:
            assert upset_orders(images) == list(range(7, -1, -1))
            assert all(((images[i] - images[j]) & -(images[i] - images[j]))
                       == ((i - j) & -(i - j)) for i in range(8) for j in range(i))
    assert accepted == 128
    result.append({"length": 8, "all_permutations": 40320, "accepted": accepted})
    return result


def generated_classification():
    # All 2^(8+3-1) classified maps, including nonlinear noninjective ones.
    count, image_sizes = 0, Counter()
    for diagonal in product((0, 1), repeat=3):
        for truth_bits in product((0, 1), repeat=7):
            tables = [truth_bits[:1], truth_bits[1:3], truth_bits[3:7]]
            images = triangular_map(diagonal, tables)
            assert decode_triangular_map(images) == (list(diagonal), [list(x) for x in tables])
            assert preserves_flag(images)
            if all(diagonal):
                assert upset_orders(images) == list(range(7, -1, -1))
            image_sizes[len(set(images))] += 1
            count += 1
    assert count == 1024
    return {"length": 8, "classified_maps_checked": count,
            "image_size_counts": dict(sorted(image_sizes.items()))}


def check_kernel(images, full_fock=False):
    seed = filtered_pascal_kernel(images)
    length, matrix, pascal = seed["n"], seed["matrix"], seed["pascal"]
    first = [row[:length] for row in matrix[:length]]
    second = [row[length:] for row in matrix[length:]]
    target = [[polynomial_power(3, 2 * length) if value else 0 for value in row] for row in pascal]
    assert matrix_product(first, second) == target
    assert matrix_product(matrix, inverse_kernel(seed)) == [
        [int(i == j) for j in range(2 * length)] for i in range(2 * length)]
    for i, row in enumerate(second):
        for j, value in enumerate(row):
            order = (value & -value).bit_length() - 1 if value else None
            assert order is None or order >= 2 * (length - 1 - (i | j))
            if all(seed["diagonal"]):
                assert order == 2 * (length - 1 - (i | j))
    heights = natural_heights(seed)
    values, shifts, rows, columns = clifford_data(seed, 0)
    row, column = [heights[0][i] for i in rows], [heights[0][i] for i in columns]
    actual = exact_weighted_profile(values, shifts, row + [-x for x in row],
                                    column + [-x for x in column])
    expected = sorted((seed["zero_energies"] + [-x for x in seed["zero_energies"]]) * 2)
    assert actual == expected
    # The actual infinity profile is retained without imposing an energy-list claim.
    infinity = exact_weighted_profile(matrix, [[0] * (2 * length) for _ in matrix],
                                     heights[1][2 * length:], heights[1][:2 * length], True)
    record = {"length": length, "images": list(images), "diagonal": seed["diagonal"],
              "exponents": seed["exponents"], "zero_clifford_orders": actual,
              "actual_infinity_orders_at_reference_gauge": infinity,
              "original_fock_orders_checked": 0}
    if full_fock:
        tensor = exterior_tensor(seed)
        full_zero = original_fock_profile(seed, tensor, heights[0], 0)
        full_infinity = original_fock_profile(seed, tensor, heights[1], 1)
        assert full_zero == subset_sums(seed["zero_energies"] * 2)
        assert full_infinity == subset_sums(infinity)
        record["original_fock_orders_checked"] = len(full_zero) + len(full_infinity)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    classification = complete_small_classification()
    generated = generated_classification()
    samples = [([0, 0], True), ([0, 1], True), ([1, 0], True), ([1, 1], True),
               ([0, 0, 2, 2], True), ([0, 3, 2, 1], True),
               ([0, 1, 2, 7, 0, 1, 2, 7], False)]
    rng = random.Random(20261101)
    for dimension in (3, 3, 4, 4, 5):
        diagonal = [rng.randrange(2) for _ in range(dimension)]
        tables = [[rng.randrange(2) for _ in range(1 << j)] for j in range(dimension)]
        samples.append((triangular_map(diagonal, tables), False))
    kernels = [check_kernel(images, full_fock) for images, full_fock in samples]
    assert kernels[2]["actual_infinity_orders_at_reference_gauge"] == [-3, -3, 3, 3]
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), Path(__file__).with_name("filtered_pascal_exponents.py"),
               Path(__file__).with_name("check_weighted_pascal_flags.py"),
               Path(__file__).with_name("weighted_pascal_flags.py"),
               Path(__file__).with_name("tensor_pascal_box_crossings.py"),
               Path(__file__).with_name("cube_symmetric_exterior_blocks.py"),
               root / "research/filtered-pascal-exponents.md"]
    report = {"status": "pass", "complete_small_classification": classification,
              "generated_classification": generated, "kernel_cases": kernels,
              "factoring_score": "1/2",
              "source_sha256": {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest()
                                for path in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "pass", "all_maps_checked": 260,
                      "all_permutations_checked": 40320, "classified_maps_checked": 1024,
                      "kernel_cases": len(kernels), "original_fock_orders_checked":
                      sum(case["original_fock_orders_checked"] for case in kernels)}, indent=2))


if __name__ == "__main__":
    main()
