#!/usr/bin/env python3
"""Dense transform controls and readback of explicit signed projections."""

from fractions import Fraction
import json
from pathlib import Path
import random

from bell_projection import project, quadratic_transform, split


def main():
    rng = random.Random(118709)
    dense_checks = 0
    for seed in range(20261031, 20261039):
        v, d, _ = split(9, 4, seed, 8)
        for vectors in (v, d):
            size = len(vectors)
            root = 1 << ((size.bit_length() - 1) // 2)
            gaussian = sum((-1) ** (x.bit_count() // 2 % 2) for x in vectors)
            values = tuple(rng.randrange(3) for _ in vectors)
            dense = tuple(gaussian // root * pow(root, -1, 3) * sum(
                (-1) ** ((x ^ y).bit_count() // 2 % 2) * z
                for y, z in zip(vectors, values)) % 3 for x in vectors)
            assert quadratic_transform(values, vectors) == dense
            assert quadratic_transform(dense, vectors) == values
            assert quadratic_transform((1,) * size, vectors) == (1,) * size
            dense_checks += 3
    root = Path(__file__).resolve().parents[1]
    receipt = json.loads((root / "evidence/bell-signed-projection.json").read_text())
    for example in receipt["examples"]:
        rebuilt = project(example["physical_bits"], example["boundary_bits"],
                          example["seed"], example["mixes"], example["boundary"],
                          example.get("weights"))
        for key, actual in rebuilt.items():
            assert json.loads(json.dumps(actual)) == example[key], key
        assert sum(example["valuation_histogram"]) == 1 << (2 * example["physical_bits"])
        common = example["common_divisor"]
        a = example["removed_initial_order"]
        assert example["degree"] == example["ambient_degree"] - a - (len(common) - 1)
        assert sum(common) % 3 == 1
        x = Fraction(*example["scale"])
        h = example["degree"]
        exact = x ** h * sum(Fraction(count) / x ** i
                            for i, count in enumerate(example["valuation_histogram"])) ** 2
        assert exact == Fraction(*example["pressure"])
        assert example["strict_target"] == (exact ** 577 < 2 ** (4 * example["physical_bits"] * 408))
    assert receipt["examples"][0]["removed_initial_order"] == 1
    assert receipt["examples"][0]["degree"] == receipt["examples"][0]["ambient_degree"] - 2
    print(f"Passed {dense_checks} independent dense transform controls and two full signed-projection receipts.")


if __name__ == "__main__":
    main()
