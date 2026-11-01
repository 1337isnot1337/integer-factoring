#!/usr/bin/env python3
"""Finite controls for the explicit-rational ordinary-interval selector bound."""

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random

from rational_weight_valuation_partition import polynomial_valuation_partition, polynomial_value

ROOT = Path(__file__).resolve().parents[1]


def multiply_linear(coefficients, root):
    result = [0] * (len(coefficients) + 1)
    for i, value in enumerate(coefficients):
        result[i] -= root * value
        result[i + 1] += value
    return result


def controls():
    newton_patterns = step_patterns = degree_budget_controls = 0
    for half_size in (1, 2, 4, 8):
        for degree in range(9):
            for mask in range(1 << (degree + 1)):
                values = [sum(((mask >> j) & 1) * comb(x, j)
                              for j in range(min(x, degree) + 1)) % 2
                          for x in range(2 * half_size)]
                if values == [0] * half_size + [1] * half_size:
                    assert degree >= half_size
                    step_patterns += 1
                newton_patterns += 1
    rng = random.Random(2026110118)
    obstruction_examples = []
    for _ in range(150):
        bits = rng.randrange(1, 9)
        coefficients = [rng.randrange(-8, 9) for _ in range(rng.randrange(1, 7))]
        if not any(coefficients):
            continue
        domain = ((1, 1 << bits),)
        try:
            partition = polynomial_valuation_partition(domain, coefficients)
        except ValueError:
            continue
        degree = partition["active_degree_bound"]
        depth = min(leaf["bits"] for leaf in partition["leaves"])
        if depth <= bits:
            assert (1 << max(0, depth - 2)) <= degree
        else:
            assert (1 << (bits - 1)) <= degree
        degree_budget_controls += 1
        witnesses = [leaf for leaf in partition["leaves"]
                     if leaf["bits"] <= bits and leaf["valuation"] < (1 << (bits - leaf["bits"]))]
        if witnesses and len(obstruction_examples) < 8:
            leaf = min(witnesses, key=lambda value: value["bits"])
            obstruction_examples.append({"selector_bit": bits, "denominator": coefficients,
                "residue": leaf["residue"], "residue_bits": leaf["bits"],
                "denominator_valuation": leaf["valuation"],
                "required_local_polynomial_degree": 1 << (bits - leaf["bits"]),
                "scope": "No integer-polynomial numerator can give the upper-half parity indicator with this denominator."})
    sharp_regular_examples, selector_values = [], 0
    for bits in range(1, 9):
        degree = 1 << (bits - 1)
        denominator_exponent = degree + (degree - degree.bit_count())
        numerator = [1]
        for j in range(degree):
            numerator = multiply_linear(numerator, 2 * j + 1)
        denominator = 1 << denominator_exponent
        partition = polynomial_valuation_partition(((1, 1 << bits),), [denominator])
        for a in range(1, 1 << (bits + 1), 2):
            value = polynomial_value(numerator, a)
            assert value % denominator == 0
            assert (value // denominator) % 2 == int(a >= (1 << bits))
            selector_values += 1
        assert (1 << bits) <= 4 * partition["active_degree_bound"] * partition["value_bit_bound"]
        sharp_regular_examples.append({"selector_bit": bits, "numerator_degree": degree,
            "denominator_degree": 0, "denominator_power_of_two": denominator_exponent,
            "denominator_value_bit_bound": partition["value_bit_bound"],
            "leaf_residue_bits": partition["leaves"][0]["bits"],
            "leaf_valuation": partition["leaves"][0]["valuation"]})
    return {"integer_valued_polynomial_parity_patterns": newton_patterns,
            "actual_step_patterns": step_patterns,
            "primitive_degree_budget_depth_controls": degree_budget_controls,
            "explicit_upper_half_selector_values": selector_values,
            "upper_half_selector_families": sharp_regular_examples,
            "denominator_only_obstruction_examples": obstruction_examples}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "An ordinary upper-half parity selector requires 2^K <=4 max(1,degree P) H in an explicit rational representation",
              "controls": controls(), "factoring_lower_bound": False,
              "source_sha256": {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, "rational_weight_valuation_partition.py")},
              "scope": "A bound on a single explicitly written rational weight, not on factoring or arbitrary oracle algorithms."}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report["controls"].items() if isinstance(value, int)}, indent=2))


if __name__ == "__main__":
    main()
