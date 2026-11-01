#!/usr/bin/env python3
"""Exact valuation-tree, pointwise compiler, and rational inverse controls."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random
from time import perf_counter

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_rational_weights import inverse_rational_weight_contraction
from rational_weight_valuation_partition import (
    polynomial_valuation_partition, polynomial_value, compile_partition_weight,
)

ROOT = Path(__file__).resolve().parents[1]


def evaluate_compiled(compiled, a):
    result = Quad()
    for leaf in compiled["leaf_weights"]:
        pole = Quad(*(Fraction(value) for value in leaf["pole"]))
        reciprocal = 1 / (a - pole)
        power, value = Quad(Fraction(1)), Quad()
        for coefficient in leaf["coefficients"]:
            power *= reciprocal
            value += Quad(*(Fraction(x) for x in coefficient)) * power
        result += value / (1 << leaf["denominator_bits"])
    return result


def controls():
    rng = random.Random(2026110117)
    partitions = identities = zeros = pointwise_checks = 0
    examples = []
    for index in range(500):
        polynomial = [rng.randrange(-12, 13) for _ in range(rng.randrange(1, 7))]
        if not any(polynomial):
            continue
        nodes = list(range(1, rng.randrange(5, 31, 2), 2))
        domain = ((1, len(nodes)),)
        true_zeros = [a for a in nodes if polynomial_value(polynomial, a) == 0]
        try:
            partition = polynomial_valuation_partition(domain, polynomial)
        except ValueError:
            assert true_zeros
            zeros += 1
            continue
        assert not true_zeros
        for a in nodes:
            leaves = [leaf for leaf in partition["leaves"]
                      if (a - leaf["residue"]) % (1 << leaf["bits"]) == 0]
            assert len(leaves) == 1
            leaf = leaves[0]
            value = polynomial_value(polynomial, a)
            assert (value & -value).bit_length() - 1 == leaf["valuation"]
            parameter = (a - leaf["residue"]) // (1 << leaf["bits"])
            assert value == (1 << leaf["valuation"]) * polynomial_value(leaf["unit_polynomial"], parameter)
            assert leaf["unit_polynomial"][0] % 2
            assert all(coefficient % (1 << i) == 0 for i, coefficient in enumerate(leaf["unit_polynomial"][1:], 1))
            identities += 1
        partitions += 1
        if index < 65:
            numerator = [rng.randrange(-4, 5) for _ in range(rng.randrange(1, 4))]
            precision = 3 + partition["maximum_denominator_valuation"]
            compiled = compile_partition_weight(partition, numerator, precision)
            for a in nodes[:4]:
                expected = Quad(Fraction(polynomial_value(numerator, a), polynomial_value(polynomial, a)))
                actual = evaluate_compiled(compiled, a)
                difference = actual - expected
                assert not difference or difference.valuation() >= 3
                pointwise_checks += 1
        if len(examples) < 8:
            examples.append({"denominator": polynomial, "node_count": len(nodes),
                             "value_bit_bound": partition["value_bit_bound"],
                             "processed_primitive_nodes": partition["processed_primitive_nodes"],
                             "leaves": len(partition["leaves"]),
                             "maximum_active_degree_budget": partition["maximum_active_degree_budget"]})
    regression = polynomial_valuation_partition(((1, 10),), [3, -2, 1])
    assert len(regression["leaves"]) == 2
    assert [leaf["bits"] for leaf in regression["leaves"]] == [2, 2]
    assert [leaf["unit_polynomial"] for leaf in regression["leaves"]] == [[1, 0, 8], [3, 8, 8]]
    dense_controls = 0
    denominator_catalogue = ([1], [2], [3, -2, 1], [11, -6, 1], [1, 0, 1], [8, -4, 1])
    for index in range(12):
        size = 1 + index % 3
        left, right = ((1, size),), ((3, size),)
        aa, bb = expand(left), expand(right)
        n, k = 15 + 2 * (index % 4), 6
        row_denominator = denominator_catalogue[index % len(denominator_catalogue)]
        column_denominator = denominator_catalogue[(index + 2) % len(denominator_catalogue)]
        row_numerator, column_numerator = [1, 1], [2, -1]
        normalized = index % 2 == 0
        actual = inverse_rational_weight_contraction(n, left, right, k,
            row_numerator, row_denominator, column_numerator, column_denominator, 2,
            normalize_inverse=normalized)
        dense = inverse(matrix_at(n, aa, bb, k))
        expected = sum((dense[j][i]
                        * Fraction(polynomial_value(row_numerator, a), polynomial_value(row_denominator, a))
                        * Fraction(polynomial_value(column_numerator, b), polynomial_value(column_denominator, b))
                        for i, a in enumerate(aa) for j, b in enumerate(bb)), Quad())
        expected *= 1 << actual["inverse_normalization_exponent"]
        assert tuple(actual["residue"]) == mod_pair(expected * (1 << (-actual["scale_exponent"])),
                                                    actual["working_precision"])
        assert actual["working_precision"] + actual["scale_exponent"] == 2
        dense_controls += 1
    return {"nonzero_denominator_partitions": partitions, "literal_leaf_identities": identities,
            "actual_domain_zero_rejections": zeros, "pointwise_rational_weight_checks": pointwise_checks,
            "unit_leaf_refinement_regression": regression,
            "dense_rational_inverse_contractions": dense_controls, "examples": examples}


def implicit_controls():
    size = (1 << 200) + 13
    deep = polynomial_valuation_partition(((1, size),), [(1 << 140) + 1, -2, 1])
    assert deep["domain_size"] == size
    assert deep["maximum_denominator_valuation"] == 141
    assert deep["maximum_active_degree_budget"] <= 2
    n, left, right = (1 << 401) + 45, ((1, size),), ((3, size),)
    k = max(abs(n - 3), abs(n - (2 * size - 1) * (2 * size + 1))).bit_length() + 12
    started = perf_counter()
    lower = inverse_rational_weight_contraction(n, left, right, k,
        [1, 1], [3, -2, 1], [2, -1], [11, -6, 1], 2, normalize_inverse=True)
    elapsed = perf_counter() - started
    higher = inverse_rational_weight_contraction(n, left, right, k,
        [1, 1], [3, -2, 1], [2, -1], [11, -6, 1], 4, normalize_inverse=True)
    common = min(lower["scale_exponent"], higher["scale_exponent"])
    modulus = 1 << (2 - common)
    assert [(x << (lower["scale_exponent"] - common)) % modulus for x in lower["residue"]] == [
            (x << (higher["scale_exponent"] - common)) % modulus for x in higher["residue"]]
    return {"deep_denominator_partition": deep, "input_bits": n.bit_length(), "n": n,
            "domain_size": size, "k": k, "lower_result": lower, "higher_result": higher,
            "lower_precision_seconds": elapsed, "enumerated_matrix_nodes": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "Complete implicit inverse contractions for arbitrary explicit separated rational weights",
              "controls": controls(), "factoring_exponent_improved": False}
    if not args.small_only:
        report["implicit_controls"] = implicit_controls()
    report["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "rational_weight_valuation_partition.py",
                     "pole_cauchy_rational_weights.py", "pole_cauchy_pole_polynomials.py")}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report["controls"].items() if isinstance(value, int)}, indent=2))


if __name__ == "__main__":
    main()
