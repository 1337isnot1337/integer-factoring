#!/usr/bin/env python3
"""Independent pointwise and dense checks of binary product-power weights."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random
from time import perf_counter

from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_product_power_weights import (
    product_power_partition, compile_normalized_product_power, inverse_product_power_contraction,
)
from rational_weight_valuation_partition import polynomial_value

ROOT = Path(__file__).resolve().parents[1]


def valuation(value):
    assert value
    return (value & -value).bit_length() - 1


def normalized_value(factors, node, minimum, precision):
    """Direct scalar evaluation; no partition or local power-series code."""
    order, units = 0, []
    modulus = 1 << precision
    for polynomial, exponent in factors:
        if exponent == 0:
            continue
        value = polynomial_value(polynomial, node)
        s = valuation(value)
        order += s * exponent
        units.append((value >> s, exponent))
    gap = order - minimum
    assert gap >= 0
    if gap >= precision:
        return 0
    result = 1 << gap
    for unit, exponent in units:
        result = result * pow(unit, exponent, modulus) % modulus
    return result


def compiled_value(compiled, node):
    value = Quad()
    for leaf in compiled["leaf_weights"]:
        pole = Quad(*(Fraction(x) for x in leaf["pole"]))
        reciprocal = 1 / (node - pole)
        power = Quad(Fraction(1))
        for coefficient in leaf["coefficients"]:
            power *= reciprocal
            value += Quad(*(Fraction(x) for x in coefficient)) * power
    return value


def exact_weight(factors, node):
    value = Fraction(1)
    for polynomial, exponent in factors:
        value *= Fraction(polynomial_value(polynomial, node)) ** exponent
    return value


def controls():
    rng = random.Random(2026110123)
    catalog = ([1], [2], [3], [-2], [1, 1], [3, 1], [5, -2, 1],
               [13, -6, 1], [3, -2, 1], [1, 0, 1], [8, -4, 1])
    partitions = identities = pointwise = dense = exact_scaled = rejected = extremizers = 0
    examples = []
    for index in range(80):
        first = rng.randrange(1, 8, 2)
        size = rng.randrange(1, 10)
        domain = ((first, size),) if index % 7 else ((first, 2), (first + 10, 3))
        factors = [(list(rng.choice(catalog)), rng.randrange(-4, 5))
                   for _ in range(rng.randrange(0, 5))]
        if index % 8 == 0:
            factors.append(([5, -2, 1], -((1 << (128 + index)) + 3)))
        partition = product_power_partition(domain, factors)
        literal_orders = []
        for node in expand(domain):
            order = sum(exponent * valuation(polynomial_value(polynomial, node))
                        for polynomial, exponent in factors if exponent)
            literal_orders.append(order)
            leaves = [leaf for leaf in partition["leaves"]
                      if (node - leaf["residue"]) % (1 << leaf["bits"]) == 0]
            assert len(leaves) == 1
            leaf = leaves[0]
            assert leaf["valuation"] == order
            parameter = (node - leaf["residue"]) // (1 << leaf["bits"])
            for record, unit, s in zip(partition["factors"], leaf["unit_polynomials"], leaf["factor_valuations"]):
                assert polynomial_value(record["polynomial"], node) == (1 << s) * polynomial_value(unit, parameter)
                identities += 1
        assert partition["minimum_valuation"] == min(literal_orders)
        assert partition["maximum_valuation"] == max(literal_orders)
        assert len(partition["leaves"]) <= partition["source_leaf_bound"]
        partitions += 1
        if index < 30:
            for precision in (1, 2, 4):
                compiled = compile_normalized_product_power(partition, precision)
                for node in expand(domain):
                    actual = compiled_value(compiled, node)
                    assert mod_pair(actual, precision) == (
                        normalized_value(factors, node, partition["minimum_valuation"], precision), 0)
                    pointwise += 1
        if len(examples) < 5:
            examples.append({"factors": factors, "domain": domain,
                             "refined_leaves": len(partition["leaves"]),
                             "source_leaf_bound": partition["source_leaf_bound"],
                             "minimum_valuation": partition["minimum_valuation"]})
    for index in range(30):
        size = 1 + index % 4
        left, right = ((1, size),), ((3, size),)
        if index % 10 == 9:
            left, right = ((1, 1), (7, 1)), ((3, 1), (11, 1))
        aa, bb = expand(left), expand(right)
        row = [(list(catalog[(2 * index + 6) % len(catalog)]), -1 - index % 3),
               ([2], index % 3 - 1)]
        column = [(list(catalog[(index + 7) % len(catalog)]), 1 + index % 4)]
        if index % 6 == 0:
            row.append(([5, -2, 1], -((1 << (129 + index)) + 5)))
            column.append(([13, -6, 1], (1 << (131 + index)) + 7))
        n, k, precision = 15 + 2 * (index % 4), 8, 1 + index % 4
        actual = inverse_product_power_contraction(n, left, right, k, row, column, precision)
        inv = inverse(matrix_at(n, aa, bb, k))
        modulus = 1 << precision
        expected = [0, 0]
        loss = actual["inverse_normalization_exponent"]
        for i, a in enumerate(aa):
            fa = normalized_value(row, a, actual["row_normalization_exponent"], precision)
            for j, b in enumerate(bb):
                gb = normalized_value(column, b, actual["column_normalization_exponent"], precision)
                entry = mod_pair(inv[j][i] * (1 << loss), precision)
                expected = [(x + fa * gb * y) % modulus for x, y in zip(expected, entry)]
        assert actual["residue"] == expected
        assert actual["maximum_pole_order"] <= 6 * precision - 3
        dense += 1
        if index % 6:
            original = sum((inv[j][i] * exact_weight(row, a) * exact_weight(column, b)
                            for i, a in enumerate(aa) for j, b in enumerate(bb)), Quad())
            scale = actual["scale_exponent"]
            normalized = original / (1 << scale) if scale >= 0 else original * (1 << -scale)
            assert mod_pair(normalized, precision) == tuple(actual["residue"])
            exact_scaled += 1
    for polynomial in catalog:
        nodes = tuple(range(1, 32, 2))
        maximum = max(valuation(polynomial_value(polynomial, a)) for a in nodes)
        for precision in range(1, 7):
            exponent = -(1 << (precision - 1))
            for a in nodes:
                actual = normalized_value([(polynomial, exponent)], a, maximum * exponent, precision)
                assert actual == int(valuation(polynomial_value(polynomial, a)) == maximum)
                extremizers += 1
    for factors in ([([-3, 1], -1)], [([-3, 1], 2)], [([0], 1)]):
        try:
            product_power_partition(((1, 4),), factors)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("A forbidden domain zero was accepted")
    for precision in (0, -1):
        try:
            compile_normalized_product_power(product_power_partition(((1, 2),), []), precision)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Nonpositive precision was accepted")
    return {"product_power_partitions": partitions, "literal_factor_leaf_identities": identities,
            "pointwise_normalized_checks": pointwise, "dense_normalized_inverse_contractions": dense,
            "exact_original_scaled_contractions": exact_scaled, "invalid_contract_rejections": rejected,
            "maximum_valuation_selector_checks": extremizers,
            "examples": examples}


def implicit_controls():
    size = (1 << 200) + 13
    n, left, right, k = (1 << 401) + 45, ((1, size),), ((3, size),), 414
    row = [([5, -2, 1], -((1 << 512) + 3))]
    column = [([13, -6, 1], (1 << 383) + 5)]
    started = perf_counter()
    lower = inverse_product_power_contraction(n, left, right, k, row, column, 5)
    elapsed = perf_counter() - started
    higher = inverse_product_power_contraction(n, left, right, k, row, column, 7)
    assert lower["scale_exponent"] == higher["scale_exponent"]
    assert lower["residue"] == [value % 32 for value in higher["residue"]]
    assert lower["residue"] != [0, 0]
    assert lower["row_skipped_valuation_leaves"] == 1
    assert lower["column_skipped_valuation_leaves"] == 2
    assert lower["pole_polynomial_contractions"] == 2
    return {"input_bits": n.bit_length(), "domain_size": size, "k": k,
            "row_power_bits": 513, "column_power_bits": 384,
            "lower_precision_seconds": elapsed, "lower_result": lower, "higher_result": higher,
            "enumerated_matrix_nodes": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"claim": "Scaled inverse contractions for binary product-power weights with nonunit bases",
              "controls": controls(), "factoring_exponent_improved": False,
              "base_commit": "c404b544a538454b9cca8b810ebb224b2c74a0ec"}
    if not args.small_only:
        report["implicit_controls"] = implicit_controls()
    report["source_sha256"] = {name: sha256((ROOT / "tools" / name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, "pole_cauchy_product_power_weights.py",
                     "rational_weight_valuation_partition.py", "pole_cauchy_pole_polynomials.py")}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report["controls"].items() if isinstance(value, int)}, indent=2))
    if not args.small_only:
        case = report["implicit_controls"]
        print(json.dumps({"implicit_residue": case["lower_result"]["residue"],
                          "higher_residue": case["higher_result"]["residue"],
                          "scale_encoding_bits": abs(case["lower_result"]["scale_exponent"]).bit_length(),
                          "lower_seconds": case["lower_precision_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
