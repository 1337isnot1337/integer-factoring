#!/usr/bin/env python3
"""Exact and huge-domain checks for complete univariate pole traces."""

from __future__ import annotations

import argparse
import json
import random
from fractions import Fraction
from pathlib import Path

from circuit_box_bezout import Builder, valuation_two
from circuit_box_pole_trace import (
    dyadic_valuation_partition, normalized_circuit_pole_trace,
    pointwise_circuit_polynomial,
)


def exact_value(circuit, point):
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = Fraction(point)
        elif operation == "const":
            value = Fraction(instruction[1])
        elif operation == "add":
            value = values[instruction[1]] + values[instruction[2]]
        elif operation == "mul":
            value = values[instruction[1]] * values[instruction[2]]
        else:
            assert operation == "inv"
            value = 1 / values[instruction[1]]
        values.append(value)
    return values[circuit["output"]]


def residue(value, precision):
    value = Fraction(value)
    assert value.denominator & 1
    return value.numerator * pow(value.denominator, -1, 1 << precision) % (1 << precision)


def polynomial_value(coefficients, point):
    value = 0
    for coefficient in reversed(coefficients):
        value = value * point + coefficient
    return value


def literal_trace(circuit, interval, target, scale_exponent, precision):
    scale = 1 << scale_exponent
    answer = [Fraction(0), Fraction(0)]
    roots = 0
    for point in range(*interval):
        value = exact_value(circuit, point) - target
        norm = value * value + scale * value + scale * scale
        answer[0] += -scale * (value + scale) / norm
        answer[1] += -scale * scale / norm
        roots += residue(value, scale_exponent) == 0 if scale_exponent else 1
    return [residue(value, precision) for value in answer], roots


def examples():
    result = [{"nodes": [["const", 0]], "output": 0},
              {"nodes": [["const", 12]], "output": 0},
              {"nodes": [["var", 0]], "output": 0},
              {"nodes": [["var", 0], ["mul", 0, 0]], "output": 1}]
    generator = random.Random(20261101)
    for _ in range(8):
        builder = Builder()
        indices = [0, builder.constant(generator.randrange(-7, 8))]
        for _ in range(5):
            first, second = generator.choice(indices), generator.choice(indices)
            indices.append(builder.operation(generator.choice(("add", "mul")), first, second))
        result.append(builder.finish(indices[-1]))
    for original in result[4:]:
        nodes = [list(node) for node in original["nodes"]]
        twice = len(nodes)
        nodes.append(["const", 2])
        product = len(nodes)
        nodes.append(["mul", original["output"], twice])
        one = len(nodes)
        nodes.append(["const", 1])
        denominator = len(nodes)
        nodes.append(["add", product, one])
        output = len(nodes)
        nodes.append(["inv", denominator])
        result.append({"nodes": nodes, "output": output})
    return result


def controls():
    circuits = examples()
    compression_checks = 0
    partition_checks = 0
    trace_checks = 0
    root_count_checks = 0
    for circuit in circuits:
        for precision in (1, 3, 9):
            compressed = pointwise_circuit_polynomial(circuit, 3, precision)
            for point in range(-7, 13):
                expected = residue(exact_value(circuit, point) - 3, precision)
                assert polynomial_value(compressed["coefficients"], point) % (1 << precision) == expected
                compression_checks += 1
        compressed = pointwise_circuit_polynomial(circuit, 3, 9)
        partition = dyadic_valuation_partition(compressed["coefficients"], 5, 9, (-7, 13))
        for point in range(-7, 13):
            containing = [leaf for leaf in partition["leaves"]
                          if (point - leaf["residue"]) % (1 << leaf["residue_exponent"]) == 0]
            assert len(containing) == 1
            leaf = containing[0]
            local = (point - leaf["residue"]) >> leaf["residue_exponent"]
            reconstructed = polynomial_value(leaf["coefficients"], local) << leaf["valuation"]
            actual = residue(exact_value(circuit, point) - 3, 9)
            assert reconstructed % 512 == actual
            actual_valuation = valuation_two(actual) if actual else 9
            assert min(actual_valuation, 5) == leaf["valuation"]
            partition_checks += 1
        for interval in ((-7, 11), (3, 15), (5, 5)):
            for scale_exponent in (0, 2, 7):
                for precision in (1, 5):
                    computed = normalized_circuit_pole_trace(circuit, interval, 3,
                                                              scale_exponent, precision)
                    expected, roots = literal_trace(circuit, interval, 3, scale_exponent, precision)
                    assert computed["normalized_trace"] == expected
                    assert computed["congruence_root_count"] == roots
                    trace_checks += 1
                    root_count_checks += 1

    # The K+p pointwise precision is necessary in the worst case. Values
    # 0 and 2^(K+p-1) agree one bit earlier but their scaled poles differ.
    sharp_controls = []
    for scale_exponent in (1, 7, 40):
        for precision in (1, 3, 8):
            first = {"nodes": [["const", 0]], "output": 0}
            second = {"nodes": [["const", 1 << (scale_exponent + precision - 1)]], "output": 0}
            a = normalized_circuit_pole_trace(first, (0, 1), 0, scale_exponent, precision)
            b = normalized_circuit_pole_trace(second, (0, 1), 0, scale_exponent, precision)
            assert a["normalized_trace"] != b["normalized_trace"]
            sharp_controls.append({"scale_exponent": scale_exponent, "precision": precision,
                                   "first_trace": a["normalized_trace"], "second_trace": b["normalized_trace"]})

    identity = {"nodes": [["var", 0]], "output": 0}
    linear = normalized_circuit_pole_trace(identity, (0, (1 << 160) + 13), 0, 80, 12)
    assert linear["normalized_trace"] == [4095, 4095]
    assert linear["congruence_root_count"] == (1 << 80) + 1

    builder = Builder()
    output = 0
    for _ in range(500):
        output = builder.operation("mul", output, output)
    enormous = builder.finish(output)
    interval = 7, 7 + 3 * (1 << 300) + 17
    even_count = (interval[1] + 1) // 2 - (interval[0] + 1) // 2
    large = normalized_circuit_pole_trace(enormous, interval, 0, 80, 12)
    assert large["normalized_trace"] == [(-even_count) % 4096] * 2
    assert large["congruence_root_count"] == even_count

    nodes = [list(node) for node in enormous["nodes"]]
    two = len(nodes); nodes.append(["const", 2])
    doubled = len(nodes); nodes.append(["mul", output, two])
    one = len(nodes); nodes.append(["const", 1])
    denominator = len(nodes); nodes.append(["add", one, doubled])
    inverse = len(nodes); nodes.append(["inv", denominator])
    rational = normalized_circuit_pole_trace({"nodes": nodes, "output": inverse}, interval, 1, 80, 12)
    assert rational["normalized_trace"] == [(-even_count) % 4096] * 2
    assert rational["congruence_root_count"] == even_count

    # A bounded one-variable candidate curve yields an exact zero count.
    builder = Builder()
    first = builder.affine(0, 2, 1)
    second = builder.affine(0, 1, 2)
    product = builder.finish(builder.operation("mul", first, second))
    curve = normalized_circuit_pole_trace(product, (0, 21), 77, 15, 5)
    assert curve["normalized_trace"] == [31, 31]
    assert curve["congruence_root_count"] == 1
    assert curve["root_balls"][0]["first_index"] == 5
    curve["verified_factors"] = [7, 11]
    return {"pointwise_compression_checks": compression_checks,
            "literal_valuation_partition_checks": partition_checks,
            "literal_complete_trace_checks": trace_checks,
            "literal_congruence_root_count_checks": root_count_checks,
            "sharp_precision_controls": sharp_controls,
            "large_linear_control": linear,
            "large_unrestricted_degree_control": large,
            "large_regular_rational_control": rational,
            "bounded_candidate_curve_control": curve,
            "scope": "Complete one-variable singular normalized trace; arbitrary two-variable factor rectangles remain outside the theorem."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    receipt = controls()
    payload = json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
