#!/usr/bin/env python3
"""Independent normalized-average checks and exact quotient thresholds."""

from __future__ import annotations

import argparse
import json
from itertools import product
from math import factorial, prod
from pathlib import Path
from random import Random

from circuit_box_averages import (
    average_order, normalized_box_average, universal_average_degree,
    valuation_two,
)


def circuit_value(circuit, point, modulus):
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = point[instruction[1]]
        elif operation == "const":
            value = instruction[1]
        elif operation == "add":
            value = values[instruction[1]] + values[instruction[2]]
        elif operation == "inv":
            value = pow(values[instruction[1]], -1, modulus)
        else:
            value = values[instruction[1]] * values[instruction[2]]
        values.append(value % modulus)
    return values[circuit["output"]]


def literal_average(circuit, intervals, precision):
    dimensions = len(intervals)
    cardinality = prod(upper - lower for lower, upper in intervals)
    loss = valuation_two(cardinality)
    modulus = 1 << (precision + loss)
    total = sum(circuit_value(circuit, point, modulus)
                for point in product(*(range(lower, upper) for lower, upper in intervals)))
    numerator = (1 << dimensions) * total
    assert numerator % (1 << loss) == 0
    return (numerator >> loss) * pow(cardinality >> loss, -1, 1 << precision) % (1 << precision)


def falling_circuit(degree, extra_power=0):
    nodes = [["var", 0], ["const", 1]]
    output = 1
    for root in range(degree):
        nodes.append(["const", -root])
        nodes.append(["add", 0, len(nodes) - 1])
        nodes.append(["mul", output, len(nodes) - 1])
        output = len(nodes) - 1
    for _ in range(extra_power):
        nodes.append(["mul", output, 0])
        output = len(nodes) - 1
    return {"nodes": nodes, "output": output}


def controls():
    rng = Random(913286)
    thresholds = []
    minimality_checks = 0
    ideal_checks = 0
    for precision in range(1, 49):
        degree = universal_average_degree(precision)
        assert degree <= 8 * precision
        assert average_order(degree - 1) < precision
        assert all(average_order(index) >= precision for index in range(degree, 16 * precision))
        witness = falling_circuit(degree - 1)
        actual = normalized_box_average(witness, ((0, degree),), precision)["scaled_average_residue"]
        expected_valuation = average_order(degree - 1)
        assert actual and valuation_two(actual) == expected_valuation
        assert actual == literal_average(witness, ((0, degree),), precision)
        minimality_checks += 1
        thresholds.append({"precision": precision, "degree": degree,
                           "minimality_witness_order": expected_valuation})
        if precision <= 12:
            for extra in (0, 1, 5, 13):
                circuit = falling_circuit(degree, extra)
                interval = (0, 2 * degree + 3)
                assert literal_average(circuit, (interval,), precision) == 0
                assert normalized_box_average(circuit, (interval,), precision)["scaled_average_residue"] == 0
                ideal_checks += 1

    literal_checks = 0
    rational_checks = 0
    for dimensions, repetitions, max_precision in ((1, 40, 9), (2, 32, 7), (3, 16, 4)):
        for repeat in range(repetitions):
            nodes = [["var", index] for index in range(dimensions)]
            nodes += [["const", rng.randrange(-8, 9)] for _ in range(2)]
            output = 0
            for axis in range(1, dimensions):
                nodes.append(["add", output, axis])
                output = len(nodes) - 1
            for _ in range(5 + repeat % 7):
                nodes.append(["mul", output, output])
                output = len(nodes) - 1
            nodes.append(["add", output, rng.randrange(dimensions + 2)])
            circuit = {"nodes": nodes, "output": len(nodes) - 1}
            intervals = []
            for _ in range(dimensions):
                lower = rng.randrange(-4, 5)
                intervals.append((lower, lower + rng.randrange(1, 6)))
            precision = 1 + repeat % max_precision
            result = normalized_box_average(circuit, intervals, precision)
            assert result["scaled_average_residue"] == literal_average(circuit, intervals, precision)
            literal_checks += 1
            # 1+2*q is globally odd even when q is a coupled circuit.
            nodes.append(["const", 2])
            nodes.append(["mul", circuit["output"], len(nodes) - 1])
            nodes.append(["const", 1])
            nodes.append(["add", len(nodes) - 2, len(nodes) - 1])
            nodes.append(["inv", len(nodes) - 1])
            rational = {"nodes": nodes, "output": len(nodes) - 1}
            rational_result = normalized_box_average(rational, intervals, precision)
            assert rational_result["scaled_average_residue"] == literal_average(rational, intervals, precision)
            rational_checks += 1

    # A unit need not be coefficientwise one modulo two. This polynomial
    # is odd at every integer but has nilpotents in the quotient residue.
    nilpotent = {"nodes": [["var", 0], ["mul", 0, 0], ["add", 0, 1],
                             ["const", 1], ["add", 2, 3], ["inv", 4]], "output": 5}
    for precision in range(1, 13):
        result = normalized_box_average(nilpotent, ((-5, 19),), precision)
        assert result["scaled_average_residue"] == literal_average(nilpotent, ((-5, 19),), precision)
        rational_checks += 1
    rejected = False
    try:
        normalized_box_average({"nodes": [["var", 0], ["inv", 0]], "output": 1}, ((1, 2),), 4)
    except ValueError:
        rejected = True
    assert rejected

    # Every odd x has x^(2^500)=1 mod2^502, and every even x gives
    # zero to far greater precision. The interval has equally many
    # odd and even integers, so the scaled average is exactly 1 mod2^16.
    nodes = [["var", 0]]
    for _ in range(500):
        nodes.append(["mul", len(nodes) - 1, len(nodes) - 1])
    circuit = {"nodes": nodes, "output": len(nodes) - 1}
    length = 3 << 400
    large = normalized_box_average(circuit, ((7, 7 + length),), 16)
    assert large["scaled_average_residue"] == 1
    assert large["discarded_numerical_valuation"] == 400
    nodes.append(["const", 2])
    nodes.append(["mul", circuit["output"], len(nodes) - 1])
    nodes.append(["const", 1])
    nodes.append(["add", len(nodes) - 2, len(nodes) - 1])
    nodes.append(["inv", len(nodes) - 1])
    rational_circuit = {"nodes": nodes, "output": len(nodes) - 1}
    large_rational = normalized_box_average(rational_circuit, ((7, 7 + length),), 16)
    assert large_rational["scaled_average_residue"] == 4 * pow(3, -1, 1 << 16) % (1 << 16)

    # This coupled two-variable example removes 250 count-valuation
    # bits using an eight-bit quotient. Half the pairs have odd sum.
    nodes = [["var", 0], ["var", 1], ["add", 0, 1]]
    for _ in range(300):
        nodes.append(["mul", len(nodes) - 1, len(nodes) - 1])
    circuit = {"nodes": nodes, "output": len(nodes) - 1}
    first_length, second_length = 3 << 100, 5 << 150
    coupled = normalized_box_average(circuit,
                                      ((-7, -7 + first_length), (5, 5 + second_length)), 8)
    assert coupled["scaled_average_residue"] == 2
    assert coupled["discarded_numerical_valuation"] == 250
    return {"exact_universal_thresholds": thresholds,
            "sharp_degree_witnesses": minimality_checks,
            "literal_ideal_annihilation_checks": ideal_checks,
            "literal_normalized_circuit_checks": literal_checks,
            "literal_regular_rational_checks": rational_checks,
            "nonglobal_unit_rejected": rejected,
            "large_single_variable_control": large,
            "large_regular_rational_control": large_rational,
            "large_coupled_control": coupled,
            "scope": "Normalized regular-circuit averages with sharp universal falling-factorial degrees; no singular factor-count evaluator."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    arguments = parser.parse_args()
    result = controls()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        Path(arguments.output).write_text(text)
    print(text, end="")
