#!/usr/bin/env python3
"""Exact-inverse checks for polynomial-circuit inverse moments."""

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, factorial
from pathlib import Path
from random import Random

from pole_cauchy_circuit_moments import (compressed_inverse_context,
                                               inverse_circuit_moment,
                                               monomial_circuit)


def qadd(x, y):
    return x[0] + y[0], x[1] + y[1]


def qnegative(x):
    return -x[0], -x[1]


def qmultiply(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def qinverse(x):
    a, b = x
    norm = a * a - a * b + b * b
    return (a - b) / norm, -b / norm


def exact_matrix_inverse(matrix):
    size = len(matrix)
    work = [row[:] + [(Fraction(int(i == j)), Fraction(0))
                      for j in range(size)] for i, row in enumerate(matrix)]
    for j in range(size):
        pivot = next(i for i in range(j, size) if work[i][j] != (0, 0))
        work[j], work[pivot] = work[pivot], work[j]
        inverse = qinverse(work[j][j])
        work[j] = [qmultiply(inverse, value) for value in work[j]]
        for i in range(size):
            if i != j:
                coefficient = work[i][j]
                work[i] = [qadd(a, qnegative(qmultiply(coefficient, b)))
                           for a, b in zip(work[i], work[j])]
    return [row[size:] for row in work]


def node_list(domain):
    return [first + 2 * i for first, length in domain for i in range(length)]


def normalized_inverse(number, left, right, displacement):
    rows, columns = node_list(left), node_list(right)
    scale = 1 << displacement
    matrix = []
    for a in rows:
        row = []
        for b in columns:
            delta = number - a * b
            norm = delta * delta - delta * scale + scale * scale
            row.append((Fraction(scale * scale, norm), Fraction(scale * delta, norm)))
        matrix.append(row)
    inverse = exact_matrix_inverse(matrix)
    return [[(scale * a, scale * b) for a, b in row] for row in inverse]


def v2(value):
    assert value > 0
    return (value & -value).bit_length() - 1


def circuit_at_nodes(circuit, row, column, modulus):
    values = []
    for node in circuit["nodes"]:
        operation = node[0]
        if operation == "row":
            value = row
        elif operation == "column":
            value = column
        elif operation == "const":
            value = node[1]
        elif operation == "add":
            value = values[node[1]] + values[node[2]]
        else:
            assert operation == "mul"
            value = values[node[1]] * values[node[2]]
        values.append(value % modulus)
    return values[circuit["output"]]


def literal_moment(inverse, left, right, circuit, precision):
    rows, columns = node_list(left), node_list(right)
    denominator_bits = max(v2(coordinate.denominator)
                           for row in inverse for pair in row for coordinate in pair)
    modulus = 1 << (precision + denominator_bits)
    total = [0, 0]
    for i, b in enumerate(columns):
        for j, a in enumerate(rows):
            weight = circuit_at_nodes(circuit, a, b, modulus)
            for coordinate in (0, 1):
                value = inverse[i][j][coordinate] * (1 << denominator_bits)
                assert value.denominator % 2
                residue = value.numerator * pow(value.denominator, -1, modulus) % modulus
                total[coordinate] = (total[coordinate] + residue * weight) % modulus
    assert all(value % (1 << denominator_bits) == 0 for value in total)
    return [value >> denominator_bits for value in total]


def coupled_power_circuit(exponent):
    # (a*b+3*a-2*b+5)^exponent, written as ordinary circuit operations.
    nodes = [["row"], ["column"], ["const", 3], ["const", -2], ["const", 5],
             ["mul", 0, 1], ["mul", 2, 0], ["mul", 3, 1], ["add", 5, 6],
             ["add", 8, 7], ["add", 9, 4], ["const", 1]]
    base, result = 10, 11
    while exponent:
        if exponent & 1:
            nodes.append(["mul", result, base])
            result = len(nodes) - 1
        exponent >>= 1
        if exponent:
            nodes.append(["mul", base, base])
            base = len(nodes) - 1
    return {"nodes": nodes, "output": result}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def check():
    rng = Random(20261101)
    comparisons = compressed_cases = 0
    for index in range(48):
        size = rng.randrange(1, 9)
        if index % 2:
            left = tuple((x, 1) for x in sorted(rng.sample(range(1, 42, 2), size)))
            right = tuple((x, 1) for x in sorted(rng.sample(range(1, 42, 2), size)))
        else:
            left = ((2 * rng.randrange(6) + 1, size),)
            right = ((2 * rng.randrange(6) + 1, size),)
        number = 2 * rng.randrange(1, 256) + 1
        differences = [abs(number - a * b) for a in node_list(left) for b in node_list(right)]
        displacement = max(1, max(differences).bit_length()) + 1
        precision = rng.randrange(1, 9)
        context = compressed_inverse_context(number, left, right, displacement, precision)
        compressed_cases += context["quotient_degree"] < size
        inverse = normalized_inverse(number, left, right, displacement)
        circuits = [monomial_circuit(index % 17, (3 * index) % 19),
                    monomial_circuit((1 << 80) + 13, (1 << 67) + 29),
                    coupled_power_circuit((1 << 13) + index + 1)]
        for circuit in circuits:
            actual = inverse_circuit_moment(context, circuit)["residue"]
            expected = literal_moment(inverse, left, right, circuit, precision)
            assert actual == expected, (index, precision, actual, expected)
            comparisons += 1

    # Absolute-precision example with a quotient smaller than the matrix.
    number, left, right, displacement, precision = 105, ((1, 16),), ((3, 16),), 10, 14
    context = compressed_inverse_context(number, left, right, displacement, precision)
    assert context["quotient_degree"] == 13 < 16
    circuit = coupled_power_circuit((1 << 31) + 7)
    actual = inverse_circuit_moment(context, circuit)
    inverse = normalized_inverse(number, left, right, displacement)
    assert actual["residue"] == literal_moment(inverse, left, right, circuit, precision)
    absolute_control = {"number": number, "matrix_dimension": 16,
                        "precision": precision, "absolute_error_bits": precision - displacement,
                        "quotient_degree": 13, "circuit_sha256": digest(circuit),
                        "result": actual}

    # Complete odd residue cycles cancel from the origin/infinity ratios.
    # Since k>=P, multiplication by z after inversion permutes odd integer
    # residues modulo 2^P. The remaining three nodes can be checked by a
    # literal 3x3 inverse. No long interval is enumerated.
    number = (1 << 201) + 123
    dimension = (1 << 100) + 3
    left, right = ((1, dimension),), ((5, dimension),)
    displacement, precision = 211, 16
    context = compressed_inverse_context(number, left, right, displacement, precision)
    circuit = monomial_circuit((1 << 128) + 17, (1 << 144) + 51)
    actual = inverse_circuit_moment(context, circuit)
    small_left, small_right = ((1, 3),), ((5, 3),)
    inverse = normalized_inverse(number, small_left, small_right, displacement)
    expected = literal_moment(inverse, small_left, small_right, circuit, precision)
    assert actual["residue"] == expected
    large_control = {"number_bits": number.bit_length(), "matrix_dimension": dimension,
                     "row_exponent_bits": 129, "column_exponent_bits": 145,
                     "context_sha256": digest(context), "circuit_sha256": digest(circuit),
                     "independent_remainder_dimension": 3, "result": actual}

    # Sharp denominator bounds, including literal cardinal polynomials.
    cardinal_checks = 0
    denominator_rows = []
    for s in range(0, 10):
        length = 1 << s
        degree = length - 1
        bound = degree + degree - degree.bit_count()
        assert bound == 2 * length - 2 - s
        for selected in range(length):
            assert comb(degree, selected) % 2
            denominator = (1 << degree) * factorial(selected) * factorial(degree - selected)
            assert v2(denominator) == bound
            if length <= 16:
                coefficients = [Fraction(1)]
                chosen = 3 + 2 * selected
                signed_denominator = 1
                for j in range(length):
                    if j == selected:
                        continue
                    node = 3 + 2 * j
                    following = [Fraction(0)] * (len(coefficients) + 1)
                    for power, value in enumerate(coefficients):
                        following[power] -= node * value
                        following[power + 1] += value
                    coefficients = following
                    signed_denominator *= chosen - node
                coefficients = [value / signed_denominator for value in coefficients]
                assert max(v2(value.denominator) for value in coefficients) == bound
                for j in range(length):
                    value = sum(coefficient * (3 + 2 * j) ** power
                                for power, coefficient in enumerate(coefficients))
                    assert value == int(j == selected)
                    cardinal_checks += 1
        denominator_rows.append({"odd_row_count": length, "minimum_denominator_bits": bound})

    difference_checks = 0
    for _ in range(128):
        degree = rng.randrange(25)
        order = rng.randrange(1, 17)
        start = rng.randrange(-8, 9)
        coefficients = [rng.randrange(-100, 101) for _ in range(degree + 1)]
        values = [sum(coefficient * (start + 2 * j) ** power
                      for power, coefficient in enumerate(coefficients))
                  for j in range(order + 1)]
        difference = sum((-1) ** (order - j) * comb(order, j) * values[j]
                         for j in range(order + 1))
        assert difference % ((1 << order) * factorial(order)) == 0
        difference_checks += 1

    huge_length = 1 << 200
    huge_bound = 2 * huge_length - 202
    assert huge_bound == (huge_length - 1) * 2 - (huge_length - 1).bit_count()
    return {"schema": 1, "scope": "Polynomial-circuit inverse moments; no factoring-bound improvement",
            "exact_matrix_instances": 48, "genuinely_compressed_instances": compressed_cases,
            "literal_inverse_circuit_comparisons": comparisons,
            "absolute_precision_control": absolute_control,
            "implicit_complete_cycle_control": large_control,
            "cardinal_polynomial_value_checks": cardinal_checks,
            "integer_step_two_difference_checks": difference_checks,
            "sharp_denominator_rows": denominator_rows,
            "binary_large_denominator_control": {"odd_row_count": huge_length,
                                                 "minimum_denominator_bits": huge_bound}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(encoded)
    print(encoded)
