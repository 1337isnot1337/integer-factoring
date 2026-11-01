#!/usr/bin/env python3
"""Literal rational and implicit-domain controls for circuit Bezout pairs."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

from circuit_box_bezout import (
    Builder, IsometricContext, circuit_isometric_bezout, circuit_value,
    odd_factorial_unit, product_jet, valuation_two,
)


def exact_value(circuit, point):
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = point
        elif operation == "const":
            value = instruction[1]
        elif operation == "add":
            value = values[instruction[1]] + values[instruction[2]]
        else:
            value = values[instruction[1]] * values[instruction[2]]
        values.append(value)
    return values[circuit["output"]]


def polynomial_add(first, second):
    result = [Fraction(0)] * max(len(first), len(second))
    for i, value in enumerate(first):
        result[i] += value
    for i, value in enumerate(second):
        result[i] += value
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result


def polynomial_multiply(first, second):
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result


def polynomial_divide(numerator, denominator):
    numerator = list(numerator)
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        offset = len(numerator) - len(denominator)
        value = numerator[-1] / denominator[-1]
        quotient[offset] = value
        for index, coefficient in enumerate(denominator):
            numerator[offset + index] -= value * coefficient
        while len(numerator) > 1 and not numerator[-1]:
            numerator.pop()
    assert not any(numerator)
    return quotient


def fraction_residue(value, precision):
    modulus = 1 << precision
    assert value.denominator % 2
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def divided_difference(nodes, values):
    work = [Fraction(value) for value in values]
    for step in range(1, len(nodes)):
        work = [(work[i + 1] - work[i]) / (nodes[i + step] - nodes[i])
                for i in range(len(work) - 1)]
    return work[0]


def literal_pair(context):
    """Independent rational Lagrange interpolation, followed by division."""
    count = context.count
    first_nodes = [context.first_a + context.step_a *
                   (i + 2 * exact_value(context.circuit_a, i)) for i in range(count)]
    second_nodes = [context.first_b + context.step_b *
                    (i + 2 * exact_value(context.circuit_b, i)) for i in range(count)]
    shifted = []
    for first, second in zip(first_nodes, second_nodes):
        shifted += [Fraction(first - 1), Fraction(context.number, second) - 1]
    interpolant = [Fraction(0)]
    for index in range(1, 2 * count, 2):
        numerator = [Fraction(1)]
        denominator = Fraction(1)
        for other, root in enumerate(shifted):
            if other != index:
                numerator = polynomial_multiply(numerator, [-root, Fraction(1)])
                denominator *= shifted[index] - root
        interpolant = polynomial_add(interpolant, [value / denominator for value in numerator])
    p = [Fraction(1)]
    q = [Fraction(1)]
    for first, second in zip(first_nodes, second_nodes):
        p = polynomial_multiply(p, [Fraction(1 - first), Fraction(1)])
        q = polynomial_multiply(q, [Fraction(context.number - second), Fraction(-second)])
    u = polynomial_divide(interpolant, p)
    v = polynomial_divide(polynomial_add([Fraction(1)], [-value for value in interpolant]), q)
    identity = polynomial_add(polynomial_multiply(u, p), polynomial_multiply(v, q))
    assert identity == [Fraction(1)]
    return u, v, [value / (1 << context.height) for value in shifted]


def circuit_examples():
    zero = {"nodes": [["const", 0]], "output": 0}
    variable = {"nodes": [["var", 0]], "output": 0}
    square = {"nodes": [["var", 0], ["mul", 0, 0]], "output": 1}
    cube_plus = {"nodes": [["var", 0], ["mul", 0, 0], ["mul", 0, 1],
                            ["const", 1], ["add", 2, 3]], "output": 4}
    negative_square = {"nodes": [["var", 0], ["mul", 0, 0], ["const", -1],
                                  ["mul", 1, 2]], "output": 3}
    return zero, variable, square, cube_plus, negative_square


def literal_product_jet(circuit, interval, precision, degree):
    result = [1] + [0] * degree
    modulus = 1 << precision
    for index in range(*interval):
        value = exact_value(circuit, index)
        result = [(result[j] - (value * result[j - 1] if j else 0)) % modulus
                  for j in range(degree + 1)]
    return result


def residue_bucket_factorial(number, precision):
    """Independent finite-precision reference; enumerates 2^(p-1) residues."""
    modulus = 1 << precision
    period = modulus // 2
    result = 1
    while number:
        count = (number + 1) // 2
        for index in range(period):
            multiplicity = (count - 1 - index) // period + 1 if index < count else 0
            result = result * pow(2 * index + 1, multiplicity, modulus) % modulus
        number //= 2
    return result


class QuadraticRing:
    def __init__(self, precision):
        self.modulus = 1 << precision

    def add(self, first, second):
        return tuple((a + b) % self.modulus for a, b in zip(first, second))

    def scale(self, value, scalar):
        return tuple(scalar * a % self.modulus for a in value)

    def multiply(self, first, second):
        a, b = first
        c, d = second
        return (a * c - b * d) % self.modulus, (a * d + b * c - b * d) % self.modulus

    def power(self, value, exponent):
        result = (1, 0)
        while exponent:
            if exponent & 1:
                result = self.multiply(result, value)
            exponent >>= 1
            if exponent:
                value = self.multiply(value, value)
        return result


def implicit_identity_checks(context, result, points):
    precision = result["parameters"]["precision"]
    modulus = 1 << precision
    quadratic = QuadraticRing(precision)
    count = context.count
    evaluated = []
    for point in points:
        p_value = q_value = (1, 0)
        for index in range(modulus):
            multiplicity = (count - 1 - index) // modulus + 1 if index < count else 0
            first = circuit_value(context.a_circuit, index, modulus)
            second = circuit_value(context.b_circuit, index, modulus)
            p_factor = quadratic.add(point, (-first, 0))
            q_factor = quadratic.add((context.number, 0), quadratic.scale(point, -second))
            p_value = quadratic.multiply(p_value, quadratic.power(p_factor, multiplicity))
            q_value = quadratic.multiply(q_value, quadratic.power(q_factor, multiplicity))
        values = []
        for label in ("u", "v"):
            value = (0, 0)
            shifted = quadratic.add(point, (-1, 0))
            for offset, coefficient in enumerate(result[label + "_reversed_coefficients"]):
                value = quadratic.add(value, quadratic.scale(
                    quadratic.power(shifted, count - 1 - offset), coefficient))
            values.append(value)
        total = quadratic.add(quadratic.multiply(values[0], p_value), quadratic.multiply(values[1], q_value))
        assert total == (pow(2, result["bezout_exponent"], modulus), 0)
        assert all(a % 2 or b % 2 for a, b in values)
        evaluated.append({"point": point, "u": values[0], "v": values[1],
                          "p": p_value, "q": q_value})
    return evaluated


def controls():
    examples = circuit_examples()
    jet_checks = 0
    for circuit in examples:
        for precision in (1, 3, 6):
            for degree in (0, 2, 7):
                computed = product_jet(circuit, (1, 8), precision, degree)
                assert computed == literal_product_jet(circuit, (1, 8), precision, degree)
                jet_checks += 1

    factorial_checks = 0
    for number in range(0, 48, 3):
        exact = factorial(number)
        unit = exact >> valuation_two(exact)
        for precision in (1, 3, 7):
            assert odd_factorial_unit(number, precision) == unit % (1 << precision)
            factorial_checks += 1

    rational_pairs = 0
    complete_pair_checks = 0
    coefficient_checks = 0
    difference_checks = 0
    isometry_checks = 0
    exact_scale_checks = 0
    pairings = ((0, 0), (1, 2), (2, 1), (2, 3), (3, 4), (4, 2))
    for height in (1, 2):
        for count in range(1, 6):
            for first_index, second_index in pairings:
                first, second = examples[first_index], examples[second_index]
                number = 1 + (1 << height) + (1 << (height + 1))
                first_a, first_b = 1, 1 + (1 << (height + 1))
                step_a, step_b = 1 << (height + 1), 3 << (height + 1)
                context = IsometricContext(number, first_a, step_a, first_b, step_b,
                                            count, height, first, second)
                u, v, z_nodes = literal_pair(context)
                rational_pairs += 1
                for i in range(len(z_nodes)):
                    for j in range(i):
                        difference = z_nodes[i] - z_nodes[j]
                        actual = valuation_two(difference.numerator) - valuation_two(difference.denominator)
                        assert actual == valuation_two(i - j)
                        isometry_checks += 1
                exponent = height * (2 * count - 1) - (count - 1).bit_count()
                actual_scale = max(0, max(valuation_two(value.denominator)
                                          - valuation_two(value.numerator)
                                          for value in u + v if value))
                assert exponent == actual_scale
                exact_scale_checks += 1
                for precision in (1, 3, 7):
                    result = circuit_isometric_bezout(number, first_a, step_a, first_b, step_b,
                                                       count, height, precision, first, second)
                    assert result["bezout_exponent"] == exponent
                    for label, polynomial in (("u", u), ("v", v)):
                        computed = result[label + "_reversed_coefficients"]
                        padded = computed + [0] * (count - len(computed))
                        expected = [fraction_residue(value * (1 << exponent), precision)
                                    for value in reversed(polynomial)]
                        assert padded == expected
                        coefficient_checks += count
                    for certificate in result["retained_newton_coefficients"]:
                        order = certificate["order"]
                        exact_difference = divided_difference(z_nodes[:order + 1],
                                                               [i % 2 for i in range(order + 1)])
                        normalized = exact_difference * Fraction((-1) ** (order + 1) * factorial(order),
                                                                  1 << (order - 1))
                        assert normalized.numerator % 2 and normalized.denominator % 2
                        assert certificate["unit"] == fraction_residue(normalized, certificate["precision"])
                        difference_checks += 1
                    complete_pair_checks += 1

    # A high-degree map on more than 2^80 nodes is actually compiled.
    # Only a 64-residue reference is used to verify the complete polynomial
    # identity at unit points; it never enumerates the original domains.
    builder = Builder()
    high_output = builder.power(0, (1 << 200) + 1)
    high = builder.finish(high_output)
    count = (1 << 80) + 17
    parameters = (7, 1, 4, 5, 12, count, 1)
    result = circuit_isometric_bezout(*parameters, 6, high, high)
    lower_precision = circuit_isometric_bezout(*parameters, 4, high, high)
    for label in ("u", "v"):
        first = result[label + "_reversed_coefficients"]
        second = lower_precision[label + "_reversed_coefficients"]
        assert [value % 16 for value in first] == second + [0] * (len(first) - len(second))
    factorial_reference = residue_bucket_factorial(2 * count - 1, 6)
    assert odd_factorial_unit(2 * count - 1, 6) == factorial_reference
    context = IsometricContext(*parameters, high, high)
    modulus = 64
    evaluated = implicit_identity_checks(context, result,
                                         ((a, b) for a in range(4) for b in (1, 3)))
    assert result["maximum_working_precision"] < 200
    # A count at an exact power of two forces the large factorial-ratio
    # precision loss. This directly exercises the cost that the first
    # large control's convenient low bits do not require.
    guard_parameters = (7, 1, 4, 5, 12, 1 << 80, 1)
    guard_result = circuit_isometric_bezout(*guard_parameters, 6, high, high)
    guard_lower = circuit_isometric_bezout(*guard_parameters, 4, high, high)
    for label in ("u", "v"):
        first = guard_result[label + "_reversed_coefficients"]
        second = guard_lower[label + "_reversed_coefficients"]
        assert [value % 16 for value in first] == second + [0] * (len(first) - len(second))
    guard_context = IsometricContext(*guard_parameters, high, high)
    guard_values = implicit_identity_checks(guard_context, guard_result, ((0, 1), (1, 1)))
    assert guard_result["maximum_working_precision"] > 80
    return {"reverse_product_jet_checks": jet_checks,
            "literal_factorial_unit_checks": factorial_checks,
            "literal_rational_pairs": rational_pairs,
            "literal_isometry_distances": isometry_checks,
            "exact_optimal_scale_checks": exact_scale_checks,
            "complete_normalized_pair_checks": complete_pair_checks,
            "normalized_coefficient_checks": coefficient_checks,
            "literal_normalized_divided_difference_checks": difference_checks,
            "large_control": {"result": result,
                              "input_circuit_sha256": hashlib.sha256(json.dumps(high, sort_keys=True).encode()).hexdigest(),
                              "lower_precision_consistency": True,
                              "independent_factorial_unit": factorial_reference,
                              "unit_point_global_identity_checks": evaluated,
                              "reference_residue_classes": modulus},
            "large_guard_control": {"result": guard_result,
                                    "lower_precision_consistency": True,
                                    "unit_point_global_identity_checks": guard_values},
            "scope": "Complete normalized Bezout compilation on separated circuit-isometric domains; the unrestricted factor-root oracle remains open."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    arguments = parser.parse_args()
    result = controls()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        Path(arguments.output).write_text(encoded)
    print(encoded, end="")
