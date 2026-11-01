#!/usr/bin/env python3
"""Independent exact checks for count-independent Bezout precision."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import factorial
from pathlib import Path
from types import SimpleNamespace

from circuit_box_bezout import (
    Builder, IsometricContext, circuit_isometric_bezout,
    normalized_parity_difference,
)
from circuit_box_bezout_averages import (
    circuit_isometric_bezout_by_average, normalized_moment_differences,
    normalized_parity_difference_by_average,
)
from check_circuit_box_bezout import (
    circuit_examples, divided_difference, exact_value, fraction_residue,
    implicit_identity_checks, literal_pair,
)


def exact_node(context, index):
    point = index // 2
    if index % 2 == 0:
        value = (context.first_a + context.step_a
                 * (point + 2 * exact_value(context.circuit_a, point)))
    else:
        denominator = (context.first_b + context.step_b
                       * (point + 2 * exact_value(context.circuit_b, point)))
        value = Fraction(context.number, denominator)
    return (Fraction(value) - 1) / (1 << context.height)


def literal_weighted_ratio(context, order, degree):
    """2^j H_(k,j)/((k+1)...(k+j)) by exact rational series."""
    numerator = [Fraction(1)] + [Fraction(0)] * degree
    for index in range(order + degree):
        numerator = [numerator[j] - (index * numerator[j - 1] if j else 0)
                     for j in range(degree + 1)]
    inverse = [Fraction(1)] + [Fraction(0)] * degree
    for index in range(order + 1):
        node = exact_node(context, index)
        following = [inverse[0]]
        for j in range(1, degree + 1):
            following.append(inverse[j] + node * following[-1])
        inverse = following
    coefficient = sum(numerator[j] * inverse[degree - j] for j in range(degree + 1))
    denominator = 1
    for index in range(1, degree + 1):
        denominator *= order + index
    return (1 << degree) * coefficient / denominator


def controls():
    examples = circuit_examples()
    pairings = ((0, 0), (1, 2), (2, 3), (3, 4))
    moment_checks = 0
    ratio_checks = 0
    difference_checks = 0
    pair_checks = 0
    coefficient_checks = 0
    for first_index, second_index in pairings:
        context = IsometricContext(7, 1, 4, 5, 12, 5, 1,
                                    examples[first_index], examples[second_index])
        for length in range(1, 13):
            nodes = [exact_node(context, index) for index in range(length)]
            for precision in (3, 7):
                moments = normalized_moment_differences(context, length, 6, precision)
                for degree in range(7):
                    expected = sum((node ** degree - index ** degree
                                    for index, node in enumerate(nodes)), Fraction(0)) / length
                    assert expected.denominator % 2 == 1
                    assert moments[degree] == fraction_residue(expected, precision)
                    moment_checks += 1
        for order in range(1, 10):
            nodes = [exact_node(context, index) for index in range(order + 1)]
            difference = divided_difference(nodes, [index % 2 for index in range(order + 1)])
            exact_unit = (-1) ** (order + 1) * factorial(order) * difference / (1 << (order - 1))
            for precision in (1, 3, 7):
                result = normalized_parity_difference_by_average(context, order, precision)
                assert result["unit"] == fraction_residue(exact_unit, precision)
                reference = normalized_parity_difference(context, order, precision)
                assert result["unit"] == reference["unit"]
                difference_checks += 1
                if precision == 7:
                    for certificate in result["short_ratio_certificates"]:
                        weighted = literal_weighted_ratio(context, order, certificate["degree"])
                        assert certificate["weighted_contribution"] == fraction_residue(weighted, precision)
                        ratio_checks += 1
    for height in (1, 2):
        for count in (1, 3, 5):
            for first_index, second_index in pairings:
                parameters = (1 + (1 << height) + (1 << (height + 1)),
                              1, 1 << (height + 1),
                              1 + (1 << (height + 1)), 3 << (height + 1), count, height)
                first, second = examples[first_index], examples[second_index]
                context = IsometricContext(*parameters, first, second)
                exact_u, exact_v, _ = literal_pair(context)
                for precision in (3, 7):
                    result = circuit_isometric_bezout_by_average(*parameters, precision, first, second)
                    scale = 1 << result["bezout_exponent"]
                    for label, coefficients in (("u", exact_u), ("v", exact_v)):
                        literal = [fraction_residue(scale * value, precision)
                                   for value in reversed(coefficients)]
                        computed = result[label + "_reversed_coefficients"]
                        while len(literal) < count:
                            literal.insert(0, 0)
                        assert computed == literal[:len(computed)]
                        assert not any(literal[len(computed):])
                        coefficient_checks += len(literal)
                    pair_checks += 1

    builder = Builder()
    root = 0
    for _ in range(200):
        root = builder.operation("mul", root, root)
    perturbation = builder.finish(builder.operation("mul", root, 0))
    large_controls = []
    for count in ((1 << 80) + 17, 1 << 80):
        arguments = 7, 1, 4, 5, 12, count, 1, 6, perturbation, perturbation
        computed = circuit_isometric_bezout_by_average(*arguments)
        reference = circuit_isometric_bezout(*arguments)
        for label in ("u", "v"):
            assert computed[label + "_reversed_coefficients"] == reference[label + "_reversed_coefficients"]
        context = IsometricContext(7, 1, 4, 5, 12, count, 1, perturbation, perturbation)
        computed["independent_unit_point_checks"] = implicit_identity_checks(context, computed, [(0, 1), (1, 1)])
        computed["previous_working_precision"] = reference["maximum_working_precision"]
        large_controls.append(computed)

    # A closed-form huge-order reference has z_i = 5i+6. Scaling all
    # pairwise differences gives D_k=5^(-k), without any root expansion.
    even = {"nodes": [["var", 0], ["const", 10], ["mul", 0, 1],
                       ["const", 6], ["add", 2, 3]], "output": 4}
    odd = {"nodes": [["var", 0], ["const", 10], ["mul", 0, 1],
                      ["const", 11], ["add", 2, 3]], "output": 4}
    affine = SimpleNamespace(_regular_z_cache=(even, odd))
    order = (1 << 4096) - 1
    huge = normalized_parity_difference_by_average(affine, order, 8)
    assert huge["unit"] == pow(pow(5, -1, 256), order, 256)
    assert huge["working_precision"] == 13
    old_guard = max(order.bit_count() - (order + j).bit_count() for j in range(8))
    huge["previous_working_precision"] = 8 + old_guard
    huge["order_bit_length"] = order.bit_length()
    huge["closed_form_reference"] = "5^(-k) modulo 256 for z_i=5i+6"
    return {"normalized_moment_checks": moment_checks,
            "literal_weighted_ratio_checks": ratio_checks,
            "literal_divided_difference_checks": difference_checks,
            "literal_complete_pair_checks": pair_checks,
            "literal_pair_coefficients_checked": coefficient_checks,
            "large_complete_pair_controls": large_controls,
            "huge_order_affine_control": huge,
            "scope": "Count-independent coefficient precision on certified isometric domains; no arbitrary factor-root oracle."}


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
    print("Receipt SHA256:", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
