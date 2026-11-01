#!/usr/bin/env python3
"""Independent finite controls for unrestricted-degree circuit box products."""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from math import comb, prod
from pathlib import Path
from random import Random

from circuit_box_products import (
    IntegerFunctionQuotient, _elementaries_mod, _weighted_packet_mod, circuit_box_moments,
    circuit_box_product, circuit_box_unit_trace, evaluate_product_circuit,
    factorial_valuation,
)


def literal_value(circuit, point, modulus):
    values = []
    for instruction in circuit["nodes"]:
        operation = instruction[0]
        if operation == "var":
            value = point[instruction[1]]
        elif operation == "const":
            value = instruction[1]
        elif operation == "add":
            value = values[instruction[1]] + values[instruction[2]]
        else:
            assert operation == "mul"
            value = values[instruction[1]] * values[instruction[2]]
        values.append(value % modulus)
    return values[circuit["output"]]


def power_node(nodes, index, exponent):
    nodes.append(["const", 1])
    result = len(nodes) - 1
    while exponent:
        if exponent & 1:
            nodes.append(["mul", result, index])
            result = len(nodes) - 1
        exponent >>= 1
        if exponent:
            nodes.append(["mul", index, index])
            index = len(nodes) - 1
    return result


class PolynomialRing:
    def __init__(self, precision):
        self.precision = precision
        self.modulus = 1 << precision

    def clean(self, values):
        values = [value % self.modulus for value in values]
        while len(values) > 1 and values[-1] == 0:
            values.pop()
        return tuple(values)

    def constant(self, value):
        return (value % self.modulus,)

    def add(self, first, second):
        values = [0] * max(len(first), len(second))
        for i, value in enumerate(first):
            values[i] += value
        for i, value in enumerate(second):
            values[i] += value
        return self.clean(values)

    def multiply(self, first, second):
        values = [0] * (len(first) + len(second) - 1)
        for i, a in enumerate(first):
            for j, b in enumerate(second):
                values[i + j] += a * b
        return self.clean(values)


class QuadraticRing:
    """Exact residues in (Z/2^p)[w]/(w^2+w+1)."""

    def __init__(self, precision):
        self.precision = precision
        self.modulus = 1 << precision

    def constant(self, value):
        return value % self.modulus, 0

    def add(self, first, second):
        return tuple((a + b) % self.modulus for a, b in zip(first, second))

    def scale(self, value, coefficient):
        return tuple(a * coefficient % self.modulus for a in value)

    def multiply(self, first, second):
        a, b = first
        c, d = second
        return (a * c - b * d) % self.modulus, (a * d + b * c - b * d) % self.modulus

    def inverse(self, value):
        a, b = value
        norm = (a * a - a * b + b * b) % self.modulus
        if not norm & 1:
            raise ValueError("The requested quadratic element is not a unit")
        inverse = pow(norm, -1, self.modulus)
        return ((a - b) * inverse) % self.modulus, (-b * inverse) % self.modulus


def exact_elementaries(values, degree):
    result = [1] + [0] * degree
    for value in values:
        for j in range(degree, 0, -1):
            result[j] += value * result[j - 1]
    return result


def controls():
    rng = Random(5162763)
    null_checks = 0
    for precision in range(1, 21):
        ring = IntegerFunctionQuotient(1, precision)
        assert factorial_valuation(ring.degree) >= precision
        assert factorial_valuation(ring.degree - 1) < precision
        for point in range(-35, 65):
            value = sum(coefficient * point ** exponent
                        for exponent, coefficient in enumerate(ring.null_polynomial))
            assert value % ring.modulus == 0
            null_checks += 1

    quotient_product_checks = 0
    for variables in (1, 2, 3):
        for precision in range(1, 7):
            ring = IntegerFunctionQuotient(variables, precision)
            first = {tuple(rng.randrange(ring.degree) for _ in range(variables)):
                     rng.randrange(ring.modulus) for _ in range(12)}
            second = {tuple(rng.randrange(ring.degree) for _ in range(variables)):
                      rng.randrange(ring.modulus) for _ in range(12)}
            reduced = ring.multiply(first, second)
            assert all(all(exponent < ring.degree for exponent in key) for key in reduced)
            for _ in range(12):
                point = tuple(rng.randrange(-8, 9) for _ in range(variables))

                def value(polynomial):
                    return sum(coefficient * prod(x ** exponent for x, exponent in zip(point, key))
                               for key, coefficient in polynomial.items()) % ring.modulus

                assert value(reduced) == value(first) * value(second) % ring.modulus
                quotient_product_checks += 1

    newton_checks = 0
    weighted_newton_checks = 0
    for precision in range(1, 13):
        for count in range(18):
            degree = min(count, precision - 1)
            working = precision + factorial_valuation(degree)
            roots = [rng.randrange(-150, 151) for _ in range(count)]
            moments = [sum(value ** j for value in roots) % (1 << working)
                       for j in range(degree + 1)]
            computed = _elementaries_mod(moments, degree, working, precision)
            exact = exact_elementaries(roots, degree)
            assert computed == tuple(value % (1 << precision) for value in exact)
            newton_checks += 1
            reduced_precision = max(1, precision - 1)
            reduced_moments = [sum(value ** j for value in roots) % (1 << reduced_precision)
                               for j in range(degree + 1)]
            weighted = _weighted_packet_mod(reduced_moments, degree, reduced_precision, precision)
            assert weighted == tuple(((-2) ** j * value) % (1 << precision)
                                      for j, value in enumerate(exact))
            weighted_newton_checks += 1

    dense_checks = 0
    moment_checks = 0
    trace_checks = 0
    examples = []
    for variables, maximum_precision, repetitions in ((1, 8, 36), (2, 6, 28), (3, 3, 12)):
        for repeat in range(repetitions):
            precision = 1 + repeat % maximum_precision
            nodes = [["var", i] for i in range(variables)]
            nodes += [["const", rng.randrange(-9, 10)] for _ in range(3)]
            for _ in range(8):
                operation = rng.choice(("add", "mul"))
                nodes.append([operation, rng.randrange(len(nodes)), rng.randrange(len(nodes))])
            random_output = len(nodes) - 1
            # Retain every input variable in a positive nonlinear branch,
            # rather than letting random dead nodes make most cases linear.
            base = 0
            for axis in range(1, variables):
                nodes.append(["add", base, axis])
                base = len(nodes) - 1
            nonlinear = power_node(nodes, base, 5 + 2 * (repeat % 4))
            nodes.append(["add", random_output, nonlinear])
            # A large formal exponent is still independently evaluated by
            # modular powering at the finitely many literal test points.
            output = power_node(nodes, len(nodes) - 1, (1 << (12 + repeat % 7)) + 3)
            nodes.append(["add", output, repeat % variables])
            circuit = {"nodes": nodes, "output": len(nodes) - 1}
            intervals = []
            for _ in range(variables):
                lower = rng.randrange(3)
                intervals.append((lower, lower + rng.randrange(1, 4)))
            if repeat == 0:
                intervals[0] = (0, 0)
            emitted = circuit_box_product(circuit, intervals, precision)
            poly_ring = PolynomialRing(precision)
            actual = evaluate_product_circuit(emitted, poly_ring, (0, 1))
            literal = (1,)
            points = tuple(product(*(range(a, b) for a, b in intervals)))
            for point in points:
                root = literal_value(circuit, point, 1 << precision)
                literal = poly_ring.multiply(literal, ((-root) % (1 << precision), 1))
            assert actual == literal
            assert len(emitted["nodes"]) <= 10 * precision + 8 * max(1, len(points).bit_length()) + 20
            dense_checks += 1

            working = emitted["working_precision"]
            degree = min(len(points), precision - 1)
            moment_data = circuit_box_moments(circuit, intervals, degree, working)
            direct = [[0] * (degree + 1), [0] * (degree + 1)]
            counts = [0, 0]
            for point in points:
                value = literal_value(circuit, point, 1 << (working + 1))
                parity = value % 2
                normalized = (value - parity) // 2
                counts[parity] += 1
                for exponent in range(degree + 1):
                    direct[parity][exponent] = (direct[parity][exponent]
                                               + pow(normalized, exponent, 1 << working)) % (1 << working)
            assert moment_data["counts"] == counts and moment_data["moments"] == direct
            moment_checks += sum(len(row) for row in direct)

            quadratic = QuadraticRing(precision)
            trace = circuit_box_unit_trace(circuit, intervals, quadratic, (0, 1))
            direct_trace = (0, 0)
            for point in points:
                root = literal_value(circuit, point, 1 << precision)
                direct_trace = quadratic.add(direct_trace, quadratic.inverse((-root % quadratic.modulus, 1)))
            assert trace == direct_trace
            trace_checks += 1
            if repeat == repetitions - 1:
                examples.append({"variables": variables, "precision": precision,
                                 "formal_degree_bound": emitted["formal_degree_bound"],
                                 "input_nodes": len(nodes), "output_nodes": len(emitted["nodes"]),
                                 "product_degree": emitted["input_degree"]})

    # Both numerical box lengths are over 2^100. The independent reference
    # visits only residue classes modulo a small, explicit power of two.
    nodes = [["var", 0], ["var", 1], ["const", 3], ["mul", 1, 2], ["add", 0, 3]]
    first = power_node(nodes, 4, (1 << 81) + 17)
    nodes.append(["mul", 0, 1])
    second = power_node(nodes, len(nodes) - 1, (1 << 65) + 3)
    nodes.append(["add", first, second])
    circuit = {"nodes": nodes, "output": len(nodes) - 1}
    intervals = ((3, (1 << 101) + 10), (5, (1 << 103) + 18))
    precision = 4
    emitted = circuit_box_product(circuit, intervals, precision)
    working = emitted["working_precision"]
    period = 1 << (working + 1)
    histogram = [0] * period

    def residue_count(interval, residue):
        lower, upper = interval
        return (upper - 1 - residue) // period - (lower - 1 - residue) // period

    for first_residue in range(period):
        first_count = residue_count(intervals[0], first_residue)
        for second_residue in range(period):
            count = first_count * residue_count(intervals[1], second_residue)
            value = literal_value(circuit, (first_residue, second_residue), period)
            histogram[value] += count
    assert sum(histogram) == emitted["input_degree"]
    for packet in emitted["packets"]:
        parity = packet["parity"]
        degree = packet["retained_degree"]
        elementaries = [1] + [0] * degree
        count = 0
        for value in range(parity, period, 2):
            multiplicity = histogram[value]
            count += multiplicity
            normalized = (value - parity) // 2
            factor = [comb(multiplicity, j) * normalized ** j % (1 << precision)
                      for j in range(degree + 1)]
            elementaries = [sum(elementaries[j - k] * factor[k] for k in range(j + 1))
                           % (1 << precision) for j in range(degree + 1)]
        expected = [((-2) ** j * value) % (1 << precision)
                    for j, value in enumerate(elementaries)]
        assert packet["count"] == count and packet["coefficients"] == expected

    # The unit evaluator must reject a singular scalar center. A small
    # modular product does not authorize division by its lost content.
    singular_rejected = False
    try:
        circuit_box_unit_trace({"nodes": [["var", 0]], "output": 0},
                               ((0, 2),), QuadraticRing(4), (1, 0))
    except ValueError:
        singular_rejected = True
    assert singular_rejected

    return {"null_polynomial_integer_checks": null_checks,
            "quotient_product_point_checks": quotient_product_checks,
            "modular_newton_checks": newton_checks,
            "weighted_packet_newton_checks": weighted_newton_checks,
            "full_polynomial_coefficient_checks": dense_checks,
            "independent_normalized_moment_checks": moment_checks,
            "literal_quadratic_unit_trace_checks": trace_checks,
            "small_examples": examples,
            "large_box_control": {
                "intervals": [list(interval) for interval in intervals],
                "precision": precision, "working_precision": working,
                "input_degree": emitted["input_degree"],
                "formal_degree_bound": emitted["formal_degree_bound"],
                "input_circuit_nodes": len(circuit["nodes"]),
                "output_circuit_nodes": len(emitted["nodes"]),
                "quotient_degree": emitted["quotient_degree"],
                "quotient_dimension": emitted["quotient_dimension"],
                "reference_residue_pairs_enumerated": period ** 2,
                "packets": emitted["packets"],
                "circuit_sha256": hashlib.sha256(json.dumps(circuit, sort_keys=True).encode()).hexdigest(),
            },
            "singular_center_rejected": singular_rejected,
            "claim_scope": "Global modular products and unit traces; no singular normalization or factoring exponent improvement."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    arguments = parser.parse_args()
    result = controls()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        Path(arguments.output).write_text(text)
    print(text, end="")
