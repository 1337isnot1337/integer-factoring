#!/usr/bin/env python3
"""Independent exact controls for global bilinear modular packets."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb, factorial
from pathlib import Path
from random import Random

from bilinear_pole_packets import (bilinear_parity_moments,
                                   bilinear_product_circuit, unit_pole_trace)
from interval_packet_circuits import evaluate_circuit
from interval_unit_jets import QuadraticJets


class DensePolynomials:
    def __init__(self, precision, degree):
        self.modulus, self.size = 1 << precision, degree + 1

    def constant(self, value):
        return (value % self.modulus,) + (0,) * (self.size - 1)

    def argument(self):
        return (0, 1) + (0,) * (self.size - 2) if self.size > 1 else (0,)

    def add(self, left, right):
        return tuple((a + b) % self.modulus for a, b in zip(left, right))

    def multiply(self, left, right):
        result = [0] * self.size
        nz_left = [(i, a) for i, a in enumerate(left) if a]
        nz_right = [(j, b) for j, b in enumerate(right) if b]
        for i, a in nz_left:
            for j, b in nz_right:
                if i + j < self.size:
                    result[i + j] += a * b
        return tuple(value % self.modulus for value in result)


def roots_in_box(coefficients, interval_u, interval_v):
    a, b, c, d = coefficients
    return [a * u * v + b * u + c * v + d
            for u in range(*interval_u) for v in range(*interval_v)]


def dense_product(roots, modulus):
    result = [1]
    for root in roots:
        following = [0] * (len(result) + 1)
        for i, value in enumerate(result):
            following[i] -= root * value
            following[i + 1] += value
        result = [value % modulus for value in following]
    return tuple(result)


def unit_trace_from_moments(moments, ring, argument):
    result = ring.constant(0)
    for parity in (0, 1):
        if not moments[parity][0]:
            continue
        inverse = ring.inverse(ring.add(argument, ring.constant(-parity)))
        for j, moment in enumerate(moments[parity]):
            result = ring.add(result, ring.scale(ring.power(inverse, j + 1),
                                                 (1 << j) * moment))
    return result


def logarithmic_circuit_jet(circuit, degree):
    """Differentiate the polynomial via one extra Taylor coefficient."""
    large = QuadraticJets(circuit["precision"], degree + 1)
    values = evaluate_circuit(circuit, large, large.argument())
    ring = QuadraticJets(circuit["precision"], degree)
    derivative = tuple(((j + 1) * values[j + 1][0] % ring.modulus,
                        (j + 1) * values[j + 1][1] % ring.modulus)
                       for j in range(degree + 1))
    return ring.multiply(derivative, ring.inverse(values[:degree + 1]))


def power_sums_stirling(interval, degree):
    """Independent exact basis: sum t^k=sum S(k,j)j!binomial(t,j+1)."""
    lower, upper = interval
    row = [1]
    result = []
    for k in range(degree + 1):
        result.append(sum(row[j] * factorial(j) *
                          (comb(upper, j + 1) - comb(lower, j + 1))
                          for j in range(k + 1)))
        following = [0] * (k + 2)
        for j in range(k + 2):
            following[j] = ((j * row[j] if j <= k else 0)
                            + (row[j - 1] if j else 0))
        row = following
    return tuple(result)


def odd_power_sums_stirling(interval, degree):
    lower, upper = interval
    # Odd t=1+2x has lower ceil((lower-1)/2) and upper ceil((upper-1)/2).
    powers = power_sums_stirling((lower // 2, upper // 2), degree)
    return tuple(sum(comb(k, j) * (1 << j) * powers[j]
                     for j in range(k + 1)) for k in range(degree + 1))


def product_parity_moments_separable(interval_u, interval_v, degree):
    all_u, all_v = (power_sums_stirling(interval, degree)
                    for interval in (interval_u, interval_v))
    odd_u, odd_v = (odd_power_sums_stirling(interval, degree)
                    for interval in (interval_u, interval_v))
    total = [a * b for a, b in zip(all_u, all_v)]
    odd = [a * b for a, b in zip(odd_u, odd_v)]
    even_moments, odd_moments = [], []
    for j in range(degree + 1):
        numerator_even = total[j] - odd[j]
        numerator_odd = sum(comb(j, k) * (-1) ** (j - k) * odd[k]
                            for k in range(j + 1))
        assert numerator_even % (1 << j) == 0
        assert numerator_odd % (1 << j) == 0
        even_moments.append(numerator_even // (1 << j))
        odd_moments.append(numerator_odd // (1 << j))
    return tuple(even_moments), tuple(odd_moments)


def circuit_digest(circuit):
    encoded = json.dumps(circuit, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def check():
    rng = Random(20261101)
    cases = [((0, 0, 0, 0), (0, 4), (0, 3)),
             ((0, 0, 0, -7), (1, 4), (2, 5)),
             ((1, 0, 0, 0), (0, 0), (1, 7)),
             ((1, 0, 0, 0), (1, 7), (3, 3))]
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                for d in (0, 1):
                    cases.append(((a, b, c, d), (0, 5), (1, 5)))
    for _ in range(280):
        coefficients = tuple(rng.randrange(-5, 6) for _ in range(4))
        lo_u, lo_v = rng.randrange(6), rng.randrange(6)
        cases.append((coefficients, (lo_u, lo_u + rng.randrange(6)),
                      (lo_v, lo_v + rng.randrange(6))))

    moment_checks = coefficient_checks = jet_checks = 0
    for index, (coefficients, interval_u, interval_v) in enumerate(cases):
        roots = roots_in_box(coefficients, interval_u, interval_v)
        moments = bilinear_parity_moments(coefficients, interval_u, interval_v, 9)
        for parity in (0, 1):
            reduced = [(root - parity) // 2 for root in roots if root % 2 == parity]
            expected = tuple(sum(value ** j for value in reduced) for j in range(10))
            assert moments[parity] == expected
            moment_checks += 1
        for precision in (1, 3, 6, 9):
            circuit = bilinear_product_circuit(coefficients, interval_u, interval_v,
                                               precision)
            ring = DensePolynomials(precision, len(roots))
            actual = evaluate_circuit(circuit, ring, ring.argument())
            assert actual == dense_product(roots, ring.modulus)
            assert all(instruction[0] in ("x", "const", "add", "mul")
                       for instruction in circuit["nodes"])
            assert all(0 <= instruction[1] < ring.modulus
                       for instruction in circuit["nodes"] if instruction[0] == "const")
            coefficient_checks += 1
        if index < 64:
            for precision in (1, 3, 8, 12):
                circuit = bilinear_product_circuit(coefficients, interval_u, interval_v,
                                                   precision)
                for degree in (0, 2):
                    ring = QuadraticJets(precision, degree)
                    argument = ring.argument()
                    expected = ring.constant(0)
                    for root in roots:
                        expected = ring.add(expected, ring.inverse(
                            ring.add(argument, ring.constant(-root))))
                    actual = unit_pole_trace(coefficients, interval_u, interval_v,
                                             ring, argument)
                    assert actual == expected
                    assert logarithmic_circuit_jet(circuit, degree) == expected
                    jet_checks += 1

    # An ordinary integer argument cannot invert both residue classes.
    ring = QuadraticJets(8, 0)
    try:
        unit_pole_trace((1, 0, 0, 0), (1, 4), (1, 4), ring, ring.constant(1))
    except ValueError:
        rejected_singular_argument = True
    else:
        raise AssertionError("The singular pole argument was accepted as a unit")

    large_cases = [((1, (1 << 200) + 12), (1, (1 << 180) + 7), 32),
                   (((1 << 120) + 3, (1 << 120) + (1 << 80) + 5),
                    ((1 << 170) + 17, (1 << 170) + (1 << 90) + 1), 40),
                   ((0, (1 << 144) + 1), ((1 << 99) + 6, (1 << 100) + 9), 48)]
    large_receipts = []
    for interval_u, interval_v, precision in large_cases:
        moments = bilinear_parity_moments((1, 0, 0, 0), interval_u, interval_v,
                                         precision - 1)
        independent = product_parity_moments_separable(interval_u, interval_v,
                                                      precision - 1)
        assert moments == independent
        circuit = bilinear_product_circuit((1, 0, 0, 0), interval_u, interval_v,
                                           precision)
        ring = QuadraticJets(precision, 3)
        expected = unit_trace_from_moments(independent, ring, ring.argument())
        assert logarithmic_circuit_jet(circuit, 3) == expected
        large_receipts.append({"interval_u": list(interval_u),
                               "interval_v": list(interval_v),
                               "precision": precision,
                               "input_degree": circuit["input_degree"],
                               "nodes": len(circuit["nodes"]),
                               "circuit_sha256": circuit_digest(circuit),
                               "unit_trace_jet": [list(pair) for pair in expected]})

    # General bilinear large input: q=3(u+alpha)(v+beta)+gamma.
    interval_u = (3, (1 << 128) + 11)
    interval_v = (5, (1 << 96) + 3)
    alpha, beta, gamma, scale = (1 << 32) + 7, 5, -17, 3
    coefficients = (scale, scale * beta, scale * alpha, scale * alpha * beta + gamma)
    precision = 32
    circuit = bilinear_product_circuit(coefficients, interval_u, interval_v, precision)
    shifted_u = tuple(value + alpha for value in interval_u)
    shifted_v = tuple(value + beta for value in interval_v)
    product_circuit = bilinear_product_circuit((1, 0, 0, 0), shifted_u, shifted_v,
                                              precision)
    ring = QuadraticJets(precision, 3)
    transformed_argument = ring.scale(ring.add(ring.argument(), ring.constant(-gamma)),
                                     pow(scale, -1, ring.modulus))
    transformed = evaluate_circuit(product_circuit, ring, transformed_argument)
    transformed = ring.scale(transformed, pow(scale, circuit["input_degree"], ring.modulus))
    assert evaluate_circuit(circuit, ring, ring.argument()) == transformed
    transformed_receipt = {"coefficients": list(coefficients),
                           "interval_u": list(interval_u), "interval_v": list(interval_v),
                           "precision": precision, "input_degree": circuit["input_degree"],
                           "nodes": len(circuit["nodes"]),
                           "circuit_sha256": circuit_digest(circuit)}

    # Same zero raw reduction, different exact normalized pole counts.
    singular_controls = []
    interval = (1, 17)
    roots = roots_in_box((1, 0, 0, 0), interval, interval)
    precision = 32
    circuit = bilinear_product_circuit((1, 0, 0, 0), interval, interval, precision)
    ring = QuadraticJets(precision, 0)
    for number, expected_count in ((289, 0), (225, 1)):
        differences = [number - root for root in roots]
        count = differences.count(0)
        assert count == expected_count
        k = precision + max(abs(value) for value in differences).bit_length()
        valuation = count * k
        for difference in differences:
            if difference:
                value = abs(difference)
                valuation += (value & -value).bit_length() - 1
        assert valuation >= precision
        assert evaluate_circuit(circuit, ring, ring.constant(number)) == ring.constant(0)
        singular_controls.append({"number": number, "factor_pairs": count,
                                  "filter_displacement_exponent": k,
                                  "exact_product_valuation": valuation,
                                  "raw_modulus_precision": precision})

    return {"schema": 1, "status": "exact controls passed; no factoring-bound improvement",
            "small_bilinear_cases": len(cases), "parity_moment_comparisons": moment_checks,
            "full_coefficient_comparisons": coefficient_checks,
            "unit_trace_and_log_derivative_jet_comparisons": jet_checks,
            "singular_argument_rejected": rejected_singular_argument,
            "large_separable_controls": large_receipts,
            "large_affine_product_control": transformed_receipt,
            "raw_content_controls": singular_controls}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(encoded)
    else:
        print(encoded, end="")
    print("Bilinear pole packets: exact controls passed.")
