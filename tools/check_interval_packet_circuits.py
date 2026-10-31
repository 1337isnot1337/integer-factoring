#!/usr/bin/env python3
"""Coefficient and large-jet controls for global interval-product circuits."""

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random

from interval_orbit_gcd import interval_product
from interval_packet_circuits import (evaluate_circuit, interval_product_circuit,
                                      paired_elementaries)
from interval_unit_jets import QuadraticJets, interval_elementaries, interval_product_jet


class DensePolynomials:
    def __init__(self, precision):
        self.modulus = 1 << precision

    def constant(self, value):
        return (value % self.modulus,)

    def add(self, a, b):
        values = [0] * max(len(a), len(b))
        for polynomial in (a, b):
            for i, coefficient in enumerate(polynomial):
                values[i] = (values[i] + coefficient) % self.modulus
        while len(values) > 1 and not values[-1]:
            values.pop()
        return tuple(values)

    def multiply(self, a, b):
        values = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                values[i + j] = (values[i + j] + x * y) % self.modulus
        while len(values) > 1 and not values[-1]:
            values.pop()
        return tuple(values)


def check():
    rng = random.Random(20261035)
    elementary_checks = coefficient_checks = spurious_factor_checks = 0
    for _ in range(256):
        lower, pairs, degree = rng.randrange(1, 20), rng.randrange(17), rng.randrange(20)
        expected = [1] + [0] * degree
        for t in range(pairs):
            value = (lower + t) * (lower + 2 * pairs - 1 - t) // 2
            for k in range(degree, 0, -1):
                expected[k] += value * expected[k - 1]
        assert paired_elementaries(lower, pairs, degree) == tuple(expected)
        elementary_checks += 1
    for precision in range(1, 17):
        ring = DensePolynomials(precision)
        for _ in range(32):
            lower = rng.randrange(1, 24)
            interval = (lower, lower + rng.randrange(34))
            for number in (None, rng.randrange(100)):
                circuit = interval_product_circuit(interval, precision, number)
                actual = evaluate_circuit(circuit, ring, (0, 1))
                expected = ring.add(tuple(a % ring.modulus
                                          for a in interval_product(interval, number)), (0,))
                assert actual == expected
                coefficient_checks += 1
    # Rational gcd one, but the modular inputs have a large common factor.
    for pairs_x in range(1, 9):
        for pairs_y in range(1, 9):
            center_x, center_y = 2 * pairs_x + 1, 2 * pairs_y + 1
            number = center_x * center_y
            assert (center_x - 1) * (center_y - 1) < number
            for precision in range(1, min(pairs_x, pairs_y) + 1):
                exponent = min(pairs_x, pairs_y) - precision + 1
                for interval, selected in (((1, center_x), None), ((1, center_y), number)):
                    polynomial = interval_product(interval, selected)
                    for k in range(exponent):
                        shifted = sum(comb(i, k) * a * center_x ** (i - k)
                                      for i, a in enumerate(polynomial) if i >= k)
                        assert shifted % (1 << precision) == 0
                spurious_factor_checks += 1
    large = []
    for precision, degree, endpoint_bits, length_bits in ((32, 6, 120, 110),
                                                          (48, 8, 240, 200)):
        ring = QuadraticJets(precision, degree)
        argument = ring.argument()
        lower = (1 << endpoint_bits) + 7
        interval = (lower, lower + (1 << length_bits) + 12)
        number = (1 << (endpoint_bits + 7)) - 1
        for selected in (None, number):
            circuit = interval_product_circuit(interval, precision, selected)
            assert evaluate_circuit(circuit, ring, argument) == interval_product_jet(
                ring, interval, argument, selected)
            # At zero, Q's first coefficients have a separate direct Newton formula.
            if selected is not None:
                zero_argument = list(ring.constant())
                zero_argument[1] = (1, 0)
                actual = evaluate_circuit(circuit, ring, tuple(zero_argument))
                elementary = interval_elementaries(*interval, degree)
                expected = tuple(((-1) ** k * elementary[k]
                                  * pow(number, interval[1] - interval[0] - k, ring.modulus)
                                  % ring.modulus, 0) for k in range(degree + 1))
                assert actual == expected
            large.append({'precision': precision, 'jet_degree': degree,
                          'interval': [str(a) for a in interval],
                          'number': None if selected is None else str(selected),
                          'nodes': len(circuit['nodes']),
                          'input_degree': str(circuit['input_degree']),
                          'circuit_sha256': sha256(json.dumps(circuit, sort_keys=True).encode()).hexdigest()})
    return {
        'scope': 'Global polynomial circuits modulo 2^p, not just pointwise function equality; no global gcd algorithm.',
        'seed': 20261035,
        'paired_symmetric_controls': elementary_checks,
        'complete_coefficient_controls': coefficient_checks,
        'spurious_modular_common_factor_controls': spurious_factor_checks,
        'large_unit_jet_cross_checks': 4,
        'large_zero_jet_cross_checks': 2,
        'large_circuits': large,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    result = check()
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
