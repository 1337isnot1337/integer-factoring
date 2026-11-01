#!/usr/bin/env python3
"""Exact controls for the sharp precision cost of regular residue selectors."""

import argparse
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import random

from pole_cauchy_analytic_weights import compile_regular_odd_circuit

ROOT = Path(__file__).resolve().parents[1]


def valuation(value):
    value = abs(value)
    return (value & -value).bit_length() - 1 if value else None


def evaluate(circuit, x, modulus):
    values = []
    for gate in circuit:
        op = gate[0]
        if op == 'var': result = x
        elif op == 'const': result = gate[1]
        elif op == 'add': result = values[gate[1]] + values[gate[2]]
        elif op == 'sub': result = values[gate[1]] - values[gate[2]]
        elif op == 'mul': result = values[gate[1]] * values[gate[2]]
        elif op == 'pow': result = pow(values[gate[1]], gate[2], modulus)
        elif op == 'inverse': result = pow(values[gate[1]], -1, modulus)
        else: raise AssertionError(op)
        values.append(result % modulus)
    return values[-1]


def difference(values):
    degree = len(values) - 1
    return sum((-1 if (degree - j) & 1 else 1) * comb(degree, j) * value
               for j, value in enumerate(values))


def random_circuits():
    rng = random.Random(202611013)
    circuits = []
    for _ in range(40):
        nodes = [('var',), ('const', rng.randrange(-8, 9))]
        parity = [1, nodes[-1][1] & 1]
        for _ in range(rng.randrange(3, 13)):
            i, j = rng.randrange(len(nodes)), rng.randrange(len(nodes))
            op = rng.choice(('add', 'sub', 'mul', 'pow', 'inverse'))
            if op == 'inverse' and parity[i]:
                nodes.append((op, i)); parity.append(1)
            elif op == 'pow':
                exponent = rng.randrange(-4 if parity[i] else 0, 10)
                nodes.append((op, i, exponent)); parity.append(1 if not exponent else parity[i])
            else:
                if op == 'inverse': op = 'sub'
                nodes.append((op, i, j))
                parity.append(parity[i] * parity[j] if op == 'mul' else parity[i] ^ parity[j])
        circuits.append(nodes)
    return circuits


def sharp_selector_circuit(h):
    degree = (1 << h) - 1
    odd_factorial = factorial(degree)
    factorial_valuation = valuation(odd_factorial)
    odd_factorial >>= factorial_valuation
    nodes = [('var',), ('const', 1)]
    product = 1
    for j in range(1, degree + 1):
        nodes.append(('const', 2 * j - 1))
        nodes.append(('add', 0, len(nodes) - 1))
        nodes.append(('mul', product, len(nodes) - 1))
        product = len(nodes) - 1
    nodes.append(('const', odd_factorial))
    nodes.append(('inverse', len(nodes) - 1))
    nodes.append(('mul', product, len(nodes) - 1))
    return nodes, degree + factorial_valuation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    difference_controls = 0
    for circuit in random_circuits():
        for degree in range(1, 13):
            bound = degree + valuation(factorial(degree))
            modulus = 1 << (bound + 3)
            for offset in (-5, 0, 3, 12):
                values = [evaluate(circuit, 1 + 2 * (offset + j), modulus)
                          for j in range(degree + 1)]
                assert difference(values) % (1 << bound) == 0
                difference_controls += 1
    sharp = []
    total_selector_values = 0
    for h in range(1, 9):
        circuit, denominator_bits = sharp_selector_circuit(h)
        degree, period = (1 << h) - 1, 1 << h
        assert denominator_bits == (1 << (h + 1)) - h - 2
        precision = denominator_bits + 1
        modulus = 1 << precision
        compiled = compile_regular_odd_circuit(circuit, precision)
        residues = []
        for x in range(-period, 2 * period):
            odd = 1 + 2 * x
            actual = evaluate(circuit, odd, modulus)
            polynomial = sum(c * pow(odd, i, modulus)
                             for i, c in enumerate(compiled['coefficients'])) % modulus
            assert actual == polynomial
            assert actual % (1 << denominator_bits) == 0
            assert (actual >> denominator_bits) == int(x % period == 0)
            if 0 <= x <= degree:
                residues.append(actual)
            total_selector_values += 1
        assert valuation(difference(residues) % modulus) == denominator_bits
        sharp.append({'residue_bits_in_centered_variable': h, 'odd_input_modulus': 1 << (h + 1),
                      'degree': degree, 'sharp_denominator_bits': denominator_bits,
                      'compiled_degree': compiled['degree'], 'precision_for_parity': precision,
                      'tested_integer_arguments': 3 * period, 'circuit_gates': len(circuit)})
    report = {'claim': 'Sharp dyadic denominator cost for a regular analytic residue-class parity selector',
              'finite_difference_controls': difference_controls,
              'sharp_selector_values': total_selector_values,
              'sharp_selectors': sharp,
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'pole_cauchy_analytic_weights.py')},
              'factoring_lower_bound': False,
              'scope': 'Regular odd-input circuits divided by a charged external power of two; pointwise parity selection.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
