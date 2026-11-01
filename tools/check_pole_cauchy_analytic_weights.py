#!/usr/bin/env python3
"""Exact checks for analytic weight compilation and complete inverse moments."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_pole_cauchy_conditioning import inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_analytic_weights import (
    compile_regular_odd_circuit, inverse_analytic_bilinear_form,
)


ROOT = Path(__file__).resolve().parents[1]


def evaluate(circuit, argument, modulus):
    values = []
    for node in circuit:
        op = node[0]
        if op == 'var':
            value = argument
        elif op == 'const':
            value = node[1]
        elif op == 'add':
            value = values[node[1]] + values[node[2]]
        elif op == 'sub':
            value = values[node[1]] - values[node[2]]
        elif op == 'mul':
            value = values[node[1]] * values[node[2]]
        elif op == 'pow':
            value = pow(values[node[1]], node[2], modulus)
        elif op == 'inverse':
            value = pow(values[node[1]], -1, modulus)
        else:
            raise AssertionError(op)
        values.append(value % modulus)
    return values[-1]


def shifted_power(shift, exponent):
    return [('var',), ('const', 1 << shift), ('mul', 0, 1),
            ('const', 1), ('add', 2, 3), ('pow', 4, exponent)]


def controls():
    rng = random.Random(202611011)
    circuits = [[('var',)], [('var',), ('inverse', 0)],
                [('var',), ('mul', 0, 0), ('add', 0, 1),
                 ('const', 1), ('add', 2, 3), ('inverse', 4)],
                shifted_power(1, (1 << 512) + (1 << 257) + 123),
                shifted_power(3, -((1 << 383) + (1 << 193) + 7))]
    for _ in range(35):
        nodes, parity = [('var',), ('const', rng.randrange(-7, 8))], [1]
        parity.append(nodes[1][1] & 1)
        for _ in range(rng.randrange(3, 12)):
            first, second = rng.randrange(len(nodes)), rng.randrange(len(nodes))
            operation = rng.choice(('add', 'sub', 'mul', 'pow', 'inverse'))
            if operation == 'inverse' and parity[first]:
                nodes.append(('inverse', first)); parity.append(1)
            elif operation == 'pow':
                exponent = rng.randrange(-5 if parity[first] else 0, 9)
                nodes.append(('pow', first, exponent))
                parity.append(1 if exponent == 0 else parity[first])
            else:
                if operation == 'inverse':
                    operation = 'add'
                nodes.append((operation, first, second))
                parity.append(parity[first] * parity[second] if operation == 'mul'
                              else parity[first] ^ parity[second])
        circuits.append(nodes)
    values_checked = compilations = 0
    for circuit in circuits:
        for precision in range(1, 11):
            compiled = compile_regular_odd_circuit(circuit, precision)
            modulus = 1 << precision
            assert compiled['degree'] < precision
            for argument in range(1, modulus, 2):
                polynomial = 0
                for coefficient in reversed(compiled['coefficients']):
                    polynomial = (polynomial * argument + coefficient) % modulus
                assert polynomial == evaluate(circuit, argument, modulus)
                values_checked += 1
            compilations += 1
    rejected = 0
    for circuit, precision in (([('var',), ('const', 1), ('sub', 0, 1), ('inverse', 2)], 8),
                               ([('const', 2), ('pow', 0, -1)], 8),
                               ([('add', 0, 0)], 8), ([('var',)], 0)):
        try:
            compile_regular_odd_circuit(circuit, precision)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('An invalid circuit was accepted.')

    matrix_controls = entry_controls = bilinear_controls = 0
    for _ in range(40):
        size = rng.randrange(1, 7)
        first_a, first_b = rng.randrange(1, 10, 2), rng.randrange(1, 10, 2)
        left, right = ((first_a, size),), ((first_b, size),)
        rows = [first_a + 2 * i for i in range(size)]
        columns = [first_b + 2 * j for j in range(size)]
        n = rng.randrange(31, 250, 2)
        k = max(abs(n - a * b) for a in rows for b in columns).bit_length() + 2
        dense = inverse(matrix_at(n, rows, columns, k))
        normalized = [[value * (1 << k) for value in row] for row in dense]
        assert all(value.valuation() >= 0 for row in normalized for value in row)
        entry_controls += size * size
        for _ in range(2):
            row_circuit, column_circuit = rng.choice(circuits), rng.choice(circuits)
            precision = 11
            modulus = 1 << precision
            aa = [evaluate(row_circuit, a, modulus) for a in rows]
            bb = [evaluate(column_circuit, b, modulus) for b in columns]
            expected = sum((normalized[j][i] * aa[i] * bb[j]
                            for i in range(size) for j in range(size)), Quad())
            result = inverse_analytic_bilinear_form(
                n, left, right, k, row_circuit, column_circuit, precision,
                row_denominator_bits=2, column_denominator_bits=3)
            assert tuple(result['residue']) == mod_pair(expected, precision)
            assert result['scale_exponent'] == -k - 5
            assert result['absolute_precision'] == precision - k - 5
            bilinear_controls += 1
        matrix_controls += 1
    try:
        inverse_analytic_bilinear_form(101, ((1, 1), (9, 1)), ((3, 2),), 10,
                                      [('var',)], [('var',)], 8)
    except ValueError:
        rejected += 1
    else:
        raise AssertionError('An unproved domain-conditioning contract was accepted.')
    return {'regular_circuits': len(circuits), 'compiled_precision_pairs': compilations,
            'exhaustive_odd_residue_values': values_checked,
            'dense_matrix_controls': matrix_controls,
            'integral_normalized_inverse_entries': entry_controls,
            'complete_dense_bilinear_comparisons': bilinear_controls,
            'invalid_contract_rejections': rejected}


def implicit_control():
    n, size = (1 << 401) + 45, (1 << 200) + 13
    left, right = ((1, size),), ((3, size),)
    bound = max(abs(n - 3), abs(n - (2 * size - 1) * (2 * size + 1)))
    k = bound.bit_length() + 12
    row = shifted_power(64, (1 << 512) + (1 << 257) + 123)
    column = shifted_power(80, -((1 << 383) + (1 << 193) + 7))
    result = inverse_analytic_bilinear_form(n, left, right, k, row, column, k + 12)
    stronger = inverse_analytic_bilinear_form(n, left, right, k, row, column, k + 16)
    assert result['residue'] == [value % (1 << (k + 12)) for value in stronger['residue']]
    return {'integer': n, 'input_bits': n.bit_length(), 'domain_size': size,
            'left': left, 'right': right, 'row_circuit': row, 'column_circuit': column,
            'k': k, 'result': result, 'higher_precision_result': stronger,
            'symbolic_row_power_bits': row[-1][2].bit_length(),
            'symbolic_column_power_bits': abs(column[-1][2]).bit_length(),
            'enumerated_matrix_nodes': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'claim': 'Compressed regular analytic weights for implicit complete Cauchy inverse moments',
              'controls': controls(), 'implicit_control': implicit_control(),
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'pole_cauchy_analytic_weights.py',
                                             'pole_cauchy_polynomial_solve.py', 'modular_interval_jets.py')},
              'factoring_exponent_improved': False,
              'scope': 'Regular odd-input weight circuits; no singular factor-count trace is computed.'}
    encoded = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({'controls': report['controls'],
                      'input_bits': report['implicit_control']['input_bits'],
                      'domain_size': report['implicit_control']['domain_size'],
                      'absolute_precision': report['implicit_control']['result']['absolute_precision'],
                      'compiled_degrees': [report['implicit_control']['result'][name]['degree']
                                           for name in ('row_compiler', 'column_compiler')]}, indent=2))


if __name__ == '__main__':
    main()

