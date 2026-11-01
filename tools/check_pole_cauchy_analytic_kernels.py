#!/usr/bin/env python3
"""Independent exact controls for complete bivariate inverse contractions."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_pole_cauchy_conditioning import inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_analytic_kernels import compile_regular_odd_kernel
from pole_cauchy_analytic_kernels import inverse_analytic_kernel_contraction
from pole_cauchy_analytic_kernels import inverse_polynomial_kernel_contraction
from pole_cauchy_polynomial_solve import inverse_polynomial_bilinear_form

ROOT = Path(__file__).resolve().parents[1]


def evaluate(circuit, a, b, modulus):
    values = []
    for node in circuit:
        op = node[0]
        if op == 'var': value = (a, b)[node[1]]
        elif op == 'const': value = node[1]
        elif op == 'add': value = values[node[1]] + values[node[2]]
        elif op == 'sub': value = values[node[1]] - values[node[2]]
        elif op == 'mul': value = values[node[1]] * values[node[2]]
        elif op == 'pow': value = pow(values[node[1]], node[2], modulus)
        elif op == 'inverse': value = pow(values[node[1]], -1, modulus)
        else: raise AssertionError(op)
        values.append(value % modulus)
    return values[-1]


def large_circuit():
    return [('var', 0), ('var', 1), ('mul', 0, 1),
            ('const', 1 << 64), ('mul', 2, 3), ('const', 1),
            ('add', 4, 5), ('pow', 6, -((1 << 512) + (1 << 257) + 123)),
            ('add', 0, 1), ('const', 1 << 80), ('mul', 8, 9),
            ('add', 10, 5), ('pow', 11, -((1 << 383) + (1 << 193) + 7)),
            ('mul', 7, 12)]


def catalogue(rng):
    circuits = [[('var', 0)], [('var', 1)],
                [('var', 0), ('var', 1), ('add', 0, 1), ('const', 1),
                 ('add', 2, 3), ('inverse', 4)],
                [('var', 0), ('var', 1), ('mul', 0, 1), ('const', 2),
                 ('mul', 2, 3), ('const', 1), ('add', 4, 5),
                 ('pow', 6, -((1 << 120) + 17))], large_circuit()]
    for _ in range(15):
        nodes = [('var', 0), ('var', 1), ('const', rng.randrange(-4, 5))]
        parity = [1, 1, nodes[-1][1] & 1]
        for _ in range(rng.randrange(2, 9)):
            i, j = rng.randrange(len(nodes)), rng.randrange(len(nodes))
            operation = rng.choice(('add', 'sub', 'mul', 'pow', 'inverse'))
            if operation == 'inverse' and parity[i]:
                nodes.append(('inverse', i)); parity.append(1)
            elif operation == 'pow':
                exponent = rng.randrange(-3 if parity[i] else 0, 6)
                nodes.append(('pow', i, exponent)); parity.append(1 if not exponent else parity[i])
            else:
                if operation == 'inverse': operation = 'add'
                nodes.append((operation, i, j))
                parity.append(parity[i] * parity[j] if operation == 'mul' else parity[i] ^ parity[j])
        circuits.append(nodes)
    return circuits


def controls():
    rng = random.Random(202611012)
    circuits = catalogue(rng)
    pair_controls = compiled_controls = 0
    for circuit in circuits:
        for precision in range(1, 8):
            compiled = compile_regular_odd_kernel(circuit, precision)
            assert compiled['total_degree'] < precision
            modulus = 1 << precision
            for a in range(1, modulus, 2):
                for b in range(1, modulus, 2):
                    polynomial = sum(value * pow(a, d, modulus) * pow(b, e, modulus)
                                     for d, e, value in compiled['terms']) % modulus
                    assert polynomial == evaluate(circuit, a, b, modulus)
                    pair_controls += 1
            compiled_controls += 1
    dense_controls = polynomial_controls = 0
    for _ in range(30):
        size, first_a, first_b = rng.randrange(1, 6), rng.randrange(1, 10, 2), rng.randrange(1, 10, 2)
        left, right = ((first_a, size),), ((first_b, size),)
        rows, columns = [first_a + 2 * i for i in range(size)], [first_b + 2 * i for i in range(size)]
        n = rng.randrange(31, 250, 2)
        k = max(abs(n - a * b) for a in rows for b in columns).bit_length() + 2
        precision, modulus = 10, 1 << 10
        dense = inverse(matrix_at(n, rows, columns, k))
        for _ in range(2):
            circuit = rng.choice(circuits)
            expected = sum((dense[j][i] * (1 << k) * evaluate(circuit, a, b, modulus)
                            for i, a in enumerate(rows) for j, b in enumerate(columns)), Quad())
            actual = inverse_analytic_kernel_contraction(n, left, right, k, circuit, precision,
                                                          denominator_bits=3)
            assert tuple(actual['residue']) == mod_pair(expected, precision)
            assert actual['scale_exponent'] == -k - 3
            dense_controls += 1
        f, g = tuple(rng.randrange(-5, 6) for _ in range(8)), tuple(rng.randrange(-5, 6) for _ in range(7))
        old = inverse_polynomial_bilinear_form(n, left, right, k, f, g, precision)
        terms = [(d, e, a * b) for d, a in enumerate(f) for e, b in enumerate(g)]
        new = inverse_polynomial_kernel_contraction(n, left, right, k, terms, precision)
        assert old['residue'] == new['residue']
        assert new['moment_recurrence_cells'] == len(f) * len(g)
        polynomial_controls += 1
    rejected = 0
    for circuit in ([('var', 0), ('var', 1), ('mul', 0, 1), ('const', 101),
                     ('sub', 3, 2), ('inverse', 4)], [('var', 2)], [('inverse', 0)]):
        try: compile_regular_odd_kernel(circuit, 12)
        except ValueError: rejected += 1
        else: raise AssertionError('An invalid/singular circuit was accepted.')
    return {'circuits': len(circuits), 'compiled_precision_pairs': compiled_controls,
            'exhaustive_odd_pair_values': pair_controls, 'dense_joint_kernel_contractions': dense_controls,
            'quadratic_recurrence_comparisons': polynomial_controls, 'invalid_or_singular_rejections': rejected}


def implicit_control():
    n, size = (1 << 401) + 45, (1 << 200) + 13
    left, right = ((1, size),), ((3, size),)
    k = max(abs(n - 3), abs(n - (2 * size - 1) * (2 * size + 1))).bit_length() + 12
    circuit = large_circuit()
    lower = inverse_analytic_kernel_contraction(n, left, right, k, circuit, k + 12)
    upper = inverse_analytic_kernel_contraction(n, left, right, k, circuit, k + 16)
    assert lower['residue'] == [value % (1 << (k + 12)) for value in upper['residue']]
    return {'integer': n, 'input_bits': n.bit_length(), 'domain_size': size,
            'left': left, 'right': right, 'k': k, 'circuit': circuit,
            'result': lower, 'higher_precision_result': upper, 'enumerated_matrix_nodes': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'claim': 'Complete implicit inverse contractions with regular joint analytic kernels',
              'controls': controls(), 'implicit_control': implicit_control(),
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'pole_cauchy_analytic_kernels.py',
                                             'pole_cauchy_polynomial_solve.py', 'modular_interval_jets.py')},
              'factoring_exponent_improved': False,
              'scope': 'Unit-denominator joint circuits; the singular factor-count kernel remains excluded.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'controls': report['controls'], 'input_bits': report['implicit_control']['input_bits'],
                      'domain_size': report['implicit_control']['domain_size'],
                      'compiled_total_degree': report['implicit_control']['result']['compiler']['total_degree'],
                      'moment_recurrence_cells': report['implicit_control']['result']['moment_recurrence_cells'],
                      'absolute_precision': report['implicit_control']['result']['absolute_precision']}, indent=2))


if __name__ == '__main__':
    main()
