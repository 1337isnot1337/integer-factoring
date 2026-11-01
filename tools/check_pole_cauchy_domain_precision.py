#!/usr/bin/env python3
"""Independent exact controls for implicit domain precision certificates."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_pole_cauchy_analytic_kernels import evaluate, large_circuit
from check_pole_cauchy_conditioning import expand, inverse, mod_pair
from check_pole_smith_stabilization import Quad, matrix_at
from pole_cauchy_conditioning import balanced_padding
from pole_cauchy_domain_precision import inverse_certified_analytic_kernel
from pole_cauchy_domain_precision import inverse_precision_certificate, residue_count_range

ROOT = Path(__file__).resolve().parents[1]


def controls():
    rng = random.Random(202611015)
    range_controls = 0
    for _ in range(120):
        domain, first = [], rng.randrange(1, 8, 2)
        for _ in range(rng.randrange(1, 6)):
            length = rng.randrange(1, 15)
            domain.append((first, length))
            first += 2 * (length + rng.randrange(1, 12))
        nodes = expand(tuple(domain))
        for level in range(1, 10):
            modulus = 1 << level
            counts = [sum(a % modulus == c for a in nodes) for c in range(1, modulus, 2)]
            actual = residue_count_range(domain, level)
            assert (actual['minimum'], actual['maximum']) == (min(counts), max(counts))
            assert actual['event_positions'] <= 4 * len(domain) + 2
            range_controls += 1
    dense_matrices = entries = contractions = balanced_cases = unbalanced_cases = 0
    examples = []
    for index in range(40):
        n = rng.randrange(51, 180, 2)
        if index % 2 == 0:
            left, right, _ = balanced_padding(n, (1, rng.randrange(1, 6)), (3, rng.randrange(1, 6)))
        else:
            size = rng.randrange(2, 5)
            left = tuple((a, 1) for a in sorted(rng.sample(list(range(1, 40, 2)), size)))
            right = tuple((a, 1) for a in sorted(rng.sample(list(range(1, 40, 2)), size)))
        rows, columns = expand(left), expand(right)
        k = max(abs(n - a * b) for a in rows for b in columns).bit_length() + 2
        certificate = inverse_precision_certificate(n, left, right, k)
        dense = inverse(matrix_at(n, rows, columns, k))
        loss = certificate['inverse_loss_bound']
        actual_loss = max(0, -min(entry.valuation() for row in dense for entry in row))
        assert actual_loss <= loss
        entries += len(rows) ** 2
        if certificate['balanced_through_displacement']: balanced_cases += 1
        else: unbalanced_cases += 1
        circuit = [('var', 0), ('var', 1), ('add', 0, 1), ('const', 1),
                   ('add', 2, 3), ('inverse', 4)]
        result = inverse_certified_analytic_kernel(n, left, right, k, circuit, 3,
                                                   denominator_bits=2, maximum_compilation_precision=1000)
        modulus = 1 << result['compilation_precision']
        expected = sum((dense[j][i] * (1 << k) * evaluate(circuit, a, b, modulus)
                        for i, a in enumerate(rows) for j, b in enumerate(columns)), Quad())
        assert result['residue'] == list(mod_pair(expected, result['moment_precision']))
        assert result['scale_exponent'] + result['moment_precision'] == 3
        contractions += 1
        dense_matrices += 1
        if len(examples) < 6:
            examples.append({'n': n, 'left': left, 'right': right, 'k': k,
                             'certified_loss': loss, 'exact_entry_loss': actual_loss,
                             'balanced': certificate['balanced_through_displacement']})
    rejected = 0
    for args in ((101, ((1, 2),), ((3, 3),), 10),
                 (101, ((1, 3), (3, 2)), ((1, 5),), 10),
                 (101, ((2, 2),), ((3, 2),), 10)):
        try: inverse_precision_certificate(*args)
        except ValueError: rejected += 1
        else: raise AssertionError('An invalid domain was accepted.')
    try:
        inverse_certified_analytic_kernel(101, ((1, 2),), ((3, 2),), 10,
                                          [('const', 1)], 10, maximum_compilation_precision=2)
    except ValueError: rejected += 1
    else: raise AssertionError('The supplied precision cap was ignored.')
    return {'residue_range_controls': range_controls, 'dense_matrices': dense_matrices,
            'exact_inverse_entries': entries, 'complete_contractions': contractions,
            'balanced_cases': balanced_cases, 'unbalanced_cases': unbalanced_cases,
            'invalid_contracts_rejected': rejected, 'examples': examples}


def large_control():
    n, size = (1 << 401) + 45, 1 << 200
    left, right, padding_modulus = balanced_padding(n, (1, (1 << 180) + 7), (3, size))
    max_a, max_b = left[-1][0] + 2 * (left[-1][1] - 1), right[-1][0] + 2 * (right[-1][1] - 1)
    k = max(abs(n - left[0][0] * right[0][0]), abs(n - max_a * max_b)).bit_length() + 12
    result = inverse_certified_analytic_kernel(n, left, right, k, large_circuit(), 12)
    higher = inverse_certified_analytic_kernel(n, left, right, k, large_circuit(), 16)
    assert result['residue'] == [value % (1 << result['moment_precision']) for value in higher['residue']]
    assert result['precision_certificate']['balanced_through_displacement']
    assert result['precision_certificate']['inverse_loss_bound'] == k - 201
    return {'n': n, 'input_bits': n.bit_length(), 'left': left, 'right': right,
            'padding_modulus': padding_modulus, 'domain_size': size, 'k': k,
            'result': result, 'higher_precision_result': higher, 'enumerated_matrix_nodes': 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'claim': 'Implicit residue-range certificates with fully charged analytic precision',
              'controls': controls(), 'large_control': large_control(),
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'pole_cauchy_domain_precision.py',
                                             'pole_cauchy_analytic_kernels.py')},
              'factoring_exponent_improved': False}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    large = report['large_control']
    print(json.dumps({'controls': report['controls'], 'input_bits': large['input_bits'],
                      'domain_size': large['domain_size'], 'k': large['k'],
                      'inverse_loss': large['result']['precision_certificate']['inverse_loss_bound'],
                      'compilation_precision': large['result']['compilation_precision'],
                      'moment_precision': large['result']['moment_precision'],
                      'compiled_degree': large['result']['compiler']['total_degree']}, indent=2))


if __name__ == '__main__':
    main()
