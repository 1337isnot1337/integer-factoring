#!/usr/bin/env python3
"""Independent controls for the local Smith engine and a complete nine-letter probe."""
import argparse
from fractions import Fraction
import json
from math import exp, log
from pathlib import Path
import random
import subprocess
import tempfile

from check_quadratic_crossings import ROOT, binary_rank, matrices, orders
from check_quadratic_inverse_paths import compose, construct, inverse, kernel_flag


def native(executable, cases, precision):
    payload = ''.join(f'{len(matrix)} {len(matrix[0])} {precision}\n' +
                      '\n'.join(' '.join(map(str, row)) for row in matrix) + '\n'
                      for matrix in cases)
    output = subprocess.check_output([str(executable)], input=payload, text=True)
    return [json.loads(line) for line in output.splitlines()]


def prefix_rank(matrix, precision):
    columns = len(matrix[0])
    rows = []
    for output in range(precision):
        for row in matrix:
            packed = 0
            for input_degree in range(output+1):
                for j, entry in enumerate(row):
                    if entry >> (output-input_degree) & 1:
                        packed |= 1 << (columns*input_degree+j)
            rows.append(packed)
    return binary_rank(rows)


def nine_letter():
    q, n = 9, 81
    unit = [(1 << (q*b+a)) ^ ((1 << (q*(b-1)+a)) if b else 0)
            for a in range(q) for b in range(q)]
    nilpotent = [column ^ 1 << i for i, column in enumerate(unit)]
    basis, levels = kernel_flag(nilpotent)
    inverse_basis = inverse(basis)
    triangular = compose(compose(inverse_basis, nilpotent), basis)
    physical, _ = construct(q, basis, levels, triangular, inverse_basis, (32,))
    degree = max(value.bit_length()-1 for row in physical for value in row)
    row = [[physical[q*x+u][q*v+y] for u in range(q) for v in range(q)]
           for x in range(q) for y in range(q)]
    # Column input transposition is an invertible constant permutation, so it
    # can be omitted for Smith orders of the reversed operator matrix.
    reverse = lambda p: sum((p >> i & 1) << (degree-i) for i in range(degree+1))
    column = [[reverse(entry) for entry in row] for row in physical]
    assert degree == 61 and len(row) == len(column) == n
    return degree, row, column


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nine-letter', action='store_true')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='factoring-local-smith-') as directory:
        executable = Path(directory)/'local-smith'
        subprocess.run(['clang++', '-O3', '-std=c++17', '-Wall', '-Wextra', '-Werror',
                        str(ROOT/'tools/polynomial_local_smith.cpp'), '-o', str(executable)], check=True)
        rng = random.Random(8403)
        cases = [matrix for degree in (2, 3, 4) for _ in range(24)
                 for matrix in matrices(rng.getrandbits(16*(degree-1)), degree)]
        cases += [matrix + [[a ^ b for a, b in zip(matrix[0], matrix[1])], list(matrix[2])]
                  for matrix in cases[:24]]
        results = native(executable, cases, 32)
        assert len(results) == len(cases)
        for matrix, result in zip(cases, results):
            assert result['complete'] and tuple(result['orders']) == orders(matrix)
        short = native(executable, cases, 4)
        for full, partial in zip(results, short):
            assert partial['orders'] == [value for value in full['orders'] if value < 4]
        singular = native(executable, [[[0]*4 for _ in range(4)], [[4, 0], [0, 0]]], 8)
        assert singular[0]['orders'] == [] and not singular[0]['complete']
        assert singular[1]['orders'] == [2] and not singular[1]['complete']
        report = {'independent_complete_small_matrices': len(cases),
                  'truncated_precision_controls': len(cases), 'singular_controls': 2}
        if args.nine_letter:
            degree, row, column = nine_letter()
            precision = 128
            while True:
                profiles = native(executable, [row, column], precision)
                if all(result['complete'] for result in profiles):
                    break
                precision *= 2
                assert precision <= 8192
            mu, nu = [result['orders'] for result in profiles]
            assert sum(nu) == 81*degree
            assert sum(mu) <= 81*degree
            assert sum(mu[:9])+sum(nu[:9]) <= 9*degree
            assert mu[1] <= 2*(max(nu)-degree)
            for matrix, profile in zip((row, column), (mu, nu)):
                for cutoff in (1, 2, 3, 4):
                    assert prefix_rank(matrix, cutoff) == sum(max(0, cutoff-value) for value in profile)
            def log_pressure(t):
                return degree*t+log(sum(exp(-value*t) for value in mu))+log(sum(exp(-value*t) for value in nu))
            lo, hi = 0, 2
            for _ in range(100):
                a, b = (2*lo+hi)/3, (lo+2*hi)/3
                if log_pressure(a) < log_pressure(b):
                    hi = b
                else:
                    lo = a
            t = (lo+hi)/2
            scale = Fraction(round(exp(t)*10**6), 10**6)
            pressure = scale**degree * sum(scale**(-value) for value in mu) * sum(scale**(-value) for value in nu)
            channel_count = mu.count(2)
            cutoff_scale = Fraction(6, 5)
            small_floor = 81*channel_count/cutoff_scale**4
            large_floor = channel_count*cutoff_scale**(degree-2)
            assert max(nu) <= degree+2 and min(nu) == 0
            assert small_floor**7 > 9**20 and large_floor**7 > 9**20
            report['nine_letter'] = {'degree': degree, 'precision': precision,
                                     'infinity_orders': mu, 'zero_orders': nu,
                                     'independent_prefix_checks': 8,
                                     'approximate_optimal_scale': exp(t),
                                     'approximate_rho': log_pressure(t)/(4*log(9)),
                                     'exact_evaluation_scale': str(scale),
                                     'exact_pressure_numerator_bits': pressure.numerator.bit_length(),
                                     'exact_pressure_denominator_bits': pressure.denominator.bit_length(),
                                     'all_scale_rho_lower_bound': '5/7',
                                     'all_scale_split': str(cutoff_scale),
                                     'small_scale_pressure_lower_bound': str(small_floor),
                                     'large_scale_pressure_lower_bound': str(large_floor),
                                     'scope': 'Complete exact profiles; numerical optimization is exploratory'}
        print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
