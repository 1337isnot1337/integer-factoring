#!/usr/bin/env python3
"""Exact canonical binary-packet visibility bounds at lengths two and four.

The length-four search uses an optional C++17 meet-in-the-middle engine.
The projections are necessary conditions, never probabilistic conclusions.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import permutations
import json
from math import comb
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import random

from cauchy_packet_polynomials import cauchy_polynomials, multiply, shift


PRECISION = 32
MODULUS = 1 << PRECISION


def coefficient_matrices(length):
    degree = length * (2 * length - 1)
    count = 1 << length
    result = []
    kernels = cauchy_polynomials(length)
    for center in range(-length + 1, length, 2):
        local = {difference: [value.numerator * pow(value.denominator, -1, MODULUS) % MODULUS
                              for value in shift(polynomial, center, 2 * length)]
                 for difference, polynomial in kernels.items()}
        block = [[[0] * count for _ in range(count)] for _ in range(degree + 1)]
        for column in range(count):
            images = [[sum(local[y - x][power] for y in range(length) if column >> y & 1) % MODULUS
                       for power in range(2 * length)] for x in range(length)]
            for row in range(count):
                polynomial = [1]
                for x, image in enumerate(images):
                    factor = image[:]
                    if not row >> x & 1:
                        factor = [-value for value in factor]
                        factor[0] += 1
                    polynomial = multiply(polynomial, factor, MODULUS)
                for power, coefficient in enumerate(polynomial):
                    block[power][row][column] = coefficient
        result.extend(block)
    return result


def content(matrices, row_mask, column_mask):
    count = len(matrices[0])
    aggregate = 0
    for matrix in matrices:
        value = sum(matrix[row][column] for row in range(count) if row_mask >> row & 1
                    for column in range(count) if column_mask >> column & 1) % MODULUS
        aggregate |= value
    return (aggregate & -aggregate).bit_length() - 1 if aggregate else PRECISION


def check():
    small = coefficient_matrices(2)
    maximum = -1
    pairs = 0
    for row in range(16):
        for column in range(16):
            if (row & column).bit_count() % 2:
                pairs += 1
                maximum = max(maximum, content(small, row, column))
    assert pairs == 120 and maximum == 2

    matrices = coefficient_matrices(4)
    assert len(matrices) == 116
    assert content(matrices, 1 << 5, 1 << 5) == 11
    adaptive_checks = 0
    moment_checks = 0
    rng = random.Random(67213)
    for length, all_matrices in ((2, small), (4, matrices)):
        count = 1 << length
        block_size = length * (2 * length - 1) + 1
        positive_endpoint = all_matrices[(length - 1) * block_size]
        for _ in range(64):
            row_precision = rng.randrange(1, 5)
            column_precision = rng.randrange(1, 5)
            row = [rng.randrange(1 << row_precision) for _ in range(count)]
            column = [rng.randrange(1 << column_precision) for _ in range(count)]
            if not any(row):
                row[0] = 1
            if not any(column):
                column[0] = 1
            maximum_state = max(range(count), key=row.__getitem__)
            translated_row = [row[state ^ maximum_state] for state in range(count)]
            translated_column = [column[state ^ maximum_state] for state in range(count)]
            assert sum(a * b for a, b in zip(row, column)) == sum(
                a * b for a, b in zip(translated_row, translated_column))
            expected = translated_row[0] * sum(translated_column) + comb(2 * length - 1, length) * (
                translated_row[0] - translated_row[1]) * sum(
                    value for state, value in enumerate(translated_column) if state >> (length - 1) & 1)
            actual = sum(translated_row[y] * positive_endpoint[y][x] * translated_column[x]
                         for y in range(count) for x in range(count)) % MODULUS
            assert actual == expected and 0 < expected < 1 << (row_precision + column_precision + 3 * length - 1)
            adaptive_checks += 1
        for order_limit in range(length - 1):
            active = order_limit + 2
            for permutation in permutations(range(length)):
                permuted = [sum(((state >> bit) & 1) << permutation[bit] for bit in range(active))
                            for state in range(1 << active)]
                for translation in range(count):
                    even = [permuted[state] ^ translation for state in range(1 << active)
                            if state.bit_count() % 2 == 0]
                    odd = [permuted[state] ^ translation for state in range(1 << active)
                           if state.bit_count() % 2 == 1]
                    assert translation in even and translation not in odd
                    for endpoint in (0, (length - 1) * block_size):
                        for order in range(order_limit + 1):
                            matrix = all_matrices[endpoint + order]
                            assert (sum(matrix[translation][state] for state in even)
                                    - sum(matrix[translation][state] for state in odd)) % MODULUS == 0
                            moment_checks += 1
    source = Path(__file__).with_name('canonical_packet_exhaustive.cpp')
    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional exhaustive packet checker requires C++17.')
    with tempfile.TemporaryDirectory(prefix='integer-factoring-canonical-') as directory:
        temporary = Path(directory)
        binary = temporary / 'exhaustive'
        data = temporary / 'coefficients.txt'
        run = subprocess.run(compiler + ['-O3', '-std=c++17', str(source), '-o', str(binary)],
                             text=True, capture_output=True)
        if run.returncode:
            raise RuntimeError(run.stderr[-12000:])
        with data.open('w') as output:
            output.write(str(len(matrices)) + '\n')
            for matrix in matrices:
                for row in matrix:
                    output.write(' '.join(map(str, row)) + '\n')
        records = []
        for precision, seed in ((12, 123456789123), (12, 808005192913), (11, 123456789123)):
            run = subprocess.run([str(binary), str(data), str(precision), str(seed)],
                                 text=True, capture_output=True)
            if run.returncode:
                raise RuntimeError(run.stderr[-12000:])
            result = json.loads(run.stdout)
            result['projection_seed'] = seed
            if precision == 12:
                assert result['found'] is False
                assert result['column_tensors_tested'] == 65535
                assert result['right_subset_queries'] == 65535 * 256
                assert result['odd_original_pairs_covered'] == 65535 * 32768
            else:
                assert result['found'] is True
                row, column = result['row_tensor_mask'], result['column_tensor_mask']
                assert (row & column).bit_count() % 2 == 1
                assert content(matrices, row, column) == 11
            records.append(result)

    root = Path(__file__).resolve().parent.parent
    sources = ['tools/canonical_packet_exhaustive.cpp', 'tools/check_canonical_packet_exhaustive.py',
               'tools/cauchy_packet_polynomials.py', 'research/canonical-packet-upper-frontier.md']
    return {'status': 'pass', 'creation_precision': 1,
            'length_two_odd_pairs': pairs, 'length_two_maximum_content': maximum,
            'length_four_maximum_content': 11, 'length_four_coefficient_matrices': len(matrices),
            'adaptive_positive_endpoint_checks': adaptive_checks,
            'translated_endpoint_moment_checks': moment_checks,
            'length_four_all_binary_pairs': 1 << 32, 'runs': records,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in sources},
            'scope': 'Exact finite upper and lower bounds on odd-center visibility for '
                     'canonical binary tensors. This does not prove an all-size upper bound '
                     'or construct a decoder of original parity from the local data.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.write_text(json.dumps({'status': 'fail', 'error': repr(error)}, indent=2) + '\n')
        raise
    rendered = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
