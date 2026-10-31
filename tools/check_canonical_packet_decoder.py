#!/usr/bin/env python3
"""Exact local-data separation of canonical binary columns at L=2 and L=4."""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

from check_canonical_packet_exhaustive import coefficient_matrices, MODULUS, PRECISION


def difference_content(matrices, row_mask, column):
    count = len(column)
    aggregate = 0
    for matrix in matrices:
        value = sum(matrix[row][state] * column[state]
                    for row in range(count) if row_mask >> row & 1
                    for state in range(count)) % MODULUS
        aggregate |= value
    return (aggregate & -aggregate).bit_length() - 1 if aggregate else PRECISION


def check():
    small = coefficient_matrices(2)
    maximum_odd = maximum_nonzero = -1
    nonzero_pairs = odd_pairs = 0
    for row in range(16):
        for column in product((-1, 0, 1), repeat=4):
            original = sum(column[state] for state in range(4) if row >> state & 1)
            if not original:
                continue
            nonzero_pairs += 1
            content = difference_content(small, row, column)
            maximum_nonzero = max(maximum_nonzero, content)
            if original % 2:
                odd_pairs += 1
                maximum_odd = max(maximum_odd, content)
    assert maximum_odd == 2 and maximum_nonzero == 3

    matrices = coefficient_matrices(4)
    alternating = 5
    row = (1 << alternating) | (1 << (alternating | 8))
    column = [(-1)**state.bit_count() - 2 * int(state == alternating) for state in range(16)]
    assert all(value in (-1, 1) for value in column)
    assert sum(column[state] for state in range(16) if row >> state & 1) == -2
    assert difference_content(matrices, row, column) == 12
    high_precision_column = [0, 2161744, 1759216, 2879648, 1447568, 1, 3606864, 13792,
                             3581056, 4179184, 633632, 3252192, 2959904, 450048, 1666656, 3691904]
    assert high_precision_column[5] == 1 and max(high_precision_column) < 1 << 22
    assert difference_content(matrices, 1 << 5, high_precision_column) == 22

    compiler = shlex.split(os.environ.get('CXX', 'c++'))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError('The optional canonical decoder checker requires C++17.')
    source = Path(__file__).with_name('canonical_packet_decoder.cpp')
    runs = []
    with tempfile.TemporaryDirectory(prefix='integer-factoring-canonical-decoder-') as directory:
        temporary = Path(directory)
        binary, data = temporary / 'decoder', temporary / 'coefficients.txt'
        run = subprocess.run(compiler + ['-O3', '-std=c++17', str(source), '-o', str(binary)],
                             text=True, capture_output=True)
        if run.returncode:
            raise RuntimeError(run.stderr[-12000:])
        with data.open('w') as output:
            output.write(str(len(matrices)) + '\n')
            for matrix in matrices:
                for values in matrix:
                    output.write(' '.join(map(str, values)) + '\n')
        for precision, mode in ((12, 'odd'), (13, 'nonzero'), (12, 'nonzero'), (11, 'odd')):
            run = subprocess.run([str(binary), str(data), str(precision), mode, '123456789123'],
                                 text=True, capture_output=True)
            if run.returncode:
                raise RuntimeError(run.stderr[-12000:])
            result = json.loads(run.stdout)
            assert result['precision'] == precision and result['mode'] == mode
            expected_found = (precision, mode) in ((12, 'nonzero'), (11, 'odd'))
            assert result['found'] == expected_found
            if expected_found:
                tensor = result['row_tensor_mask']
                difference = result['column_difference']
                original = sum(difference[state] for state in range(16) if tensor >> state & 1)
                assert original == result['original_difference']
                assert (original % 2 != 0) if mode == 'odd' else original != 0
                assert difference_content(matrices, tensor, difference) == precision
            else:
                assert result['row_tensors_tested'] == 65535
                assert result['right_half_queries'] == 65535 * 3**8
            runs.append(result)

    root = Path(__file__).resolve().parent.parent
    sources = ['tools/canonical_packet_decoder.cpp', 'tools/check_canonical_packet_decoder.py',
               'tools/check_canonical_packet_exhaustive.py', 'tools/cauchy_packet_polynomials.py',
               'research/canonical-packet-decoder-frontier.md']
    return {'status': 'pass', 'creation_precision': 1,
            'length_two_binary_row_ternary_column_pairs': 6**4,
            'length_two_nonzero_original_pairs': nonzero_pairs,
            'length_two_odd_original_pairs': odd_pairs,
            'length_two_exact_parity_precision': maximum_odd + 1,
            'length_two_exact_count_precision': maximum_nonzero + 1,
            'length_four_binary_row_ternary_column_pairs': 6**16,
            'length_four_exact_parity_precision': 12,
            'length_four_exact_count_precision': 13,
            'length_four_odd_pairs_covered': (6**16 - 2**16) // 2,
            'explicit_count_collision': {'row_tensor_mask': row, 'column_difference': column,
                                         'original_difference': -2, 'minimum_content': 12},
            'higher_creation_precision_control': {'row_tensor_mask': 1 << 5,
                'column_creation_precision': 22, 'column': high_precision_column,
                'original_value': 1, 'minimum_content': 22},
            'runs': runs,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in sources},
            'scope': 'Exact finite parity and full-count separation for canonical binary '
                     'columns with the binary row tensor fixed and known. The proof allows '
                     'an exhaustive decoder; no efficient all-size decoder or factoring '
                     'exponent improvement is established.'}


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
