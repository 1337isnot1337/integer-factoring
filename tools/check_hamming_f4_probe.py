#!/usr/bin/env python3
"""Replay an extension-field counterexample to the prime-field row-rank lemma."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_extension_descent import multiply, power
from check_hamming_joint_lower import ELL, source_basis

ROOT = Path(__file__).resolve().parents[1]


def rank(matrix, modulus):
    a = [row[:] for row in matrix]
    rows, columns = len(a), len(a[0])
    result = 0
    q = 1 << (modulus.bit_length()-1)
    for column in range(columns):
        pivot = next((i for i in range(result, rows) if a[i][column]), None)
        if pivot is None:
            continue
        a[result], a[pivot] = a[pivot], a[result]
        inverse = power(a[result][column], q-2, modulus)
        a[result] = [multiply(value, inverse, modulus) for value in a[result]]
        for i in range(rows):
            if i != result and a[i][column]:
                scale = a[i][column]
                a[i] = [x ^ multiply(scale, y, modulus) for x, y in zip(a[i], a[result])]
        result += 1
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = json.loads((ROOT/'evidence/hamming-joint-f4-probe.json').read_text())
    generators = source_basis()
    coefficient = int(result['coefficient_mask'], 16)
    values = [(coefficient >> (2*j)) & 3 for j in range(len(generators))]
    matrix = [[0]*16 for _ in range(16)]
    separator = 0
    for j, (generator, value) in enumerate(zip(generators, values)):
        if (ELL >> j) & 1:
            separator ^= value
        for row in range(16):
            for column in range(16):
                if (generator >> (16*row+column)) & 1:
                    matrix[row][column] ^= value
    other = [[matrix[4*x+y][4*u+v] for y in range(4) for v in range(4)]
             for x in range(4) for u in range(4)]
    actual = rank(matrix, result['field_modulus']), rank(other, result['field_modulus'])
    assert separator == result['separator_value'] == 1
    assert actual == (result['row_rank'], result['column_rank']) == (4, 10)
    assert actual[0]*actual[1] == result['rank_product'] == 40
    paths = [Path(__file__), ROOT/'tools/check_extension_descent.py',
             ROOT/'tools/check_hamming_joint_lower.py']
    result['sha256'] = {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'sha256'}))


if __name__ == '__main__':
    main()
