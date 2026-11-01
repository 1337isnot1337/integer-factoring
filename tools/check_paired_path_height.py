#!/usr/bin/env python3
"""Exact fixed-height walk counts and higher q-splitting determinant controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from bounded_spin_crossing import suffix_counts as label_counts
from bounded_spin_support import suffix_counts as support_counts
from paired_path_crossing import MULTIPLY as M,power

ROOT = Path(__file__).resolve().parents[1]


def walks(vertices,start,length):
    values = [int(i == start) for i in range(vertices)]
    for _ in range(length):
        values = [(values[i-1] if i else 0)+(values[i+1] if i+1 < vertices else 0)
                  for i in range(vertices)]
    return values[start]


def adjacent_matrix_power(vertices,start,length):
    # Independent integer matrix powering, not the one-step walk recurrence.
    a = [[int(abs(i-j) == 1) for j in range(vertices)] for i in range(vertices)]
    result = [[int(i == j) for j in range(vertices)] for i in range(vertices)]

    def multiply(x,y):
        return [[sum(x[i][k]*y[k][j] for k in range(vertices)) for j in range(vertices)]
                for i in range(vertices)]

    while length:
        if length&1:
            result = multiply(result,a)
        a = multiply(a,a)
        length >>= 1
    return result[start][start]


def check():
    assert all(power(2,k) != 1 for k in range(1,15)) and power(2,15) == 1
    denominator = 2 ^ power(2,-1)

    def integer(k):
        return M[power(2,k) ^ power(2,-k)][power(denominator,-1)]

    blocks = 0
    for h in range(1,14):
        for m in range(1,h+1):
            determinant = M[power(2,m)][integer(h-m+1)] ^ M[power(2,-h+m-1)][integer(m)]
            assert determinant == integer(h+1) != 0
            blocks += 1
    assert integer(15) == 0
    matrix_controls = 0
    records = []
    for h in range(2,14):
        labels = [walks(h+1,0,n) for n in range(2,33,2)]
        support = [walks(2*h+1,h,n) for n in range(2,33,2)]
        assert all(s >= l for l,s in zip(labels,support))
        for n in (0,2,6,10):
            assert walks(h+1,0,n) == adjacent_matrix_power(h+1,0,n)
            assert walks(2*h+1,h,n) == adjacent_matrix_power(2*h+1,h,n)
            matrix_controls += 2
        if h == 3:
            assert labels == [label_counts(n)[n][0] for n in range(2,33,2)]
            assert support == [support_counts(n)[n][3] for n in range(2,33,2)]
        records.append({'height':h,'even_lengths':list(range(2,33,2)),
                        'label_counts':labels,'support_counts':support})
    files = ('tools/check_paired_path_height.py','tools/paired_path_crossing.py',
             'tools/bounded_spin_crossing.py','tools/bounded_spin_support.py',
             'tools/check_fourier_multiplication_gates.py')
    return {'status':'pass','local_F16_determinants':blocks,
            'q_order':15,'last_invertible_height':13,'first_forbidden_q_integer':15,
            'independent_adjacency_power_controls':matrix_controls,'walk_counts':records,
            'scope':'Finite exact height/formula controls; all-height entropy and rank rates are proved in the note.',
            'source_sha256':{f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
