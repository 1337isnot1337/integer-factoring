#!/usr/bin/env python3
"""Independent exhaustive-minor replay for entry-dependent Laurent shifts.

The reference clears one global Laurent power and enumerates determinants
over the exact polynomial ring. It does not use the native elimination.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path
import random

from native_matrix_smith_f4 import MatrixOracle
from weighted_smith_oracle_f4 import polynomial_multiply


def planes(value):
    low = high = 0
    for degree in range((value.bit_length() + 1) // 2):
        coefficient = value >> (2 * degree) & 3
        low |= (coefficient & 1) << degree
        high |= (coefficient >> 1) << degree
    return low, high


def determinant(matrix):
    result = 0
    for permutation in permutations(range(len(matrix))):
        term = 1
        for row, column in enumerate(permutation):
            term = polynomial_multiply(term, matrix[row][column])
        result ^= term
    return result


def genuinely_entry_dependent(shifts, matrix):
    n = len(matrix)
    return any(matrix[i][j] and matrix[i][l] and matrix[k][j] and matrix[k][l]
               and shifts[i][j] + shifts[k][l] != shifts[i][l] + shifts[k][j]
               for i, k in combinations(range(n), 2)
               for j, l in combinations(range(n), 2))


def reference_minors(matrix, shifts):
    n = len(matrix)
    floor = min(shifts[i][j] for i in range(n) for j in range(n) if matrix[i][j])
    cleared = [[value << (2 * (shifts[i][j] - floor))
                for j, value in enumerate(row)] for i, row in enumerate(matrix)]
    result = []
    examined = 0
    for size in range(1, n + 1):
        nonzero = {}
        for rows in combinations(range(n), size):
            for columns in combinations(range(n), size):
                value = determinant([[cleared[i][j] for j in columns] for i in rows])
                examined += 1
                if value:
                    zero = ((value & -value).bit_length() - 1) // 2 + size * floor
                    infinity = -((value.bit_length() - 1) // 2) - size * floor
                    nonzero[(rows, columns)] = (zero, infinity)
        result.append(nonzero)
    rank = max(size for size, values in enumerate(result, 1) if values)
    return result, rank, examined


def check():
    rng = random.Random(91073)
    cases = []
    for n in (2, 3, 4):
        for _ in range(8):
            while True:
                matrix = [[rng.randrange(1, 1 << rng.randrange(2, 9, 2))
                           for _ in range(n)] for _ in range(n)]
                shifts = [[rng.randrange(-7, 8) for _ in range(n)] for _ in range(n)]
                if genuinely_entry_dependent(shifts, matrix):
                    break
            cases.append((matrix, shifts))

    # This polynomial matrix is the product of a 3-by-2 and a 2-by-3
    # matrix, and its 2-by-2 upper-left minor is nonzero. Its rank is two.
    # Removing each entry's own leading power gives nonseparable shifts.
    for n in (3, 4):
        raw = [[4, 5, 5], [21, 17, 0], [1, 0, 17]]
        if n == 4:
            raw = [row + [0] for row in raw] + [[0, 0, 0, 4]]
        matrix = []
        shifts = []
        for row in raw:
            normalized = []
            bases = []
            for value in row:
                exponent = ((value & -value).bit_length() - 1) // 2 if value else 0
                normalized.append(value >> (2 * exponent))
                bases.append(exponent - 3)
            matrix.append(normalized)
            shifts.append(bases)
        cases.append((matrix, shifts))

    largest_error = 0.0
    minor_count = 0
    prefix_count = 0
    records = []
    for index, (matrix, shifts) in enumerate(cases):
        n = len(matrix)
        assert genuinely_entry_dependent(shifts, matrix)
        minors, rank, examined = reference_minors(matrix, shifts)
        minor_count += examined
        oracle = MatrixOracle([[planes(value) for value in row] for row in matrix], rank, shifts)
        rows = [Fraction(rng.randrange(-17, 18), 19) for _ in range(n)]
        columns = [Fraction(rng.randrange(-17, 18), 23) for _ in range(n)]
        scale = Fraction(3, 11)
        for infinity in (False, True):
            prefixes = [Fraction(0)]
            for values in minors[:rank]:
                prefixes.append(min(scale * valuations[int(infinity)]
                                    + sum(rows[i] for i in rr)
                                    + sum(columns[j] for j in cc)
                                    for (rr, cc), valuations in values.items()))
            expected = [right - left for left, right in zip(prefixes, prefixes[1:])]
            actual = oracle.profile(rows, columns, scale, infinity)
            error = max(abs(float(a) - b) for a, b in zip(expected, actual['orders']))
            largest_error = max(largest_error, error)
            assert error < 1e-10, (index, infinity, expected, actual)
            for size in range(1, rank + 1):
                rr = tuple(sorted(actual['rows'][:size]))
                cc = tuple(sorted(actual['cols'][:size]))
                expected_base = minors[size - 1][(rr, cc)][int(infinity)]
                assert actual['base_prefix'][size - 1] == expected_base
                prefix_count += 1
        records.append({'order': n, 'rank': rank, 'entry_dependent_shifts': True})

    assert records[-2]['rank'] == 2 and records[-1]['rank'] == 3
    root = Path(__file__).resolve().parent.parent
    source_names = ['tools/native_matrix_smith_f4.py', 'tools/native_matrix_smith_f4.cpp',
                    'tools/native_rectangular_smith_f4.cpp',
                    'tools/check_matrix_f4_entry_shifts.py']
    return {'status': 'pass', 'random_seed': 91073, 'matrix_cases': len(cases),
            'both_cuts_profiles': 2 * len(cases), 'all_minors_examined': minor_count,
            'attaining_prefix_minors_checked': prefix_count,
            'maximum_floating_error': largest_error,
            'rank_deficient_entry_dependent_cases': 2,
            'records': records,
            'source_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in source_names},
            'scope': 'Exact all-minor controls at order at most four for genuinely '
                     'entry-dependent Laurent shifts, at zero and infinity. Rank is '
                     'proved by exhaustive exact minors; floating weights remain exploratory.'}


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
