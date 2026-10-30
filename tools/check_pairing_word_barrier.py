#!/usr/bin/env python3
"""Exact local blocks and common-functionals for pairing type recombinations."""

import argparse
from hashlib import sha256
from math import comb
import json
from pathlib import Path

from check_modular_recombination import crossing, insert, represent
from check_word_symmetry_barrier import log_bounds


ZERO = set()
ONE = {(0, 0, 0)}
R, H, C = {(1, 0, 0)}, {(0, 1, 0)}, {(0, 0, 1)}
D = R | H | C


def multiply(left, right):
    result = set()
    for a in left:
        for b in right:
            monomial = tuple(x+y for x, y in zip(a, b))
            result.symmetric_difference_update((monomial,))
    return result


def power(polynomial, exponent):
    result = ONE
    while exponent:
        if exponent & 1:
            result = multiply(result, polynomial)
        polynomial = multiply(polynomial, polynomial)
        exponent >>= 1
    return result


def row_sum(rows):
    entries = [ONE]
    for row in rows:
        entries = [multiply(a, b) for a in entries for b in row]
    result = set()
    for entry in entries:
        result.symmetric_difference_update(entry)
    return result


def local_blocks(d):
    n = d*d
    diagonal = sum(1 << (i*d+i) for i in range(d))
    r = [diagonal if i//d == i%d else 0 for i in range(n)]
    e, identity = crossing(d), [1 << i for i in range(n)]
    h = [a ^ b ^ c for a, b, c in zip(identity, e, r)]
    basis, blocks, independent = [], [], {}
    def add(vector):
        assert insert(independent, vector, 1 << len(basis))
        basis.append(vector)
    if d % 2:
        add(diagonal)
        blocks.append(('special_odd', 1))
        for i in range(d-1):
            add((1 << (i*d+i)) ^ (1 << (n-1)))
            blocks.append(('scalar', 1))
    else:
        add(diagonal)
        add(1 << (n-1))
        blocks.append(('special_even', 2))
        for i in range(d-1):
            vector = (1 << (i*d+i)) ^ (1 << (n-1))
            if represent(independent, vector) is None:
                add(vector)
                blocks.append(('scalar', 1))
    for i in range(d):
        for j in range(i+1, d):
            add((1 << (i*d+j)) ^ (1 << (j*d+i)))
            add(1 << (i*d+j))
            blocks.append(('jordan', 2))
    assert len(basis) == n
    checked = 0
    for cut, generators in [('row', (r, h, e)), ('column', (e, h, r))]:
        for letter, matrix in zip(('r', 'h', 'c'), generators):
            columns = [represent(independent, sum(((row & v).bit_count()%2) << i
                                                  for i, row in enumerate(matrix)))
                       for v in basis]
            assert all(v is not None for v in columns)
            actual = [sum(((columns[j] >> i) & 1) << j for j in range(n)) for i in range(n)]
            expected, start = [0]*n, 0
            for kind, length in blocks:
                if kind == 'special_odd':
                    block = [1]
                elif kind == 'scalar':
                    block = [int(letter == ('c' if cut == 'row' else 'r'))]
                elif kind == 'special_even':
                    block = [1, 2] if letter == ('c' if cut == 'row' else 'r') else [2, 0]
                else:
                    full = 'c' if cut == 'row' else 'r'
                    block = [3, 2] if letter == full else [2, 0] if letter == 'h' else [0, 0]
                for j, value in enumerate(block):
                    expected[start+j] = value << start
                start += length
            assert actual == expected
            checked += n*n
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    block_entries = sum(local_blocks(d) for d in range(2, 11))
    functional_cases = 0
    for t in range(1, 9):
        for a in range(t+1):
            m = t-a
            target = multiply(power(D, a), power(H, m))
            for odd in (True, False):
                row_special = [D] if odd else [C, R | H]
                col_special = [D] if odd else [R, C | H]
                row = row_sum([row_special]*a+[[C, H | C]]*m)
                column = row_sum([col_special]*a+[[R, H | R]]*m)
                assert row == column == target and len(target)%2 == 1
                functional_cases += 1
    multiplicity_cases = 0
    for d in range(2, 17):
        l = comb(d, 2)
        for t in range(1, 65):
            weights = [comb(t, a)*l**(t-a) for a in range(t+1)]
            assert sum(weights) == (l+1)**t
            assert max(weights)*(t+1) >= (l+1)**t
            multiplicity_cases += 1
    intervals = []
    for d in range(2, 9):
        l = comb(d, 2)
        lo, hi = log_bounds(l+1)
        dlo, dhi = log_bounds(d)
        lower, upper = lo/(2*dhi), hi/(2*dlo)
        assert (2*lower*lower > 1) if d >= 5 else (2*upper*upper < 1)
        scale = 10**12
        a = lower.numerator*scale//lower.denominator
        b = -((-upper.numerator*scale)//upper.denominator)
        intervals.append({'seed_dimension': d, 'rho_lower': f'{a//scale}.{a%scale:012d}',
                          'rho_upper': f'{b//scale}.{b%scale:012d}',
                          'half_threshold_excluded_by_this_bound': d >= 5})
    for d in range(6, 1001):
        assert (d*d-d+2)**2-4*d**3 == d**3*(d-6)+5*d*d-4*d+4 > 0
    result = {'status': 'all-size characteristic-two pairing-type obstruction; no factoring lower bound',
              'local_block_entries_checked': block_entries,
              'common_functional_cases': functional_cases,
              'multiplicity_cases': multiplicity_cases, 'rho_intervals': intervals,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
