#!/usr/bin/env python3
"""A binary qutrit Weyl pencil and its exact constant-coefficient obstruction."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path


def field4_multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def polynomial_multiply(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def binary_rank(rows, low=False):
    pivots = {}
    for value in rows:
        while value:
            pivot = ((value & -value).bit_length() - 1
                     if low else value.bit_length() - 1)
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def weyl_entry(a, b, row, column):
    return (1, 2, 3)[(b * column) % 3] if (row - column - a) % 3 == 0 else 0


def local_pencil():
    terms, comparisons = [], 0
    for x, y, u, v in product(range(3), repeat=4):
        matrix = [[0] * 3 for _ in range(3)]
        for r, s in product(range(3), repeat=2):
            literal = 0
            for a, b in product(range(3), repeat=2):
                physical = field4_multiply(weyl_entry(a, b, x, y),
                                           weyl_entry(a, b, u, v))
                literal ^= field4_multiply(physical, weyl_entry(a, b, r, s))
            closed = int((x - y - u + v) % 3 == 0
                         and (r - s - x + y) % 3 == 0
                         and (s + y + v) % 3 == 0)
            assert literal == closed
            crossing = int(x == v and y == u and r == s)
            matrix[r][s] = closed | ((closed ^ crossing) << 1)
            assert ((matrix[r][s] & 1) ^ (matrix[r][s] >> 1)) == crossing
            comparisons += 1
        if any(any(row) for row in matrix):
            terms.append(((x, y, u, v), matrix))
    return terms, comparisons


def index(word):
    result = 0
    for digit in word:
        result = 3 * result + digit
    return result


def constant_rows(t):
    d = 3 ** t
    rows, columns = [0] * (d * d), [0] * (d * d)
    words = tuple(product(range(3), repeat=t))
    for x, y in product(words, repeat=2):
        differences = [(a - b) % 3 for a, b in zip(x, y)]
        if sum(differences) % 3:
            continue
        total = 0
        u, v = [], []
        for a, yi in zip(differences, y):
            total = (total + a) % 3
            vi = (total - yi) % 3
            v.append(vi)
            u.append((vi + a) % 3)
        rows[index(x) * d + index(y)] = 1 << (index(u) * d + index(v))
    for x, u in product(words, repeat=2):
        if (x[0] + u[0]) % 3:
            continue
        states = [0] + [(-x[i] - u[i]) % 3 for i in range(1, t)] + [0]
        y, v = [], []
        for i, (xi, ui) in enumerate(zip(x, u)):
            difference = (states[i] - states[i + 1]) % 3
            y.append((xi - difference) % 3)
            v.append((ui - difference) % 3)
        columns[index(x) * d + index(u)] = 1 << (index(y) * d + index(v))
    return rows, columns


def profile(terms, t):
    d = 3 ** t
    width = d * d
    rows = [[0] * width for _ in range(t + 1)]
    columns = [[0] * width for _ in range(t + 1)]
    for factors in product(terms, repeat=t):
        state, coordinates = [1, 0, 0], [0] * 4
        for physical, matrix in factors:
            coordinates = [3 * a + b for a, b in zip(coordinates, physical)]
            state = [polynomial_multiply(state[0], matrix[0][s])
                     ^ polynomial_multiply(state[1], matrix[1][s])
                     ^ polynomial_multiply(state[2], matrix[2][s])
                     for s in range(3)]
        value = state[0]
        x, y, u, v = coordinates
        for k in range(t + 1):
            if value >> k & 1:
                rows[k][x * d + y] ^= 1 << (u * d + v)
                columns[k][x * d + u] ^= 1 << (y * d + v)
    endpoint = [0] * width
    for matrix in rows:
        endpoint = [a ^ b for a, b in zip(endpoint, matrix)]
    assert endpoint == [1 << ((i % d) * d + i // d) for i in range(width)]
    assert (rows[0], columns[0]) == constant_rows(t)
    ranks = [[binary_rank(row), binary_rank(column)]
             for row, column in zip(rows, columns)]
    assert ranks == [[binary_rank(row, low=True), binary_rank(column, low=True)]
                     for row, column in zip(rows, columns)]
    assert ranks[0] == [3 ** (2 * t - 1)] * 2
    return {'length': t, 'labels': d, 'coefficient_rank_pairs': ranks,
            'constant_coefficient_rank': 3 ** (2 * t - 1),
            'coefficient_endpoint_entries_checked': width * width,
            'constant_closed_form_entries_checked_both_cuts': 2 * width * width}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    terms, entries = local_pencil()
    profiles = [profile(terms, t) for t in (1, 2, 3)]
    result = {
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'physical_alphabet': 3, 'auxiliary_dimension': 3,
        'scalar_field': 'F2',
        'weyl_derivation_field': 'F4 = F2[w]/(w^2+w+1)',
        'literal_weyl_sum_entries_checked': entries,
        'nonzero_physical_local_tuples': len(terms),
        'boundary': [1, 0, 0],
        'profiles': profiles,
        'all_length_constant_rank': '3^(2*t-1) in both cuts',
        'all_length_raw_rank_cost_lower_bound': 'd^2/3, d=3^t',
        'excludes_arbitrary_coefficient_recombinations': False,
        'scope': 'This exact linear Weyl pencil with fixed e0 boundaries and raw monomial coefficient splitting only. No obstruction to different Weyl paths, boundary conditions, or factoring algorithms is claimed.',
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
