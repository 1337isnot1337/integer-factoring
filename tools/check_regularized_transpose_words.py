#!/usr/bin/env python3
"""Independent packed-polynomial construction controls for transpose sandwiches."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random
import time

from check_rectangular_coordinate_support import rank_f4
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import F4_PRODUCT, flatten
from regularized_transpose_words import INVERSE, construct
from weighted_smith_oracle_f4 import WeightedSmithOracle, scalar_multiply
from windowed_f4_smith import pack


def shear_matrix(q, shears):
    a = [[int(i == j) for j in range(q)] for i in range(q)]
    for i, j, c in shears:
        for row in a:
            row[j] ^= F4_PRODUCT[c][row[i]]
    return a


def invert(a):
    n = len(a)
    rows = [row[:]+[int(i == j) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        i = next(i for i in range(j, n) if rows[i][j])
        rows[i], rows[j] = rows[j], rows[i]
        rows[j] = [F4_PRODUCT[INVERSE[rows[j][j]]][v] for v in rows[j]]
        for i in range(n):
            if i != j:
                c = rows[i][j]
                rows[i] = [v ^ F4_PRODUCT[c][w] for v, w in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def apply_factor(matrix, a, inverse, power):
    q = len(a)
    temporary, twisted = [[0]*q for _ in range(q)], [[0]*q for _ in range(q)]
    for i, j, k in product(range(q), repeat=3):
        temporary[i][j] ^= scalar_multiply(matrix[j][k], a[i][k])
    for i, j, k in product(range(q), repeat=3):
        twisted[i][j] ^= scalar_multiply(temporary[i][k], inverse[j][k])
    return [[matrix[i][j] ^ twisted[i][j] ^ (twisted[i][j] << (2*power))
             for j in range(q)] for i in range(q)]


def independent_columns(seed, metadata, columns):
    q = metadata['q']
    diagonal_matrices = [[[c if i == j else 0 for j, c in enumerate(diagonal)]
                          for i in range(q)] for diagonal in metadata['diagonals']]
    gate_pair = lambda gate: (shear_matrix(q, gate['shears']), gate.get('power', 1)) if isinstance(gate, dict) else (shear_matrix(q, gate), 1)
    word = ([gate_pair(s) for s in reversed(metadata['prepend'])]
            + list(zip(diagonal_matrices, metadata['base_powers']))
            + [gate_pair(s) for s in metadata['append']])
    factors = [(a, invert(a), power) for a, power in reversed(word)]
    checked = 0
    for v, u in columns:
        matrix = [[int(i == v and j == u) for j in range(q)] for i in range(q)]
        for a, inverse, power in factors:
            matrix = apply_factor(matrix, a, inverse, power)
        for x, y in product(range(q), repeat=2):
            assert matrix[x][y] == pack(seed['values'][x+q*(y+q*(u+q*v))])
            checked += 1
    return checked


def random_gates(rng, q, count, shears):
    result = []
    for _ in range(count):
        gate = []
        for _ in range(shears):
            i, j = rng.sample(range(q), 2)
            gate.append((i, j, rng.randrange(1, 4)))
        result.append(gate)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, begun = Random(121220261103), time.monotonic()
    records, entries, orders, endpoints = [], 0, 0, 0
    for q, count, shears in ((3, 2, 5), (4, 2, 7), (12, 8, 48)):
        pre, post = random_gates(rng, q, count, shears), random_gates(rng, q, count, shears)
        base_powers = (1, 1, 5) if q == 12 else (2, 2, 2) if q == 4 else (1, 2, 2)
        seed, metadata = construct(q, pre, post, base_powers=base_powers)
        columns = list(product(range(q), repeat=2)) if q <= 4 else [(0, 0), (0, 1), (3, 9), (11, 7)]
        entries += independent_columns(seed, metadata, columns)
        endpoints += q**4
        parameters = [0.]*(8*(q-1))+[1.]
        native = RectangularOracle(seed['values'], q, q).prefixes(parameters, 1)['orders']
        assert all(v == int(v) for v in native)
        n, d = q*q, metadata['degree']
        assert sum(native[n:]) == d*n
        assert min(native[:n]) == -d and native[:n].count(-d) == 1
        assert int(native[n:].count(0) > 0) == metadata['zero_R_rank']
        if metadata['zero_R_rank']:
            assert native[n:].count(0) == 1
        independent = []
        if q <= 4:
            for cut in (1, 0):
                matrix, _, _ = flatten(seed, cut)
                values, _ = WeightedSmithOracle([[pack(v) for v in row] for row in matrix], infinity=bool(cut)).profile([0]*n, [0]*n)
                independent.extend(values)
            assert native == independent
            orders += len(independent)
        for cut, power, key in ((0, 0, 'zero_R'), (1, d, 'top_C')):
            matrix, _, _ = flatten(seed, cut)
            coefficient = [[((low >> power) & 1) | (((high >> power) & 1) << 1)
                            for low, high in row] for row in matrix]
            assert rank_f4(coefficient) == metadata[key+'_rank']
            pivot = next(((i, j) for i in range(n) for j in range(n) if coefficient[i][j]), None)
            vector_ranks = [0, 0]
            if pivot:
                i, j = pivot
                column = [coefficient[k][j] for k in range(n)]
                row = coefficient[i]
                vector_ranks = [rank_f4([v[k*q:(k+1)*q] for k in range(q)]) for v in (column, row)]
            assert vector_ranks == metadata[key+'_vector_matrix_ranks']
        # Verify vector-matrix ranks with independent restriction of scalars.
        for value in range(256):
            # Small deterministic field-rank controls include rank zero and full rank.
            a = [[(value >> (2*(i+2*j))) & 3 for j in range(2)] for i in range(2)]
            from regularized_transpose_words import field_rank
            assert field_rank(a) == rank_f4(a)
        records.append({'metadata': metadata, 'complete_native_orders': native,
                        'independent_small_orders': independent,
                        'selected_input_columns': columns})
    names = ['regularized_transpose_words.py', 'check_regularized_transpose_words.py',
             'check_rectangular_coordinate_support.py', 'physical_split_crossings.py',
             'weighted_smith_oracle_f4.py', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp']
    report = {'status': 'pass', 'independent_packed_polynomial_entries': entries,
              'complete_endpoints': endpoints, 'independent_untruncated_small_orders': orders,
              'complete_native_orders': sum(len(v['complete_native_orders']) for v in records),
              'independent_leading_matrix_rank_checks': 6,
              'all_two_by_two_F4_rank_checks': 768,
              'even_degree_odd_factor_control': {'q': 4, 'degree': 10, 'factors': 7},
              'records': records, 'seconds': time.monotonic()-begun,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Exact finite construction, endpoints and profile controls. The q12 native '
                       'profile has selected independent polynomial-entry controls, not a full second '
                       '144-channel untruncated replay. No favorable pressure or score is claimed.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'independent_packed_polynomial_entries',
                     'independent_untruncated_small_orders', 'complete_native_orders', 'seconds')}))


if __name__ == '__main__':
    main()
