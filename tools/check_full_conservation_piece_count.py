#!/usr/bin/env python3
"""Checked reduction of the odd full-conservation kernel lane.

For a q-square scalar kernel K over F_p, put

    T_K(x,y,u,v) = [x+y=u+v] K(x+y, x+u).

The two standard cut ranks are the row and column convolution-rank sums.
This script checks those identities against dense tensors, verifies the
per-piece floor r*c >= q^2, and evaluates the natural two-piece catalogue
at the odd alphabets used in the Issue #52 search.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from hashlib import sha256
from pathlib import Path


def rank_mod(matrix, modulus):
    a = [[int(value) % modulus for value in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((row for row in range(rank, rows) if a[row][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][col], -1, modulus)
        a[rank] = [value * inverse % modulus for value in a[rank]]
        for row in range(rows):
            if row != rank and a[row][col]:
                factor = a[row][col]
                a[row] = [
                    (x - factor * y) % modulus for x, y in zip(a[row], a[rank])
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def convolution_rank(function, q, modulus):
    return rank_mod(
        [[function[(a + b) % q] for b in range(q)] for a in range(q)],
        modulus,
    )


def cut_ranks(kernel, q, modulus):
    row_rank = sum(
        convolution_rank([kernel[S * q + w] for w in range(q)], q, modulus)
        for S in range(q)
    )
    column_rank = sum(
        convolution_rank([kernel[S * q + w] for S in range(q)], q, modulus)
        for w in range(q)
    )
    return row_rank, column_rank


def dense_cut_ranks(kernel, q, modulus):
    zero = 0
    row_cut = [
        [zero for _ in range(q * q)]
        for _ in range(q * q)
    ]
    column_cut = [
        [zero for _ in range(q * q)]
        for _ in range(q * q)
    ]
    for x, y, u, v in itertools.product(range(q), repeat=4):
        S = (x + y) % q
        if S != (u + v) % q:
            continue
        value = kernel[S * q + (x + u) % q]
        row_cut[x * q + y][u * q + v] = value
        column_cut[x * q + u][y * q + v] = value
    return rank_mod(row_cut, modulus), rank_mod(column_cut, modulus)


def target(q):
    return q ** math.sqrt(2)


def piece_count_table():
    table = []
    for q in (5, 7, 11, 13, 17, 19, 23, 29):
        threshold = target(q)
        maximum = 0
        while (maximum + 1) * q < threshold:
            maximum += 1
        table.append({
            'q': q,
            'target': threshold,
            'two_q': 2 * q,
            'three_q': 3 * q,
            'maximum_nonzero_pieces_if_floor_q': maximum,
        })
    return table


def primitive_root_of_order(q, modulus):
    for candidate in range(2, modulus):
        if pow(candidate, q, modulus) != 1:
            continue
        if all(pow(candidate, q // p, modulus) != 1
               for p in range(2, q + 1) if q % p == 0):
            return candidate
    raise ArithmeticError(f'no primitive {q}-th root modulo {modulus}')


def catalogue(q, modulus, seed):
    rng = random.Random(seed)
    identity = [1 if i // q == i % q else 0 for i in range(q * q)]
    kernels = []

    def add(name, values):
        kernels.append({'name': name, 'values': [value % modulus for value in values]})

    for i, j in itertools.product(range(q), repeat=2):
        values = [0] * (q * q)
        values[i * q + j] = 1
        add(f'entry_{i}_{j}', values)

    root = primitive_root_of_order(q, modulus)
    for row, frequency in itertools.product(range(q), repeat=2):
        add(
            f'row_character_{row}_{frequency}',
            [pow(root, frequency * w, modulus) if S == row else 0
             for S in range(q) for w in range(q)],
        )
    for column, frequency in itertools.product(range(q), repeat=2):
        add(
            f'column_character_{column}_{frequency}',
            [pow(root, frequency * S, modulus) if w == column else 0
             for S in range(q) for w in range(q)],
        )
    for row_frequency, column_frequency in itertools.product(range(q), repeat=2):
        add(
            f'character_product_{row_frequency}_{column_frequency}',
            [pow(root, row_frequency * S + column_frequency * w, modulus)
             for S in range(q) for w in range(q)],
        )
    for frequency in range(1, q):
        add(
            f'chirp_{frequency}',
            [pow(root, frequency * S * w, modulus)
             for S in range(q) for w in range(q)],
        )
    for shift in range(q):
        add(
            f'cyclic_permutation_{shift}',
            [1 if w == (S + shift) % q else 0
             for S in range(q) for w in range(q)],
        )

    masks = set()
    masks.add((1 << q) - 1)
    masks.add(0)
    for width in range(1, q):
        for start in range(q):
            masks.add(sum(1 << ((start + i) % q) for i in range(width)))
    for _ in range(96):
        masks.add(rng.randrange(1 << q))
    for mask in sorted(masks):
        add(
            f'diagonal_mask_{mask}',
            [mask >> S & 1 if S == w else 0 for S in range(q) for w in range(q)],
        )

    results = []
    minimum = None
    for item in kernels:
        kernel = item['values']
        complement = [
            (identity[i] - kernel[i]) % modulus for i in range(len(kernel))
        ]
        left_row, left_column = cut_ranks(kernel, q, modulus)
        right_row, right_column = cut_ranks(complement, q, modulus)
        total = math.sqrt(left_row * left_column) + math.sqrt(
            right_row * right_column
        )
        record = {
            'name': item['name'],
            'left': [left_row, left_column],
            'right': [right_row, right_column],
            'total': total,
        }
        results.append(record)
        if minimum is None or total < minimum['total']:
            minimum = record
    return {
        'q': q,
        'count': len(results),
        'minimum': minimum,
        'target': target(q),
        'passes_target': minimum['total'] < target(q),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    identity_checks = []
    for q in (3, 5, 7):
        modulus = 1000003
        rng = random.Random(100 + q)
        for trial in range(4 if q == 7 else 8):
            kernel = [rng.randrange(modulus) for _ in range(q * q)]
            expected = cut_ranks(kernel, q, modulus)
            actual = dense_cut_ranks(kernel, q, modulus)
            if expected != actual:
                raise AssertionError((q, trial, expected, actual))
            identity_checks.append({
                'q': q,
                'trial': trial,
                'expected': expected,
                'actual': actual,
            })

    floor_checks = []
    for q, modulus, trials in (
        (3, 2, None),
        (3, 5, None),
        (5, 11, 400),
        (7, 29, 200),
        (13, 53, 80),
    ):
        if trials is None:
            kernels = []
            for index in range(1, modulus ** (q * q)):
                value = index
                kernel = []
                for _ in range(q * q):
                    kernel.append(value % modulus)
                    value //= modulus
                kernels.append(kernel)
        else:
            rng = random.Random(200 + q)
            kernels = [
                [rng.randrange(modulus) for _ in range(q * q)]
                for _ in range(trials)
            ]
        for kernel in kernels:
            row_rank, column_rank = cut_ranks(kernel, q, modulus)
            if row_rank and row_rank * column_rank < q * q:
                raise AssertionError(('piece floor failed', q, row_rank, column_rank))
        floor_checks.append({
            'q': q,
            'modulus': modulus,
            'kernels': len(kernels),
            'floor': q * q,
        })

    counts = piece_count_table()
    expected_counts = {5: 1, 7: 2, 11: 2, 13: 2, 17: 3, 19: 3, 23: 3, 29: 4}
    assert {row['q']: row['maximum_nonzero_pieces_if_floor_q'] for row in counts} == expected_counts

    catalogues = []
    for q, modulus in ((7, 29), (11, 23), (13, 53), (17, 103), (19, 191)):
        record = catalogue(q, modulus, 300 + q)
        if record['passes_target']:
            raise AssertionError(('catalogue crossed the target', q, record))
        catalogues.append(record)

    payload = {
        'scope': (
            'Odd-alphabet full-conservation kernel lane. The per-piece floor '
            'r*c >= q^2 is general; the two-piece catalogue is finite and does '
            'not exclude arbitrary two-piece decompositions. No factoring '
            'exponent or complete-bound change.'
        ),
        'cut_rank_identity_controls': identity_checks,
        'piece_floor_controls': floor_checks,
        'piece_count_reduction': counts,
        'two_piece_catalogue': catalogues,
        'status': 'PASS',
    }
    raw = json.dumps(payload, separators=(',', ':'), sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + '\n')
    print(json.dumps({
        'status': 'PASS',
        'identity_controls': len(identity_checks),
        'floor_controls': sum(row['kernels'] for row in floor_checks),
        'catalogue_minima': {
            row['q']: row['minimum']['total'] for row in catalogues
        },
        'sha256': sha256(raw.encode()).hexdigest(),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
