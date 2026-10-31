#!/usr/bin/env python3
"""Exact finite controls for the all-field degree-one Smith obstruction."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import subprocess
import tempfile


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            column = row.bit_length() - 1
            if column in pivots:
                row ^= pivots[column]
            else:
                pivots[column] = row
                break
    return len(pivots)


def reference_orders(matrix):
    """Recover orders from ranks over F_2[t]/t^k, without polynomial minors."""
    cumulative, previous = [], 0
    for precision in range(1, 6):
        rows = []
        for degree in range(precision):
            for i in range(4):
                row = 0
                for input_degree in range(degree + 1):
                    for j in range(4):
                        if matrix[i][j] >> (degree - input_degree) & 1:
                            row |= 1 << (4 * input_degree + j)
                rows.append(row)
        rank = binary_rank(rows)
        cumulative.append(rank - previous)
        previous = rank
    assert cumulative[-1] == 4
    result, previous = [], 0
    for order, count in enumerate(cumulative):
        result.extend([order] * (count - previous))
        previous = count
    return tuple(result)


def scalar_pencil(mask):
    infinity = [[0] * 4 for _ in range(4)]
    zero = [[0] * 4 for _ in range(4)]
    for x, y, u, v in product(range(2), repeat=4):
        index = ((2 * x + y) * 2 + u) * 2 + v
        a, e = mask >> index & 1, int(x == v and y == u)
        b = a ^ e
        infinity[2 * x + y][2 * u + v] = b | a << 1
        zero[2 * x + u][2 * y + v] = a | b << 1
    return infinity, zero


def profile(mask):
    return ':'.join(''.join(map(str, reference_orders(matrix))) for matrix in scalar_pencil(mask))


def prime_rank(matrix, prime):
    work = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = pow(work[rank][column], -1, prime)
        work[rank] = [value * scale % prime for value in work[rank]]
        for i in range(rank + 1, len(work)):
            scale = work[i][column]
            work[i] = [(a - scale * b) % prime for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


def check_small_rank_lemma(prime):
    matrices = list(product(range(prime), repeat=4))
    counts = {}
    for X, Y in product(matrices, repeat=2):
        matrix = [[int(x == v and y == u) - X[2 * x + u] * Y[2 * y + v]
                   for u, v in product(range(2), repeat=2)]
                  for x, y in product(range(2), repeat=2)]
        rank = prime_rank(matrix, prime)
        assert rank in (1, 3, 4)
        inverse_pair = all(sum(Y[2 * i + k] * X[2 * k + j] for k in range(2)) % prime == int(i == j)
                           for i, j in product(range(2), repeat=2))
        assert (rank == 1) == inverse_pair
        counts[str(rank)] = counts.get(str(rank), 0) + 1
    return {"prime": prime, "matrix_pairs": len(matrices) ** 2, "rank_counts": counts}


def partitions(total, maximum=None):
    if total == 0:
        yield ()
        return
    if maximum is None:
        maximum = total
    for first in range(min(maximum, total), 0, -1):
        for rest in partitions(total - first, first):
            yield (first,) + rest


def check_integral_budget():
    checks = 0
    for q in (2, 3, 4):
        n = q * q
        for total in range(n + 1):
            for positive in partitions(total):
                zeros = n - len(positive)
                for x in map(Fraction, ('1', '1/2', '1/3', '2/3', '1/10')):
                    actual = zeros + sum(x ** value for value in positive)
                    lower = n * x + zeros * (1 - x) ** 2
                    assert actual >= lower
                    checks += 1
    return checks


def equality_masks():
    result = []
    for a, b, c, d in product(range(2), repeat=4):
        if (a * d - b * c) % 2 != 1:
            continue
        S, inverse = (a, b, c, d), (d, b, c, a)
        mask = sum((S[2 * x + u] * inverse[2 * y + v]) << (((2 * x + y) * 2 + u) * 2 + v)
                   for x, y, u, v in product(range(2), repeat=4))
        result.append(mask)
    return sorted(result)


def run():
    source = Path(__file__).with_name('enumerate_binary_linear_crossings.cpp')
    with tempfile.TemporaryDirectory(prefix='linear-smith-') as directory:
        executable = str(Path(directory) / 'enumerator')
        subprocess.run(['c++', '-O3', '-std=c++17', str(source), '-o', executable], check=True)
        catalogue = json.loads(subprocess.check_output([executable], text=True))
    assert catalogue['pencils'] == sum(row['count'] for row in catalogue['profiles']) == 65536
    assert len(catalogue['profiles']) == 60
    assert catalogue['minimizers'] == equality_masks()
    assert len(catalogue['minimizers']) == 6
    records = {row['profile']: row for row in catalogue['profiles']}
    assert records['0112:0112']['count'] == 6
    reference_checks = 0
    for row in catalogue['profiles']:
        assert profile(row['witness']) == row['profile']
        reference_checks += 1
        first, second = (tuple(map(int, part)) for part in row['profile'].split(':'))
        assert first.count(0) + second.count(0) >= 2
        assert sum(first) <= 4 and sum(second) <= 4
        if row['profile'] != '0112:0112':
            assert first.count(0) != 1 or second.count(0) != 1
    randomizer = random.Random(20261030)
    for mask in sorted(set(randomizer.randrange(65536) for _ in range(1024)) | set(equality_masks())):
        assert profile(mask) in records
        reference_checks += 1
    assert all(profile(mask) == '0112:0112' for mask in equality_masks())
    x = Fraction(1, 3)
    assert (1 + x) ** 4 / x == Fraction(256, 27)
    assert (1 + 3 * x) ** 2 / x == 12
    return {
        'status': 'pass', 'catalogue': catalogue,
        'independent_toeplitz_profile_checks': reference_checks,
        'integral_budget_checks': check_integral_budget(),
        'rank_lemma_controls': [check_small_rank_lemma(p) for p in (2, 3)],
        'binary_sharp_B': '256/27', 'odd_characteristic_scalar_sharp_B': '12',
        'degree_one_all_alphabet_lower_B': 'q^3',
        'enumerator_sha256': sha256(source.read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Bare degree-one Smith certificates; no factoring lower bound or improved exponent.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = run()
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
