#!/usr/bin/env python3
"""Check that the ten-site coefficient table does not lift to larger primes.

The rank=trace identity over F5 would give exact cut-module ranks if the
coefficient table were the reduction of an integral table over a prime
larger than the module dimension.  This checker rebuilds the ten-site
entry-functional table over several primes with full-rank witnesses.
The characteristic-five table is commutative and unital; the larger-prime
tables fail the same algebra identities.  This is a negative diagnostic,
not a rank bound.
"""
from __future__ import annotations

import argparse
import json
import random
import time
from math import comb
from pathlib import Path

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry

LENGTH = 10
BASIS = [0, 2, 3, 4, 5, 6, 7, 10]
PRIMES = [5, 7, 11, 10007]


def original_form(four_paths, prime):
    n = len(four_paths[0])
    result = [0] * (n + 1)
    for (_, b, c), value in monomial_entry(four_paths, prime=prime).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % prime
    return result


def rank_mod(matrix, prime):
    a = [[int(x) % prime for x in row] for row in matrix]
    rows = len(a)
    columns = len(a[0]) if rows else 0
    r = 0
    for column in range(columns):
        pivot = next((i for i in range(r, rows) if a[i][column]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inverse = pow(a[r][column], -1, prime)
        a[r] = [x * inverse % prime for x in a[r]]
        for i in range(rows):
            if i != r and a[i][column]:
                factor = a[i][column]
                a[i] = [(x - factor * y) % prime for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


def solve_mod(matrix, rhs, prime):
    a = [[int(x) % prime for x in row] for row in matrix]
    b = [int(x) % prime for x in rhs]
    rows, columns = len(a), len(a[0])
    augmented = [row + [value] for row, value in zip(a, b)]
    r = 0
    pivots = []
    for column in range(columns):
        pivot = next((i for i in range(r, rows) if augmented[i][column]), None)
        if pivot is None:
            continue
        augmented[r], augmented[pivot] = augmented[pivot], augmented[r]
        inverse = pow(augmented[r][column], -1, prime)
        augmented[r] = [x * inverse % prime for x in augmented[r]]
        for i in range(rows):
            if i != r and augmented[i][column]:
                factor = augmented[i][column]
                augmented[i] = [(x - factor * y) % prime
                                for x, y in zip(augmented[i], augmented[r])]
        pivots.append(column)
        r += 1
        if r == rows:
            break
    for i in range(r, rows):
        if not any(augmented[i][:columns]) and augmented[i][-1]:
            raise ValueError('inconsistent system')
    solution = [0] * columns
    for i, column in enumerate(pivots):
        solution[column] = augmented[i][-1] % prime
    return solution


def witness(ps, prime, seed):
    q = len(ps)
    rng = random.Random(seed)
    quads = []
    while len(quads) < len(BASIS):
        ids = tuple(rng.randrange(q) for _ in range(4))
        vector = original_form(tuple(ps[i] for i in ids), prime)
        current = [original_form(tuple(ps[i] for i in t), prime) for t in quads]
        if rank_mod(current + [vector], prime) > len(quads):
            quads.append(ids)
    return quads


def table(ps, prime, quads):
    q = len(ps)
    entry = [[original_form(tuple(ps[i] for i in ids), prime)[k]
              for k in BASIS] for ids in quads]
    rows, columns = [], []
    for ids in quads:
        rows.append([original_form((ps[ids[0]], ps[ids[1]], ps[a], ps[b]),
                                   prime) for a in range(q) for b in range(q)])
        columns.append([original_form((ps[a], ps[b], ps[ids[2]], ps[ids[3]]),
                                      prime) for a in range(q) for b in range(q)])
    unit_rhs = [int((ids[0] * q + ids[1]) == (ids[2] * q + ids[3]))
                for ids in quads]
    unit = solve_mod(entry, unit_rhs, prime)
    result = [[[0] * len(BASIS) for _ in BASIS] for _ in BASIS]
    for i in range(len(BASIS)):
        for j in range(len(BASIS)):
            rhs = [sum(rows[l][s][BASIS[i]] * columns[l][s][BASIS[j]]
                       for s in range(q * q)) % prime for l in range(len(BASIS))]
            result[i][j] = solve_mod(entry, rhs, prime)
    commutative = all(result[i][j] == result[j][i]
                      for i in range(len(BASIS)) for j in range(len(BASIS)))
    unit_ok = True
    for i in range(len(BASIS)):
        product = [sum(result[i][k][j] * unit[k] for k in range(len(BASIS))) % prime
                   for j in range(len(BASIS))]
        if product != [int(j == i) for j in range(len(BASIS))]:
            unit_ok = False
    return {'prime': prime, 'witness_rank': rank_mod(entry, prime),
            'commutative': commutative, 'unit_ok': unit_ok,
            'unit': unit, 'table': result}


def check():
    started = time.time()
    ps = paths(LENGTH)
    results = []
    for prime in PRIMES:
        quads = witness(ps, prime, 7)
        results.append(table(ps, prime, quads))
    by_prime = {item['prime']: item for item in results}
    assert by_prime[5]['witness_rank'] == len(BASIS)
    assert by_prime[5]['commutative'] and by_prime[5]['unit_ok']
    for prime in PRIMES[1:]:
        assert by_prime[prime]['witness_rank'] == len(BASIS)
        assert not (by_prime[prime]['commutative'] and by_prime[prime]['unit_ok'])
    return {
        'status': 'pass',
        'length': LENGTH,
        'basis_indices': BASIS,
        'primes': PRIMES,
        'results': [{'prime': item['prime'],
                     'witness_rank': item['witness_rank'],
                     'commutative': item['commutative'],
                     'unit_ok': item['unit_ok']} for item in results],
        'seconds': round(time.time() - started, 1),
        'scope': ('Negative diagnostic: the characteristic-five ten-site coefficient '
                  'table is not the reduction of an integral table over larger primes '
                  'with the same basis; the larger-prime tables fail commutativity or '
                  'unitality. No cut-rank bound or score change is claimed.'),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
