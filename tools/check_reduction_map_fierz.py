#!/usr/bin/env python3
"""Exact reduction-map Fierz identity and diagonal rank-one word transfer."""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import product
import json
from math import prod
from pathlib import Path
import time

from check_boundary_stable_catalytic import (add, apply, endpoint_reachable,
    identity, left_right, multiply, rank, rank_factor, unit, zero)


def data(d):
    b = d*d
    eta = sum(1 << (d*i+i) for i in range(d))
    W = tuple((1 << i) ^ (eta if eta >> i & 1 else 0) for i in range(b))
    assert multiply(W, W) == identity(b)
    assert apply(W, eta) == eta
    tables = [{}, {}]
    for x, y, u, v in product(range(d), repeat=4):
        q = x, y, u, v
        tables[0][q] = unit(b, d*x+v, d*y+u)
        tables[1][q] = multiply(multiply(W, unit(b, d*x+v, d*u+y)), W)
    return eta, W, tables


def extended_rows(table, d, column):
    rows = [0]*(d**4)
    for (x, y, u, v), M in table.items():
        left = x*d+(u if column else y)
        right = (y if column else u)*d+v
        for output, row in enumerate(M):
            a, c = divmod(output, d)
            while row:
                bit = row & -row
                b, e = divmod(bit.bit_length()-1, d)
                rows[left*d*d+a*d+b] |= 1 << (right*d*d+c*d+e)
                row ^= bit
    return rows


def transition(d, i, j, y, u, x, v):
    aa = int(y == x and u == v)
    bb = int(u == x and y == v)
    ee = int(y == u and x == v)
    return ((aa, aa ^ ee), (bb ^ ee, bb))[i][j]


def transition_ranks(d):
    answer = [[None]*2 for _ in range(2)]
    for i, j in product(range(2), repeat=2):
        row = [sum(transition(d, i, j, y, u, x, v) << (u*d+v)
                   for u, v in product(range(d), repeat=2))
               for y, x in product(range(d), repeat=2)]
        col = [sum(transition(d, i, j, y, u, x, v) << (y*d+v)
                   for y, v in product(range(d), repeat=2))
               for u, x in product(range(d), repeat=2)]
        answer[i][j] = rank(row), rank(col)
    return answer


def predicted(d, word, ranks):
    boundary = d*(d-1 if word[0] else 1)
    r = boundary*prod(ranks[i][j][0] for i, j in zip(word, word[1:]))
    c = boundary*prod(ranks[i][j][1] for i, j in zip(word, word[1:]))
    return r, c


def full_word_ranks(d, word, eta, tables):
    k = len(word)
    labels = list(product(range(d), repeat=k))
    n = len(labels)
    rows, columns = [0]*(n*n), [0]*(n*n)
    for ix, iy, iu, iv in product(range(n), repeat=4):
        x, y, u, v = labels[ix], labels[iy], labels[iu], labels[iv]
        state = eta
        for position in range(k-1, -1, -1):
            q = x[position], y[position], u[position], v[position]
            state = apply(tables[word[position]][q], state)
        if state & 1:
            rows[ix*n+iy] |= 1 << (iu*n+iv)
            columns[ix*n+iu] |= 1 << (iy*n+iv)
    return rank(rows), rank(columns)


def pair_add(a, b):
    return a[0]+b[0], a[1]+b[1]


def pair_multiply(a, b, q):
    return a[0]*b[0]+q*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def pair_power(a, exponent, q):
    answer = 1, 0
    for _ in range(exponent):
        answer = pair_multiply(answer, a, q)
    return answer


def check_markov(d, ranks):
    q = d*(d-1)//2
    checks = 0
    for k in range(1, 8):
        classes = {}
        total = 0, 0
        for word in product(range(2), repeat=k):
            counts = tuple(sum((a, b) == (i, j) for a, b in zip(word, word[1:]))
                           for i, j in product(range(2), repeat=2))
            kind = word[0], word[-1], counts
            profile = predicted(d, word, ranks)
            if kind in classes:
                assert classes[kind] == profile
            else:
                classes[kind] = profile
            switches = counts[1]+counts[2]
            coefficient = d**k*(d-1 if word[0] else 1)*q**(switches//2)
            cost = (0, coefficient) if switches % 2 else (coefficient, 0)
            # The formal pair is a+b sqrt(q); its square is the rank product.
            square = pair_multiply(cost, cost, q)
            assert square[1] == 0 and square[0] == profile[0]*profile[1]
            total = pair_add(total, cost)
        target = pair_multiply((d*d, 0), pair_power((d, d), k-1, q), q)
        assert total == target
        assert len(classes) <= 2*(k+1)**4
        checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    endpoint_checks = factor_rows = transition_checks = word_profiles = markov_checks = 0
    summaries = []
    for d in (2, 4):
        eta, W, tables = data(d)
        endpoint = {q: add(tables[0][q], tables[1][q]) for q in tables[0]}
        values = {q: int(q[0] == q[3] and q[1] == q[2]) for q in endpoint}
        for q, M in endpoint.items():
            assert apply(M, eta) == values[q]*eta
            endpoint_checks += 1
        accepted, reach_dimension = endpoint_reachable(endpoint, eta, 1, values)
        assert accepted and reach_dimension == 1
        qminus = d*(d-1)//2
        extended = []
        for table in tables:
            new_table = {q: multiply(multiply(W, M), W) for q, M in table.items()}
            profile = []
            for column, current in ((False, table), (True, new_table)):
                rows = extended_rows(current, d, column)
                U, V = rank_factor(rows)
                for coordinate, row in zip(U, rows):
                    reconstructed = 0
                    for j, basis_row in enumerate(V):
                        if coordinate >> j & 1:
                            reconstructed ^= basis_row
                    assert reconstructed == row
                    factor_rows += 1
                profile.append(len(V))
            extended.append(profile)
        assert extended == [[1, d*d*qminus], [d*d*qminus, 1]]
        tranks = transition_ranks(d)
        assert tranks == [[(1, d*d), (d*d, qminus)], [(qminus, d*d), (d*d, 1)]]
        for i, j, y, u, x, v in product(range(2), range(2), *([range(d)]*4)):
            ket = 1 << (d*x+v)
            if j:
                ket = apply(W, ket)
            covector = 1 << (d*(u if i else y)+(y if i else u))
            if i:
                covector = apply(W, covector)  # W is symmetric.
            actual = (covector & ket).bit_count() & 1
            assert actual == transition(d, i, j, y, u, x, v)
            transition_checks += 1
        for i, y, u, x, v in product(range(2), *([range(d)]*4)):
            assert transition(d, i, 0, y, u, x, v) ^ transition(d, i, 1, y, u, x, v) == (y == u and x == v)
        markov_checks += check_markov(d, tranks)
        if d == 2:
            for k in (1, 2, 3):
                for word in product(range(2), repeat=k):
                    assert full_word_ranks(d, word, eta, tables) == predicted(d, word, tranks)
                    word_profiles += 1
        summaries.append({'external_base': d, 'extended_cut_profiles': extended,
                          'transition_rank_profiles': tranks,
                          'raw_base': f'{2*d} sqrt({qminus})',
                          'markov_base': f'{d}+{d} sqrt({qminus})'})
    ranks4 = transition_ranks(4)
    a, b = predicted(4, (0, 0, 1), ranks4), predicted(4, (0, 1, 0), ranks4)
    assert a == (64, 384) and b == (384, 384)
    assert a[0]*b[1] != b[0]*a[1]
    combined_square = (a[0]+b[0])*(a[1]+b[1])
    p, q = a[0]*a[1], b[0]*b[1]
    assert combined_square-p-q > 0
    assert (combined_square-p-q)**2 > 4*p*q
    output = {'status': 'passed', 'matrix_vector_endpoint_checks': endpoint_checks,
              'exact_extended_factor_rows': factor_rows,
              'auxiliary_transition_identity_checks': transition_checks,
              'full_individual_word_rank_profiles': word_profiles,
              'exact_markov_path_and_type_checks': markov_checks,
              'symbol_count_only_negative_control': {'AAB': a, 'ABA': b, 'strict_cauchy_gap_verified': True},
              'profiles': summaries,
              'scope': 'Exact concrete crossing and sharper word accounting; its common terminal projection reduces its base to an ordinary transition-row seed.',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'seconds': time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
