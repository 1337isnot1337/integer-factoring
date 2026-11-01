#!/usr/bin/env python3
"""Exact finite controls for observable transition-table realization."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
from random import Random
import time

from check_binary_two_piece_f4_boundary import Field, QUADS, cuts, matmul, rank


def inverse(a, field):
    n = len(a)
    a = [row[:]+[int(i == j) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        scale = field.inv(a[j][j])
        a[j] = [field.mul(scale, x) for x in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                scale = a[i][j]
                a[i] = [field.sub(x, field.mul(scale, y)) for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


def matvec(a, v, field):
    return [field.dot(row, v) for row in a]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    field = Field(2)
    w = [[1, 1, 1], [1, 1, 0], [0, 1, 1]]
    wi = inverse(w, field)
    boundary = [1, 1, 1]
    assert [row[1] for row in w] == wi[1] == boundary
    table = {}
    basis = []
    for q in QUADS:
        x, y, u, v = q
        e, r, c = int(x == v and y == u), int(x == y and u == v), int(x == u and y == v)
        t0 = [[0, r, 0], [0, e, c], [0, 0, 0]]
        table[q] = matmul(matmul(w, t0, field), wi, field)
        assert table[q] == [[e ^ r ^ c, e ^ r ^ c, e ^ r],
                            [e ^ r ^ c, e ^ r ^ c, e ^ r],
                            [e ^ c, e ^ c, e]]
        basis.append([e, r, c])
    assert rank(basis, field) == 3

    eigen_constraints = []
    left_eigen_constraints = []
    for q, t in table.items():
        e = int(q[0] == q[3] and q[1] == q[2])
        a = [[t[i][j] ^ (e if i == j else 0) for j in range(3)] for i in range(3)]
        eigen_constraints.extend(a)
        left_eigen_constraints.extend(zip(*a))
    assert rank(eigen_constraints, field) == rank(left_eigen_constraints, field) == 3

    reachable = [boundary+[1]]
    for vector in reachable:
        for q, t in table.items():
            e = int(q[0] == q[3] and q[1] == q[2])
            candidate = matvec(t, vector[:3], field)+[e*vector[3]]
            if rank(reachable+[candidate], field) > len(reachable):
                reachable.append(candidate)
    assert len(reachable) == 2
    assert all(field.dot(boundary+[1], v) == 0 for v in reachable)

    word_checks = 0
    for length in range(4):
        for words in product(QUADS, repeat=length):
            vector = boundary[:]
            expected = 1
            for q in reversed(words):
                vector = matvec(table[q], vector, field)
                expected *= int(q[0] == q[3] and q[1] == q[2])
            assert field.dot(boundary, vector) == expected
            word_checks += 1

    profiles = [[cuts({q: table[q][i][j] for q in QUADS}, field)
                 for j in range(3)] for i in range(3)]
    assert profiles == [[(2, 2), (2, 2), (4, 1)],
                        [(2, 2), (2, 2), (4, 1)],
                        [(1, 4), (1, 4), (4, 4)]]
    g = [[2, 2, 2], [2, 2, 2], [2, 2, 4]]
    # Twice the rational positive vector (1,1,3/2).
    assert all(sum(a*b for a, b in zip(row, [2, 2, 3])) <= 7*h
               for row, h in zip(g, [2, 2, 3]))
    assert matvec(g, [1, -1, 0], Field(3)) == [0, 0, 0]

    eta = [1]+[0]*12
    phi = [1]+[int(x == v) for _ in range(3) for x, v in product(range(2), repeat=2)]
    kets = {(i, x, v): [int(j == 1+4*i+2*x+v) for j in range(13)]
            for i in range(3) for x, v in product(range(2), repeat=2)}
    bras = {(i, y, u): [int(y == u)]+[table[x, y, u, v][i][j]
             for j in range(3) for x, v in product(range(2), repeat=2)]
            for i in range(3) for y, u in product(range(2), repeat=2)}
    lift_transition_checks = 0
    for x, y, u, v in QUADS:
        for i, j in product(range(3), repeat=2):
            assert field.dot(bras[i, y, u], kets[j, x, v]) == table[x, y, u, v][i][j]
            lift_transition_checks += 1
    initial = [rank([[field.dot(phi, kets[i, x, v]) for v in range(2)]
                     for x in range(2)], field) for i in range(3)]
    terminal = [rank([[field.dot(bras[i, y, u], eta) for u in range(2)]
                      for y in range(2)], field) for i in range(3)]
    assert initial == terminal == [2, 2, 2]

    rng = Random(20261030)
    lift_word_checks = 0
    for length in range(1, 5):
        words_to_check = list(product(QUADS, repeat=length)) if length <= 2 else [
            tuple(rng.choice(QUADS) for _ in range(length)) for _ in range(128)]
        for words in words_to_check:
            vector = eta[:]
            expected = 1
            for x, y, u, v in reversed(words):
                following = [0]*13
                for i in range(3):
                    scalar = field.dot(bras[i, y, u], vector)
                    if scalar:
                        following[1+4*i+2*x+v] ^= 1
                vector = following
                expected *= int(x == v and y == u)
            assert field.dot(phi, vector) == expected
            lift_word_checks += 1
    # The chosen lift boundary is not a common right eigenvector.
    x = y = u = v = 0
    actual = [0]*13
    for i in range(3):
        actual[1+4*i] ^= field.dot(bras[i, y, u], eta)
    assert actual != eta and field.dot(phi, actual) == 1

    protected_units = 0
    for i in range(3):
        for j in range(i, 3):
            unit = [[int(a == i and b == j) for b in range(3)] for a in range(3)]
            conjugate = matmul(matmul(w, unit, field), wi, field)
            assert conjugate[2][2] == int(i == j == 1)
            protected_units += 1
    result = {
        'field': 'F2', 'states': 3, 'lift_auxiliary_dimension': 13,
        'right_E_eigenvector_constraint_rank': rank(eigen_constraints, field),
        'left_E_eigenvector_constraint_rank': rank(left_eigen_constraints, field),
        'augmented_reachable_dimension': len(reachable),
        'all_matrix_words_through_length3': word_checks,
        'entry_cut_rank_pairs': profiles,
        'transition_cost_matrix': g,
        'exact_spectral_base': '4+2sqrt(2)',
        'rational_upper_certificate': {'beta': 7, 'twice_positive_vector': [2, 2, 3]},
        'lift_transition_entry_checks': lift_transition_checks,
        'lift_word_endpoint_checks': lift_word_checks,
        'initial_boundary_ranks': initial, 'terminal_boundary_ranks': terminal,
        'protected_diagonal_matrix_units_checked': protected_units,
        'scope': 'Exact diagnostic construction outside common E eigenvectors; not an improving seed.',
        'seconds': time.monotonic()-started,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
