#!/usr/bin/env python3
"""Exact controls for the all-field binary lemma and F4 boundary pencil.

The universal statements have written proofs. This checker uses small finite
fields, sampled free tensors, exact elimination, and integer radical tests.
It does not interpret numerical experiments as a proof over all fields.
"""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import isqrt
from pathlib import Path
from random import Random
import time


class Field:
    def __init__(self, q):
        assert q in (2, 3, 4)
        self.q = q

    def add(self, a, b):
        return a ^ b if self.q == 4 else (a+b) % self.q

    def neg(self, a):
        return a if self.q == 4 else (-a) % self.q

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        if self.q != 4:
            return a*b % self.q
        answer = 0
        while b:
            if b & 1:
                answer ^= a
            a <<= 1
            b >>= 1
            if a & 4:
                a ^= 7
        return answer

    def inv(self, a):
        assert a
        return next(b for b in range(1, self.q) if self.mul(a, b) == 1)

    def dot(self, a, b):
        result = 0
        for x, y in zip(a, b):
            result = self.add(result, self.mul(x, y))
        return result


def rank(a, field):
    a = [list(row) for row in a]
    used = 0
    for column in range(len(a[0])):
        pivot = next((i for i in range(used, len(a)) if a[i][column]), None)
        if pivot is None:
            continue
        a[used], a[pivot] = a[pivot], a[used]
        inverse = field.inv(a[used][column])
        a[used] = [field.mul(inverse, value) for value in a[used]]
        for i in range(len(a)):
            if i != used and a[i][column]:
                coefficient = a[i][column]
                a[i] = [field.sub(x, field.mul(coefficient, y))
                        for x, y in zip(a[i], a[used])]
        used += 1
    return used


QUADS = tuple(product(range(2), repeat=4))  # ordinary x,y,u,v order


def cuts(tensor, field):
    return (
        rank([[tensor[x, y, u, v] for u, v in product(range(2), repeat=2)]
              for x, y in product(range(2), repeat=2)], field),
        rank([[tensor[x, y, u, v] for y, v in product(range(2), repeat=2)]
              for x, u in product(range(2), repeat=2)], field))


def matmul(a, b, field):
    return [[field.dot(row, column) for column in zip(*b)] for row in a]


def kron(a, b, field):
    return [[field.mul(a[i][j], b[k][l])
             for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def two_roots_below_four(p, q):
    return p+q < 16 and 4*p*q < (16-p-q)**2


def pair_cost_at_least(p, q, target):
    if p+q >= target*target:
        return True
    return 4*p*q >= (target*target-p-q)**2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('binary-two-piece-f4-boundary.json'))
    args = parser.parse_args()
    started = time.monotonic()
    rng = Random(20261030)

    exceptional = []
    for r, c, s, t in product(range(1, 5), repeat=4):
        if r+s < 4 or c+t < 4:
            continue
        if any(a == 1 and b == 3 for a, b in [(r, c), (c, r), (s, t), (t, s)]):
            continue
        if two_roots_below_four(r*c, s*t):
            exceptional.append([r, c, s, t])
    assert exceptional == [[1, 2, 3, 2], [2, 1, 2, 3],
                           [2, 3, 2, 1], [3, 2, 1, 2]]

    swap = [[int(i//2 == j % 2 and i % 2 == j//2) for j in range(4)]
            for i in range(4)]
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    rank_one_controls = 0
    arbitrary_two_piece_controls = 0
    for q in (2, 3, 4):
        field = Field(q)
        for _ in range(128):
            u = [[rng.randrange(q) for _ in range(2)] for _ in range(2)]
            v = [[rng.randrange(q) for _ in range(2)] for _ in range(2)]
            a = {(x, y, z, w): field.mul(u[x][y], v[z][w]) for x, y, z, w in QUADS}
            r, c = cuts(a, field)
            assert c == rank(u, field)*rank(v, field)
            if r == 1:
                assert c in (1, 2, 4)
            t = matmul(swap, kron(u, v, field), field)
            square = matmul(t, t, field)
            assert square == kron(matmul(v, u, field), matmul(u, v, field), field)
            complement = [[field.sub(identity[i][j], t[i][j]) for j in range(4)]
                          for i in range(4)]
            if rank(t, field) == rank(complement, field) == 2:
                assert square == t
            if rank(u, field)*rank(v, field) == 2:
                assert rank(square, field) <= 1
                assert rank(complement, field) != 2
            rank_one_controls += 1

            a = {quad: rng.randrange(q) for quad in QUADS}
            b = {quad: field.sub(int(quad[0] == quad[3] and quad[1] == quad[2]), a[quad])
                 for quad in QUADS}
            ra, ca = cuts(a, field)
            rb, cb = cuts(b, field)
            assert pair_cost_at_least(ra*ca, rb*cb, 4)
            arbitrary_two_piece_controls += 1

    field = Field(4)
    omega = 2
    assert field.add(field.add(field.mul(omega, omega), omega), 1) == 0
    eta = [1, 0, 0, 0, 0]
    initial = [rank([[int(x == v) if i < 2 else 0 for v in range(2)]
                     for x in range(2)], field) for i in range(3)]
    assert initial == [2, 2, 0]
    terminal = [1, 1, 2]
    h = [1, 1, 2]
    endpoints = transitions = lower_rows = 0
    zero_attainment = None
    # Include the zero tensors and all scalar-multiple roots before arbitrary
    # free tensors; this exercises the zero-transition branches explicitly.
    families = []
    for scalar in range(4):
        fs = []
        for i in range(3):
            fs.append({(x, y, u, v): field.mul(scalar,
                       (int(y == u == 0) if i == 0 else int(y == u == 1) if i == 1
                        else int(y == u == 0) ^ (2 if y == u == 1 else 0)))
                       if x == v else 0 for x, y, u, v in QUADS})
        families.append(fs)
    for _ in range(64):
        families.append([{quad: rng.randrange(4) for quad in QUADS} for _ in range(3)])

    for number, fs in enumerate(families):
        tensors = [[{} for _ in range(3)] for _ in range(3)]
        for x, y, u, v in QUADS:
            delta = int(x == v)
            beta = [int(y == u == 0), int(y == u == 1),
                    int(y == u == 0) ^ (2 if y == u == 1 else 0)]
            pair = [int(2*x+v == j) for j in range(4)]
            ket = [[delta]+pair, [delta]+[field.mul(omega, a) for a in pair], [0]+pair]
            # b_i is evaluated at y,u; all pair coordinates are independent.
            bra = [[beta[i]]+[fs[i][xx, y, u, vv] for xx, vv in product(range(2), repeat=2)]
                   for i in range(3)]
            endpoint = [0]*5
            for i in range(3):
                for coordinate in range(5):
                    endpoint[coordinate] = field.add(endpoint[coordinate],
                        field.mul(ket[i][coordinate], beta[i]))
            assert endpoint == [int(x == v and y == u)]+[0]*4
            endpoints += 1
            for i in range(3):
                d = field.mul(delta, beta[i])
                f = fs[i][x, y, u, v]
                expected = [field.add(f, d), field.add(field.mul(omega, f), d), f]
                for j in range(3):
                    value = field.dot(bra[i], ket[j])
                    assert value == expected[j]
                    tensors[i][j][x, y, u, v] = value
                    transitions += 1
        products = []
        for i in range(3):
            row = []
            for j in range(3):
                r, c = cuts(tensors[i][j], field)
                row.append(r*c)
            products.append(row)
            # All possible equality cases here have integer square roots.
            # A floor gives a rigorous lower bound for every other value.
            precision = 32
            lower = sum(isqrt(p << (2*precision))*weight for p, weight in zip(row, h))
            assert lower >= 4*h[i]*(1 << precision)
            if i == 2:
                assert pair_cost_at_least(row[0], row[2], 4)
                assert pair_cost_at_least(row[1], row[2], 4)
            lower_rows += 1
        if number == 0:
            assert products == [[4, 4, 0], [4, 4, 0], [16, 16, 0]]
            zero_attainment = {'rank_products': products, 'live_states': [0, 1],
                              'live_matrix': [[2, 2], [2, 2]], 'exact_base': 4}

    result = {
        'exceptional_rank_profiles_before_projection_argument': exceptional,
        'fields_in_finite_controls': [2, 3, 4],
        'rank_one_Kronecker_and_projection_controls': rank_one_controls,
        'arbitrary_two_piece_controls': arbitrary_two_piece_controls,
        'F4_free_tensor_families_checked': len(families),
        'vector_endpoint_checks': endpoints,
        'transition_entry_checks': transitions,
        'exact_lower_supervector_rows': lower_rows,
        'lower_vector': h,
        'initial_boundary_ranks': initial,
        'terminal_boundary_ranks': terminal,
        'attaining_construction': zero_attainment,
        'scope': 'Finite exact controls. The all-field lemma and universal F4 base>=4 are proved in the companion note.',
        'seconds': time.monotonic()-started,
    }
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
