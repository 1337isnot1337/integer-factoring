#!/usr/bin/env python3
"""Check Smith channels, Toeplitz ranks, and virtual ring factor queries.

Dense reference algorithms establish the finite examples independently.
The bounded oracle samples recompute the existing weighted-Gram rank and
prefix selection queries without memoization. Python memory is not a
measurement of the polylogarithmic workspace theorem.
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
if not (HERE / "check_polylog_space_rank.py").exists():
    sys.path.insert(0, str(HERE.parent / "integer-factoring" / "tools"))

from check_polylog_space_rank import Counters, rank_query, requested_coordinate


def ring_add(a, b, p):
    return tuple((x+y) % p for x, y in zip(a, b))


def ring_neg(a, p):
    return tuple(-x % p for x in a)


def ring_mul(a, b, p):
    return tuple(sum(a[j]*b[i-j] for j in range(i+1)) % p for i in range(len(a)))


def valuation(a):
    return next((i for i, x in enumerate(a) if x), len(a))


def divide_by_pivot(a, pivot, p):
    """Return one quotient when valuation(a)>=valuation(pivot)."""
    h, exponent = len(a), valuation(pivot)
    assert exponent < h and valuation(a) >= exponent
    answer = [0]*h
    inverse_lead = pow(pivot[exponent], -1, p)
    for i in range(h-exponent):
        residual = a[i+exponent]-sum(pivot[exponent+j]*answer[i-j] for j in range(1, i+1))
        answer[i] = residual*inverse_lead % p
    answer = tuple(answer)
    assert ring_mul(pivot, answer, p) == a
    return answer


def smith_valuations(matrix, p):
    """Independent row/column elimination over the principal chain ring."""
    a = [[tuple(x) for x in row] for row in matrix]
    m, n, h = len(a), len(a[0]), len(a[0][0])
    exponents = []
    for k in range(min(m, n)):
        exponent, i, j = min((valuation(a[i][j]), i, j)
                            for i in range(k, m) for j in range(k, n))
        if exponent == h:
            break
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        pivot = a[k][k]
        for i in range(k+1, m):
            multiplier = divide_by_pivot(a[i][k], pivot, p)
            for j in range(k, n):
                a[i][j] = ring_add(a[i][j], ring_neg(ring_mul(multiplier, a[k][j], p), p), p)
        for j in range(k+1, n):
            multiplier = divide_by_pivot(a[k][j], pivot, p)
            for i in range(m):
                a[i][j] = ring_add(a[i][j], ring_neg(ring_mul(multiplier, a[i][k], p), p), p)
        exponents.append(exponent)
    assert exponents == sorted(exponents)
    return exponents


def field_rank(matrix, p):
    a = [list(row) for row in matrix]
    if not a:
        return 0
    m, n, rank = len(a), len(a[0]), 0
    for j in range(n):
        pivot = next((i for i in range(rank, m) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][j], -1, p)
        a[rank] = [x*inverse % p for x in a[rank]]
        for i in range(rank+1, m):
            multiplier = a[i][j]
            a[i] = [(x-multiplier*y) % p for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == m:
            break
    return rank


def field_solution(matrix, rhs, p):
    m, n = len(matrix), len(matrix[0])
    a = [list(row)+[value] for row, value in zip(matrix, rhs)]
    pivots = []
    for j in range(n):
        pivot = next((i for i in range(len(pivots), m) if a[i][j]), None)
        if pivot is None:
            continue
        rank = len(pivots)
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][j], -1, p)
        a[rank] = [x*inverse % p for x in a[rank]]
        for i in range(m):
            if i != rank:
                multiplier = a[i][j]
                a[i] = [(x-multiplier*y) % p for x, y in zip(a[i], a[rank])]
        pivots.append(j)
    assert all(any(row[:n]) or not row[-1] for row in a)
    answer = [0]*n
    for i, j in enumerate(pivots):
        answer[j] = a[i][-1]
    return answer


def toeplitz(matrix, s, shift=0):
    m, n = len(matrix), len(matrix[0])
    return [[matrix[i][j][a-b-shift] if a >= b+shift else 0
             for j in range(n) for b in range(s)]
            for i in range(m) for a in range(s)]


def module_prefix(matrix, j):
    h = len(matrix[0][0])
    background = toeplitz(matrix, h, 1)
    return [row+[matrix[i][column][a] for column in range(j)]
            for i in range(len(matrix)) for a, row in enumerate(background[i*h:(i+1)*h])]


def dense_factorization(matrix, p):
    m, n, h = len(matrix), len(matrix[0]), len(matrix[0][0])
    previous = field_rank(module_prefix(matrix, 0), p)
    selected = []
    for j in range(n):
        current = field_rank(module_prefix(matrix, j+1), p)
        if current > previous:
            selected.append(j)
        previous = current
    u = [[row[j] for j in selected] for row in matrix]
    r = len(selected)
    v = [[None]*n for _ in range(r)]
    if r:
        linearized = toeplitz(u, h)
        for j in range(n):
            rhs = [matrix[i][j][a] for i in range(m) for a in range(h)]
            solution = field_solution(linearized, rhs, p)
            for a in range(r):
                v[a][j] = tuple(solution[a*h:(a+1)*h])
    return selected, u, v


class VirtualFactor:
    """No stored pivot list; every U entry reselects its original column."""

    def __init__(self, matrix, p, counters):
        self.matrix, self.p, self.counters = matrix, p, counters
        self.m, self.n, self.h = len(matrix), len(matrix[0]), len(matrix[0][0])
        self.base_rank = rank_query(self.prefix_entry, self.m*self.h, self.n*self.h, p, counters)
        total = rank_query(self.prefix_entry, self.m*self.h, self.n*self.h+self.n, p, counters)
        self.r = total-self.base_rank

    def prefix_entry(self, row, column):
        i, a = divmod(row, self.h)
        if column >= self.n*self.h:
            return self.matrix[i][column-self.n*self.h][a]
        j, b = divmod(column, self.h)
        return self.matrix[i][j][a-b-1] if a > b else 0

    def selected_index(self, ordinal):
        low, high = 1, self.n
        while low < high:
            middle = (low+high)//2
            value = rank_query(self.prefix_entry, self.m*self.h, self.n*self.h+middle,
                               self.p, self.counters)-self.base_rank
            if value >= ordinal+1:
                high = middle
            else:
                low = middle+1
        return low-1

    def u_coefficient(self, i, j, a):
        return self.matrix[i][self.selected_index(j)][a]

    def linearized_u_entry(self, row, column):
        i, a = divmod(row, self.h)
        j, b = divmod(column, self.h)
        return self.u_coefficient(i, j, a-b) if a >= b else 0

    def v_coefficient(self, i, j, a):
        rhs = lambda row: self.matrix[row//self.h][j][row % self.h]
        return requested_coordinate(self.linearized_u_entry, rhs, lambda _: 0,
                                    self.m*self.h, self.r*self.h, self.p,
                                    i*self.h+a, self.counters)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    rng = random.Random(20261030)
    ring4 = list(product(range(2), repeat=2))
    cases = [(2, [[values[0], values[1]], [values[2], values[3]]])
             for values in product(ring4, repeat=4)]
    cases += [(2, [[a, b]]) for a, b in product(product(range(2), repeat=3), repeat=2)]
    for p in (2, 3, 5):
        for h in (1, 2, 3, 4):
            for m, n in ((1, 3), (3, 1), (2, 3), (3, 2), (3, 3)):
                cases.append((p, [[tuple(rng.randrange(p) for _ in range(h))
                                   for _ in range(n)] for _ in range(m)]))
    rank_identities = factor_entries = generator_checks = 0
    distributions = {}
    for p, matrix in cases:
        m, n, h = len(matrix), len(matrix[0]), len(matrix[0][0])
        exponents = smith_valuations(matrix, p)
        ranks = [0]
        for s in range(1, h+1):
            rank = field_rank(toeplitz(matrix, s), p)
            assert rank == sum(max(s-e, 0) for e in exponents)
            ranks.append(rank)
            rank_identities += 1
        assert ranks[h]-ranks[h-1] == len(exponents)
        selected, u, v = dense_factorization(matrix, p)
        assert len(selected) == len(exponents)
        generator_checks += 1
        zero = (0,)*h
        for i, j in product(range(m), range(n)):
            actual = zero
            for k in range(len(selected)):
                actual = ring_add(actual, ring_mul(u[i][k], v[k][j], p), p)
            assert actual == matrix[i][j]
            factor_entries += 1
        if p == 2 and h == 2 and m == n == 2:
            key = str(tuple(exponents))
            distributions[key] = distributions.get(key, 0)+1

    # Exhaustively certify minimality through one channel over F2[t]/t².
    one_channel = set()
    for u, v in product(product(ring4, repeat=2), repeat=2):
        one_channel.add(tuple(ring_mul(u[i], v[j], 2)
                              for i, j in product(range(2), repeat=2)))
    minimum_checks = 0
    for values in product(ring4, repeat=4):
        rank = len(smith_valuations([list(values[:2]), list(values[2:])], 2))
        assert (values in one_channel) == (rank <= 1)
        minimum_checks += 1

    # Actual recomputed weighted-Gram rank queries on small nonreduced examples.
    counters = Counters()
    samples = [([[(0, 1), (1, 1)]], 2),
               ([[(0, 1), (0, 0)]], 2),
               ([[(0, 0), (0, 0)]], 2),
               ([[(0, 0, 1), (0, 1, 1)]], 2)]
    oracle_generator_checks = oracle_factor_checks = oracle_coordinate_checks = 0
    for matrix, p in samples:
        virtual = VirtualFactor(matrix, p, counters)
        selected, u, v = dense_factorization(matrix, p)
        assert virtual.r == len(selected)
        assert [virtual.selected_index(i) for i in range(virtual.r)] == selected
        oracle_generator_checks += 1
        for i, j, a in product(range(virtual.m), range(virtual.r), range(virtual.h)):
            assert virtual.u_coefficient(i, j, a) == u[i][j][a]
            oracle_factor_checks += 1
        # Nested coefficient recovery is checked on the order-two example;
        # larger examples already have the independent dense factor checks.
        if matrix == samples[0][0]:
            for i, j, a in product(range(virtual.r), range(virtual.n), range(virtual.h)):
                assert virtual.v_coefficient(i, j, a) == v[i][j][a]
                oracle_coordinate_checks += 1

    nilpotent_diagonal = [[(0, 1), (0, 0)], [(0, 0), (0, 1)]]
    assert smith_valuations(nilpotent_diagonal, 2) == [1, 1]
    assert field_rank(toeplitz(nilpotent_diagonal, 1), 2) == 0
    assert field_rank(toeplitz(nilpotent_diagonal, 2), 2) == 2
    nilpotent_scalar = [[(0, 1)]]
    scalar_linearized = toeplitz(nilpotent_scalar, 2)
    impossible_rhs = [row+[int(i == 0)] for i, row in enumerate(scalar_linearized)]
    possible_rhs = [row+[int(i == 1)] for i, row in enumerate(scalar_linearized)]
    assert len(smith_valuations(nilpotent_scalar, 2)) == 1
    assert len(smith_valuations([[(0, 1), (1, 0)]], 2)) == 1
    assert field_rank(scalar_linearized, 2) == field_rank(possible_rhs, 2) == 1
    assert field_rank(impossible_rhs, 2) == 2
    assert dense_factorization([[(0, 1), (1, 0)]], 2)[0] == [1]
    output = {"status": "passed", "cases": len(cases),
              "toeplitz_smith_rank_identities": rank_identities,
              "generator_prefix_checks": generator_checks,
              "factor_entry_checks": factor_entries,
              "exhaustive_minimum_channel_checks": minimum_checks,
              "one_channel_matrices_over_f2_t2": len(one_channel),
              "f2_t2_2x2_smith_profiles": distributions,
              "oracle_generator_checks": oracle_generator_checks,
              "oracle_u_coefficient_checks": oracle_factor_checks,
              "oracle_v_coefficient_checks": oracle_coordinate_checks,
              "oracle_counters": vars(counters),
              "nonzero_nilpotent_rank_control": True,
              "ring_channel_count_consistency_negative_control": True,
              "full_image_radical_generator_control": True,
              "scope": "Finite exact identities; Python memory is not the workspace proof.",
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "seconds": time.monotonic()-started}
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
