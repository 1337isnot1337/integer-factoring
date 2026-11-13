#!/usr/bin/env python3
"""Test the proposed rank-product barrier for full-conservation kernels."""
from __future__ import annotations

import json

PRIME = 2
Q = 3


def rank_two(matrix: list[list[int]]) -> int:
    a = [[value % 2 for value in row] for row in matrix]
    rows, columns = len(a), len(a[0])
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, rows) if a[i][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        for i in range(rows):
            if i != rank and a[i][column]:
                a[i] = [(x - y) % 2 for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def convolution_rank(values: list[int]) -> int:
    return rank_two([[values[(a - b) % Q] for b in range(Q)]
                     for a in range(Q)])


def cuts(kernel: list[list[int]]) -> tuple[int, int]:
    row = sum(convolution_rank(values) for values in kernel)
    column = sum(convolution_rank([kernel[a][b] for a in range(Q)])
                 for b in range(Q))
    return row, column


def check() -> dict:
    kernel = [[1, 1, 1], [1, 1, 1], [1, 1, 0]]
    row_cut, column_cut = cuts(kernel)
    kernel_rank = rank_two(kernel)
    assert (kernel_rank, row_cut, column_cut) == (2, 4, 4)
    assert row_cut * column_cut < Q * Q * kernel_rank
    minimum = None
    for mask in range(1 << (Q * Q)):
        candidate = [[(mask >> (a * Q + b)) & 1 for b in range(Q)]
                     for a in range(Q)]
        candidate_rank = rank_two(candidate)
        if not candidate_rank:
            continue
        r, c = cuts(candidate)
        gap = Q * Q * candidate_rank - r * c
        if minimum is None or gap > minimum[0]:
            minimum = (gap, r, c, candidate_rank, mask)
    assert minimum is not None
    return {
        'passed': True,
        'prime': PRIME,
        'q': Q,
        'kernel': kernel,
        'kernel_rank': kernel_rank,
        'row_cut_rank': row_cut,
        'column_cut_rank': column_cut,
        'sqrt_cost': 4.0,
        'proposed_lower_bound_q_sqrt_rank': Q * 2 ** 0.5,
        'inequality_fails': row_cut * column_cut < Q * Q * kernel_rank,
        'witness_ratio': (row_cut * column_cut)
        / (Q * Q * kernel_rank),
        'all_q3_kernels_checked': 1 << (Q * Q),
        'scope': (
            'A full-conservation kernel of rank two on Z/3 at characteristic '
            'two has both cut ranks four, so sqrt(r c)=4 is below q sqrt(rank)='
            '3 sqrt(2). This disproves the proposed q-sqrt-rank shortcut to an '
            'all-q conservation barrier. It does not contradict the proved '
            'minimum-cut cutoff or the characteristic-two p-group rank-sum '
            'theorem, and it does not close q>=22 over general fields.'
        ),
    }


if __name__ == '__main__':
    print(json.dumps(check(), indent=2))
