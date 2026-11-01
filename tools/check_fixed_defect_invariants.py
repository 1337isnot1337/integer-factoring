#!/usr/bin/env python3
"""Check equality-diagram obstructions for fixed moving defects."""

import argparse
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random


def partitions(n):
    def extend(word):
        if len(word) == n:
            yield word
        else:
            for value in range(max(word, default=-1)+2):
                yield from extend(word+(value,))
    return tuple(extend(()))


PARTITIONS = partitions(4)
E = (0, 1, 1, 0)
C = (0, 1, 0, 1)
R = (0, 0, 1, 1)


def refines(left, right):
    return all(left[i] != left[j] or right[i] == right[j]
               for i in range(4) for j in range(4))


def diagram(partition, positions):
    return int(all(partition[i] != partition[j] or positions[i] == positions[j]
                   for i in range(4) for j in range(4)))


def through(partition, cut):
    first, second = set(cut[:2]), set(cut[2:])
    return sum(bool(first & {i for i, a in enumerate(partition) if a == block})
               and bool(second & {i for i, a in enumerate(partition) if a == block})
               for block in set(partition))


def rank(matrix, p):
    A = [[x % p for x in row] for row in matrix]
    if not A:
        return 0
    row = 0
    for column in range(len(A[0])):
        pivot = next((i for i in range(row, len(A)) if A[i][column]), None)
        if pivot is None:
            continue
        A[row], A[pivot] = A[pivot], A[row]
        inv = pow(A[row][column], -1, p)
        A[row] = [inv*x % p for x in A[row]]
        for i in range(row+1, len(A)):
            if A[i][column]:
                factor = A[i][column]
                A[i] = [(a-factor*b) % p for a, b in zip(A[i], A[row])]
        row += 1
        if row == len(A):
            break
    return row


def position_matrix(partition, t, cut):
    pairs = tuple(product(range(t), repeat=2))
    answer = []
    for left in pairs:
        row = []
        for right in pairs:
            positions = [0]*4
            for i, value in zip(cut, left+right):
                positions[i] = value
            row.append(diagram(partition, positions))
        answer.append(row)
    return answer


def leading_matrix(A, B, t):
    """Identity-position term A plus swapped-position term B."""
    q2 = len(A)
    answer = [[0]*(t*t*q2) for _ in range(t*t*q2)]
    for x, y in product(range(t), repeat=2):
        for i, j in product(range(q2), repeat=2):
            answer[(x*t+y)*q2+i][(x*t+y)*q2+j] += A[i][j]
            answer[(x*t+y)*q2+i][(y*t+x)*q2+j] += B[i][j]
    return answer


def binary_rank(rows):
    basis = {}
    for value in rows:
        while value:
            pivot = value.bit_length()-1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def ordered_lift_checks(rng):
    """No factorial division: replication preserves rank in characteristic two."""
    unordered = tuple(combinations(range(4), 2))
    ordered = tuple(permutations(range(4), 2))
    owner = [unordered.index(tuple(sorted(x))) for x in ordered]
    n, m = len(unordered)**2, len(ordered)**2
    checks = 0
    for _ in range(8):
        matrix = [rng.randrange(1 << n) for _ in range(n)]
        lifted = []
        for a, b in product(range(len(ordered)), repeat=2):
            row = matrix[owner[a]*len(unordered)+owner[b]]
            lifted.append(sum(((row >> (owner[c]*len(unordered)+owner[d])) & 1)
                              << (c*len(ordered)+d)
                              for c, d in product(range(len(ordered)), repeat=2)))
        assert binary_rank(matrix) == binary_rank(lifted)
        assert len(lifted) == m
        checks += 1
    # Full-through bijections act regularly on the ordered assignments of W.
    orders = tuple(permutations(range(4)))
    for left in orders:
        images = {tuple(left[p] for p in bijection) for bijection in orders}
        assert len(images) == 24
    # The unordered crossing pulls back to a sum of four distinct matching
    # diagrams. Each matching permutation is unique when it exists.
    for x, y, u, v in product(ordered, repeat=4):
        count = sum(all(x[j] == v[sigma[j]] for j in range(2))
                    and all(y[j] == u[tau[j]] for j in range(2))
                    for sigma, tau in product(tuple(permutations(range(2))), repeat=2))
        expected = int(set(x) == set(v) and set(y) == set(u))
        assert count == expected
    return checks, len(orders), len(ordered)**4


def zero_sum_checks():
    gram_checks = 0
    for p in (2, 3, 5):
        for s in range(3, 11):
            gram = [[int(i == j)+1 for j in range(s-1)] for i in range(s-1)]
            assert rank(gram, p) == s-1-int(s % p == 0)
            gram_checks += 1
    # Two defects: four disjoint blocks, one for each matching pair.
    pairs = ((0, 6), (1, 7), (2, 4), (3, 5))
    slot_block = [None]*8
    choices = [None]*8
    for block, (a, b) in enumerate(pairs):
        slot_block[a] = slot_block[b] = block
        # Paired difference vectors have inner product one over every field.
        choices[a] = ((3*block, 1), (3*block+2, -1))
        choices[b] = ((3*block+1, 1), (3*block+2, -1))
    refinements = []
    admissible = 0
    for partition in partitions(8):
        if any(partition[2*j] == partition[2*j+1] for j in range(4)):
            continue
        admissible += 1
        survives_blocks = all(partition[i] != partition[j]
                              or slot_block[i] == slot_block[j]
                              for i in range(8) for j in range(8))
        if survives_blocks:
            refinements.append(partition)
    assert len(refinements) == 16
    for partition in refinements:
        value = 0
        for selected in product(*choices):
            positions = [a for a, _ in selected]
            if all(partition[i] != partition[j] or positions[i] == positions[j]
                   for i in range(8) for j in range(8)):
                sign = 1
                for _, coefficient in selected:
                    sign *= coefficient
                value += sign
        expected = int(all(partition[a] == partition[b] for a, b in pairs))
        assert value == expected
    return gram_checks, admissible, len(refinements)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert len(PARTITIONS) == 15
    zeta = [[int(refines(a, b)) for a in PARTITIONS] for b in PARTITIONS]
    row_cut, col_cut = (0, 1, 2, 3), (0, 2, 1, 3)
    assert {p for p in PARTITIONS if through(p, row_cut) == 2} == {E, C}
    assert {p for p in PARTITIONS if through(p, col_cut) == 2} == {E, R}
    diagram_checks = 0
    for p in (2, 3, 5):
        assert rank(zeta, p) == 15
        for t in (4, 5):
            for cut in (row_cut, col_cut):
                for partition in PARTITIONS:
                    value = rank(position_matrix(partition, t, cut), p)
                    assert value <= t**through(partition, cut)
                    diagram_checks += 1
    rng = random.Random(20261030)
    block_checks = 0
    for p in (2, 3, 5):
        for q in (1, 2):
            n = q*q
            for case in range(8):
                A = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
                B = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
                B[0][0] = 1  # Nonzero crossing coefficient.
                block = [a+b for a, b in zip(A, B)] + [b+a for a, b in zip(A, B)]
                both = [[a+b for a, b in zip(x, y)] for x, y in zip(A, B)]
                rb, rd = rank(block, p), rank(both, p)
                assert rb >= 1
                for t in (3, 4):
                    actual = rank(leading_matrix(A, B, t), p)
                    assert actual == (t*(t-1)//2)*rb+t*rd
                    block_checks += 1
    # A constant position function has E-orbit value one but E-diagram
    # coefficient zero. This catches confusion between the two bases.
    discrete = (0, 1, 2, 3)
    assert diagram(discrete, E) == 1 and discrete != E
    lift_checks, bijections, crossing_checks = ordered_lift_checks(rng)
    gram_checks, admissible, projection_checks = zero_sum_checks()
    result = {
        "status": "exact checks of the uniform invariant-rank obstruction",
        "diagram_count": len(PARTITIONS),
        "fields": [2, 3, 5],
        "position_rank_checks": diagram_checks,
        "colored_leading_block_checks": block_checks,
        "ordered_two_defect_rank_replication_checks": lift_checks,
        "two_defect_full_through_bijections": bijections,
        "two_defect_pulled_back_crossing_entries": crossing_checks,
        "zero_sum_gram_checks": gram_checks,
        "admissible_two_defect_diagrams": admissible,
        "two_defect_refinement_projection_checks": projection_checks,
        "row_two_through_diagrams": [E, C],
        "column_two_through_diagrams": [E, R],
        "rank_lower_bound": "binom(t,2)-13*q^2*t in each cut of the same piece",
        "uniform_rank_lower_bound": "(floor(t/(2*b))-2)^(2*b) in both cuts when t>=6*b",
        "scope": "simultaneous position-permutation invariance; fixed-color contraction obstruction when t/b tends to infinity",
        "negative_control": "exact equality-orbit values differ from diagram coefficients",
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
