"""Exact coherent physical projections of the six-vertex matching alphabet."""

from perfect_matching_crossing import perfect_matchings, dense_seed, binary_rank
from physical_split_crossings import retract, check_split


def binary_solution(equations, right_hand_sides):
    pivots = {}
    for row, rhs in zip(equations, right_hand_sides):
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row, rhs
                break
            previous, value = pivots[pivot]
            row ^= previous
            rhs ^= value
        if not row and rhs:
            raise ValueError("inconsistent binary system")
    solution = 0
    for pivot, (row, rhs) in sorted(pivots.items()):
        solution |= (rhs ^ ((row & solution).bit_count() % 2)) << pivot
    return solution


def canonical_matching_split(include_trivial=True):
    """Isotropic duals for im(A+I+J), optionally with the trivial line.

    A connects perfect matchings that share exactly one edge. The image
    has dimension four, is isotropic, and is orthogonal to the all-ones
    vector. Both returned Gram matrices are diag(0,0,0,0,1) when the
    trivial line is included, and zero when it is omitted.
    """
    labels = tuple(map(set, perfect_matchings(tuple(range(6)))))
    n = len(labels)
    nilpotent = [sum((int(len(a & b) == 1) ^ int(i == j) ^ 1) << j
                     for j, b in enumerate(labels)) for i, a in enumerate(labels)]
    assert binary_rank(nilpotent) == 4
    for row in nilpotent:
        square = 0
        for j in range(n):
            if row >> j & 1:
                square ^= nilpotent[j]
        assert square == 0
    columns = []
    for row in nilpotent:
        if binary_rank(columns + [row]) > len(columns):
            columns.append(row)
    one = (1 << n) - 1
    assert all((a & b).bit_count() % 2 == 0 for a in columns for b in columns)
    assert all((a & one).bit_count() % 2 == 0 for a in columns)
    duals = [binary_solution(columns + [one], [int(i == j) for j in range(4)] + [0])
             for i in range(4)]
    gram = [[(a & b).bit_count() % 2 for b in duals] for a in duals]
    assert all(gram[i][i] == 0 for i in range(4))
    repaired = duals[:]
    for i in range(4):
        for j in range(i + 1, 4):
            if gram[i][j]:
                repaired[i] ^= columns[j]
    assert all((a & b).bit_count() % 2 == 0 for a in repaired for b in repaired)
    if include_trivial:
        columns.append(one)
        repaired.append(one)
    left = [[(row >> j) & 1 for j in range(n)] for row in repaired]
    right = [[(row >> j) & 1 for j in range(n)] for row in columns]
    check_split(left, right)
    return left, right


def canonical_projection(include_trivial=True):
    left, right = canonical_matching_split(include_trivial)
    seed = retract(dense_seed(3), left, right)
    seed.update({'constructor': 'isotropic_matching_projection',
                 'include_trivial_line': include_trivial,
                 'left': left, 'right': right})
    return seed
