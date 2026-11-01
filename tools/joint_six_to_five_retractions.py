#!/usr/bin/env python3
"""A six-to-five rational compression at a regraded positive-deficit shape."""
from functools import lru_cache
from itertools import product

from joint_hot_retraction_module import deflate
from native_matrix_smith_f4 import MatrixOracle
from polynomial_retraction_spaces import bounded_laurent_space, row_reduced_basis
from restrict_implicit_reflection import restricted_tensor

HEIGHTS = (0, 0, 6, 8, 9)
H = tuple(x+y for x in HEIGHTS for y in HEIGHTS)
MU_ZERO, MU_INFINITY = 10, -16


@lru_cache(None)
def source():
    letters = (0, 8, 4, 12, 2, 10)
    seed = restricted_tensor(8, letters)
    common = min((p & -p).bit_length()-1 for p, _ in seed['values'] if p)
    assert common == 8
    values = [deflate(p >> common) for p, high in seed['values']]
    matrix = [[values[r+36*(c//6)+216*(c%6)] for c in range(36)] for r in range(36)]
    selected = tuple(x+6*y for x in range(5) for y in range(5))
    excluded = tuple(i for i in range(36) if i not in selected)
    base = [[matrix[i][j] for j in selected] for i in selected]
    rows = [[matrix[i][j] for j in selected] for i in excluded]
    return dict(matrix=matrix, retained=base, source_rows=rows, selected=selected,
                excluded=excluded, letters=letters, common_w_shift=common)


@lru_cache(None)
def module():
    return row_reduced_basis(source()['source_rows'])


@lru_cache(None)
def row_space(row_height, protected=0):
    return bounded_laurent_space(module()[0], MU_ZERO+protected-row_height,
                                 [h-MU_INFINITY-protected for h in H])


def crossing(assignments=()):
    """Apply (target row, bounded-space basis index, F4 scalar) corrections."""
    base = source()['retained']
    clearing = [row_space(h)['clearing'] for h in H]
    matrix = [[(p << c, 0) for p in row] for row, c in zip(base, clearing)]
    for row, index, scalar in assignments:
        if row not in range(25) or scalar not in (1, 2, 3):
            raise ValueError('Use a target row and nonzero F4 scalar.')
        space = row_space(H[row])
        if index not in range(space['dimension']):
            raise ValueError('The basis index is outside the declared bounded row space.')
        for j, p in enumerate(space['rows'][index]):
            a, b = ((0, 0), (p, 0), (0, p), (p, p))[scalar]
            low, high = matrix[row][j]
            matrix[row][j] = low ^ a, high ^ b
    for i, j in product(range(25), repeat=2):
        low, high = matrix[i][j]
        assert (low.bit_count() % 2, high.bit_count() % 2) == (int(i == j), 0)
    return dict(matrix=matrix, row_clearing=clearing)


def profiles(code):
    matrix, clear = code['matrix'], code['row_clearing']
    shifts = [[-c]*25 for c in clear]
    opposite = [[None]*25 for _ in range(25)]
    opposite_shifts = [[0]*25 for _ in range(25)]
    for x, y, u, v in product(range(5), repeat=4):
        opposite[5*x+u][5*v+y] = matrix[5*x+y][5*v+u]
        opposite_shifts[5*x+u][5*v+y] = -clear[5*x+y]
    zero = MatrixOracle(matrix, 25, shifts=shifts).profile(H)['orders']
    infinity = MatrixOracle(opposite, 25, shifts=opposite_shifts).profile(
        [y for _ in HEIGHTS for y in HEIGHTS], [x for x in HEIGHTS for _ in HEIGHTS], infinity=True)['orders']
    assert all(v == int(v) for v in zero+infinity)
    return dict(raw_zero_orders=[int(v) for v in zero], raw_infinity_orders=[int(v) for v in infinity],
                normalized_zero_orders=[int(v)-MU_ZERO for v in zero],
                normalized_infinity_orders=[int(v)-MU_INFINITY for v in infinity],
                declared_deficit=-MU_ZERO-MU_INFINITY,
                scope='Raw orders include every Laurent shift; the normalized pressure needs the declared deficit term.')
