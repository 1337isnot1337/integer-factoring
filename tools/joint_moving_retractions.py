#!/usr/bin/env python3
"""Joint one-sided Laurent retractions of the literal five-letter operator."""
from itertools import product

from joint_five_letter_splits import HEIGHTS, coordinate_start, source
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import scalar_multiply


def coordinate_indices():
    return tuple(row.index(1) for row in coordinate_start()[0])


def moving_crossing(inside, outside, power, scalar=1, side='column'):
    """P=P0+z^power*c*e_out*e_in^T, L=L0, or its row counterpart.

    The added ambient direction is outside the coordinate image, so LP=I
    identically. This is a joint map, without separate physical leg factors.
    Returns a polynomial numerator z^clearing*C, with its clearing charged.
    """
    indices = coordinate_indices()
    if inside not in range(16) or outside not in range(25) or outside in indices:
        raise ValueError('Use a target coordinate and a discarded source coordinate.')
    if not isinstance(power, int) or scalar not in (1, 2, 3) or side not in ('column', 'row'):
        raise ValueError('Use an integer Laurent power, nonzero F4 scalar and row or column side.')
    ambient, _, _ = source()
    clearing = max(0, -power)
    matrix = [[tuple(value << clearing for value in ambient[i][j]) for j in indices] for i in indices]
    for j in range(16):
        row, col = (j, inside) if side == 'column' else (inside, j)
        value = ambient[indices[j]][outside] if side == 'column' else ambient[outside][indices[j]]
        update = scalar_multiply(value, scalar)
        matrix[row][col] = tuple(a ^ (b << (power+clearing)) for a, b in zip(matrix[row][col], update))
    for i, j in product(range(16), repeat=2):
        low, high = matrix[i][j]
        assert (low.bit_count() % 2, high.bit_count() % 2) == (int(i == j), 0)
    return dict(matrix=matrix, clearing=clearing, inside=inside, outside=outside,
                power=power, scalar=scalar, side=side)


def profiles(code):
    matrix, clearing = code['matrix'], code['clearing']
    opposite = [[None]*16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        opposite[4*x+u][4*v+y] = matrix[4*x+y][4*v+u]
    first = MatrixOracle(matrix, 16).profile([x+y for x in HEIGHTS for y in HEIGHTS], [0]*16)['orders']
    second = MatrixOracle(opposite, 16).profile([y for _ in HEIGHTS for y in HEIGHTS],
                                              [x for x in HEIGHTS for _ in HEIGHTS], infinity=True)['orders']
    assert all(value == int(value) for value in first+second)
    return dict(zero_orders=[int(v)-8-clearing for v in first],
                infinity_orders=[int(v)+8+clearing for v in second])
