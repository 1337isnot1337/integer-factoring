#!/usr/bin/env python3
"""F4 rectangular joint energy splits of the literal five-letter source."""
from functools import lru_cache
from itertools import product

from joint_energy_conjugations import model
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import F4_PRODUCT, check_split, linear_combination

HEIGHTS = (0, 1, 3, 4)


@lru_cache(None)
def source():
    data = model()
    return data['operator'], data['energies'], data['half']


def coordinate_start():
    _, _, heights = source()
    physical = [heights.index(h) for h in HEIGHTS]
    vectors = [[int(j == a+5*b) for j in range(25)] for a in physical for b in physical]
    return [row[:] for row in vectors], [row[:] for row in vectors]


def shear(embedding, retraction, a, b, coefficient):
    """P <- E_ab P; L <- L E_ab, with E_ab inverse equal to itself."""
    _, energies, _ = source()
    if a == b or energies[a] != energies[b] or coefficient not in (1, 2, 3):
        raise ValueError('Use distinct coordinates of equal joint energy and a nonzero F4 scalar.')
    p, l = [row[:] for row in embedding], [row[:] for row in retraction]
    for row in p:
        row[a] ^= F4_PRODUCT[coefficient][row[b]]
    for row in l:
        row[b] ^= F4_PRODUCT[coefficient][row[a]]
    return p, l


def crossing(embedding, retraction):
    if check_split(retraction, embedding) != 25 or len(embedding) != 16:
        raise ValueError('Use a sixteen-state split of the twenty-five-state joint space.')
    ambient, energies, _ = source()
    heights = [a+b for a in HEIGHTS for b in HEIGHTS]
    if any(energies[j] != h for vector, h in zip(embedding+retraction, heights+heights)
           for j, coefficient in enumerate(vector) if coefficient):
        raise ValueError('Both joint maps must preserve the new additive energy.')
    temporary = [[linear_combination((ambient[i][j], value) for j, value in enumerate(column) if value)
                  for column in embedding] for i in range(25)]
    operator = [[linear_combination((temporary[i][j], value) for i, value in enumerate(row) if value)
                 for j in range(16)] for row in retraction]
    for i, j in product(range(16), repeat=2):
        low, high = operator[i][j]
        assert (low.bit_count() % 2, high.bit_count() % 2) == (int(i == j), 0)
        assert all(not ((value << heights[i]) & 255) for value in (low, high))
        assert max(low.bit_length(), high.bit_length())-1 <= heights[j]+8
    return operator


def profiles(operator):
    opposite = [[None]*16 for _ in range(16)]
    for x, y, u, v in product(range(4), repeat=4):
        opposite[4*x+u][4*v+y] = operator[4*x+y][4*v+u]
    first = MatrixOracle(operator, 16).profile([a+b for a in HEIGHTS for b in HEIGHTS], [0]*16)['orders']
    second = MatrixOracle(opposite, 16).profile([b for _ in HEIGHTS for b in HEIGHTS],
                                              [a for a in HEIGHTS for _ in HEIGHTS], infinity=True)['orders']
    assert all(value == int(value) for value in first+second)
    return dict(zero_orders_t=sorted(int(value)-8 for value in first),
                infinity_orders_t=sorted(int(value)+8 for value in second),
                heights=list(HEIGHTS), completed=True, native_precision=512)
