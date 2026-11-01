#!/usr/bin/env python3
"""Selected F4-valued original Fock minors and complete occupation sectors."""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import product

from cube_f4_bivariate_smith import ZERO, add, multiply
from native_matrix_smith_f4 import MatrixOracle
from pascal_coordinate_sectors import SectorOracle
from pascal_letter_gauges import PascalLetterPressure


def coordinate_exterior_tensor(seed, alphabet):
    alphabet = list(alphabet)
    n = seed['n']
    if (n > 16 or not 2 <= len(alphabet) <= 64 or len(set(alphabet)) != len(alphabet)
            or any(type(v) is not int or not 0 <= v < 1 << n for v in alphabet)):
        raise ValueError('Use 2--64 distinct coordinate words on at most sixteen modes.')
    matrix = [[tuple(v) if isinstance(v, (tuple, list)) else (v, 0) for v in row]
              for row in seed['matrix']]

    @lru_cache(1 << 20)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return ZERO
        if not rows:
            return (1, 0)
        i = min((a for a in range(2*n) if rows >> a & 1),
                key=lambda a: sum(matrix[a][b] != ZERO for b in range(2*n) if columns >> b & 1))
        remaining, value = columns, ZERO
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            j = bit.bit_length()-1
            if matrix[i][j] != ZERO:
                value = add(value, multiply(matrix[i][j], minor(rows ^ (1 << i), columns ^ bit)))
        return value

    values = {}
    occupations = {word: word.bit_count() for word in alphabet}
    for b, x, v, y in product(alphabet, repeat=4):
        if occupations[b]+occupations[x] != occupations[v]+occupations[y]:
            continue
        value = minor(v | (y << n), b | (x << n))
        if value != ZERO:
            values[b, x, v, y] = value
    cache = minor.cache_info()
    minor.cache_clear()
    return {'alphabet': alphabet, 'values': values, 'minor_cache_entries': cache.currsize,
            'minor_cache_misses': cache.misses}


class F4SectorOracle(SectorOracle):
    def __init__(self, tensor, place, precision=512):
        alphabet = tensor['alphabet']
        self.q, self.place = len(alphabet), place
        self.legs = (0, 2, 1, 3) if place == 0 else (0, 1, 2, 3)
        index = {word: i for i, word in enumerate(alphabet)}
        rows, columns = defaultdict(list), defaultdict(list)
        for a, b in product(range(self.q), repeat=2):
            charge = (alphabet[a].bit_count()-alphabet[b].bit_count() if place == 0 else
                      alphabet[a].bit_count()+alphabet[b].bit_count())
            rows[charge].append(a+self.q*b)
            columns[-charge if place == 0 else charge].append(a+self.q*b)
        blocks, locations = {}, {}
        for charge in sorted(rows):
            assert len(rows[charge]) == len(columns[charge])
            blocks[charge] = [[ZERO]*len(columns[charge]) for _ in rows[charge]]
            locations[charge] = ({v: i for i, v in enumerate(rows[charge])},
                                 {v: i for i, v in enumerate(columns[charge])})
        for labels, value in tensor['values'].items():
            b, x, v, y = labels
            charge = b.bit_count()-v.bit_count() if place == 0 else b.bit_count()+x.bit_count()
            row = index[labels[self.legs[0]]]+self.q*index[labels[self.legs[1]]]
            column = index[labels[self.legs[2]]]+self.q*index[labels[self.legs[3]]]
            ri, ci = locations[charge]
            blocks[charge][ri[row]][ci[column]] = value
        self.blocks = [(charge, rows[charge], columns[charge],
                        MatrixOracle(blocks[charge], len(rows[charge]), precision=precision))
                       for charge in sorted(rows)]
        assert sum(len(r) for _, r, _, _ in self.blocks) == self.q**2


class F4SectorPascalLetterPressure(PascalLetterPressure):
    def __init__(self, seed, alphabet, precision=512):
        self.seed, self.n = seed, seed['n']
        tensor = coordinate_exterior_tensor(seed, alphabet)
        self.alphabet, self.q = tensor['alphabet'], len(tensor['alphabet'])
        self.size, self.queries, self.last = 4*self.q, 0, None
        self.construction = {'nonzero_original_entries': len(tensor['values']),
                             'nonbinary_original_entries': sum(bool(v[1]) for v in tensor['values'].values()),
                             'minor_cache_entries': tensor['minor_cache_entries']}
        self.cuts = [(F4SectorOracle(tensor, place, precision), legs)
                     for place, legs in enumerate(((0, 2, 1, 3), (0, 1, 2, 3)))]
