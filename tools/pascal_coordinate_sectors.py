#!/usr/bin/env python3
"""Selected original Fock entries and complete occupation-sector profiles.

No deletion closure or Gaussian profile substitution is used. Each original
minor and every nonempty sector is retained over the exact field F2.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from itertools import product

from native_matrix_smith_f4 import MatrixOracle
from pascal_letter_gauges import PascalLetterPressure
from tensor_pascal_box_crossings import polynomial_product


def coordinate_exterior_tensor(seed, alphabet):
    alphabet = list(alphabet)
    n = seed['n']
    if (len(alphabet) < 2 or len(set(alphabet)) != len(alphabet)
            or any(type(value) is not int or not 0 <= value < (1 << n) for value in alphabet)):
        raise ValueError('Supply at least two distinct ambient occupation words.')
    if n > 16 or len(alphabet) > 64:
        raise ValueError('The finite selected-entry implementation limits modes to 16 and letters to 64.')
    matrix = seed['matrix']

    @lru_cache(1 << 20)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        candidates = [i for i in range(2 * n) if rows >> i & 1]
        i = min(candidates, key=lambda k: sum(bool(matrix[k][j])
                for j in range(2 * n) if columns >> j & 1))
        pending, value = columns, 0
        while pending:
            other = pending & -pending
            j = other.bit_length() - 1
            if matrix[i][j]:
                value ^= polynomial_product(matrix[i][j], minor(rows ^ (1 << i), columns ^ other))
            pending ^= other
        return value

    values = {}
    occupations = {word: word.bit_count() for word in alphabet}
    for b, x, v, y in product(alphabet, repeat=4):
        if occupations[b] + occupations[x] != occupations[v] + occupations[y]:
            continue
        value = minor(v | (y << n), b | (x << n))
        if value:
            values[b, x, v, y] = value
    cache = minor.cache_info()
    minor.cache_clear()
    return {'alphabet': alphabet, 'values': values,
            'minor_cache_entries': cache.currsize,
            'minor_cache_misses': cache.misses}


class SectorOracle:
    def __init__(self, tensor, place, precision=512):
        alphabet = tensor['alphabet']
        self.q = len(alphabet)
        self.place = place
        index = {word: i for i, word in enumerate(alphabet)}
        self.legs = (0, 2, 1, 3) if place == 0 else (0, 1, 2, 3)
        rows, columns = defaultdict(list), defaultdict(list)
        for b, second in product(range(self.q), repeat=2):
            first_occupation = alphabet[b].bit_count()
            second_occupation = alphabet[second].bit_count()
            charge = first_occupation - second_occupation if place == 0 else first_occupation + second_occupation
            rows[charge].append(b + self.q * second)
            column_charge = -charge if place == 0 else charge
            columns[column_charge].append(b + self.q * second)
        blocks = {}
        locations = {}
        for charge in sorted(rows):
            row_indices, column_indices = rows[charge], columns[charge]
            assert len(row_indices) == len(column_indices)
            blocks[charge] = [[(0, 0) for _ in column_indices] for _ in row_indices]
            locations[charge] = ({value: i for i, value in enumerate(row_indices)},
                                 {value: i for i, value in enumerate(column_indices)})
        for labels, value in tensor['values'].items():
            b, x, v, y = labels
            charge = b.bit_count() - v.bit_count() if place == 0 else b.bit_count() + x.bit_count()
            row = index[labels[self.legs[0]]] + self.q * index[labels[self.legs[1]]]
            column = index[labels[self.legs[2]]] + self.q * index[labels[self.legs[3]]]
            ri, ci = locations[charge]
            blocks[charge][ri[row]][ci[column]] = (value, 0)
        self.blocks = [(charge, rows[charge], columns[charge],
                        MatrixOracle(blocks[charge], len(rows[charge]), precision=precision))
                       for charge in sorted(rows)]
        assert sum(len(r) for _, r, _, _ in self.blocks) == self.q ** 2

    def profile(self, row_heights, column_heights, scale, infinity=False):
        if bool(infinity) != bool(self.place):
            raise ValueError('Use the valuation place of the compiled occupation sectors.')
        records, merged = [], []
        for charge, rows, columns, oracle in self.blocks:
            result = oracle.profile([row_heights[i] for i in rows],
                                    [column_heights[j] for j in columns], scale,
                                    infinity=infinity)
            prefixes = [0] + result['base_prefix']
            for i, (order, r, c) in enumerate(zip(result['orders'], result['rows'], result['cols'])):
                merged.append((order, rows[r], columns[c], prefixes[i + 1] - prefixes[i]))
            records.append({'charge': charge, 'dimension': len(rows), 'orders': result['orders']})
        merged.sort()
        base_prefix, total = [], 0
        for _, _, _, base in merged:
            total += base
            base_prefix.append(total)
        return {'orders': [record[0] for record in merged],
                'rows': [record[1] for record in merged],
                'cols': [record[2] for record in merged], 'base_prefix': base_prefix,
                'occupation_blocks': records}


class SectorPascalLetterPressure(PascalLetterPressure):
    def __init__(self, seed, alphabet, precision=512):
        self.seed, self.n = seed, seed['n']
        tensor = coordinate_exterior_tensor(seed, alphabet)
        self.alphabet, self.q = tensor['alphabet'], len(tensor['alphabet'])
        self.size, self.queries, self.last = 4 * self.q, 0, None
        self.construction = {'nonzero_original_entries': len(tensor['values']),
                             'minor_cache_entries': tensor['minor_cache_entries']}
        self.cuts = [(SectorOracle(tensor, place, precision), legs)
                     for place, legs in enumerate(((0, 2, 1, 3), (0, 1, 2, 3)))]
