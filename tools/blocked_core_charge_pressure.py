#!/usr/bin/env python3
"""Exact charge-block acceleration for physical-core pressure searches.

The block split uses two already conserved charges: the integer Gaussian
charge, and the projected physical F3 charge. It changes no tensor,
profile, or allowed physical gauge. Every advertised off-block entry is
checked to be identically zero before any smaller oracle is constructed.
"""
from __future__ import annotations

from collections import defaultdict

from native_matrix_smith_f4 import MatrixOracle
from strict_core_charge_gauges import CoreChargePressure, gaussian_coordinates


class BlockMatrixOracle:
    """A full-rank square matrix, split into certified exact square blocks."""

    def __init__(self, matrix, shifts, row_keys, column_keys, precision=512):
        self.n = len(matrix)
        self.rank = self.n
        self.precision = precision
        if len(row_keys) != self.n or len(column_keys) != self.n:
            raise ValueError('Provide one charge key per row and column.')
        if any(len(row) != self.n for row in matrix):
            raise ValueError('The full matrix must be square.')
        row_groups, column_groups = defaultdict(list), defaultdict(list)
        for i, key in enumerate(row_keys):
            row_groups[key].append(i)
        for j, key in enumerate(column_keys):
            column_groups[key].append(j)
        if set(row_groups) != set(column_groups):
            raise ArithmeticError('The conserved charge sectors disagree.')
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                if value != (0, 0) and row_keys[i] != column_keys[j]:
                    raise ArithmeticError('An alleged off-block entry is nonzero.')
        self.blocks = []
        for key in sorted(row_groups):
            rows, columns = row_groups[key], column_groups[key]
            if len(rows) != len(columns):
                raise ArithmeticError('A promised full-rank charge block is not square.')
            small = [[matrix[i][j] for j in columns] for i in rows]
            offsets = [[shifts[i][j] for j in columns] for i in rows]
            oracle = MatrixOracle(small, len(rows), offsets, precision)
            self.blocks.append((key, rows, columns, oracle))

    @property
    def block_sizes(self):
        return [len(rows) for _, rows, _, _ in self.blocks]

    def profile(self, row_weights=None, column_weights=None, scale=1.0, infinity=False):
        row_weights = [0.0] * self.n if row_weights is None else list(row_weights)
        column_weights = [0.0] * self.n if column_weights is None else list(column_weights)
        if len(row_weights) != self.n or len(column_weights) != self.n:
            raise ValueError('Provide one weight per full row and column.')
        entries = []
        for _, rows, columns, oracle in self.blocks:
            result = oracle.profile([row_weights[i] for i in rows],
                                    [column_weights[j] for j in columns], scale, infinity)
            previous = 0
            for order, prefix, i, j in zip(result['orders'], result['base_prefix'],
                                            result['rows'], result['cols']):
                entries.append((order, rows[i], columns[j], prefix - previous))
                previous = prefix
        entries.sort()
        cumulative, prefixes = 0, []
        for _, _, _, base in entries:
            cumulative += base
            prefixes.append(cumulative)
        return {'orders': [entry[0] for entry in entries], 'base_prefix': prefixes,
                'rows': [entry[1] for entry in entries], 'cols': [entry[2] for entry in entries]}


class BlockedCoreChargePressure(CoreChargePressure):
    """The existing objective with independently checked charge block oracles."""

    def __init__(self, interval, charges, core=2, family='strict_low'):
        super().__init__(interval, charges, core, family)
        for side, record in enumerate(self.records):
            data, metadata, left, right = record
            gaussian, shifts, _, _, _ = gaussian_coordinates(data['place'], interval, core, family)

            def signed(basis):
                return [sum((1 if letter in ('x', 'u') else -1) * self.charges[labels[i]]
                            for letter, labels in basis['labels'].items()) % 3
                        for i in range(len(gaussian))]

            first, second = signed(left), signed(right)
            row_keys = list(zip(left['row_charges'], first))
            column_keys = list(zip(right['column_charges'], [(-value) % 3 for value in second]))
            matrix = [[value if (first[i] + second[j]) % 3 == 0 else (0, 0)
                       for j, value in enumerate(row)] for i, row in enumerate(gaussian)]
            self.oracles[side] = BlockMatrixOracle(matrix, shifts, row_keys, column_keys,
                                                   metadata['profile_precision'])
