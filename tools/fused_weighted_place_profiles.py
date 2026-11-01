#!/usr/bin/env python3
"""Complete weighted local profiles for literal fused polynomial crossings.

Every determinantal prefix is minimized over all nonzero block minors.
The exact API accepts integers or Fractions. A separate NumPy API supports
exploratory continuous heights and returns no certification claim.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from fused_six_vertex_crossings import construct
from check_fused_multipoint_profiles import determinant
from check_fused_closed_profiles import factor_valuation, from_canonical


def labelled_blocks(matrix, alphabet, cut):
    size = alphabet * alphabet
    entries = [[{} for _ in range(size)] for _ in range(size)]
    labels = [[None for _ in range(size)] for _ in range(size)]
    for key in product(range(alphabet), repeat=4):
        x, y, u, v = key
        r, c = ((alphabet * x + y, alphabet * u + v) if cut == 0 else
                (alphabet * x + u, alphabet * y + v))
        entries[r][c] = matrix[alphabet * y + x][alphabet * u + v]
        labels[r][c] = key
    unseen = set(range(size))
    result = []
    while unseen:
        rows, columns = {min(unseen)}, set()
        while True:
            new_columns = {j for i in rows for j in range(size) if entries[i][j]}
            new_rows = {i for j in new_columns for i in range(size) if entries[i][j]}
            if new_rows == rows and new_columns == columns:
                break
            rows, columns = new_rows, new_columns
        assert rows and rows <= unseen and len(rows) == len(columns)
        unseen -= rows
        rows, columns = sorted(rows), sorted(columns)
        result.append({'matrix': [[entries[i][j] for j in columns] for i in rows],
                       'labels': [[labels[i][j] for j in columns] for i in rows]})
    return result


def physical_counts(labels, rows, columns, alphabet, cut):
    counts = [0] * (4 * alphabet)
    row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
    for i in rows:
        key = labels[i][0]
        for leg in row_legs:
            counts[leg * alphabet + key[leg]] += 1
    for j in columns:
        key = labels[0][j]
        for leg in column_legs:
            counts[leg * alphabet + key[leg]] += 1
    assert all(sum(counts[leg * alphabet:(leg + 1) * alphabet]) == len(rows)
               for leg in range(4))
    return tuple(counts)


class FusedPlaceOracle:
    def __init__(self, m, theta, catalogue=None):
        self.source = construct(m, theta)
        assert self.source['invariant'] and self.source['endpoint_identity']
        self.alphabet, self.degree = m + 1, self.source['degree']
        if catalogue is None:
            catalogue = json.loads((Path(__file__).parent.parent / 'evidence' /
                                    'fused-closed-profiles.json').read_text())
        self.record = next(case for case in catalogue['cases']
                           if (case['m'], case['theta']) == (m, theta))
        self.points = [tuple(point['factor']) for point in self.record['closed_points']] + [None]
        self.residue_degrees = {tuple(point['factor']): point['residue_degree']
                                for point in self.record['closed_points']}
        self.residue_degrees[None] = 1
        self.blocks = [labelled_blocks(self.source['matrix'], self.alphabet, cut) for cut in (0, 1)]
        self.minors, self.groups, self.counts, self.valuations = [], [], [], []
        for cut in (0, 1):
            minors, groups = [], []
            for block in self.blocks[cut]:
                size = len(block['matrix'])
                block_groups = []
                for k in range(1, size + 1):
                    start = len(minors)
                    for rows in combinations(range(size), k):
                        for columns in combinations(range(size), k):
                            value = determinant([[block['matrix'][i][j] for j in columns] for i in rows])
                            if value:
                                counts = physical_counts(block['labels'], rows, columns, self.alphabet, cut)
                                minors.append({'polynomial': value, 'counts': counts, 'size': k})
                    assert len(minors) > start
                    block_groups.append((start, len(minors)))
                groups.append(block_groups)
            valuations = {}
            for point in self.points:
                factor = None if point is None else from_canonical(point)
                valuations[point] = tuple(
                    self.degree * minor['size'] - max(minor['polynomial']) if point is None else
                    factor_valuation(minor['polynomial'], factor) for minor in minors)
            self.minors.append(minors)
            self.groups.append(groups)
            self.counts.append(np.asarray([minor['counts'] for minor in minors], dtype=float))
            self.valuations.append(valuations)

    def profile(self, cut, point, heights, scale=1):
        """Exact restored profile. The output includes all channels."""
        if len(heights) != 4 * self.alphabet:
            raise ValueError('Supply one height for every letter at all four legs.')
        if any(not isinstance(value, (int, Fraction)) for value in heights):
            raise TypeError('The exact API accepts only integers or Fractions.')
        if not isinstance(scale, (int, Fraction)) or scale < 0:
            raise ValueError('The valuation scale must be a nonnegative integer or Fraction.')
        values = [scale * value + sum(count * height for count, height in zip(minor['counts'], heights))
                  for minor, value in zip(self.minors[cut], self.valuations[cut][point])]
        orders = []
        for block in self.groups[cut]:
            previous = 0
            block_orders = []
            for start, stop in block:
                current = min(values[start:stop])
                block_orders.append(current - previous)
                previous = current
            assert block_orders == sorted(block_orders)
            orders.extend(block_orders)
        return tuple(sorted(orders))

    def profile_float(self, cut, point, heights, scale=1.0):
        """Exploratory float profile; exact replay is needed for a certificate."""
        heights = np.asarray(heights, dtype=float)
        if heights.shape != (4 * self.alphabet,) or scale < 0:
            raise ValueError('Invalid physical height vector or valuation scale.')
        values = self.counts[cut] @ heights + scale * np.asarray(self.valuations[cut][point], dtype=float)
        orders = []
        for block in self.groups[cut]:
            previous = 0.0
            for start, stop in block:
                current = float(np.min(values[start:stop]))
                orders.append(current - previous)
                previous = current
        return np.sort(np.asarray(orders))

    def log_profile_with_gradient(self, cut, point, heights, scale=1.0):
        """Exploratory log partition sum and an active-minor subgradient."""
        heights = np.asarray(heights, dtype=float)
        if heights.shape != (4 * self.alphabet,) or scale < 0:
            raise ValueError('Invalid physical height vector or valuation scale.')
        base_values = np.asarray(self.valuations[cut][point], dtype=float)
        values = self.counts[cut] @ heights + scale * base_values
        orders, height_derivatives, scale_derivatives = [], [], []
        for block in self.groups[cut]:
            previous = 0.0
            previous_counts = np.zeros(4 * self.alphabet)
            previous_base = 0.0
            for start, stop in block:
                pivot = start + int(np.argmin(values[start:stop]))
                current = float(values[pivot])
                counts = self.counts[cut][pivot]
                base = base_values[pivot]
                orders.append(current - previous)
                height_derivatives.append(counts - previous_counts)
                scale_derivatives.append(base - previous_base)
                previous, previous_counts, previous_base = current, counts, base
        orders = np.asarray(orders)
        minimum = float(np.min(orders))
        exponentials = np.exp(-(orders - minimum))
        probabilities = exponentials / np.sum(exponentials)
        value = -minimum + float(np.log(np.sum(exponentials)))
        gradient_heights = -probabilities @ np.asarray(height_derivatives)
        gradient_scale = -float(probabilities @ np.asarray(scale_derivatives))
        return value, gradient_heights, gradient_scale, np.sort(orders)

    def linear_heights(self, slopes):
        if len(slopes) != 4:
            raise ValueError('Four physical-leg slopes are required.')
        return tuple(slope * a for slope in slopes for a in range(self.alphabet))

    def type_shift(self, heights, type_x, type_y):
        if len(type_x) != self.alphabet or len(type_y) != self.alphabet:
            raise ValueError('Types must use the actual physical alphabet.')
        q = self.alphabet
        return (sum(type_x[a] * (heights[a] + heights[3 * q + a]) for a in range(q))
                + sum(type_y[a] * (heights[q + a] + heights[2 * q + a]) for a in range(q)))
