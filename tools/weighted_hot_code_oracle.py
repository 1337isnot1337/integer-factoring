#!/usr/bin/env python3
"""Complete weighted profiles of selected hot-code compound blocks.

Full-letter diagonal heights are supported. Equality with a natural
reflection cut requires the separate Gaussian-factor integrality cone.
Every selected block retains its complete polynomial coefficients and
supplies its endpoint-proved rank to the native matrix evaluator.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix

from hot_code_charge_profiles import compound_block, hot_energy, one_particle, validate_code
from native_matrix_smith_f4 import MatrixOracle


def additions(letters):
    lookup = {value: index for index, value in enumerate(letters)}
    result = []
    for target_index, target in enumerate(letters):
        for bit in range(target.bit_length()):
            if target >> bit & 1:
                source = target ^ (1 << bit)
                if source in lookup:
                    result.append((lookup[source], target_index, bit))
    return result


class HotWeightedOracle:
    def __init__(self, xletters, yletters=None):
        xletters = tuple(xletters)
        yletters = xletters if yletters is None else tuple(yletters)
        self.modes = max(1, max((value.bit_length() for value in xletters + yletters), default=1))
        self.xletters = validate_code(xletters, self.modes)
        self.yletters = validate_code(yletters, self.modes)
        self.qx, self.qy = len(xletters), len(yletters)
        self.dimensions = (self.qx, self.qy, self.qy, self.qx)
        self.offsets = np.cumsum((0,) + self.dimensions)
        self.height_count = int(self.offsets[-1])
        self.energies_x = np.array(list(map(hot_energy, xletters)), dtype=float)
        self.energies_y = np.array(list(map(hot_energy, yletters)), dtype=float)
        self.blocks = []
        matrix = one_particle(self.modes)
        pairs = [(x, y) for x in range(self.qx) for y in range(self.qy)]
        charges = sorted({xletters[x].bit_count() + yletters[y].bit_count() for x, y in pairs})
        for charge in charges:
            selected = [(x, y) for x, y in pairs
                        if xletters[x].bit_count() + yletters[y].bit_count() == charge]
            if len(selected) > 1024:
                raise ValueError('The current weighted matrix interface supports charge blocks of order at most 1024.')
            rows = [xletters[x] | (yletters[y] << self.modes) for x, y in selected]
            columns = [yletters[y] | (xletters[x] << self.modes) for x, y in selected]
            data = compound_block(matrix, rows, columns, coefficients=True)
            coefficients = [[(value, 0) for value in row] for row in data['matrix']]
            precision = data['precision_t']
            self.blocks.append({'charge': charge, 'pairs': np.asarray(selected, dtype=int),
                                'coefficients': coefficients, 'precision': precision,
                                'oracle': MatrixOracle(coefficients, len(selected), precision=precision)})
        self.active_heights = np.array([self.offsets[leg] + symbol for leg, size in enumerate(self.dimensions)
                                       for symbol in range(1, size)], dtype=int)

    def log_profile_with_gradient(self, heights, scale):
        """Return log sum exp(-nu), height derivative, scale derivative, orders.

        Heights are already multiplied by the logarithmic local scale.
        ``scale`` multiplies valuations in the original variable w.
        """
        heights = np.asarray(heights, dtype=float)
        if heights.shape != (self.height_count,) or not np.all(np.isfinite(heights)):
            raise ValueError('Supply one finite real height per physical letter on each leg.')
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('The logarithmic local scale must be positive.')
        orders, labels, base_values = [], [], []
        for block in self.blocks:
            pairs = block['pairs']
            rows = heights[self.offsets[0] + pairs[:, 0]] + heights[self.offsets[1] + pairs[:, 1]]
            columns = heights[self.offsets[2] + pairs[:, 1]] + heights[self.offsets[3] + pairs[:, 0]]
            while True:
                try:
                    result = block['oracle'].profile(rows, columns, 2 * float(scale))
                    break
                except ArithmeticError:
                    if block['precision'] >= 8192:
                        raise
                    block['precision'] *= 2
                    block['oracle'] = MatrixOracle(block['coefficients'], len(pairs), precision=block['precision'])
            orders.extend(result['orders'])
            base_values.extend(2 * np.diff(np.r_[0.0, result['base_prefix']]))
            chosen_rows, chosen_columns = pairs[np.asarray(result['rows'])], pairs[np.asarray(result['cols'])]
            labels.extend(np.column_stack((chosen_rows[:, 0], chosen_rows[:, 1],
                                           chosen_columns[:, 1], chosen_columns[:, 0])))
        orders, labels, base_values = np.asarray(orders), np.asarray(labels), np.asarray(base_values)
        if len(orders) != self.qx * self.qy:
            raise ArithmeticError('The weighted blocks did not recover every promised channel.')
        minimum = float(np.min(orders))
        weights = np.exp(-(orders - minimum))
        probabilities = weights / np.sum(weights)
        value = -minimum + float(np.log(np.sum(weights)))
        gradient = np.zeros(self.height_count)
        for leg, size in enumerate(self.dimensions):
            gradient[self.offsets[leg]:self.offsets[leg + 1]] = -np.bincount(
                labels[:, leg], weights=probabilities, minlength=size)
        return value, gradient, -float(probabilities @ base_values), np.sort(orders)

    def cone(self, scale):
        """Return A,b for a stabilized ambient realization, with A h >= b.

        Every finite code admits a sufficiently large such interval.
        At a smaller non-stabilized interval, same-leg reflection edges
        require additional constraints not supplied by this method.
        """
        edges_x, edges_y = additions(self.xletters), additions(self.yletters)
        rows, columns, values, bounds = [], [], [], []
        for first_leg, second_leg, first_edges, second_edges in (
                (0, 1, edges_x, edges_y), (2, 3, edges_y, edges_x)):
            for source_a, target_a, bit_a in first_edges:
                for source_b, target_b, bit_b in second_edges:
                    row = len(bounds)
                    for column, value in ((self.offsets[first_leg] + target_a, 1),
                                          (self.offsets[first_leg] + source_a, -1),
                                          (self.offsets[second_leg] + target_b, 1),
                                          (self.offsets[second_leg] + source_b, -1)):
                        rows.append(row); columns.append(int(column)); values.append(value)
                    bounds.append(-4 * float(scale) * (bit_a & bit_b))
        return csr_matrix((values, (rows, columns)), shape=(len(bounds), self.height_count)), np.asarray(bounds)

    def type_mean(self, heights, type_x, type_y):
        marginals = np.concatenate((type_x, type_y, type_y, type_x))
        return float(marginals @ np.asarray(heights)), marginals
