#!/usr/bin/env python3
"""Finite exact-profile screens of binary occupation-transfer kernels.

The finite-field profiles are exact. The charge/scale optimization is
floating point and supplies search evidence, not a pressure certificate.
Every pressure includes the unchanged spectators and paired type costs.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from coupled_identity_operator_filters import (
    _geometry, identity_kernel, infinity_profiles, normalized_kernel,
    states, zero_profiles,
)
from prime_state_filter_pressure import gaussian_pressure


def encode_kernel(kernel):
    return [{'row': row, 'column': column, 'polynomial': value}
            for (row, column), value in sorted(kernel.items())]


def decode_kernel(entries):
    return {(tuple(entry['row']), tuple(entry['column'])):
            {int(degree): value for degree, value in entry['polynomial'].items()}
            for entry in entries}


def _binary_matrix(code, size):
    return tuple(tuple((code >> (size * i + j)) & 1 for j in range(size))
                 for i in range(size))


def _idempotent(matrix):
    size = len(matrix)
    return all(sum(matrix[i][k] * matrix[k][j] for k in range(size)) % 2 == matrix[i][j]
               for i in range(size) for j in range(size))


def _reverse(matrix, transpose=False):
    size = len(matrix)
    return tuple(tuple(matrix[size - 1 - (j if transpose else i)]
                             [size - 1 - (i if transpose else j)]
                       for j in range(size)) for i in range(size))


def candidates(modes, ramification, selected, family='projector'):
    """Return a deduplicated, explicitly finite kernel catalogue."""
    alphabet = states(2, modes)
    levels = sorted(set((1, ramification)))
    if modes == 2:
        middle = [state for state in alphabet if sum(state) == 1]
        block_pairs = [((middle, _binary_matrix(code, 2)),) for code in range(16)]
    elif modes == 3:
        one = [state for state in alphabet if sum(state) == 1]
        two = [state for state in alphabet if sum(state) == 2]
        zero = ((0, 0, 0),) * 3
        matrices = [_binary_matrix(code, 3) for code in range(512)
                    if family == 'all' or _idempotent(_binary_matrix(code, 3))]
        pairs = set()
        for matrix in matrices:
            pairs.add((matrix, zero))
            pairs.add((zero, matrix))
            pairs.add((matrix, _reverse(matrix)))
            pairs.add((matrix, _reverse(matrix, True)))
        block_pairs = [((one, first), (two, second)) for first, second in sorted(pairs)]
    else:
        raise ValueError('This finite catalogue is defined for two or three binary modes.')
    result = {}
    for blocks, top_shift, endpoint_flag, height in product(block_pairs, (0, 1), (0, 1), levels):
        kernel = identity_kernel(2, modes)
        if top_shift:
            top = (1,) * modes
            kernel[top, top] = {-height: 1}
        for labels, matrix in blocks:
            for i, row in enumerate(labels):
                for j, column in enumerate(labels):
                    index = (row, column)
                    value = dict(kernel.get(index, {}))
                    if endpoint_flag and sum(a * b for a, b in zip(selected, row)) > sum(
                            a * b for a, b in zip(selected, column)):
                        value[0] = value.get(0, 0) ^ 1
                    if matrix[i][j]:
                        for degree in (0, -height):
                            value[degree] = value.get(degree, 0) ^ 1
                    value = {degree: coefficient for degree, coefficient in value.items() if coefficient}
                    if value:
                        kernel[index] = value
                    else:
                        kernel.pop(index, None)
        kernel = normalized_kernel(2, modes, kernel, 'filtered', selected)
        key = tuple((a, b, tuple(sorted(value.items()))) for (a, b), value in sorted(kernel.items()))
        result[key] = kernel
    return tuple(result[key] for key in sorted(result))


class Pressure:
    """Full Gaussian pressure plus the exact active replacement."""
    def __init__(self, n, ramification, modes, geometry, blocks, infinity):
        self.n, self.ramification, self.modes = n, ramification, modes
        data = _geometry(2, n, modes, geometry)
        self.selected, self.retained = data['selected'], data['retained']
        self.zero = np.array([order / ramification for block in blocks for order in block['orders']])
        self.charges = np.array([block['charge'] for block in blocks for order in block['orders']])
        self.infinity = np.array(infinity) / ramification

    @staticmethod
    def partition(orders, scale):
        ground = np.min(orders, axis=0)
        return -ground * scale + logsumexp(-(orders - ground) * scale, axis=0)

    def grid(self, scales, charges, relative=False):
        scales, charges = np.broadcast_arrays(np.asarray(scales), np.asarray(charges))
        zero = self.zero[:, None] + self.charges[:, None] * (1 - charges)
        replacement = self.partition(zero, scales)
        replacement += self.partition(self.infinity[:, None], scales)
        for j in self.selected:
            replacement -= np.logaddexp(0, -(2 * j + 1 - charges) * scales)
            replacement -= np.logaddexp(0, -(2 * j + 1 - 2 * self.n) * scales)
        for j in self.retained:
            replacement -= np.logaddexp(0, -(2 * j + 1 + charges) * scales)
        if relative:
            return replacement
        return replacement + np.array([_gaussian(self.n, float(J), float(t))
                                       for t, J in zip(scales.flat, charges.flat)]).reshape(scales.shape)

    def __call__(self, log_scale, charge, relative=False):
        return float(self.grid([math.exp(log_scale)], [charge], relative)[0])


@lru_cache(maxsize=32768)
def _gaussian(n, charge, scale):
    return gaussian_pressure(2, n, charge, scale)


def run(n, ramification, modes, geometry, family='projector', refinements=30,
        limit=None, progress=True):
    start = time.monotonic()
    selected = _geometry(2, n, modes, geometry)['selected']
    catalogue = candidates(modes, ramification, selected, family)
    if limit is not None:
        catalogue = catalogue[:limit]
    scale_values = math.log(2 * n) / n * np.array((.25, .5, .75, 1., 1.5, 2.))
    charge_values = sorted(set((-2 * n, -n, -4, -2, -1, 0, 1, 2, 4, n / 2, n, 2 * n)))
    grid = tuple(product(scale_values, charge_values))
    scales = np.array([point[0] for point in grid])
    charges = np.array([point[1] for point in grid])
    rows, best, relative_best = [], math.inf, 0.
    for index, kernel in enumerate(catalogue):
        blocks = zero_profiles(2, n, ramification, modes, kernel, geometry, 'filtered')
        infinity = infinity_profiles(2, n, ramification, modes, kernel, geometry, 'filtered')
        pressure = Pressure(n, ramification, modes, geometry, blocks, infinity)
        values = pressure.grid(scales, charges)
        relative = pressure.grid(scales, charges, True)
        i, j = int(np.argmin(values)), int(np.argmin(relative))
        row = {'kernel': encode_kernel(kernel), 'zero_blocks': blocks, 'infinity_active': infinity,
               'critical_grid': {'pressure': float(values[i]), 'scale': float(scales[i]),
                                 'charge': float(charges[i])},
               'relative_grid': {'change': float(relative[j]), 'scale': float(scales[j]),
                                 'charge': float(charges[j])}}
        rows.append(row)
        if progress and values[i] < best - 1e-9:
            best = float(values[i])
            print(json.dumps({'progress': index + 1, 'total': len(catalogue), 'best': row['critical_grid']}), flush=True)
        if progress and relative[j] < relative_best - 1e-9:
            relative_best = float(relative[j])
            print(json.dumps({'progress': index + 1, 'relative_improvement': row['relative_grid']}), flush=True)
        if progress and (index + 1) % 100 == 0:
            print(json.dumps({'progress': index + 1, 'total': len(catalogue), 'seconds': time.monotonic() - start}), flush=True)
    rows.sort(key=lambda row: row['critical_grid']['pressure'])
    for row in rows[:refinements]:
        pressure = Pressure(n, ramification, modes, geometry, row['zero_blocks'], row['infinity_active'])
        grid_best = row['critical_grid']
        optimum = minimize(lambda point: pressure(*point),
                           [math.log(grid_best['scale']), grid_best['charge']], method='Nelder-Mead',
                           bounds=[(-math.log(n) - 4, math.log(1.5)), (-2 * n, 2 * n)],
                           options={'maxiter': 200, 'fatol': 1e-10})
        refined = {'pressure': float(optimum.fun), 'scale': math.exp(float(optimum.x[0])),
                   'charge': float(optimum.x[1]), 'optimizer_success': bool(optimum.success)}
        row['refined'] = refined if refined['pressure'] <= grid_best['pressure'] else dict(grid_best)
    return {'status': 'Finite exact-profile screen with numerical pressure optimization; not an obstruction theorem.',
            'prime': 2, 'n': n, 'ramification': ramification, 'modes': modes, 'geometry': geometry,
            'family': family, 'endpoint_policy': 'filtered', 'candidate_count': len(catalogue),
            'grid_size': len(grid), 'refinement_count': min(refinements, len(rows)),
            'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in rows),
            'best_grid_relative_change': min(row['relative_grid']['change'] for row in rows),
            'rows': rows, 'seconds': time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--P', type=int, default=1)
    parser.add_argument('--modes', type=int, choices=(2, 3), default=2)
    parser.add_argument('--geometry', choices=('initial', 'mixed', 'terminal'), default='initial')
    parser.add_argument('--family', choices=('projector', 'all'), default='projector')
    parser.add_argument('--refinements', type=int, default=30)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.n, args.P, args.modes, args.geometry, args.family, args.refinements, args.limit)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'rows'}), flush=True)


if __name__ == '__main__':
    main()
