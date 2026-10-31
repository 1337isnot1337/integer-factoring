#!/usr/bin/env python3
"""Exploratory continuous mode gauges for the direct polynomial Schur kernel.

Uses exact polynomial determinantal valuations and a separation oracle.
The optimization uses floating-point arithmetic and does not certify bounds.
Requires NumPy and SciPy. Fixed mode counts justify the direct endpoint for
arbitrary index sets; see research/fixed-mode-crossings.md.
"""
from __future__ import annotations

import argparse
import json
from math import log, log1p
from pathlib import Path
import time

import numpy as np
from scipy.optimize import LinearConstraint, brentq, minimize
from scipy.special import expit

from check_schur_downset_crossing import sparse_kernel
from weighted_smith_oracle import WeightedSmithOracle


class ModePressureSearch:
    def __init__(self, indices, ambient, orientation=1, max_rounds=80):
        self.indices = tuple(indices)
        self.ambient = ambient
        self.n = n = len(indices)
        if not n or len(set(indices)) != n or any(a < 0 or a >= ambient for a in indices):
            raise ValueError('indices must be a nonempty subset of [0, ambient)')
        if ambient < 1 or ambient & (ambient - 1):
            raise ValueError('ambient must be a power of two')
        self.max_rounds = max_rounds
        self.dimension = 24 * n + 1
        self.scale_index = 16 * n
        self.order_base = self.scale_index + 1
        self.blocks = []
        self.constraints = {}
        equalities = []
        kernel = sparse_kernel(indices, ambient)
        self.type_map = np.zeros((2 * n, self.dimension))
        for endpoint in (0, 1):
            base = 8 * n * endpoint
            for group in range(4):
                for i in range(n):
                    self.type_map[(0 if group >= 2 else 1) * n + i, base + group * n + i] = 1
                    self.type_map[(1 if group >= 2 else 0) * n + i, base + 4 * n + group * n + i] = 1
        for endpoint, cut in ((0, 1 - orientation), (1, orientation)):
            groups = ((1, 3), (2, 0)) if cut == 0 else ((1, 0), (2, 3))
            base = 8 * n * endpoint
            for direction in (0, 1):
                rows = tuple(g * n + i for g in groups[direction] for i in range(n))
                columns = tuple(g * n + i for g in groups[1 - direction] for i in range(n))
                matrix = [[kernel[i][j] for j in columns] for i in rows]
                block = {'endpoint': endpoint, 'rows': rows, 'columns': columns,
                         'oracle': WeightedSmithOracle(matrix, infinity=bool(endpoint))}
                self.blocks.append(block)
                index = len(self.blocks) - 1
                for i, row in enumerate(matrix):
                    for j, value in enumerate(row):
                        if value:
                            self.add_minor(index, 1, 1 << i, 1 << j, block['oracle'].valuation(value))
                _, prefix = block['oracle'].profile([0] * (2 * n), [0] * (2 * n), 0)
                last = prefix[-1]
                equalities.append(self.minor_vector(index, last.size, last.row_mask,
                                                     last.column_mask, last.valuation))
            for group in groups:
                for a in group:
                    for b in group:
                        for i in range(n):
                            for j in range(n):
                                value = kernel[a * n + i][b * n + j]
                                if not value:
                                    continue
                                valuation = 1 - value.bit_length() if endpoint else (value & -value).bit_length() - 1
                                row = np.zeros(self.dimension)
                                row[base + a * n + i] = 1
                                row[base + 4 * n + b * n + j] = 1
                                row[self.scale_index] = valuation
                                self.constraints['internal', endpoint, a, b, i, j] = row
        for block in range(4):
            for i in range(2 * n - 1):
                row = np.zeros(self.dimension)
                row[self.order_base + 2 * n * block + i + 1] = 1
                row[self.order_base + 2 * n * block + i] = -1
                self.constraints['sorted', block, i] = row
        row = np.zeros(self.dimension)
        row[self.scale_index] = 1
        self.constraints['nonnegative_scale'] = row
        self.equalities = np.array(equalities)

    def minor_vector(self, index, size, row_mask, column_mask, valuation):
        block = self.blocks[index]
        n = self.n
        base = 8 * n * block['endpoint']
        row = np.zeros(self.dimension)
        row[self.scale_index] = valuation
        for i in range(2 * n):
            if row_mask >> i & 1:
                row[base + block['rows'][i]] = 1
            if column_mask >> i & 1:
                row[base + 4 * n + block['columns'][i]] = 1
        row[self.order_base + 2 * n * index:self.order_base + 2 * n * index + size] = -1
        return row

    def add_minor(self, index, size, row_mask, column_mask, valuation):
        key = ('minor', index, size, row_mask, column_mask)
        if key not in self.constraints:
            self.constraints[key] = self.minor_vector(index, size, row_mask, column_mask, valuation)

    def initial_guess(self):
        n, m = self.n, self.ambient
        x = np.zeros(self.dimension)
        ell = x[self.scale_index] = 0.5 / m
        rows = [0] * n + [1 + 2 * a for a in self.indices] * 2 + [0] * n
        columns = ([2 * m] * n + [4 * m - 1 - 2 * a for a in self.indices]
                   + [2 * m - 1 - 2 * a for a in self.indices] + [0] * n)
        x[8 * n:16 * n] = ell * np.array(rows + columns)
        for index, block in enumerate(self.blocks):
            base = 8 * n * block['endpoint']
            orders, _ = block['oracle'].profile([x[base + i] for i in block['rows']],
                                                [x[base + 4 * n + j] for j in block['columns']], ell)
            x[self.order_base + 2 * n * index:self.order_base + 2 * n * (index + 1)] = orders
        return x

    def solve_fixed_type(self, probabilities, start=None):
        n = self.n
        probabilities = np.asarray(probabilities)
        if len(probabilities) != 2 * n or any(not 0 < p < 1 for p in probabilities):
            raise ValueError('one probability strictly between zero and one per species and mode is required')
        entropy = sum(-p * log(p) - (1 - p) * log1p(-p) for p in probabilities)
        linear = probabilities @ self.type_map
        def objective(x):
            gradient = linear.copy()
            gradient[self.order_base:] = -expit(-x[self.order_base:])
            value = linear @ x + np.logaddexp(0, -x[self.order_base:]).sum()
            return float(value), gradient
        x = self.initial_guess() if start is None else np.asarray(start)
        total_iterations = 0
        for round_index in range(self.max_rounds):
            matrix = np.array(list(self.constraints.values()))
            result = minimize(objective, x, jac=True, method='SLSQP',
                              constraints=[LinearConstraint(self.equalities, 0, 0),
                                           LinearConstraint(matrix, 0, np.inf)],
                              options={'maxiter': 800, 'ftol': 1e-11})
            x = result.x
            total_iterations += result.nit
            violations = 0
            worst = 0.0
            for index, block in enumerate(self.blocks):
                base = 8 * n * block['endpoint']
                _, prefixes = block['oracle'].profile(
                    [x[base + i] for i in block['rows']],
                    [x[base + 4 * n + j] for j in block['columns']],
                    max(0.0, x[self.scale_index]))
                for prefix in prefixes[:-1]:
                    row = self.minor_vector(index, prefix.size, prefix.row_mask,
                                            prefix.column_mask, prefix.valuation)
                    slack = float(row @ x)
                    worst = min(worst, slack)
                    if slack < -2e-8:
                        self.add_minor(index, prefix.size, prefix.row_mask,
                                       prefix.column_mask, prefix.valuation)
                        violations += 1
            if not violations:
                break
        return {'rho': float(result.fun / (4 * entropy)), 'probabilities': probabilities.tolist(),
                'log_s': float(x[self.scale_index]), 'optimizer_success': bool(result.success),
                'optimizer_message': result.message, 'oracle_minimum_slack': worst,
                'listed_minimum_slack': float(min(matrix @ x)),
                'determinant_equality_error': float(max(abs(self.equalities @ x))),
                'cutting_rounds': round_index + 1, 'optimizer_iterations': total_iterations,
                'constraints': len(self.constraints),
                'oracle_cache_sizes': [len(block['oracle'].cache) for block in self.blocks],
                'solution': x.tolist()}

    def optimal_independent_type(self, solution):
        x = np.asarray(solution)
        weights = self.type_map @ x
        raw_log_pressure = np.logaddexp(0, -x[self.order_base:]).sum()
        def difference(r):
            return raw_log_pressure - r * np.logaddexp(0, -weights / r).sum()
        lower, upper = 1e-8, 40.0
        if difference(lower) <= 0:
            raise ArithmeticError('near-zero pressure requires an exact audit')
        while difference(upper) > 0:
            upper *= 2
        r = brentq(difference, lower, upper)
        return expit(-weights / r), r / 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--indices', default='0,1,2,4')
    parser.add_argument('--ambient', type=int, default=8)
    parser.add_argument('--orientation', type=int, choices=(0, 1), default=1)
    parser.add_argument('--iterations', type=int, default=20)
    parser.add_argument('--cut-rounds', type=int, default=80)
    parser.add_argument('--p-a', type=float, default=0.23)
    parser.add_argument('--p-b', type=float, default=0.23)
    parser.add_argument('--fixed-type', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    indices = tuple(map(int, args.indices.split(',')))
    search = ModePressureSearch(indices, args.ambient, args.orientation, args.cut_rounds)
    probabilities = np.array([args.p_a] * len(indices) + [args.p_b] * len(indices))
    solution = None
    history = []
    start = time.monotonic()
    for iteration in range(args.iterations):
        result = search.solve_fixed_type(probabilities, solution)
        history.append(result)
        solution = result['solution']
        next_probabilities, next_rho = search.optimal_independent_type(solution)
        print(json.dumps({'iteration': iteration, 'rho': result['rho'], 'next_rho': next_rho,
                          'next_probabilities': next_probabilities.tolist(),
                          'optimizer_success': result['optimizer_success'],
                          'oracle_minimum_slack': result['oracle_minimum_slack'],
                          'cutting_rounds': result['cutting_rounds'],
                          'seconds': time.monotonic() - start}), flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            report = {'status': 'exploratory floating-point search; no certified exponent',
                      'indices': indices, 'ambient': args.ambient, 'orientation': args.orientation,
                      'history': history}
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        if args.fixed_type or (len(history) > 1 and abs(history[-2]['rho'] - result['rho']) < 1e-9):
            break
        probabilities = np.clip(next_probabilities, 1e-9, 1 - 1e-9)


if __name__ == '__main__':
    main()
