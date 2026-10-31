#!/usr/bin/env python3
"""Bounded numerical fixed-paired-type Smith capacity search for quantum Pascal.

Primal values and LP lower estimates are numerical, not exact certificates.
The separate minor-flag proof justifies the convex supporting cuts. Missing
native pivots trigger higher precision; they are never treated as zero.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb, log, sqrt
from pathlib import Path

from native_matrix_smith_f4 import MatrixOracle
from quantum_pascal_crossing import construct, cut_blocks


def binomial_type(d, probability):
    return [comb(d-1, i)*probability**i*(1-probability)**(d-1-i) for i in range(d)]


class Pressure:
    def __init__(self, d, probability_x, probability_y):
        import numpy as np
        self.np, self.d = np, d
        tensor, self.degree, _ = construct(d)
        self.types = np.array([binomial_type(d, probability_x), binomial_type(d, probability_y)])
        self.entropy = float(-np.sum(self.types*np.log(self.types)))
        self.leg_types = self.types[[0, 1, 1, 0]]
        self.blocks = []
        for cut in (0, 1):
            side = []
            for _, rows, columns, matrix in cut_blocks(tensor, d, cut):
                packed = [[(value, 0) for value in row] for row in matrix]
                side.append(dict(matrix=packed, rows=rows, columns=columns, precision=512,
                                 oracle=MatrixOracle(packed, len(rows), precision=512)))
            self.blocks.append(side)
        self.calls = 0

    def evaluate(self, parameters):
        np, d = self.np, self.d
        scale = float(parameters[0])
        energies = np.zeros((2, 4, d))
        energies[:, :, 1:] = np.array(parameters[1:]).reshape(2, 4, d-1)
        gradient = np.zeros_like(parameters)
        gauge_gradient = np.broadcast_to(self.leg_types, (2, 4, d)).copy()
        value = float(np.sum(self.leg_types*energies.sum(axis=0)))-sqrt(2)*self.entropy
        for side, blocks in enumerate(self.blocks):
            orders, descriptions = [], []
            for block in blocks:
                rows, columns = block['rows'], block['columns']
                row_weights = [energies[side, 0, a]+energies[side, 1 if side == 0 else 2, b]
                               for a, b in rows]
                column_weights = [energies[side, 2 if side == 0 else 1, a]+energies[side, 3, b]
                                  for a, b in columns]
                while True:
                    try:
                        result = block['oracle'].profile(row_weights, column_weights,
                                                         scale=scale, infinity=side == 0)
                        break
                    except ArithmeticError as error:
                        if 'status 3' not in str(error) or block['precision'] >= 8192:
                            raise
                        block['precision'] *= 2
                        block['oracle'] = MatrixOracle(block['matrix'], len(rows), precision=block['precision'])
                previous = 0
                for order, total, row, column in zip(result['orders'], result['base_prefix'],
                                                      result['rows'], result['cols']):
                    base = total-previous
                    previous = total
                    letters = (rows[row]+columns[column] if side == 0
                               else (rows[row][0], columns[column][0], rows[row][1], columns[column][1]))
                    assert abs(order-scale*base-sum(energies[side, leg, a] for leg, a in enumerate(letters))) < 1e-7
                    orders.append(order)
                    descriptions.append((base, letters))
            orders = np.array(orders)
            minimum = float(min(orders))
            weights = np.exp(-(orders-minimum))
            total = float(sum(weights))
            value += -minimum+log(total)
            weights /= total
            for weight, (base, letters) in zip(weights, descriptions):
                gradient[0] -= weight*base
                for leg, letter in enumerate(letters):
                    gauge_gradient[side, leg, letter] -= weight
        gradient[1:] = gauge_gradient[:, :, 1:].ravel()
        self.calls += 1
        return value/(2*self.entropy), gradient/(2*self.entropy)


def search(d, probability_x, probability_y, iterations):
    import numpy as np
    from scipy.optimize import linprog
    oracle = Pressure(d, probability_x, probability_y)
    size = 1+8*(d-1)
    bounds = [(.002, 1.)]+[(-2., 2.)]*(size-1)
    point = np.zeros(size)
    point[0] = 1/max(2, oracle.degree)
    cuts, samples, best, best_point, lower = [], [], float('inf'), point.copy(), float('-inf')
    target = 1/sqrt(2)
    for iteration in range(iterations):
        value, gradient = oracle.evaluate(point)
        # Check all retained supports on every new observation.
        assert all(float(np.dot(g, point))+b <= value+3e-7 for g, b in cuts)
        intercept = value-float(np.dot(gradient, point))
        assert all(v+3e-7 >= float(np.dot(gradient, x))+intercept for x, v in samples)
        samples.append((point.copy(), value))
        cuts.append((gradient, intercept))
        if value < best:
            best, best_point = value, point.copy()
        matrix = np.array([np.r_[g, -1.] for g, _ in cuts])
        vector = np.array([-b for _, b in cuts])
        result = linprog(np.r_[np.zeros(size), 1.], A_ub=matrix, b_ub=vector,
                         bounds=bounds+[(None, None)], method='highs')
        if not result.success:
            raise ArithmeticError(result.message)
        lower = float(result.fun)
        if best-lower < 2e-7:
            break
        point = result.x[:-1]
    return dict(alphabet=d, degree=oracle.degree, probability_x=probability_x, probability_y=probability_y,
                types=oracle.types.tolist(), entropy_sum=oracle.entropy,
                best_rho=target+best, numerical_domain_lower_rho=target+lower,
                numerical_gap=best-lower, best_parameters=best_point.tolist(),
                iterations=iteration+1, maximum_iterations=iterations, pressure_calls=oracle.calls,
                maximum_native_precision=max(b['precision'] for side in oracle.blocks for b in side),
                scale_bounds=[.002, 1.], physical_energy_bounds=[-2., 2.],
                letter_zero_energy=0, numerical_only=True,
                scope='Only the stated fixed paired binomial types and compact parameter box; no rational LP certification.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alphabet', type=int, default=3)
    parser.add_argument('--probability-x', type=float, default=.75)
    parser.add_argument('--probability-y', type=float, default=.75)
    parser.add_argument('--iterations', type=int, default=500)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.alphabet < 2 or args.iterations < 1 or not (0 < args.probability_x < 1 and 0 < args.probability_y < 1):
        parser.error('use an alphabet >=2, positive iteration count, and probabilities strictly between zero and one')
    result = search(args.alphabet, args.probability_x, args.probability_y, args.iterations)
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('alphabet', 'best_rho', 'numerical_domain_lower_rho', 'numerical_gap', 'iterations', 'factoring_exponent')}))


if __name__ == '__main__':
    main()
