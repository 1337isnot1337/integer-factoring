#!/usr/bin/env python3
"""Exploratory full physical-letter gauges for literal F16 fused crossings."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp, softmax

from fused_six_vertex_crossings import construct, transform, oracles


class Pressure:
    def __init__(self, matrix, q):
        self.q = q
        self.matrix = matrix
        self.precision = 512
        self.cuts = oracles(matrix, q)
        self.queries = 0
        self.last = None
        self.active = [side*4*q+leg*q+letter for side in (0, 1)
                       for leg in range(4) for letter in range(1, q)] + [8*q]

    def unpack(self, values):
        result = np.zeros(8*self.q+1)
        result[self.active] = values
        return result

    def evaluate(self, values, gradient=True):
        x = self.unpack(values)
        q = self.q
        h = x[:-1].reshape(2, 4, q)
        scale = x[-1]
        if scale <= 0:
            raise ValueError('Positive scale required.')
        derivative = np.zeros_like(x)
        # Restriction F16/F4 repeats each channel twice. Physical labels stay q.
        log_b = -2*math.log(2)
        records = []
        labels = list(product(range(q), repeat=2))
        for side, oracle in enumerate(self.cuts):
            rowlegs = (0, 1) if side == 0 else (0, 2)
            collegs = (2, 3) if side == 0 else (1, 3)
            rows = np.repeat([h[side, rowlegs[0], a]+h[side, rowlegs[1], b]
                              for a, b in labels], 2)
            cols = np.repeat([h[side, collegs[0], a]+h[side, collegs[1], b]
                              for a, b in labels], 2)
            while True:
                try:
                    result = self.cuts[side].profile(rows, cols, scale, infinity=bool(side))
                    break
                except ArithmeticError:
                    if self.precision >= 8192:
                        raise
                    self.precision *= 2
                    self.cuts = oracles(self.matrix, q, self.precision)
            orders = np.asarray(result['orders'])
            if np.max(np.abs(orders[::2]-orders[1::2])) > 2e-7:
                raise ArithmeticError('Scalar-restricted Smith multiplicities are unresolved.')
            log_b += logsumexp(-orders)
            coefficients = -softmax(-orders)
            base = np.diff([0]+result['base_prefix'])
            for coefficient, row, col, valuation in zip(coefficients, result['rows'], result['cols'], base):
                for leg, letter in zip(rowlegs, labels[row//2]):
                    derivative[side*4*q+leg*q+letter] += coefficient
                for leg, letter in zip(collegs, labels[col//2]):
                    derivative[side*4*q+leg*q+letter] += coefficient
                derivative[-1] += coefficient*valuation
            records.append(result)
        physical = [sum((h[side, leg] for side in (0, 1) for leg in legs), np.zeros(q))
                    for legs in ((0, 3), (1, 2))]
        minima = [np.min(row) for row in physical]
        centered = [row-low for row, low in zip(physical, minima)]
        constant = float(log_b+sum(minima))
        def equation(rate):
            return constant-rate*sum(logsumexp(-row/rate) for row in centered)
        if constant < 1e-10 or equation(1e-7) <= 0:
            raise ArithmeticError('Positive root unresolved near zero entropy.')
        upper = 4.
        while equation(upper) > 0:
            upper *= 2
        rate = brentq(equation, 1e-7, upper, xtol=1e-12)
        probabilities = [softmax(-row/rate) for row in centered]
        entropy = sum(logsumexp(-row/rate)+np.dot(row/rate, p)
                      for row, p in zip(centered, probabilities))
        if entropy < 1e-8:
            raise ArithmeticError('Physical entropy too small.')
        for side in (0, 1):
            for leg in range(4):
                derivative[side*4*q+leg*q:side*4*q+(leg+1)*q] += probabilities[0 if leg in (0, 3) else 1]
        self.queries += 1
        self.last = {'rho': float(rate/2), 'scale': scale, 'parameters': x.tolist(),
                     'critical_pressure': float(equation(math.sqrt(2))),
                     'physical_entropy': float(entropy), 'profiles': records,
                     'physical_energies': [row.tolist() for row in physical],
                     'physical_probabilities': [row.tolist() for row in probabilities]}
        return (rate/2, derivative[self.active]/(2*entropy)) if gradient else rate/2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spin', type=int, default=2)
    parser.add_argument('--theta', type=int, default=2)
    parser.add_argument('--starts', type=int, default=12)
    parser.add_argument('--seconds', type=float, default=900)
    parser.add_argument('--seed', type=int, default=311857)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    source = construct(args.spin, args.theta)
    assert source['invariant'] and source['endpoint_identity']
    q = args.spin+1
    # Include the natural coordinate pair and pairs involving nontrivial roots.
    pairs = [(0, None), (2, 3), (9, 15), (14, 10), (8, 10), (None, 0)]
    rng = np.random.default_rng(args.seed)
    records = []
    result = {'status': 'Finite literal F16 fusion and exploratory full-symbol gauge optimization.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source': {key: value for key, value in source.items() if key != 'matrix'},
              'requested_projective_pairs': pairs,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in (Path(__file__).name, 'fused_six_vertex_crossings.py',
                                             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')}}
    for a, b in pairs:
        if time.monotonic()-started >= args.seconds:
            break
        pressure = Pressure(transform(source['matrix'], source['degree'], a, b), q)
        for start in range(args.starts):
            if time.monotonic()-started >= args.seconds:
                break
            scale = [0.12, 0.35, 0.8, 1.6][start % 4]
            x = np.zeros(8*q+1)
            x[-1] = scale
            if start:
                weights = np.arange(q)*([1, 2, -1, -2][(start//4) % 4])*scale
                for side, leg in ((0, 0), (0, 1), (1, 2), (1, 3)):
                    x[side*4*q+leg*q:side*4*q+(leg+1)*q] = weights
                x[:-1] += rng.normal(0, 0.02*scale, len(x)-1)
                for offset in range(0, 8*q, q):
                    x[offset:offset+q] -= x[offset]
            best = None
            def objective(values):
                nonlocal best
                weighted = np.r_[values[:-1]*values[-1], values[-1]]
                value, derivative = pressure.evaluate(weighted)
                if best is None or value < best['rho']:
                    best = dict(pressure.last)
                transformed = np.r_[derivative[:-1]*values[-1],
                                    derivative[-1]+np.dot(derivative[:-1], values[:-1])]
                return value, transformed
            initial = x[pressure.active]
            initial[:-1] /= initial[-1]
            try:
                optimized = minimize(objective, initial, method='L-BFGS-B', jac=True,
                                     bounds=[(-12., 12.)]*(8*(q-1))+[(0.002, 6.)],
                                     options={'maxiter': 180, 'maxfun': 1200, 'ftol': 1e-12, 'gtol': 1e-8})
                success, message = bool(optimized.success), str(optimized.message)
            except ArithmeticError as error:
                if best is None:
                    raise
                success, message = False, 'Stopped after unresolved numerical query: '+str(error)
            records.append({'zero_image': a, 'infinity_image': b, 'start': start,
                            'optimizer_success': success,
                            'optimizer_message': message, 'queries': pressure.queries,
                            **best})
            result.update(records=records, completed_starts=len(records),
                          completed_requested_starts=len(records) == len(pairs)*args.starts,
                          elapsed_seconds=time.monotonic()-started,
                          best=min(records, key=lambda row: row['rho']))
            args.output.write_text(json.dumps(result, indent=2)+'\n')
            print(json.dumps({'pair': [a, b], 'start': start, 'rho': best['rho'],
                              'global_best': result['best']['rho'],
                              'seconds': time.monotonic()-started}), flush=True)
    print(json.dumps({'event': 'finished', 'starts': len(records), 'best': result['best']['rho'],
                      'complete': result['completed_requested_starts']}), flush=True)


if __name__ == '__main__':
    main()
