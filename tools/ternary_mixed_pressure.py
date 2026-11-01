#!/usr/bin/env python3
"""Weighted mixed-place pressure for the ternary seed.

The complete local profiles use exact polynomial-minor valuations; continuous
optimization uses floating arithmetic and supplies no certificate by itself.
"""
import argparse
from itertools import combinations, permutations, product
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, softmax, expit

from check_multipoint_hermite_crossings import mixing_seed, local_polynomial, multiply
from check_weighted_multipoint_hermite import local_gauged_profile


def determinant(matrix, prime=3):
    size = len(matrix)
    result = [0] * (2 * size + 1)
    for perm in permutations(range(size)):
        value = [1]
        for i, j in enumerate(perm):
            value = multiply(value, matrix[i][j], prime)
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(size) for j in range(i + 1, size))
        for k, coefficient in enumerate(value):
            result[k] = (result[k] + sign * coefficient) % prime
    return result


class Oracle:
    def __init__(self, point, cut):
        self.point, self.cut = point, cut
        seed = mixing_seed()
        table = {key: local_polynomial(value, point, 2, 3) for key, value in seed.items()}
        self.blocks = []
        for charge in range(3):
            if cut == 0:
                rows = [(a, b) for a, b in product(range(3), repeat=2) if (a+b) % 3 == charge]
                columns = rows
            else:
                rows = [(a, b) for a, b in product(range(3), repeat=2) if (a-b) % 3 == charge]
                columns = [(a, b) for a, b in product(range(3), repeat=2) if (b-a) % 3 == charge]
            matrix = [[table.get((a,b,c,d) if cut == 0 else (a,c,b,d), [0,0,0])
                       for c, d in columns] for a, b in rows]
            coefficients = []
            for size in range(1, 4):
                records = []
                for rs in combinations(range(3), size):
                    for cs in combinations(range(3), size):
                        polynomial = determinant([[matrix[i][j] for j in cs] for i in rs])
                        if not any(polynomial):
                            continue
                        valuation = next(i for i, x in enumerate(polynomial) if x)
                        counts = [0] * 8
                        for i in rs:
                            for leg, symbol in zip((0,1) if cut == 0 else (0,2), rows[i]):
                                if symbol:
                                    counts[2*leg + symbol-1] += 1
                        for j in cs:
                            for leg, symbol in zip((2,3) if cut == 0 else (1,3), columns[j]):
                                if symbol:
                                    counts[2*leg + symbol-1] += 1
                        records.append([valuation] + counts)
                coefficients.append(np.array(records, dtype=float))
            self.blocks.append(coefficients)

    def profiles(self, rate, weights, derivatives=False):
        values, gradients = [], []
        variable = np.concatenate(([rate], weights))
        for block in self.blocks:
            prefixes = [0.0]
            prefix_gradients = [np.zeros(9)]
            for records in block:
                candidates = records @ variable
                index = int(np.argmin(candidates))
                prefixes.append(float(candidates[index]))
                prefix_gradients.append(records[index])
            values.extend(np.diff(prefixes))
            gradients.extend(np.diff(prefix_gradients, axis=0))
        values = np.array(values)
        if derivatives:
            return values, np.array(gradients)
        return values

    def cost(self, rate, weights, px, py, derivatives=False):
        orders, gradients = self.profiles(rate, weights, True)
        means = np.concatenate((px[1:], py[1:], py[1:], px[1:]))
        z = logsumexp(-orders)
        cost = z + means @ weights
        if not derivatives:
            return cost
        grad = -softmax(-orders) @ gradients
        grad[1:] += means
        physical = weights.reshape((4,2))
        pgrad_x = np.concatenate(([0], physical[0] + physical[3]))
        pgrad_y = np.concatenate(([0], physical[1] + physical[2]))
        return cost, grad, pgrad_x, pgrad_y


ORACLES = [Oracle(2, 0), Oracle(None, 0), Oracle(0, 1)]


def unpack(vector):
    t = math.exp(vector[0])
    alpha = expit(vector[1])
    px = softmax([0, vector[2], vector[3]])
    py = softmax([0, vector[4], vector[5]])
    weights = [vector[6+8*i:14+8*i] for i in range(3)]
    return t, alpha, px, py, weights


def objective(vector, gradient=False):
    t, alpha, px, py, weights = unpack(vector)
    if min(alpha, 1-alpha) < 1e-7:
        return (100, np.zeros_like(vector)) if gradient else 100
    rates = [t/alpha, t/(1-alpha), t]
    records = [oracle.cost(rate, weight, px, py, True)
               for oracle, rate, weight in zip(ORACLES, rates, weights)]
    selected = 0 if records[0][0] >= records[1][0] else 1
    numerator = 2*t + records[selected][0] + records[2][0]
    entropy = -sum(px*np.log(px)) - sum(py*np.log(py))
    value = numerator/(2*entropy)
    if not gradient:
        return value
    grad = np.zeros_like(vector)
    physical_x = np.zeros(3)
    physical_y = np.zeros(3)
    grad[0] += 2*t
    for i in (selected, 2):
        _, local, p_x, p_y = records[i]
        grad[0] += local[0]*rates[i]
        if i == 0:
            grad[1] -= local[0]*rates[i]*(1-alpha)
        elif i == 1:
            grad[1] += local[0]*rates[i]*alpha
        grad[6+8*i:14+8*i] += local[1:]
        physical_x += p_x
        physical_y += p_y
    for start, pi, pg in ((2,px,physical_x),(4,py,physical_y)):
        jacobian = np.diag(pi)-np.outer(pi,pi)
        derivative = jacobian @ pg
        entropy_gradient = jacobian @ (-np.log(pi)-1)
        grad[start:start+2] = derivative[1:] - (numerator/entropy)*entropy_gradient[1:]
    return value, grad/(2*entropy)


def controls():
    rng = np.random.default_rng(70431)
    cases = 0
    for oracle in ORACLES:
        for _ in range(12):
            weights = rng.integers(-3,4,size=8)
            direct = local_gauged_profile(mixing_seed(),3,2,oracle.point,3,oracle.cut,
                                         [(0,int(weights[2*i]),int(weights[2*i+1])) for i in range(4)])
            assert sorted(oracle.profiles(1,weights)) == direct
            cases += 1
    for _ in range(12):
        vector = rng.normal(0,.3,size=30)
        vector[0] -= .7
        value, grad = objective(vector, True)
        numeric = np.array([(objective(vector + np.eye(30)[i]*1e-6)-objective(vector - np.eye(30)[i]*1e-6))/2e-6
                            for i in range(30)])
        assert np.max(np.abs(grad-numeric)) < 1e-5, (np.max(np.abs(grad-numeric)), grad, numeric)
    return {'weighted_profile_controls': cases, 'gradient_controls':12}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=50)
    parser.add_argument('--seed',type=int,default=41)
    parser.add_argument('--output',type=Path,default=Path('/tmp/ternary-mixed-search.json'))
    parser.add_argument('--controls',action='store_true')
    args = parser.parse_args()
    if args.controls:
        print(json.dumps(controls()), flush=True)
    rng = np.random.default_rng(args.seed)
    base = np.zeros(30)
    base[0] = math.log(.51)
    baseline = objective(base)
    best, best_x = baseline, base.copy()
    records = []
    started = time.monotonic()
    run = 0
    while time.monotonic()-started < args.seconds:
        if run == 0:
            start = base.copy()
        elif run % 3:
            start = best_x+rng.normal(0, .25 if run % 3 == 1 else .8, 30)
        else:
            start = base+rng.normal(0, .5, 30)
        result = minimize(lambda x: objective(x, True), start, jac=True, method='L-BFGS-B',
                          bounds=[(-4,1.5),(-3.5,3.5)]+[(-5,5)]*4+[(-8,8)]*24,
                          options={'maxiter':700,'maxls':40,'ftol':1e-12,'gtol':1e-7})
        records.append({'run':run,'rho':float(result.fun),'iterations':int(result.nit),
                        'success':bool(result.success)})
        if result.fun < best:
            best, best_x = float(result.fun), result.x.copy()
            t, alpha, px, py, weights = unpack(best_x)
            print(json.dumps({'run':run,'best':best,'t':t,'alpha':alpha,'px':px.tolist(),'py':py.tolist()}),flush=True)
        run += 1
    t, alpha, px, py, weights = unpack(best_x)
    report = {'baseline':baseline,'best_rho':best,'best_vector':best_x.tolist(),'t':t,'alpha':alpha,
              'px':px.tolist(),'py':py.tolist(),'scaled_heights':[w.tolist() for w in weights],
              'records':records,'elapsed_seconds':time.monotonic()-started,'seed':args.seed,
              'scope':'bounded floating search; no exact certificate or global optimality claim'}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('baseline','best_rho','elapsed_seconds')}),flush=True)


if __name__=='__main__':
    main()
