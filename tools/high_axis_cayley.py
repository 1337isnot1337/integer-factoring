#!/usr/bin/env python3
"""Specialized one-particle profiles of the explicit homogeneous Cayley kernel.

Specializing the other variables only yields optimistic generic prefix
orders. This diagnostic never promotes a favorable specialization to a
profile certificate for the multivariate kernel.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from native_matrix_smith_f4 import MatrixOracle


MUL = ((0,0,0,0), (0,1,2,3), (0,2,3,1), (0,3,1,2))


def crossblock(dimension, other_values=None):
    """Multiplication in F4[epsilon_1,...]/(epsilon_i^2), axis zero."""
    if not 2 <= dimension <= 9:
        raise ValueError('The registered diagnostic uses dimensions two through nine.')
    if other_values is None:
        other_values = [1]*(dimension-1)
    if len(other_values) != dimension-1 or any(v not in (1,2,3) for v in other_values):
        raise ValueError('Supply one nonzero F4 specialization per other variable.')
    degree, n = math.lcm(*range(1, dimension+1)), 1 << (dimension-1)
    coefficients = []
    for subset in range(n):
        order, coefficient = degree//(subset.bit_count()+1), 1
        for j, value in enumerate(other_values):
            if subset >> j & 1:
                for _ in range(order % 3):
                    coefficient = MUL[coefficient][value]
        coefficients.append(((1 << order) if coefficient & 1 else 0,
                             (1 << order) if coefficient & 2 else 0))
    matrix = [[coefficients[i ^ j] if j & ~i == 0 else (0,0)
               for j in range(n)] for i in range(n)]
    return matrix, degree


def pressure_fit(orders, degree, dimension):
    values = np.asarray(orders, dtype=float)/degree
    def objective(v):
        p, beta = 1/(1+math.exp(-float(v[0]))), math.exp(float(v[1]))
        entropy = -p*math.log(p)-(1-p)*math.log1p(-p)
        return float((np.mean(np.logaddexp(0, -beta*values))+beta*p/dimension)/entropy)
    trials = [minimize(objective, point, method='Nelder-Mead',
                       bounds=[(-12, 3), (-8, 12)], options={'maxiter': 800})
              for point in ((-2, 1), (-3, 2), (-1, 2))]
    best = min(trials, key=lambda r: r.fun)
    return {'rate': float(best.fun), 'frequency': 1/(1+math.exp(-float(best.x[0]))),
            'normalized_scale': math.exp(float(best.x[1])),
            'optimizer_success': bool(best.success), 'message': str(best.message),
            'scope': 'Numerical fit of a specialized optimistic profile, not a generic bound.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--maximum-dimension', type=int, default=9)
    args = parser.parse_args()
    report = {'status': 'running_incomplete', 'cases': [],
              'scope': 'One all-one specialization per dimension three through nine. '
                       'Specialization can increase determinantal orders; no favorable '
                       'generic profile or stronger score is certified.',
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    begun = time.monotonic()
    for dimension in range(3, args.maximum_dimension+1):
        started = time.monotonic()
        matrix, degree = crossblock(dimension)
        precision = max(128, 1 << (2*degree).bit_length())
        record = {'dimension': dimension, 'degree': degree, 'matrix_order': len(matrix),
                  'other_variable_values': [1]*(dimension-1), 'precision': precision}
        try:
            result = MatrixOracle(matrix, len(matrix), precision=precision).profile()
            orders = list(map(int, result['orders']))
            assert all(x == y for x,y in zip(orders, result['orders']))
            assert all(x+y == 2*degree for x,y in zip(orders, reversed(orders)))
            assert sum(orders) == degree*len(matrix)
            record.update(status='complete_specialized_profile', orders=orders,
                          histogram=dict(sorted(Counter(orders).items())),
                          fit=pressure_fit(orders, degree, dimension))
        except (ArithmeticError, ValueError) as error:
            record.update(status='unresolved_rejected', error=str(error))
        record['seconds'] = time.monotonic()-started
        report['cases'].append(record)
        report['seconds'] = time.monotonic()-begun
        args.output.write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({k:v for k,v in record.items() if k != 'orders'}), flush=True)
    report['status'] = 'complete_registered_specializations'
    args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
