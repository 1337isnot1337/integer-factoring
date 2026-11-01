#!/usr/bin/env python3
"""Degree-64 Cayley kernel with decomposable quadratic axis jets.

The diagnostic retains only specialized complete orders. A favorable
specialization does not certify the generic multivariate crossing profile.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from high_axis_cayley import MUL, pressure_fit
from native_matrix_smith_f4 import MatrixOracle


def power(value, exponent):
    answer = 1
    for _ in range(exponent % 3):
        answer = MUL[answer][value]
    return answer


def elementary(values, order):
    if order < 0:
        return 0
    coefficients = [1]+[0]*len(values)
    for i,value in enumerate(values):
        for j in range(i+1,0,-1):
            coefficients[j] ^= MUL[value][coefficients[j-1]]
    return coefficients[order] if order < len(coefficients) else 0


def crossblock(dimension, other_values, degree=64):
    if not 2 <= dimension <= 9 or len(other_values) != dimension-1:
        raise ValueError('Use dimensions two through nine and one value per other axis.')
    if any(v not in (1,2,3) for v in other_values):
        raise ValueError('Specialize other variables to nonzero F4 values.')
    if any(math.comb(k,degree % k) % 2 != 1 for k in range(1,dimension+1)):
        raise ValueError('The homogeneous coefficient must have all-one endpoint one.')
    n, coefficients = 1 << (dimension-1), []
    for subset in range(n):
        selected = [v for i,v in enumerate(other_values) if subset >> i & 1]
        k = len(selected)+1
        exponent, remainder = divmod(degree,k)
        multiplier = 1
        for value in selected:
            multiplier = MUL[multiplier][power(value,exponent)]
        first = MUL[multiplier][elementary(selected,remainder)]
        second = MUL[multiplier][elementary(selected,remainder-1)]
        coefficients.append(((int(bool(first & 1)) << exponent) ^ (int(bool(second & 1)) << (exponent+1)),
                             (int(bool(first & 2)) << exponent) ^ (int(bool(second & 2)) << (exponent+1))))
    matrix = [[coefficients[i ^ j] if j & ~i == 0 else (0,0)
               for j in range(n)] for i in range(n)]
    return matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    patterns = [('all-one',[1]*8), ('one-omega',[2]+[1]*7),
                ('three-symbol',[1,2,3,1,2,3,1,2])]
    report = {'status':'running_incomplete','degree':64,'dimension':9,'cases':[],
              'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope':'Exactly three registered specializations. Favorable orders cannot '
                      'establish generic profiles, rigorous pressure or a factoring improvement.'}
    begun = time.monotonic()
    for label,values in patterns:
        started = time.monotonic()
        matrix = crossblock(9,values)
        record = {'label':label,'other_values':values,'matrix_order':256,'precision':256}
        try:
            result = MatrixOracle(matrix,256,precision=256).profile()
            orders = list(map(int,result['orders']))
            assert all(a == b for a,b in zip(orders,result['orders']))
            assert all(a+b == 128 for a,b in zip(orders,reversed(orders)))
            assert sum(orders) == 64*256
            record.update(status='complete_specialized_profile',orders=orders,
                          histogram=dict(sorted(Counter(orders).items())),fit=pressure_fit(orders,64,9),
                          prefix_rows=result['rows'],prefix_columns=result['cols'],
                          determinant_prefix_orders=result['base_prefix'])
        except (ArithmeticError,ValueError) as error:
            record.update(status='unresolved_rejected',error=str(error))
        record['seconds'] = time.monotonic()-started
        report['cases'].append(record)
        report['seconds'] = time.monotonic()-begun
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in record.items() if k not in
                          ('orders','prefix_rows','prefix_columns','determinant_prefix_orders')}),flush=True)
    report['status']='complete_registered_specializations'
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
