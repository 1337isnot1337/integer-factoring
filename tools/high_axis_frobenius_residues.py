#!/usr/bin/env python3
"""A fixed Frobenius change of the degree-64 higher exterior jets."""
import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import time

from check_high_axis_plucker import toeplitz_ranks
from high_axis_cayley import MUL, pressure_fit
from high_axis_plucker import elementary, power
from native_matrix_smith_f4 import MatrixOracle


def residue(k):
    e,r = divmod(64,k)
    if not r:
        return e,0,0
    choices = [f for f in (1,2,4,8,16,32,64) if r % f == 0 and math.comb(k,r//f) % 2]
    f = max(choices)
    return e,r//f,f


def crossblock(dimension,values):
    assert 2 <= dimension <= 9 and len(values) == dimension-1 and all(v in (1,2,3) for v in values)
    n,coefficients = 1 << (dimension-1),[]
    for subset in range(n):
        selected = [v for i,v in enumerate(values) if subset >> i & 1]
        e,j,f = residue(len(selected)+1)
        factor = 1
        for v in selected:
            factor = MUL[factor][power(v,e)]
        first = MUL[factor][elementary([power(v,f) for v in selected],j)]
        second = MUL[factor][elementary([power(v,f) for v in selected],j-1)] if f else 0
        coefficients.append(((int(bool(first & 1)) << e) ^ (int(bool(second & 1)) << (e+f)),
                             (int(bool(first & 2)) << e) ^ (int(bool(second & 2)) << (e+f))))
    return [[coefficients[i ^ j] if j & ~i == 0 else (0,0) for j in range(n)] for i in range(n)]


def binary_restriction(matrix):
    n = len(matrix)
    result = [[(0,0)]*(2*n) for _ in range(2*n)]
    for i,row in enumerate(matrix):
        for j,(low,high) in enumerate(row):
            result[i][j],result[i][n+j] = (low,0),(high,0)
            result[n+i][j],result[n+i][n+j] = (high,0),(low ^ high,0)
    return result


def stopping_prefix(orders):
    for L in range(max(orders)+1):
        selected = [v for v in orders if v <= L]
        m,total = len(selected),sum(selected)
        if 3*m > 256 and 3*total <= 64*m:
            return dict(precision=L+1,channels=m,prefix_sum=total)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    patterns = [('all-one',[1]*8),('one-omega',[2]+[1]*7),('three-symbol',[1,2,3,1,2,3,1,2])]
    report = dict(scope='One new higher-jet kernel; specialized profiles are optimistic and cannot prove a generic improvement.',
                  status='running_incomplete',dimension=9,degree=64,
                  residues=[list(residue(k)) for k in range(1,10)],cases=[],
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    begun = time.monotonic()
    for label,values in patterns:
        matrix = crossblock(9,values)
        record = dict(label=label,values=values,precision=256,matrix_order=256)
        try:
            result = MatrixOracle(matrix,256,precision=256).profile()
            orders = list(map(int,result['orders']))
            assert all(a == b for a,b in zip(orders,result['orders']))
            assert all(a+b == 128 for a,b in zip(orders,reversed(orders))) and sum(orders) == 64*256
            record.update(status='complete_specialized_profile',orders=orders,
                          histogram=dict(sorted(Counter(orders).items())),fit=pressure_fit(orders,64,9))
            prefix = stopping_prefix(orders)
            if prefix:
                binary = matrix if all(v[1] == 0 for row in matrix for v in row) else binary_restriction(matrix)
                multiplier = len(binary)//256
                ranks = toeplitz_ranks(binary,prefix['precision'])
                assert all(rank == multiplier*sum(max(L-v,0) for v in orders) for L,rank in enumerate(ranks,1))
                last,previous = ranks[-1]//multiplier,ranks[-2]//multiplier if len(ranks)>1 else 0
                assert prefix['channels'] == last-previous
                assert prefix['prefix_sum'] == prefix['precision']*(last-previous)-last
                record.update(independent_stopping_prefix=prefix,binary_restriction_degree=multiplier,
                              independent_toeplitz_ranks=ranks)
        except (ArithmeticError,ValueError) as error:
            record.update(status='unresolved_rejected',error=str(error))
        report['cases'].append(record)
        report['seconds'] = time.monotonic()-begun
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in record.items() if k not in ('orders','independent_toeplitz_ranks')}),flush=True)
        if record.get('independent_stopping_prefix'):
            report['status'] = 'stopped_by_independent_all_scale_prefix'
            break
    else:
        report['status'] = 'completed_registered_specializations'
    report['seconds'] = time.monotonic()-begun
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
