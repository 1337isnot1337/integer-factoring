#!/usr/bin/env python3
"""Independent Toeplitz obstruction for the degree-64 high-axis kernel."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
import time

import numpy as np

from high_axis_plucker import MUL, crossblock
from native_matrix_smith_f4 import MatrixOracle
from regularized_transpose_words import field_rank
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import pack


def toeplitz_ranks(matrix, precision):
    """Rank of multiplication by the F2 matrix modulo x^L, all L<=precision.

    Binary rows are built literally; no polynomial elimination or native
    Smith code enters this rank calculation. Earlier rows embed unchanged
    as L grows, so one exact incremental binary elimination suffices.
    """
    n, coefficients = len(matrix), []
    for row in matrix:
        values = {}
        for j,(low,high) in enumerate(row):
            assert high == 0
            while low:
                d = (low & -low).bit_length()-1
                values[d] = values.get(d,0) ^ (1 << j)
                low &= low-1
        coefficients.append(values)
    pivots, ranks = {}, []
    for t in range(precision):
        for values in coefficients:
            row = 0
            for d,value in values.items():
                if d <= t:
                    row ^= value << (n*(t-d))
            while row:
                j = row.bit_length()-1
                if j not in pivots:
                    pivots[j] = row
                    break
                row ^= pivots[j]
        ranks.append(len(pivots))
    return ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    begun = time.monotonic()
    degree_checks = []
    for k in range(1,10):
        e,r = divmod(64,k)
        assert k*e+r == 64 and math.comb(k,r) % 2
        degree_checks.append({'subset_size':k,'exponent':e,'remainder':r,
                              'endpoint_binomial':math.comb(k,r)})
    matrix = crossblock(9,[1]*8)
    # An independent scalar subset expansion of every nonzero entry.
    entries = 0
    for i,row in enumerate(matrix):
        for j,value in enumerate(row):
            if j & ~i:
                assert value == (0,0)
                continue
            size = (i ^ j).bit_count()
            e,r = divmod(64,size+1)
            first = sum(1 for _ in combinations(range(size),r)) % 2 if r <= size else 0
            second = sum(1 for _ in combinations(range(size),r-1)) % 2 if r > 0 else 0
            assert value == ((first << e) ^ (second << (e+1)),0)
            entries += 1
    ranks = toeplitz_ranks(matrix,26)
    assert ranks[-2:] == [742,836]
    j, prefix = ranks[-1]-ranks[-2], 26*(ranks[-1]-ranks[-2])-ranks[-1]
    assert j == 94 and prefix == 1608
    assert 3*j > 256 and 3*prefix < 64*j
    # Decomposable quadratic jet at a nontrivial specialization.
    values = [2]+[1]*7
    coefficients = [0]*256
    for i,jj in combinations(range(8),2):
        # x^21=1 and x^22=x on F4^*, so this is the literal triple jet.
        coefficients[(1 << i) | (1 << jj)] = values[i] ^ values[jj]
    quadratic = np.array([[coefficients[i ^ jj] if not jj & ~i else 0
                           for jj in range(256)] for i in range(256)],dtype=np.uint8)
    assert field_rank(quadratic) == 64
    records, orders = [], 0
    for dimension in (3,4,5):
        for values in ([1]*(dimension-1), [1+i % 3 for i in range(dimension-1)]):
            small = crossblock(dimension,values)
            packed = [[pack(value) for value in row] for row in small]
            exact,_ = WeightedSmithOracle(packed).profile([0]*len(small),[0]*len(small),1)
            got = MatrixOracle(small,len(small),precision=256).profile()['orders']
            assert list(exact) == list(map(int,got))
            records.append({'dimension':dimension,'other_values':values,'orders':list(exact)})
            orders += len(exact)
    names = ['high_axis_plucker.py','check_high_axis_plucker.py','weighted_smith_oracle_f4.py']
    report = {'status':'pass','degree_and_endpoint_checks':degree_checks,
              'independent_nonzero_pattern_entries':entries,'toeplitz_ranks_precisions_1_through_26':ranks,
              'toeplitz_dimensions':[6400,6656],'prefix_channels':94,'generic_prefix_upper_bound':1608,
              'decomposable_quadratic_multiplication_rank':64,'independent_untruncated_small_orders':orders,
              'small_controls':records,'seconds':time.monotonic()-begun,
              'pressure_floor':'1/3 at every common scale and physical frequency',
              'source_sha256':{name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'scope':'Exact explicit-kernel obstruction from a binary truncated multiplication map '
                      'and a proved generic rank, not a floor for other high-axis kernels or actual '
                      'coefficient ranks. No stronger factoring bound.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','prefix_channels','generic_prefix_upper_bound',
                     'decomposable_quadratic_multiplication_rank','seconds')}))


if __name__ == '__main__':
    main()
