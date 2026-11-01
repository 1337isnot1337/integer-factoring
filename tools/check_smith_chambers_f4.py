#!/usr/bin/env python3
"""Independent complete profiles and guarded reuse of regular LDU chambers."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
from random import Random

from native_matrix_smith_f4 import MatrixOracle
from native_smith_chambers_f4 import CachedMatrixOracle, library
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_smith import matrix_product, pack, planes


def exact_comparison_controls():
    rng, records = Random(0xC4A8E2), []
    compare = library(128).smith_chamber_compare
    for i in range(512):
        exponent = rng.randint(-1000,990)
        scale = math.ldexp(rng.uniform(.5,1),exponent)
        base = rng.randrange(-1024,1025)
        ri = scale*rng.randrange(-1024,1025)
        cj = scale*rng.randrange(-1024,1025)
        rp = ri+scale*base
        cp = cj
        if i%3:
            rp = math.nextafter(rp,math.inf if i%2 else -math.inf)
        value = F(scale)*base+F(ri)-F(rp)+F(cj)-F(cp)
        sign = (value>0)-(value<0)
        assert compare(base,ri,cj,rp,cp,scale)==sign
        records.append([base,ri,cj,rp,cp,scale,sign])
    return records


def random_source(rng,n):
    units = []
    for _ in range(2):
        matrix = [[int(i==j) for j in range(n)] for i in range(n)]
        for _ in range(2*n):
            i,j = rng.sample(range(n),2)
            coefficient = rng.randrange(1,64)
            from weighted_smith_oracle_f4 import polynomial_multiply
            matrix[i] = [a^polynomial_multiply(coefficient,b) for a,b in zip(matrix[i],matrix[j])]
        units.append(matrix)
    diagonal = [[1<<(2*rng.randrange(7)) if i==j else 0 for j in range(n)] for i in range(n)]
    return matrix_product(matrix_product(units[0],diagonal),units[1])


def matrix_controls():
    rng,records,independent_orders,hits = Random(0xF4C4A8),[],0,0
    for trial in range(24):
        n = 2+trial%4
        source = random_source(rng,n)
        row_shift,col_shift = [[rng.randrange(-2,3) for _ in range(n)] for _ in range(2)]
        shifts = [[a+b for b in col_shift] for a in row_shift]
        clearing = -min(min(row) for row in shifts)
        explicit = [[source[i][j]<<(2*(shifts[i][j]+clearing)) for j in range(n)] for i in range(n)]
        matrix = [[planes(value) for value in row] for row in source]
        cached = CachedMatrixOracle(matrix,n,shifts=shifts,precision=128)
        legacy = MatrixOracle(matrix,n,shifts=shifts,precision=128)
        queries = []
        for place in (False,True):
            exact = WeightedSmithOracle(explicit,infinity=place)
            for query in range(8):
                rw,cw = [[rng.randrange(-12,13)/2 for _ in range(n)] for _ in range(2)]
                scale = F(rng.randrange(1,5),2)
                expected,_ = exact.profile(list(map(F,rw)),list(map(F,cw)),scale)
                correction = clearing*scale*(1 if place else -1)
                expected = [value+correction for value in expected]
                actual = cached.profile(rw,cw,float(scale),place)
                old = legacy.profile(rw,cw,float(scale),place)
                assert actual['orders']==old['orders']==list(map(float,expected))
                repeated = cached.profile(rw,cw,float(scale),place)
                assert repeated['orders']==actual['orders']
                assert repeated['chamber']['hit']
                independent_orders += n
                hits += 1
                queries.append({'infinity':place,'rows':rw,'columns':cw,'scale':str(scale),
                                'orders':list(map(str,expected)),'cache':actual['chamber']})
        records.append({'matrix':matrix,'shifts':shifts,'queries':queries,'stats':cached.cache_stats})
        cached.close()
    return records,independent_orders,hits


def rank_and_window_controls():
    diagonal = [[(1,0),(0,0),(0,0)],[(0,0),(2,0),(0,0)],[(0,0),(0,0),(16,0)]]
    cached = CachedMatrixOracle(diagonal,3,precision=128)
    first = cached.profile()
    changed = cached.profile([5,0,-8],[0,0,0])
    assert first['orders']==[0.,1.,4.] and changed['orders']==[-4.,1.,5.]
    assert changed['chamber']['hit'] and changed['rows']==[2,1,0]
    outside = CachedMatrixOracle([[(1,0),(2,0)],[(4,0),(9,0)]],2,precision=128)
    outside.profile()
    shifted = outside.profile([5,0],[0,0])
    reference = MatrixOracle([[(1,0),(2,0)],[(4,0),(9,0)]],2,precision=128).profile([5,0],[0,0])
    assert shifted['orders']==reference['orders']
    assert not shifted['chamber']['hit']
    singular = [[(1,0),(0,0),(0,0)],[(0,0),(4,0),(0,0)],[(0,0),(0,0),(0,0)]]
    rank_two = CachedMatrixOracle(singular,2,precision=128)
    assert rank_two.profile()['orders']==[0.,2.]
    assert rank_two.profile([4,-2,7],[0,0,0])['orders']==[0.,4.]
    matrix = [[(1|(1<<100),0),(1,0)],[(1,0),(1|(1<<100),0)]]
    too_short = CachedMatrixOracle(matrix,2,precision=128)
    try:
        too_short.profile()
    except ArithmeticError:
        unresolved = True
    else:
        raise AssertionError('A missing high pivot was silently accepted.')
    assert too_short.cache_stats['fallbacks']==1
    enough = CachedMatrixOracle(matrix,2,precision=256)
    assert enough.profile()['orders']==[0.,200.]
    assert enough.profile()['chamber']['hit']
    return {'diagonal_reordering':changed,'outside_cone_recomputed':shifted,
            'promised_rank_two':rank_two.profile(),'short_window_rejected':unresolved,
            'larger_window':enough.profile()}


def real_weight_controls():
    rng,records,comparisons = Random(0xBEA17),[],0
    for trial in range(12):
        n = 2+trial%3
        source = random_source(rng,n)
        matrix = [[planes(v) for v in row] for row in source]
        cached = CachedMatrixOracle(matrix,n,precision=256)
        row,col = [[rng.uniform(-2,2) for _ in range(n)] for _ in range(2)]
        scale = rng.uniform(.1,1.2)
        for place in (False,True):
            oracle = WeightedSmithOracle(source,infinity=place)
            for query in range(6):
                rw = [v+rng.uniform(-.001,.001) for v in row]
                cw = [v+rng.uniform(-.001,.001) for v in col]
                actual = cached.profile(rw,cw,scale,place)
                exact,_ = oracle.profile(list(map(F,rw)),list(map(F,cw)),F(scale))
                errors = [abs(v-float(w)) for v,w in zip(actual['orders'],exact)]
                assert max(errors)<1e-9
                comparisons += n
                records.append({'trial':trial,'infinity':place,'rows':rw,'columns':cw,'scale':scale,
                                'orders':actual['orders'],'error':max(errors),'cache':actual['chamber']})
        cached.close()
    return records,comparisons


def long_determinant_controls():
    rng,records,compared = Random(0x10C4A8),[],0
    n = 8
    units = []
    for _ in range(2):
        unit = [[int(i==j) for j in range(n)] for i in range(n)]
        from weighted_smith_oracle_f4 import polynomial_multiply
        for _ in range(3*n):
            i,j = rng.sample(range(n),2)
            coefficient = rng.randrange(1,4)
            unit[i] = [a^polynomial_multiply(coefficient,b) for a,b in zip(unit[i],unit[j])]
        units.append(unit)
    diagonal = [[1<<(2*15*i) if i==j else 0 for j in range(n)] for i in range(n)]
    source = matrix_product(matrix_product(units[0],diagonal),units[1])
    matrix = [[planes(v) for v in row] for row in source]
    cached = CachedMatrixOracle(matrix,n,precision=256)
    for place in (False,True):
        independent = WeightedSmithOracle(source,infinity=place)
        for trial in range(4):
            rw,cw = [[rng.randrange(-5,6) for _ in range(n)] for _ in range(2)]
            expected,_ = independent.profile(rw,cw)
            actual = cached.profile(rw,cw,infinity=place)
            assert actual['orders']==list(expected)
            assert cached.profile(rw,cw,infinity=place)['chamber']['hit']
            compared += n
            records.append({'infinity':place,'rows':rw,'columns':cw,'orders':list(expected),
                            'cache':actual['chamber']})
    assert sum(15*i for i in range(n))==420>256
    return {'source':matrix,'promised_orders_at_zero_gauge':list(range(0,106,15)),
            'determinant_order':420,'series_precision':256,'queries':records,
            'stats':cached.cache_stats},compared


def forced_prefix_reuse():
    n = 8
    matrix = [[(0,0) for _ in range(n)] for _ in range(n)]
    for i in range(6):
        matrix[i][i] = (1<<i,0)
    matrix[6][6],matrix[6][7] = (1,0),(2,0)
    matrix[7][6],matrix[7][7] = (4,0),(9,0)
    cached = CachedMatrixOracle(matrix,n,precision=128)
    cached.profile()
    rows = [10,0,0,0,0,0,5,0]
    result = cached.profile(rows,[0]*n)
    exact,_ = WeightedSmithOracle([[pack(v) for v in row] for row in matrix]).profile(rows,[0]*n)
    assert result['orders']==list(exact)==[0,1,2,3,4,5,5,10]
    assert result['chamber']['resumed_prefix_pivots']==2
    assert result['rows'][-1]==0  # A reused pivot need not stay among the smallest orders.
    return {'matrix':matrix,'rows':rows,'result':result,'independent_orders':list(exact),
            'stats':cached.cache_stats}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    signs = exact_comparison_controls()
    matrices,orders,hits = matrix_controls()
    controls = rank_and_window_controls()
    real,real_orders = real_weight_controls()
    long_det,long_orders = long_determinant_controls()
    reused = forced_prefix_reuse()
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_smith_chambers_f4.py','tools/native_smith_chambers_f4.py',
             'tools/native_smith_chambers_f4.cpp','tools/native_matrix_smith_f4.py',
             'tools/native_matrix_smith_f4.cpp','tools/native_rectangular_smith_f4.cpp',
             'tools/weighted_smith_oracle_f4.py','tools/windowed_f4_smith.py']
    result = {'status':'pass','exact_binary_rational_sign_comparisons':len(signs),
              'independent_integer_or_half_integer_orders':orders,'forced_repeat_hits':hits,
              'independent_real_binary_rational_orders':real_orders,'sign_controls':signs,
              'matrix_controls':matrices,'boundary_rank_window_controls':controls,
              'real_weight_controls':real,'source_sha256':{n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'large_determinant_order_control':long_det,
              'additional_large_order_comparisons':long_orders,
              'forced_valid_prefix_resumption':reused,
              'scope':'Adaptive regular LDU frame reuse with explicit partial-series knowledge. '
              'Every hit checks its cone; unresolved construction uses the legacy full-profile '
              'fallback. This is not an asymptotic crossing or pressure certificate.'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','exact_binary_rational_sign_comparisons',
        'independent_integer_or_half_integer_orders','forced_repeat_hits','independent_real_binary_rational_orders')}))


if __name__ == '__main__':
    main()
