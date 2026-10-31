#!/usr/bin/env python3
"""Exact characteristic-two geometric rank-budget certificate at dimension 96.

The new witnesses are literal binary Jacobians. No generic-rank assumption or
floating-point threshold is used. The preceding all-field geometric table is
recomputed from its integer allocation formula.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from random import Random
from time import monotonic

# (number of maps, input dimension, output codimension, fixed witness seed)
CASES = [(3,34,43,2026103003), (5,25,30,2026103005),
         (7,12,38,2026103007), (8,10,38,2026103008),
         (9,8,41,2026103009), (10,7,41,2026103010),
         (12,6,38,2026103012), (14,4,47,2026103014),
         (17,3,50,2026103017), (18,3,47,2026103018)]


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            position = row.bit_length()-1
            if position not in pivots:
                pivots[position] = row
                break
            row ^= pivots[position]
    return len(pivots), sorted(pivots)


def witness_matrices(d, m, s, q, seed):
    """A_i columns and D_i rows as literal binary masks."""
    rng = Random(seed)
    A = [[rng.getrandbits(d-q) for _ in range(s)] for _ in range(m)]
    D = [[rng.getrandbits(d-s) for _ in range(q)] for _ in range(m)]
    return A, D


def jacobian_rows(d, m, s, q, A, D):
    """Rows of (X,Y) -> (D_i X + Y A_i)_i, over F2.

    X has size (d-s) by s and precedes Y, of size q by (d-q), in
    row-major input order. Outputs are in (map,row,column) order.
    """
    p, ell = d-s, d-q
    assert len(A) == len(D) == m
    assert all(len(a) == s and all(0 <= x < 1 << ell for x in a) for a in A)
    assert all(len(a) == q and all(0 <= x < 1 << p for x in a) for a in D)
    rows = []
    for i in range(m):
        for u in range(q):
            spread = sum(1 << (k*s) for k in range(p) if D[i][u] >> k & 1)
            for v in range(s):
                rows.append((spread << v) | (A[i][v] << (p*s+u*ell)))
    return rows


def verify_chart_linearization(d, m, s, q, A, D):
    """Independent direct evaluation of J on several literal input masks."""
    p, ell = d-s, d-q
    rows = jacobian_rows(d,m,s,q,A,D)
    rng = Random(30102026+m)
    for _ in range(3):
        X = [rng.getrandbits(s) for _ in range(p)]
        Y = [rng.getrandbits(ell) for _ in range(q)]
        packed = sum(x << (k*s) for k,x in enumerate(X))
        packed |= sum(y << (p*s+u*ell) for u,y in enumerate(Y))
        direct = []
        for i in range(m):
            for u in range(q):
                left = 0
                for k in range(p):
                    if D[i][u] >> k & 1:
                        left ^= X[k]
                direct.extend(((left >> v) & 1) ^ ((Y[u] & A[i][v]).bit_count() & 1)
                              for v in range(s))
        assert direct == [(row & packed).bit_count() & 1 for row in rows]


@lru_cache(None)
def geometric_complement(d, m):
    if m < 2 or m >= 2*d-1:
        return 0
    best = 0
    for a in range(2,m+1):
        b = 2*m-a
        for s in range(1,(2*d-2)//b+1):
            q = min(2*(d-s)//a, (2*d-b*s)//2)
            if q > 0:
                best = max(best,s*q)
    return best


def exact_complement_table(d, improvements):
    values = [d*d,d*(d-1)//2,d*(d-1)//2]
    values.extend(max(geometric_complement(d,m),improvements.get(m,0))
                  for m in range(3,2*d))
    # Zero padding also permits the monotone envelope of any certificate list.
    for m in range(len(values)-2,-1,-1):
        values[m] = max(values[m],values[m+1])
    assert values[-1] == 0
    return values


def weighted_budgets(d, F):
    result = []
    for s in range(d+1):
        vals = [F[m]+s*s*m for m in range(2*d)]
        for m in range(2*d-2,-1,-1):
            vals[m] = min(vals[m],vals[m+1])
        result.append(vals)
    return result


def class_parameters(d,k,b,F):
    n = d*d
    a0, a1 = (k-1)**2*b+1, min(k*k*b,n)
    fb = F[b] if b < len(F) else 0
    R,C = max(n-a1,fb),n-b
    if b == 1:
        R = max(R,n-k*(k+1)//2)
        e = (k*k-3*k+4)//2
    else:
        e = max(0,a0+fb-n)
    return a0,a1,R,C,e


def class_bound(d,k,b,F,G):
    _,_,R,C,e = class_parameters(d,k,b,F)
    value = Fraction(2*d*d+e,2*k-1)
    for s in range(1,k):
        mu = Fraction(s*s*k*k-1,k**4-1)
        nu = Fraction(k*k-s*s,k**4-1)
        gb = G[s][b] if b < len(F) else s*s*b
        U = max(Fraction(gb),s*s*b+mu*R+nu*C)
        V = max(Fraction(G[s][0]),b+nu*R+mu*C)
        weight = Fraction(1,6) if s == 1 else Fraction(2,(2*s-1)*(2*s+1))
        value += weight*(U+V)
    return value


def fraction_record(x):
    return dict(numerator=x.numerator,denominator=x.denominator)


def dimension_certificate(d,F):
    G = weighted_budgets(d,F)
    # This is merely an initial class, not an assumption about the minimizer.
    best = class_bound(d,52,3,F,G)
    selected = (52,3)
    classes,full,basic = 1,0,1
    for k in range(2,d+1):
        for b in range(1,(d*d-1)//((k-1)**2)+1):
            classes += 1
            *_, e = class_parameters(d,k,b,F)
            elementary = Fraction(2*d*d+e,2*k-1)
            if elementary >= best:
                basic += 1
                continue
            value = class_bound(d,k,b,F,G)
            full += 1
            if value < best:
                best,selected = value,(k,b)
    assert classes == 1+sum((d*d-1)//((k-1)**2) for k in range(2,d+1))
    assert basic+full == classes
    assert 577**2-2*408**2 == 1
    assert best.numerator**408 > d**577*best.denominator**408
    return dict(d=d,minimizing_class=dict(k=selected[0],b=selected[1]),
                exact_lower_bound=fraction_record(best),
                exact_comparison=dict(exponent_numerator=577,exponent_denominator=408,
                                     numerator_power_exceeds_d_power_times_denominator_power=True),
                scalar_classes=classes,full_weighted_class_evaluations=full,
                classes_closed_by_nonnegative_base_term=basic,
                complement_table=F,unconditioned_weighted_budgets=[g[0] for g in G])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = monotonic();d=96;records=[];improvements={}
    for m,s,q,seed in CASES:
        A,D = witness_matrices(d,m,s,q,seed)
        rows = jacobian_rows(d,m,s,q,A,D)
        rank,pivots = binary_rank(rows)
        assert rank == m*s*q
        verify_chart_linearization(d,m,s,q,A,D)
        old,new = geometric_complement(d,m),s*q
        assert new > old
        improvements[m] = new
        records.append(dict(d=d,m=m,s=s,q=q,old_complement_bound=old,
                            new_complement_bound=new,seed=seed,
                            jacobian_rows=m*s*q,jacobian_columns=s*(d-s)+q*(d-q),
                            exact_binary_rank=rank,
                            omitted_minor_columns=sorted(set(range(s*(d-s)+q*(d-q)))-set(pivots)),
                            A_columns_hex=[[hex(x) for x in a] for a in A],
                            D_rows_hex=[[hex(x) for x in a] for a in D]))
    F = exact_complement_table(d,improvements)
    proof = dimension_certificate(d,F)
    # Independent old-table control: the old certificate genuinely falls short.
    oldF = exact_complement_table(d,{})
    old = class_bound(d,52,3,oldF,weighted_budgets(d,oldF))
    assert 1393**2-2*985**2 == -1
    assert old.numerator**985 < d**1393*old.denominator**985
    result = dict(status='PASS',
                  scope='Exact finite crossing rank-sum exclusion through 96 labels in characteristic two; no factoring exponent improvement',
                  previous_all_field_cutoff=95,new_characteristic_two_cutoff=96,
                  small_source_tensor_power_constructions_not_excluded=True,
                  witnesses=records,dimension_certificate=proof,
                  previous_table_control=dict(d=d,k=52,b=3,bound=fraction_record(old),
                    below_target_certified_by_lower_Pell_exponent=[1393,985]),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  seconds=monotonic()-started)
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('witnesses','dimension_certificate')},indent=2))
    print(json.dumps({k:v for k,v in proof.items() if k not in ('complement_table','unconditioned_weighted_budgets')},indent=2))

if __name__ == '__main__':
    main()
