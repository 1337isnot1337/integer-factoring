#!/usr/bin/env python3
"""Exact rational cut meshes: all-field cutoff96, characteristic-two cutoff97."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from time import monotonic
from check_interior_jacobian_budget import (
    binary_rank,witness_matrices,jacobian_rows,verify_chart_linearization,
    geometric_complement,exact_complement_table,weighted_budgets,
    class_parameters,class_bound,fraction_record,
)

CASES97 = [(4,28,37,2026103004),(5,23,33,2026103005),
           (6,21,28,2026103006),(7,13,36,2026103007),
           (8,10,39,2026103008),(9,9,37,2026103009),
           (10,7,42,2026104002),(11,7,37,2026103011),
           (12,6,39,2026103012),(14,4,48,2026103014),
           (17,3,51,2026103017),(18,3,48,2026103018),
           (22,3,38,2026103022)]


def mesh_weights(k,H):
    h = Q(1,H)
    base = 1/(2*k-h)
    result = []
    for j in range(H,k*H):
        s = Q(j,H)
        weight = h/(2*(2+h)) if j == H else 2*h/((2*s-h)*(2*s+h))
        assert weight > 0
        result.append((s,weight))
    assert base+sum(w for _,w in result) == Q(1,2)
    return base,result


def scalar_interpolant(t,k,H):
    base,weights = mesh_weights(k,H)
    return base*(t+1)+sum(w*(min(t,s*s)+min(s*s*t,1)) for s,w in weights)


def interpolation_controls():
    total = 0
    for H in (1,2,4):
        for k in range(2,9):
            for j in range(H,k*H+1):
                s = Q(j,H)
                assert scalar_interpolant(s*s,k,H) == s
                assert scalar_interpolant(1/(s*s),k,H) == 1/s
                total += 2
    return total


def mesh_budgets(d,F,H):
    values = []
    for j in range(H,d*H+1):
        nums = [H*H*F[m]+j*j*m for m in range(2*d)]
        for m in range(2*d-2,-1,-1):
            nums[m] = min(nums[m],nums[m+1])
        values.append(nums)
    return values


def rational_mesh_class_bound(d,k,b,F,G,H):
    a0,_,R,C,e = class_parameters(d,k,b,F)
    base,weights = mesh_weights(k,H)
    value = base*(2*d*d+e)
    for j,(s,w) in enumerate(weights,H):
        ss = s*s
        mu = (ss*k*k-1)/(k**4-1)
        nu = (k*k-ss)/(k**4-1)
        g0 = Q(G[j-H][0],H*H)
        if a0 > ss*b:
            gb = Q(G[j-H][b],H*H) if b < len(F) else ss*b
        else:
            gb = g0
        U = max(gb,min(a0,ss*b)+mu*R+nu*C)
        V = max(g0,b+nu*R+mu*C)
        value += w*(U+V)
    return value


def verify_witnesses97():
    d=97;records=[];improvements={}
    for m,s,q,seed in CASES97:
        A,D = witness_matrices(d,m,s,q,seed)
        rows = jacobian_rows(d,m,s,q,A,D)
        rank,pivots = binary_rank(rows)
        assert rank == m*s*q
        verify_chart_linearization(d,m,s,q,A,D)
        old,new = geometric_complement(d,m),s*q
        assert new > old
        improvements[m] = new
        records.append(dict(d=d,m=m,s=s,q=q,seed=seed,
            old_complement_bound=old,new_complement_bound=new,
            jacobian_rows=m*s*q,jacobian_columns=s*(d-s)+q*(d-q),
            exact_binary_rank=rank,
            omitted_minor_columns=sorted(set(range(s*(d-s)+q*(d-q)))-set(pivots)),
            A_columns_hex=[[hex(x) for x in a] for a in A],
            D_rows_hex=[[hex(x) for x in a] for a in D]))
    return exact_complement_table(d,improvements),records


def dimension_certificate(d,F,floor,H=4):
    G=weighted_budgets(d,F);MG=mesh_budgets(d,F,H)
    assert 577**2-2*408**2 == 1
    assert floor.numerator**408 > d**577*floor.denominator**408
    total,base_closed,integer_closed = 1,1,0
    refined=[];minimum=None
    for k in range(2,d+1):
        for b in range(1,(d*d-1)//((k-1)**2)+1):
            total += 1
            *_,e=class_parameters(d,k,b,F)
            if Q(2*d*d+e,2*k-1) >= floor:
                base_closed += 1
                continue
            old=class_bound(d,k,b,F,G)
            if old >= floor:
                integer_closed += 1
                continue
            new=rational_mesh_class_bound(d,k,b,F,MG,H)
            assert new >= floor,(d,k,b,float(new),float(floor))
            if minimum is None or new < minimum[0]:minimum=(new,k,b)
            refined.append(dict(k=k,b=b,old_bound=fraction_record(old),
                                rational_mesh_bound=fraction_record(new)))
    assert total == 1+sum((d*d-1)//((k-1)**2) for k in range(2,d+1))
    assert total == base_closed+integer_closed+len(refined)
    return dict(d=d,scalar_classes=total,base_term_classes=base_closed,
                integer_mesh_classes=integer_closed,quarter_mesh_classes=len(refined),
                certified_common_floor=fraction_record(floor),
                floor_exceeds_target_by_Pell_upper_exponent=[577,408],
                refined_classes=refined,
                minimum_among_refined_classes=dict(k=minimum[1],b=minimum[2],
                    bound=fraction_record(minimum[0])),
                complement_table=F)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();started=monotonic()
    controls=interpolation_controls()
    F96=exact_complement_table(96,{})
    d96=dimension_certificate(96,F96,Q(636))
    F97,witnesses=verify_witnesses97()
    d97=dimension_certificate(97,F97,Q(6453,10))
    # H=1 is exactly the preceding scalar formula, including selected-piece
    # conditioning; test representative boundary and high-ratio classes.
    for d,F in ((96,F96),(97,F97)):
        G=weighted_budgets(d,F);MG=mesh_budgets(d,F,1)
        for k,b in ((2,1),(3,2),(20,3),(52,3),(d,1)):
            assert rational_mesh_class_bound(d,k,b,F,MG,1)==class_bound(d,k,b,F,G)
    result=dict(status='PASS',scope='All-field finite-seed cutoff96 and characteristic-two finite-seed cutoff97; no factoring exponent improvement',
                all_field_cutoff=96,characteristic_two_cutoff=97,
                small_source_tensor_powers_not_excluded=True,
                interpolation_knot_controls=controls,
                all_field_dimension96=d96,characteristic_two_dimension97=d97,
                binary_dimension97_witnesses=witnesses,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                dependency_sha256=hashlib.sha256(Path(__file__).with_name('check_interior_jacobian_budget.py').read_bytes()).hexdigest(),
                seconds=monotonic()-started)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('all_field_dimension96','characteristic_two_dimension97','binary_dimension97_witnesses')},indent=2))
    for record in (d96,d97):print(json.dumps({k:v for k,v in record.items() if k not in ('refined_classes','complement_table')},indent=2))

if __name__=='__main__':main()
