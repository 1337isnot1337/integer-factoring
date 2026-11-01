#!/usr/bin/env python3
"""Shared block minors, exact Smith controls, and the all-scale certificate."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from cube_homogeneous_gaussian import (ZERO,ONE,monomial,matrix_product,
    exact_monomial_smith,minor,serializable)
from check_fixed_cube_projection_floor import partition_inequality

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/block-gaussian-two-axis-floor.json'


def transpose(matrix):return [list(row) for row in zip(*matrix)]


def inverse_binary(rows):
    n=len(rows)
    work=[row | (1 << (n+i)) for i,row in enumerate(rows)]
    for i in range(n):
        pivot=next(j for j in range(i,n) if (work[j] >> i)&1)
        work[i],work[pivot]=work[pivot],work[i]
        for j in range(n):
            if j!=i and (work[j] >> i)&1:work[j]^=work[i]
    assert all((row & ((1 << n)-1))==1 << i for i,row in enumerate(work))
    return [row >> n for row in work]


def transpose_binary(rows):
    n=len(rows)
    return [sum(((row >> j)&1) << i for i,row in enumerate(rows)) for j in range(n)]


def model(block):
    """Congruence-mixed copies of a four-mode homogeneous endpoint kernel."""
    x2=monomial(2,0,0);y2=monomial(0,2,0);xy=monomial(1,1,0)
    steps=[ZERO,y2 ^ xy,x2 ^ xy,xy]
    size=4*block
    source=[[steps[p ^ q] if i==j else ZERO
             for q in range(4) for j in range(block)]
            for p in range(4) for i in range(block)]
    upper=[(1 << i) | ((1 << (i+1)) if i+1<block else 0) for i in range(block)]
    lower=transpose_binary(upper)
    changes=[upper,lower,transpose_binary(inverse_binary(lower)),
             transpose_binary(inverse_binary(upper))]
    basis=[[ONE if p==q and ((changes[p][i] >> j)&1) else ZERO
            for q in range(4) for j in range(block)]
           for p in range(4) for i in range(block)]
    result=matrix_product(matrix_product(basis,source),transpose(basis))
    assert len(result)==size
    return result,changes


def block_controls():
    records=[];endpoint_entries=0;smith_controls=0
    for block in (1,2,3,4):
        matrix,changes=model(block)
        groups=[list(range(a*block,(a+1)*block)) for a in range(4)]
        for i,row in enumerate(matrix):
            for j,value in enumerate(row):
                assert value==matrix[j][i]
                if i==j:assert not value
                assert all(sum(e)==2 for e in value)
                assert len(value)%2==int(((i//block) ^ (j//block))==3 and i%block==j%block)
                endpoint_entries+=1
        shared=minor(matrix,groups[0],groups[3])
        assert shared==monomial(block,block,0)
        profiles=[]
        for axis,(first,second) in enumerate(((groups[0]+groups[1],groups[2]+groups[3]),
                                              (groups[0]+groups[2],groups[1]+groups[3]))):
            cross=[[matrix[i][j] for j in second] for i in first]
            exact=exact_monomial_smith(cross,axis)
            assert exact['orders']==[1]*block+[3]*block
            assert sum(exact['orders'][:block])==block
            assert sum(exact['orders'])==4*block
            smith_controls+=1
            profiles.append(dict(axis=axis,orders=exact['orders'],
                exact_transformation_sha256=sha256(json.dumps(serializable(exact),sort_keys=True).encode()).hexdigest(),
                finite_operations=len(exact['operations'])))
        records.append(dict(modes_per_parent=block,degree_cap=2,
            constant_basis_changes=changes,profiles=profiles,shared_minor=shared,
            shared_prefix_sum=2*block,shared_prefix_cap=2*block,
            total_order_sum=8*block,total_order_cap=8*block))
    return dict(records=records,endpoint_matrix_entries=endpoint_entries,
                exact_Smith_transformations=smith_controls)


def report():
    sources=['tools/check_block_gaussian_profile_floor.py',
             'tools/check_fixed_cube_projection_floor.py','tools/cube_homogeneous_gaussian.py']
    return serializable(dict(
        scope='Two-axis common-scale, common unweighted type Gaussian profile floor rho>=3/4 for every equal block multiplicity; no actual-rank or factoring lower bound.',
        controls=block_controls(),partition_certificate=partition_inequality(),
        floor=[3,4],normalization='4*L*H(p); shared L-minor bounds joint first-L prefixes',
        source_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=report()
    if args.write:RECEIPT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:assert result==json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=result['scope'],floor=result['floor'],
        blocks=len(result['controls']['records']),
        endpoint_entries=result['controls']['endpoint_matrix_entries'],
        Smith_transformations=result['controls']['exact_Smith_transformations'])))


if __name__=='__main__':main()
