#!/usr/bin/env python3
"""Retain all extreme states using an invertible nonbinary quartic residual."""
from __future__ import annotations
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,combinations_with_replacement,permutations,product
import json
from pathlib import Path
from check_rectangular_f4_repair import (ZERO,add,scale,multiply,times_binary,binary_mul,binary_power,divide_binary,c_order,ramify_binary,evaluate_one,fmul,finv,one_particle,matrix_multiply)

CORE=(1,0,2,3)
HEIGHT=(0,0,1,1)

def residual(gamma):
    if gamma not in range(4):raise ValueError('gamma must be an F4 coefficient.')
    coefficients={0:1}
    for monomial,coefficient in [(m,1) for m in (65,130,20,40,96,6)]+[(153,gamma)]:
        old=dict(coefficients)
        for mask,value in old.items():
            if not mask&monomial:
                target=mask|monomial;coefficients[target]=coefficients.get(target,0)^fmul(value,coefficient)
    return [[coefficients.get(CORE[x]|(CORE[y]<<2)|(CORE[u]<<4)|(CORE[v]<<6),0)
             for u in range(4) for x in range(4)] for y in range(4) for v in range(4)]

def inverse(matrix):
    n=len(matrix);a=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(matrix)];det=1;rank=0
    for column in range(n):
        pivot=next((i for i in range(rank,n) if a[i][column]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank];det=fmul(det,a[rank][column]);factor=finv(a[rank][column]);a[rank]=[fmul(factor,v) for v in a[rank]]
        for i in range(n):
            if i!=rank and a[i][column]:
                factor=a[i][column];a[i]=[v^fmul(factor,w) for v,w in zip(a[i],a[rank])]
        rank+=1
    if rank<n:return None,rank,0
    result=[row[n:] for row in a];identity=[[int(i==j) for j in range(n)] for i in range(n)]
    assert matrix_multiply(result,matrix)==identity and matrix_multiply(matrix,result)==identity
    return result,rank,det

def partial_repair(gamma):
    """Block-lower correction leaving a strictly height-lowering endpoint."""
    if gamma==1:raise ValueError('The full residual is singular at gamma=1.')
    diagonal=finv(1^gamma);off_diagonal=fmul(gamma,diagonal)
    C=[[int(i==j) for j in range(16)] for i in range(16)]
    for i in (10,11,14,15):C[i][i]=diagonal
    for i,j in ((10,0),(10,7),(10,13),(11,12),(14,3)):C[i][j]=off_diagonal
    U=matrix_multiply(C,residual(gamma))
    assert all(value==int(i==j) or HEIGHT[i%4]+HEIGHT[i//4]<HEIGHT[j%4]+HEIGHT[j//4]
               for i,row in enumerate(U) for j,value in enumerate(row))
    return C,U

def symbolic_partial_identity():
    """Check delta*C*F over F2[gamma], rather than only its F4 points."""
    F0=residual(0);F1=residual(1)
    F=[[a|((a^b)<<1) for a,b in zip(row,other)] for row,other in zip(F0,F1)]
    correction=[[3*int(i==j) for j in range(16)] for i in range(16)]
    for i in (10,11,14,15):correction[i][i]=1
    for i,j in ((10,0),(10,7),(10,13),(11,12),(14,3)):correction[i][j]=2
    product_matrix=[[0]*16 for _ in range(16)]
    for i,j,k in product(range(16),repeat=3):
        product_matrix[i][j]^=binary_mul(correction[i][k],F[k][j])
    for i,row in enumerate(product_matrix):
        for j,value in enumerate(row):
            if i==j:assert value==3
            elif HEIGHT[i%4]+HEIGHT[i//4]>=HEIGHT[j%4]+HEIGHT[j//4]:assert value==0
    return {'ring':'F2[gamma]','matrix_entries_checked':256,
            'cleared_endpoint_matrix':product_matrix}

def build(indices=(0,1),ambient=2,gamma=2,check_types=True,independent_exterior=False,repair_kind='partial'):
    indices=tuple(indices);matrix=one_particle(indices,ambient);n=len(indices);q=1<<n;bulk=1<<(n-2);top=1<<(n-1)
    F=residual(gamma);C,rank,det=inverse(F)
    if C is None:raise ValueError('This quartic residual is singular; full-alphabet inversion is unavailable.')
    if repair_kind=='partial':C,endpoint_target=partial_repair(gamma)
    elif repair_kind=='inverse':endpoint_target=[[int(i==j) for j in range(16)] for i in range(16)]
    else:raise ValueError('Use partial or inverse repair.')
    masks=[((a%bulk)<<1)|(1,0,top,top|1)[a//bulk] for a in range(q)]
    physical_to_logical={mask:i for i,mask in enumerate(masks)}
    def position(x,y,u,v):return x+q*(y+q*(u+q*v))
    @lru_cache(None)
    def determinant_minor(rows,columns):
        if rows.bit_count()!=columns.bit_count():return 0
        if not rows:return 1
        bit=rows&-rows;i=bit.bit_length()-1;answer=0;remaining=columns
        while remaining:
            column=remaining&-remaining;remaining^=column;j=column.bit_length()-1
            if matrix[i][j]:answer^=binary_mul(matrix[i][j],determinant_minor(rows^bit,columns^column))
        return answer
    raw=[ZERO]*(q**4)
    for x,y,u,v in product(range(q),repeat=4):
        xx,yy,uu,vv=(masks[a] for a in (x,y,u,v))
        value=(binary_mul(5,determinant_minor(vv|(yy<<n),uu|(xx<<n))),0)
        if xx&1 and uu&1 and yy&top and vv&top:
            quartic=determinant_minor((vv^top)|((yy^top)<<n),(uu^1)|((xx^1)<<n))
            value=add(value,scale((quartic,0),gamma))
        raw[position(x,y,u,v)]=value
    exterior_checks=0
    if independent_exterior:
        support={0:1}
        for i,row in enumerate(matrix):
            output=3*n+i if i<n else i
            for j,coefficient in enumerate(row):
                if not coefficient:continue
                input_index=2*n+j if j<n else j-n;monomial=(1<<output)|(1<<input_index)
                for mask,value in list(support.items()):
                    if mask&monomial:continue
                    target=mask|monomial;support[target]=support.get(target,0)^binary_mul(value,coefficient)
                    if not support[target]:del support[target]
        quartic=1|(1<<(2*n-1))|(1<<(2*n))|(1<<(4*n-1))
        exterior={mask:(binary_mul(5,value),0) for mask,value in support.items()}
        for mask,value in support.items():
            if not mask&quartic:exterior[mask|quartic]=add(exterior.get(mask|quartic,ZERO),scale((value,0),gamma))
        for x,y,u,v in product(range(q),repeat=4):
            mask=masks[x]|(masks[y]<<n)|(masks[u]<<(2*n))|(masks[v]<<(3*n))
            assert raw[position(x,y,u,v)]==exterior.get(mask,ZERO);exterior_checks+=1
    values=[ZERO]*len(raw)
    gates=[]
    for row in range(16):
        h=HEIGHT[row%4]+HEIGHT[row//4]
        gates.append([(column%4,column//4,coefficient,binary_power(3,2+HEIGHT[column%4]+HEIGHT[column//4]-h)) for column,coefficient in enumerate(C[row]) if coefficient])
    for x,y,u,v in product(range(q),repeat=4):
        value=ZERO
        for vs,ys,coefficient,factor in gates[v//bulk+4*(y//bulk)]:
            value=add(value,scale(times_binary(raw[position(x,ys*bulk+y%bulk,u,vs*bulk+v%bulk)],factor),coefficient))
        values[position(x,y,u,v)]=value
    cancellation=min(4,min(c_order(v) for v in values if v!=ZERO));denominator_degree=4-cancellation
    if cancellation:values=[divide_binary(v,binary_power(3,cancellation)) for v in values]
    heights=[sum(indices[j] for j in range(n) if mask>>j&1) for mask in masks]
    endpoint_checks=0
    for x,y,u,v in product(range(q),repeat=4):
        polynomial=tuple(ramify_binary(plane,indices[-1]) for plane in values[position(x,y,u,v)])
        exponent=heights[y]+heights[v]-heights[x]-heights[u]-denominator_degree*indices[-1]
        if exponent<0:
            assert not ((polynomial[0]|polynomial[1])&((1<<(-exponent))-1))
            polynomial=tuple(plane>>(-exponent) for plane in polynomial)
        else:polynomial=tuple(plane<<exponent for plane in polynomial)
        expected=endpoint_target[v//bulk+4*(y//bulk)][x//bulk+4*(u//bulk)] if x%bulk==v%bulk and y%bulk==u%bulk else 0
        assert ((polynomial[0]&1)|((polynomial[1]&1)<<1))==expected
        endpoint_checks+=1
    typed_checks=0
    if check_types:
        types=[sorted(set(permutations(pair))) for pair in combinations_with_replacement(range(q),2)]
        divisor=binary_power(3,2*denominator_degree)
        for X,Y in product(types,repeat=2):
            for x,v,y,u in product(X,X,Y,Y):
                polynomial=divide_binary(multiply(values[position(x[0],y[0],u[0],v[0])],values[position(x[1],y[1],u[1],v[1])]),divisor)
                assert evaluate_one(polynomial)==int(x==v and y==u);typed_checks+=1
    vacuum=physical_to_logical[0]
    vacuum_value=values[position(vacuum,vacuum,vacuum,vacuum)]
    assert vacuum_value==(binary_power(3,denominator_degree),0)
    def tensor(rows,columns):
        vv,yy=rows&((1<<n)-1),rows>>n;uu,xx=columns&((1<<n)-1),columns>>n
        return values[position(physical_to_logical[xx],physical_to_logical[yy],physical_to_logical[uu],physical_to_logical[vv])]
    witness=None;wick_checks=0
    for rows,columns in product(combinations(range(2*n),2),repeat=2):
        r0,r1=rows;c0,c1=columns
        left=multiply(tensor((1<<r0)|(1<<r1),(1<<c0)|(1<<c1)),vacuum_value)
        right=add(multiply(tensor(1<<r0,1<<c0),tensor(1<<r1,1<<c1)),multiply(tensor(1<<r0,1<<c1),tensor(1<<r1,1<<c0)))
        wick_checks+=1
        if left!=right and witness is None:witness={'rows':rows,'columns':columns,'nonzero_difference':[hex(v) for v in add(left,right)]}
    return {'field':4,'field_polynomial':'omega^2+omega+1','indices':list(indices),'ambient':ambient,'gamma':gamma,'dimensions':[q]*4,'q_x':q,'q_u':q,'bulk_size':bulk,'repair_kind':'full_alphabet_quartic_'+repair_kind,'denominator_degree':denominator_degree,'cancelled_c_order':cancellation,'strict_endpoint':False,'typed_endpoint':'independent matching full-letter types on X/V and Y/U','residual_rank':rank,'residual_determinant':det,'residual_matrix':F,'repair_matrix':C,'logical_heights_x_v':heights,'logical_heights_u_y':heights,'degree':max(max(v).bit_length()-1 for v in values),'ramified_endpoint_checks':endpoint_checks,'two_copy_type_checks':typed_checks,'independent_exterior_checks':exterior_checks,'wick_minor_checks':wick_checks,'non_gaussian_witness':witness,'coefficient_sha256':sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest(),'values':values}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);parser.add_argument('--seed-directory',type=Path);parser.add_argument('--full',action='store_true');parser.add_argument('--repair',choices=('partial','inverse','both'),default='both');args=parser.parse_args()
    F0=residual(0);inverse0=inverse(F0)[0]
    quartic=[[a^b for a,b in zip(row,base)] for row,base in zip(residual(1),F0)]
    projector=matrix_multiply(inverse0,quartic)
    assert matrix_multiply(projector,projector)==projector and inverse(projector)[1]==4
    residuals=[]
    for gamma in range(4):
        _,rank,det=inverse(residual(gamma));expected=1
        for _ in range(4):expected=fmul(expected,1^gamma)
        assert det==expected and rank==(12 if gamma==1 else 16)
        residuals.append({'gamma':gamma,'rank':rank,'determinant':det})
    records=[]
    modes=[((0,1),2),((0,1,2),4),((0,2,5),8)]+([((0,1,2,3),4)] if args.full else [])
    repairs=('partial','inverse') if args.repair=='both' else (args.repair,)
    for (indices,ambient),repair_kind,gamma in product(modes,repairs,(0,2,3)):
            seed=build(indices,ambient,gamma,check_types=True,independent_exterior=True,repair_kind=repair_kind)
            assert (seed['non_gaussian_witness'] is None)==(gamma==0)
            if args.seed_directory:
                args.seed_directory.mkdir(parents=True,exist_ok=True);mode_name='-'.join(map(str,indices));(args.seed_directory/f'full-quartic-{repair_kind}-gamma{gamma}-S{mode_name}.json').write_text(json.dumps(seed,separators=(',',':'))+'\n')
            records.append({k:v for k,v in seed.items() if k not in ('values','residual_matrix','repair_matrix')})
            print(json.dumps({'indices':indices,'gamma':gamma,'repair':repair_kind,'denominator_degree':seed['denominator_degree'],'degree':seed['degree'],'checks':seed['two_copy_type_checks'],'non_gaussian':seed['non_gaussian_witness'] is not None}),flush=True)
    result={'status':'exact full-alphabet quartic repairs; no pressure bound','symbolic_partial_identity':symbolic_partial_identity(),'residual_cases':residuals,'residual_projector_binary_row_masks':[sum(v<<j for j,v in enumerate(row)) for row in projector],'projector_rank':4,'cases':records,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
