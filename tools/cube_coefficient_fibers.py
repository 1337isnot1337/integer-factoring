#!/usr/bin/env python3
"""Exact homogeneous eight-leg coefficient fibers in the fixed Gaussian frames.

All coefficient variables are free subject to the chosen symmetry, parity,
homogeneous entry degree 3*occupation, and literal crossing endpoint.
These are fixed-frame constraints, not an exhaustive Smith-spectrum search.
"""
from __future__ import annotations
import itertools
import argparse
import hashlib
import json
from pathlib import Path
import cube_homogeneous_gaussian as c

def actions(mode='full'):
    out=[]
    for p in (itertools.permutations(range(3)) if mode in ('full','axes') else ((0,1,2),)):
        for t in (range(8) if mode in ('full','translations') else (0,)):
            vertices=[t ^ sum(((v>>i)&1)<<p[i] for i in range(3)) for v in range(8)]
            masks=[sum(((m>>i)&1)<<vertices[i] for i in range(8)) for m in range(256)]
            out.append((p,masks))
    return out

def act(pair, action):
    mask,e=pair; p,masks=action
    f=[0]*3
    for i in range(3): f[p[i]]=e[i]
    return masks[mask],tuple(f)

def basis(mode='full'):
    group=actions(mode); where={}; orbits=[]
    for mask in range(256):
        k=mask.bit_count()
        if k&1: continue
        d=3*k
        for a in range(d+1):
            for b in range(d-a+1):
                pair=(mask,(a,b,d-a-b))
                if pair in where: continue
                orbit=frozenset(act(pair,g) for g in group)
                at=len(orbits); orbits.append(orbit)
                for term in orbit: where[term]=at
    anchors={where[(mask,(mask.bit_count(),)*3)] for mask in range(256) if mask.bit_count()%2==0}
    directions=[]; tensor=[{} for _ in range(256)]
    for i,orbit in enumerate(orbits):
        if i in anchors: continue
        mask,e=min(orbit); k=mask.bit_count()
        masks={m for m,_ in orbit}
        count=len(orbit)//len(masks)
        assert all(sum(m==target for m,_ in orbit)==count for target in masks)
        direction=orbit ^ (orbits[where[(mask,(k,k,k))]] if count&1 else frozenset())
        assert all(sum(m==target for m,_ in direction)%2==0 for target in masks)
        bit=1<<len(directions); directions.append(direction)
        for m,e in direction: tensor[m][e]=tensor[m].get(e,0)^bit
    return directions,tensor,orbits,group

def frames():
    source=c.kernel(); seed=c.gaussian_tensor(source)
    U,_,V,cross=c.gaussian_cut_factors(source,0)
    assert c.matrix_product(U,U)==c.identity(16)
    assert c.matrix_product(V,V)==c.identity(16)
    frame=c.exact_monomial_smith(cross,0)
    L=c.matrix_product(c.compound(frame['left']),U)
    R=c.matrix_product(V,c.compound(frame['right']))
    diagonal=c.matrix_product(c.matrix_product(L,c.axis_matrix(seed,0)),R)
    assert all(not diagonal[i][j] for i in range(16) for j in range(16) if i!=j)
    assert all(len(diagonal[i][i])==1 for i in range(16))
    diag=[next(iter(diagonal[i][i])) for i in range(16)]
    assert all(e[0]>=0 for row in L+R for poly in row for e in poly)
    return seed,L,R,diag,frame

def equations(tensor,L,R,diag,nvars,window):
    cap=[d[0]+window for d in diag]
    table={}; terms=0; first,second=c.axis_parts(0)
    leftcols=[[(i,L[i][a]) for i in range(16) if L[i][a]] for a in range(16)]
    rightrows=[[(j,R[b][j]) for j in range(16) if R[b][j]] for b in range(16)]
    for mask,symbolic in enumerate(tensor):
        if not symbolic: continue
        a=sum(((mask>>position)&1)<<j for j,position in enumerate(first))
        b=sum(((mask>>position)&1)<<j for j,position in enumerate(second))
        for i,lp in leftcols[a]:
            for j,rp in rightrows[b]:
                pol=c.multiply(lp,rp)
                total=3*(i.bit_count()+j.bit_count())
                for e,variables in symbolic.items():
                    for h in pol:
                        x=e[0]+h[0]
                        if x>=cap[i]: continue
                        y=e[1]+h[1]; z=e[2]+h[2]
                        assert x>=0 and x+y+z==total
                        key=i,x,j,y
                        table[key]=table.get(key,0)^variables
                        terms+=1
    for i,e in enumerate(diag):
        key=i,e[0],i,e[1]
        table[key]=table.get(key,0)^(1<<nvars)
    result=[[set() for _ in range(cap[i])] for i in range(16)]
    for (i,x,j,y),row in table.items():
        if row: result[i][x].add(row)
    return [[tuple(sorted(e)) for e in row] for row in result],cap,terms,len(table)

def insert(pivots,row,mask):
    while row&mask:
        p=(row&mask).bit_length()-1
        if p not in pivots:
            pivots[p]=row; return True
        row^=pivots[p]
    return row==0

def solve(eq,targets,nvars):
    pivots={}; mask=(1<<nvars)-1
    for i,target in enumerate(targets):
        for power in range(target):
            for row in eq[i][power]:
                if not insert(pivots,row,mask): return None
    return pivots

def literal(seed,directions,value):
    out=[set(p) for p in seed]
    for i,direction in enumerate(directions):
        if not(value>>i)&1: continue
        for m,e in direction:
            if e in out[m]:out[m].remove(e)
            else:out[m].add(e)
    return tuple(frozenset(p) for p in out)

def validate(seed,L,R,tensor,targets,group):
    for mask,p in enumerate(tensor):
        assert len(p)%2==len(seed[mask])%2
        assert all(min(e)>=0 and sum(e)==3*mask.bit_count() for e in p)
    for g in group:
        for m,p in enumerate(tensor):
            for e in p:
                nm,ne=act((m,e),g)
                assert ne in tensor[nm]
    product=c.matrix_product(c.matrix_product(L,c.axis_matrix(tensor,0)),R)
    assert all(all(e[0]>=targets[i] for e in poly) for i,row in enumerate(product) for poly in row)
    return sum(len(p) for p in tensor),sum(len(p) for row in product for p in row)


def rotate_symbolic(symbolic, axis):
    permutation=list(range(3))
    permutation[0],permutation[axis]=permutation[axis],permutation[0]
    action=next(g for g in actions() if list(g[0])==permutation and g[1][1]==1)
    result=[{} for _ in range(256)]
    for mask,polynomial in enumerate(symbolic):
        for exponent,variables in polynomial.items():
            target,shifted=act((mask,exponent),action)
            result[target][shifted]=variables
    return result


def three_cut_equations(symbolic,left,right,diagonal,variables,window=1):
    result=[]
    for axis in range(3):
        equation,_,_,_=equations(rotate_symbolic(symbolic,axis),left,right,
                                 diagonal,variables,window)
        result.extend(equation)
    return result


def assignment_from_free(pivots,value,variables):
    if value<0 or value>>variables:
        raise ValueError('Assignment outside the coefficient space.')
    if any(value>>p&1 for p in pivots):
        raise ValueError('Free assignment uses a pivot coordinate.')
    for pivot,row in sorted(pivots.items()):
        value|=(((row>>variables)&1)^((row&value).bit_count()&1))<<pivot
    return value


def specialized_vector(tensor):
    result=[]
    for axis in range(3):
        matrix=c.axis_matrix(tensor,axis)
        for parity in range(2):
            labels=[i for i in range(16) if i.bit_count()%2==parity]
            for row in labels:
                for column in labels:
                    value=0
                    for exponent in matrix[row][column]:
                        value^=1<<exponent[axis]
                    result.append(value)
    return result


def floor_four_input(directions,seed,pivots,variables):
    free=[i for i in range(variables) if i not in pivots]
    assignments=[assignment_from_free(pivots,1<<i,variables) for i in free]
    zero=tuple(c.ZERO for _ in range(256))
    vectors=[specialized_vector(literal(zero,directions,value)) for value in assignments]
    image={}
    for vector in vectors:
        value=sum(entry<<(32*i) for i,entry in enumerate(vector))
        while value:
            pivot=value.bit_length()-1
            if pivot not in image:
                image[pivot]=value
                break
            value^=image[pivot]
    return dict(variables=variables,free_coordinates=free,
                free_variables=len(free),specialization_image_dimension=len(image),
                base=specialized_vector(seed),directions=vectors,
                assignments=[hex(value) for value in assignments])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    directions,symbolic,orbits,_=basis('none')
    seed,left,right,diagonal,_=frames()
    equations_all=three_cut_equations(symbolic,left,right,diagonal,len(directions))
    targets=tuple(max(0,exponent[0]-4) for exponent in diagonal)*3
    pivots=solve(equations_all,targets,len(directions))
    assert pivots is not None and not any(row>>len(directions) for row in pivots.values())
    record=floor_four_input(directions,seed,pivots,len(directions))
    record.update(scope='Complete binary coefficient fiber with all fixed Gaussian row floors lowered by four',
                  reference_orders=[exponent[0] for exponent in diagonal],row_targets=targets,
                  source_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in ('cube_coefficient_fibers.py','cube_homogeneous_gaussian.py')})
    args.output.write_text(json.dumps(record,separators=(',',':'))+'\n')
    print(json.dumps({key:record[key] for key in ('variables','free_variables','specialization_image_dimension')}))
