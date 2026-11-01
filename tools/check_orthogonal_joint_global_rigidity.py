#!/usr/bin/env python3
"""Exact algebra controls for global orthogonal joint-energy rigidity."""
from itertools import permutations, product
from pathlib import Path
import json
from joint_energy_conjugations import model, binary_product, kernel, Space, IDENTITY

# Ordinary polynomial arithmetic over F2, not Boolean arithmetic.
VARS=8
ZERO=frozenset()
ONE=frozenset({(0,)*VARS})
def var(i):
    v=[0]*VARS;v[i]=1;return frozenset({tuple(v)})
def add(*args):
    r=ZERO
    for a in args:r=r^a
    return r
def mul(a,b):
    r=set()
    for x in a:
        for y in b:
            z=tuple(i+j for i,j in zip(x,y))
            if z in r:r.remove(z)
            else:r.add(z)
    return frozenset(r)
def power(a,n):
    r=ONE
    for _ in range(n):r=mul(r,a)
    return r
def subst(a,values):
    r=ZERO
    for x in a:
        term=ONE
        for v,n in zip(values,x):term=mul(term,power(v,n))
        r=add(r,term)
    return r
def mm(A,B):return [[add(*(mul(a,b) for a,b in zip(row,col))) for col in zip(*B)] for row in A]
def tr(A):return list(map(list,zip(*A)))
def det(A):
    r=ZERO
    for perm in permutations(range(len(A))):
        term=ONE
        for i,j in enumerate(perm):term=mul(term,A[i][j])
        r=add(r,term)
    return r
def char(A,t):return det([[add(x,t if i==j else ZERO) for j,x in enumerate(row)]for i,row in enumerate(A)])
def remainder(a,rules):
    a=set(a)
    while True:
        found=False
        for mon in sorted(a,reverse=True):
            for lead,tail in rules:
                if all(x>=y for x,y in zip(mon,lead)):
                    a.remove(mon);factor=frozenset({tuple(x-y for x,y in zip(mon,lead))})
                    a.symmetric_difference_update(mul(factor,tail));found=True;break
            if found:break
        if not found:return frozenset(a)
def transpose_bits(a):return sum(((a>>(25*i+j))&1)<<(25*j+i) for i in range(25) for j in range(25))

def main():
    D=model();C=D['weighted'][-8];L=D['leading'].basis
    assert C==transpose_bits(C)
    assert binary_product(binary_product(C,C),binary_product(C,C))==IDENTITY
    order=sorted(range(5),key=lambda i:D['half'][i]);h=[D['half'][i] for i in order]
    assert h==[0,1,3,4,5]
    X=[[1,0,0,0,0],[0,1,1,0,1],[0,1,0,0,1],[0,0,0,1,0],[0,1,1,0,0]]
    for x,y,u,v in product(range(5),repeat=4):
        assert ((C>>(25*(order[x]+5*order[y])+order[v]+5*order[u]))&1)==X[x][u]*X[v][y]
    mixed=[]
    for e in range(11):
        vals=[X[x][y] for x,y in product(range(5),repeat=2) if h[x]+h[y]==e]
        if len(set(vals))>1:mixed.append(e)
    assert mixed==[4,6]
    assert all(a==transpose_bits(a) for a in L)
    assert all(binary_product(a,b)==binary_product(b,a) for a,b in product(L,repeat=2))
    squares=[binary_product(a,a) for a in L]
    K=kernel([sum(((a>>i)&1)<<j for j,a in enumerate(squares)) for i in range(625)],len(L))
    assert K==[6,24,64,384]
    radical=[D['leading'].matrix((v,0))[0] for v in K]
    assert all(binary_product(a,b)==0 for a,b in product(radical,repeat=2))
    B=D['weighted'][-7]
    comm=[binary_product(a,B)^binary_product(B,a) for a in radical]
    assert len(Space(comm).basis)==2 and comm[0]==comm[3] and comm[1]==comm[2]
    assert all(binary_product(a,b)==binary_product(b,a) for a in [radical[0]^radical[3],radical[1]^radical[2]] for b in D['weighted'].values())
    # Symbolic checks for the global step before restricting to L.
    r,s,u,v,p,q,w,t=[var(i) for i in range(VARS)]
    A=[[ONE,ZERO,ZERO,r,ZERO],[ZERO,ONE,u,ZERO,p],[ZERO,v,w,ZERO,ONE],[s,ZERO,ZERO,ONE,ZERO],[ZERO,q,ONE,ZERO,ZERO]]
    even=[[ONE,r],[s,ONE]];E=mm(even,tr(even));E2=mm(E,E)
    assert E2[0][1]==power(add(r,s),3)
    same=[r,r,u,v,p,q,w,t]
    assert subst(add(E2[0][0],ONE),same)==power(r,4)
    simplified=[ZERO,ZERO,u,u,p,q,add(p,q),t]
    def block(e):
        ix=[(x,y) for y in range(5) for x in range(5) if h[x]+h[y]==e]
        return [[mul(A[x][b],A[a][y])for a,b in ix]for x,y in ix]
    c4=subst(char(block(4),t),simplified)
    # Source charpoly is (t+1)^4. Equality gives (u^2+w)^2=1.
    w0=add(p,q)
    assert c4==mul(add(power(t,2),ONE),add(power(t,2),power(u,4),power(w0,2)))
    c6=subst(char(block(6),t),simplified)
    assert c6==add(power(t,3),mul(t,power(add(mul(p,q),power(w0,2)),2)),power(mul(w0,add(mul(p,q),u)),2))
    B3=[[ONE,u,p],[u,w,ONE],[q,ONE,ZERO]];BB=mm(B3,tr(B3));BB2=mm(BB,BB)
    e=subst(add(BB2[2][2],ONE),[r,s,u,v,p,q,add(p,q),t])
    assert e==add(power(p,2),power(q,4),mul(power(q,2),power(u,2)),power(u,2))
    # In the nontrivial cubic-root case: u=w, p=w+q,
    # q^2=w*q+w and w^2=w+1. This entry becomes 1+w*q.
    expr=subst(e,[r,s,w,v,add(w,q),q,w,t])
    rules=[(next(iter(power(q,2))),add(mul(w,q),w)),(next(iter(power(w,2))),add(w,ONE))]
    assert remainder(expr,rules)==add(ONE,mul(w,q))
    # Then q=w^2=w+1 contradicts q^2+w*q+w=0.
    contradiction=subst(add(power(q,2),mul(w,q),w),[r,s,u,v,p,add(w,ONE),w,t])
    assert remainder(contradiction,rules)==ONE
    for pp,qq in [(ZERO,ONE),(ONE,ZERO)]:
        values=[r,s,ZERO,v,pp,qq,ONE,t]
        assert any(subst(add(BB2[i][j],ONE if i==j else ZERO),values) for i,j in product(range(3),repeat=2))
    result={'completed':True,'source_leading_factor_entries':625,'leading_centralizer_dimension':10,'square_zero_dimension':4,'zero_radical_products':16,'tensor_stabilizer_dimension':2,'effective_parameters':2,'mixed_diagonal_energy_blocks':mixed,'ordinary_polynomial_identities_checked':True,'scope':'Global constant orthogonal energy-preserving gauges retaining leading realignment rank one. No all-field Smith-profile dominance claim.'}
    Path('evidence/orthogonal-joint-global-rigidity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
