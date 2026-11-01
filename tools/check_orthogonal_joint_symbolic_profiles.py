#!/usr/bin/env python3
"""All-field natural-profile classification via ordinary polynomial minors."""
import argparse
from hashlib import sha256
from itertools import accumulate
import json
from pathlib import Path
import joint_energy_conjugations as joint
import search_f4_recentered_crossings as seed_source
from joint_energy_conjugations import model,binary_product
from sparse_f2_trivariate import Poly,ZERO,ONE,A,B,Z,MASK

def source():
    data=model();L=data['leading'];P=L.matrix((6,0))[0];Q=L.matrix((24,0))[0]
    C=[[ZERO for _ in range(25)]for _ in range(25)]
    for k,V in data['weighted'].items():
        PV=binary_product(P,V);VP=binary_product(V,P)
        QV=binary_product(Q,V);VQ=binary_product(V,Q)
        terms=[(V,ONE),(PV^VP,A),(QV^VQ,B),(binary_product(PV,P),A*A),
               (binary_product(QV,Q),B*B),(binary_product(PV,Q)^binary_product(QV,P),A*B)]
        for x in range(5):
            for y in range(5):
                for u in range(5):
                    for v in range(5):
                        pos=25*(x+5*y)+v+5*u
                        for mask,c in terms:
                            if mask>>pos&1:C[5*x+u][5*y+v]+=c*Z**(k+8)
    return C

def valuation(p):return min(m&MASK for m in p.terms)
def lead(p):
    v=valuation(p);return Poly(frozenset(m-v for m in p.terms if (m&MASK)==v))
def factors(p):
    out=[]
    for f in [A,B,A+B]:
        n=0
        while p!=ONE:
            q,r=p.div(f)
            if r:break
            p=q;n+=1
        out.append(n)
    assert p==ONE, 'A pivot has an unclassified exceptional locus'
    return out

def certificate(C,stratum):
    if stratum=='a0':C=[[p.compose(A,ZERO) for p in row]for row in C]
    elif stratum=='b0':C=[[p.compose(B,ZERO) for p in row]for row in C]
    elif stratum=='equal':C=[[p.compose(B,A) for p in row]for row in C]
    elif stratum=='origin':C=[[p.compose(A,ZERO).compose(B,ZERO) for p in row]for row in C]
    rows=list(range(25));cols=list(range(25));prev=ONE;prevval=0;orders=[];cert=[]
    for step in range(25):
        _,_,i,j=min((valuation(C[i][j]),len(C[i][j]),i,j)for i in rows for j in cols if C[i][j])
        pivot=C[i][j];val=valuation(pivot);orders.append(val-prevval-8)
        powers=factors(lead(pivot))
        # Every retained leading coefficient is nonzero throughout its stratum.
        if stratum=='origin':assert powers==[0,0,0]
        if stratum=='a0':assert powers[0]==powers[2]==0
        if stratum in ('b0','equal'):assert powers[1]==powers[2]==0
        cert.append({'row':i,'column':j,'valuation':val,'leading_factor_powers_a_b_a_plus_b':powers,
                     'terms':len(pivot),'polynomial_sha256':sha256(json.dumps(sorted(pivot.terms)).encode()).hexdigest()})
        rows.remove(i);cols.remove(j)
        if step<24:
            updated={}
            for x in rows:
                for y in cols:
                    numerator=pivot*C[x][y]+C[x][j]*C[i][y]
                    updated[x,y]=numerator.exquo(prev) if numerator else ZERO
            for (x,y),p in updated.items():C[x][y]=p
        prev=pivot;prevval=val
    assert orders==sorted(orders)
    return {'stratum':stratum,'orders':orders,'prefix_minors':cert}

def arithmetic_controls():
    import random
    rng=random.Random(17731)
    for _ in range(80):
        def rand():return Poly(frozenset((rng.randrange(4)<<40)+(rng.randrange(4)<<20)+rng.randrange(8) for _ in range(12)))
        x,y,z=rand(),rand(),rand()
        assert x*(y+z)==x*y+x*z
        assert (x*y).exquo(y)==x
        assert (x+y)**2==x*x+y*y
    assert (A+B)*(A+B)==A*A+B*B
    assert A*A!=A  # No Boolean reduction.

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('evidence/orthogonal-joint-symbolic-profiles.json'))
    args=parser.parse_args();arithmetic_controls();C=source();results=[]
    baseline=[-8,-7,-5,-4,-3,-2,-1,0,0,1,1,2,2,3,3,4,4,4,5,5,6,7,7,8,12]
    for mode in ['origin','a0','b0','equal','generic']:
        r=certificate(C,mode)
        assert all(x<=y for x,y in zip(accumulate(r['orders']),accumulate(baseline)))
        if mode=='origin':assert r['orders']==baseline
        results.append(r);print(json.dumps({'stratum':mode,'orders':r['orders'],'prefix_dominated':True}),flush=True)
    root=Path(__file__).resolve().parent
    paths={n:root/n for n in ['check_orthogonal_joint_symbolic_profiles.py','sparse_f2_trivariate.py']}
    paths['joint_energy_conjugations.py']=Path(joint.__file__)
    paths['search_f4_recentered_crossings.py']=Path(seed_source.__file__)
    out={'completed':True,'ordinary_characteristic_two_polynomials':True,'strata':results,
         'source_sha256':{n:sha256(p.read_bytes()).hexdigest() for n,p in paths.items()},
         'scope':'Every constant orthogonal energy-preserving gauge retaining rank-one leading realignment, over every characteristic-two field, at the stated natural heights. No arbitrary-weight or nonorthogonal exclusion.',
         'factoring_exponent':'1/2; unchanged'}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
