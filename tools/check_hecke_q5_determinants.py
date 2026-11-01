#!/usr/bin/env python3
"""Exact Hecke determinants, irreducible factors, and two-closed-place budgets."""
import itertools,json,random
from hashlib import sha256
from pathlib import Path
import numpy as np
from f8_polynomial_determinant import determinant
from f8_polynomial_factor import add,multiply,trim,factor,irreducible
from three_row_hecke_paths import M
from three_row_hecke_crossing import coefficient_tensor,example_paths
from binary_extension_smith import BinaryExtension

def small_controls():
    count=0
    for a,b,c,d in itertools.product(range(8),repeat=4):
        C=np.array([[[a,b],[c,d]]],dtype=np.uint8)
        assert determinant(C)==trim([M[a][d]^M[b][c]]);count+=1
    rng=random.Random(445619)
    for n in range(1,5):
        for _ in range(12):
            C=np.array([[[rng.randrange(8)for _ in range(n)]for _ in range(n)]for _ in range(4)],dtype=np.uint8);expected=[]
            for perm in itertools.permutations(range(n)):
                term=[1]
                for i,j in enumerate(perm):term=multiply(term,list(map(int,C[:,i,j])))
                expected=add(expected,term)
            assert determinant(C)==expected;count+=1
    return count

def main():
    controls=small_controls();T,construction=coefficient_tensor(example_paths(5));n=25;D=len(T)-1
    known=json.loads(Path('evidence/hecke-q5-all-algebraic-places.json').read_text());assert known['coefficient_sha256']==construction['coefficient_sha256']
    field=BinaryExtension(67);mt=np.array([[field.multiply(a,b)for b in range(64)]for a in range(64)],dtype=np.uint8)
    root=next(a for a in range(2,64)if field.power(a,3)^a^1==0);embedding=np.zeros(8,dtype=np.uint8)
    for a in range(8):
        for k in range(3):
            if a>>k&1:embedding[a]^=field.power(root,k)
    assert all(int(embedding[M[a][b]])==field.multiply(int(embedding[a]),int(embedding[b]))for a,b in itertools.product(range(8),repeat=2))
    points=random.Random(846203).sample([a for a in range(64)if a not in embedding],18)
    def field_det(A):
        A=A.copy();d=1
        for j in range(len(A)):
            pivot=next((i for i in range(j,len(A))if A[i,j]),None)
            if pivot is None:return 0
            A[[j,pivot]]=A[[pivot,j]];c=int(A[j,j]);d=field.multiply(d,c);A[j]=mt[field.power(c,62),A[j]]
            for i in range(j+1,len(A)):
                if A[i,j]:A[i]^=mt[A[i,j],A[j]]
        return d
    records=[];geometric=[];closed=[]
    for cut in(0,1):
        C=T.reshape(D+1,n,n)if cut==0 else T.transpose(0,1,3,2,4).reshape(D+1,n,n)
        polynomial=determinant(C);assert polynomial
        endpoint=0
        for c in polynomial:endpoint^=c
        assert endpoint==1
        factors=factor(polynomial);assert len({tuple(f)for f,m in factors})==len(factors)
        product=[polynomial[-1]]
        for f,m in factors:
            assert irreducible(f)
            for _ in range(m):product=multiply(product,f)
        assert product==polynomial
        linear={f[0]:m for f,m in factors if len(f)==2}
        for place in known['places']:
            if place['point']is not None:assert linear.get(place['point'],0)==sum(place['profiles'][cut])
        infinity=n*D-(len(polynomial)-1);assert infinity==sum(known['places'][-1]['profiles'][cut])
        for a in points:
            value=np.zeros((n,n),dtype=np.uint8)
            for layer in C[::-1]:value=mt[a,value]^embedding[layer]
            evaluated=0
            for c in polynomial[::-1]:evaluated=field.multiply(a,evaluated)^int(embedding[c])
            assert field_det(value)==evaluated
        geometric.append(max([infinity]+[m for f,m in factors]));closed.append(max([infinity]+[(len(f)-1)*m for f,m in factors]))
        records.append({'cut':cut,'determinant_coefficients':polynomial,'determinant_degree':len(polynomial)-1,'leading_coefficient':polynomial[-1],'factors':[{'coefficients':f,'degree':len(f)-1,'multiplicity':m}for f,m in factors],'infinity_order_sum':infinity,'exact_geometric_order_maximum':geometric[-1],'exact_degree_weighted_closed_place_maximum':closed[-1]})
    assert geometric==[86,332]and closed==[200,332]and n*D-sum(closed)==368
    names=['check_hecke_q5_determinants.py','f8_polynomial_determinant.py','f8_polynomial_determinant.cpp','f8_polynomial_factor.py','coherent_hot_compounds_f8.cpp','selected_gaussian_coefficients.cpp','native_rectangular_smith_f4.cpp','three_row_hecke_paths.py','three_row_hecke_crossing.py','binary_extension_smith.py']
    out={'completed':True,'matrix_order':n,'tensor_degree':D,'tensor_sha256':construction['coefficient_sha256'],'records':records,'all_factor_products_checked':True,'all_factors_rabin_irreducible':True,'small_determinant_controls':controls,'independent_F64_determinant_evaluations':36,'F8_subfield_root_in_F64':root,'evaluation_points':points,'degree_weighted_two_closed_place_margin':n*D-sum(closed),'source_sha256':{s:sha256(Path('tools',s).read_bytes()).hexdigest()for s in names},'scope':'Exact q5 Hecke determinant factors. Excludes the unweighted uniform-type supporting criterion with one irreducible closed place per cut, including its full conjugate orbit and degree weight. General mixed-place allocations and physical gauges remain outside scope.'}
    Path('evidence/hecke-q5-polynomial-determinants.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'completed':True,'determinant_degrees':[r['determinant_degree']for r in records],'irreducible_factor_counts':[len(r['factors'])for r in records],'closed_place_order_maxima':closed,'margin':n*D-sum(closed),'small_controls':controls,'F64_evaluations':36}))
if __name__=='__main__':main()
