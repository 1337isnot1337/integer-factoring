#!/usr/bin/env python3
"""Exact carry-based swap circuits and cheap local-profile exclusions."""
import argparse,json,math
from fractions import Fraction
from pathlib import Path
from check_nilpotent_flag_crossings import (
    column_cut,rank,toeplitz_rank,reverse_at_degree,det_poly,polynomial_gcd,
)


def permutations(q):
    A=[((u+v)%q)*q+v for u in range(q) for v in range(q)]
    B=[u*q+(v-u)%q for u in range(q) for v in range(q)]
    C=[u*q+(-v)%q for u in range(q) for v in range(q)]
    return A,B,C


def permute_rows(K,P):
    out=[None]*len(P)
    for i,j in enumerate(P):out[j]=K[i]
    return out


def linear_factor(K,P):
    PK=permute_rows(K,P)
    return [[a^((a^b)<<1) for a,b in zip(row,prow)]for row,prow in zip(K,PK)]


def inverse_linear_factor(K,P,order):
    out=[row[:]for row in K];term=K
    for _ in range(1,order):
        pt=permute_rows(term,P)
        term=[[(a^b)<<1 for a,b in zip(row,prow)]for row,prow in zip(term,pt)]
        out=[[a^b for a,b in zip(row,trow)]for row,trow in zip(out,term)]
    pt=permute_rows(term,P)
    assert all(a==b for row,prow in zip(term,pt)for a,b in zip(row,prow))
    # Verify the inverse-factor identity independently of the series sum.
    assert linear_factor(out,P)==K
    return out


def inverse_permutation(P):
    out=[0]*len(P)
    for i,j in enumerate(P):out[j]=i
    return out


def construction(q,path):
    assert q>=2 and q&(q-1)==0
    n=q*q;K=[[int(i==j)for j in range(n)]for i in range(n)]
    A,B,C=permutations(q)
    if path=='direct':
        for P in [A,B,A,C]:K=linear_factor(K,P)
    elif path=='inverse':
        for P,order in [(C,2),(A,q),(B,q),(A,q)]:K=inverse_linear_factor(K,P,order)
    elif path=='binary_compressed':
        for P in [A,B,A]:
            Q=inverse_permutation(P)
            for _ in range(q.bit_length()-1):
                K=linear_factor(K,Q);Q=[Q[j]for j in Q]
            assert Q==list(range(n))
        K=linear_factor(K,C)
    else:raise ValueError(path)
    assert all(v.bit_count()%2==int(j==(i%q)*q+i//q)
               for i,row in enumerate(K)for j,v in enumerate(row))
    return K


def coefficient_rank(K,d):
    return rank([sum(((x>>d)&1)<<j for j,x in enumerate(row))for row in K])


def jet_orders(K,depth=3):
    ts=[0]+[toeplitz_rank(K,k)for k in range(1,depth+1)]
    counts=[ts[1]]+[ts[j+1]-2*ts[j]+ts[j-1]for j in range(1,depth)]
    assert all(x>=0 for x in counts)
    return counts,ts[1:]


def gate_certificate(q,L,r0,counts,used,bound):
    """AM-GM row bound and a rational upper bound for the crossing scale."""
    n=q*q;A,B=bound
    assert all(0<=a<=b for a,b in zip(used,counts))
    assert len(used)<=L  # Each retained column order is strictly below L.
    assert r0*A**L>=n*B**L
    lower=n*sum(Fraction(mu*B**j,A**j)for j,mu in enumerate(used))
    assert lower.numerator**35>q**99*lower.denominator**35
    return {'row_leading_rank':r0,'column_counts_used':used,
            'crossing_scale_upper_bound':[A,B],
            'formal_B_lower_bound':[lower.numerator,lower.denominator],
            'comparison':'B_lower_bound > q^(99/35) > q^(2*sqrt(2))'}


def verify(witness_path):
    witnesses=json.loads(witness_path.read_text())
    configurations={
        (4,'direct'):([1,3,1],(8,5)),
        (4,'inverse'):([1,3,1],(3,2)),
        (4,'binary_compressed'):([1,3,1],(8,5)),
        (8,'inverse'):([1,7],(4,3)),
        (8,'binary_compressed'):([1,7,2],(7,4)),
        (16,'inverse'):([1,15],(6,5)),
        (16,'binary_compressed'):([1,15],(5,3)),
    }
    records=[]
    for q in [4,8,16]:
        for path in ['direct','inverse','binary_compressed']:
            K=construction(q,path);C=column_cut(K,q)
            L=max(v.bit_length()-1 for row in K for v in row)
            r0,c0=coefficient_rank(K,L),coefficient_rank(C,L)
            counts,toeplitz=jet_orders(C)
            row_counts,row_toeplitz=jet_orders(reverse_at_degree(K,L),2)
            rec={'q':q,'path':path,'degree':L,'top_coefficient_ranks':[r0,c0],
                 'column_zero_order_counts_0_1_2':counts,'column_toeplitz_ranks':toeplitz,
                 'row_infinity_order_counts_0_1':row_counts,'row_toeplitz_ranks':row_toeplitz}
            if path=='direct' and q in [8,16]:
                assert (r0*c0)**35>q**99
                rec['exclusion']={'kind':'actual highest coefficient at every tensor power',
                                  'rank_product':r0*c0}
            else:
                used,bound=configurations[q,path]
                rec['exclusion']=gate_certificate(q,L,r0,counts,used,bound)
            if path=='inverse' and q in [8,16]:
                witness=next(item for item in witnesses if item['q']==q)
                g=0
                for item in witness['selected_minors']:
                    rows,cols=item['rows'],item['columns']
                    D=det_poly([[C[i][j]for j in cols]for i in rows])
                    assert D==item['determinant'];g=polynomial_gcd(g,D)
                assert g==witness['gcd'] and g and g&(g-1)==0
                k=witness['minor_size'];assert (q*q*k)**35>q**99
                rec['all_nonzero_finite_basepoints']={
                    'column_rank_at_least':k,'minor_gcd_monomial_degree':g.bit_length()-1,
                    'protected_constant_rank_product_at_least':q*q*k}
            records.append(rec)
    return {'instances':records,
            'scope':'all-size row-unimodular identities; finite q4/q8/q16 exact local-profile exclusions; no exponent improvement'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--witnesses',type=Path,required=True)
    p.add_argument('--output',type=Path);a=p.parse_args()
    result=verify(a.witnesses);rendered=json.dumps(result,indent=2)+'\n'
    if a.output:a.output.write_text(rendered)
    print(rendered,end='')
