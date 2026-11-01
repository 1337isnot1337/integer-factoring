#!/usr/bin/env python3
"""All-field finite crossing cutoff97 from compatible integer cores.

Three denominator3/4/5 families are checked by exact integer determinants.
Every maximal-ratio class is covered, with rational cut thresholds and,
where needed, exact rational subintervals. Downward rounding is integral.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd
from pathlib import Path
import hashlib
import json
from time import monotonic
from check_interior_jacobian_budget import (
    geometric_complement,weighted_budgets,class_parameters,class_bound,fraction_record,
)
from check_rational_cut_mesh_budget import mesh_budgets,rational_mesh_class_bound


def jacobian(A,D):
    m=len(A);ell=len(A[0]);s=len(A[0][0]);q=len(D[0]);p=len(D[0][0])
    rows=[]
    for Ai,Di in zip(A,D):
        for u in range(q):
            for v in range(s):
                row=[0]*(p*s+q*ell)
                for k in range(p):row[k*s+v]=Di[u][k]
                for k in range(ell):row[p*s+u*ell+k]=-Ai[k][v]
                rows.append(row)
    return rows


def determinant(matrix):
    a=[row[:] for row in matrix];n=len(a);previous=1;sign=1
    assert all(len(row)==n for row in a)
    for j in range(n-1):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];sign=-sign
        scale=a[j][j]
        for i in range(j+1,n):
            factor=a[i][j]
            for k in range(j+1,n):
                numerator=a[i][k]*scale-factor*a[j][k]
                assert numerator%previous==0
                a[i][k]=numerator//previous
            a[i][j]=0
        previous=scale
    return sign*a[-1][-1]


def block_diagonal(matrices):
    height=sum(len(a) for a in matrices);width=sum(len(a[0]) for a in matrices)
    out=[[0]*width for _ in range(height)];u=v=0
    for a in matrices:
        for i,row in enumerate(a):out[u+i][v:v+len(row)]=row
        u+=len(a);v+=len(a[0])
    return out


def assembly_control(A,D,T):
    """Check all four cross blocks in a nontrivial direct-sum Jacobian."""
    atypes=(T,1);dtypes=(T,2)
    aa=[block_diagonal([A[t][i] for t in atypes]) for i in range(3)]
    dd=[block_diagonal([D[t][i] for t in dtypes]) for i in range(3)]
    J=jacobian(aa,dd);s=sum(atypes);q=sum(dtypes)
    ell=len(aa[0]);p=len(dd[0][0]);a_in=a_out=0;covered_rows=[];covered_cols=[]
    for at in atypes:
        d_in=d_out=0
        la=len(A[at][0])
        for dt in dtypes:
            pd=len(D[dt][0][0])
            rr=[i*q*s+(d_out+u)*s+a_in+v for i in range(3) for u in range(dt) for v in range(at)]
            cc=[(d_in+k)*s+a_in+v for k in range(pd) for v in range(at)]
            cc += [p*s+(d_out+u)*ell+a_out+k for u in range(dt) for k in range(la)]
            local=jacobian(A[at],D[dt])
            assert [[J[r][c] for c in cc] for r in rr]==local
            used=set(cc)
            assert all(J[r][c]==0 for r in rr for c in range(len(J[0])) if c not in used)
            covered_rows.extend(rr);covered_cols.extend(cc)
            d_in+=pd;d_out+=dt
        a_in+=at;a_out+=la
    assert sorted(covered_rows)==list(range(len(J)))
    assert sorted(covered_cols)==list(range(len(J[0])))
    return dict(input_dimension=s,output_codimension=q,rows=len(J),columns=len(J[0]),cross_blocks=4)


def verify_cores(path):
    data=json.loads(path.read_text());records=[];descriptors=[]
    for family in data['families']:
        T=family['denominator'];p=family['p'];ell=family['ell']
        assert family['m']==3 and p==T+1 and ell==2*T-1
        A={int(k):v for k,v in family['A'].items()};D={int(k):v for k,v in family['D'].items()}
        assert sorted(A)==sorted(D)==list(range(1,T+1))
        for t in range(1,T+1):
            assert len(A[t])==len(D[t])==3
            assert all(len(a)==(ell*t+T-1)//T and all(len(row)==t for row in a) for a in A[t])
            assert all(len(a)==t and all(len(row)==(p*t+T-1)//T for row in a) for a in D[t])
        pairs=[]
        for at in range(1,T+1):
            for dt in range(1,T+1):
                J=jacobian(A[at],D[dt]);r=len(J);c=len(J[0]);cert=family['certificates'][f'{at},{dt}']
                assert r==3*at*dt and cert['rows']==r and cert['columns']==c
                determinants=[]
                for minor in cert['minors']:
                    cols=minor['columns']
                    assert len(cols)==len(set(cols))==r and all(0<=j<c for j in cols)
                    value=determinant([[row[j] for j in cols] for row in J])
                    assert value==minor['determinant']
                    determinants.append(value)
                g=0
                for value in determinants:g=gcd(g,abs(value))
                assert g==1
                pairs.append(dict(A_type=at,D_type=dt,rows=r,columns=c,determinants=determinants,gcd=1))
        records.append(dict(denominator=T,p=p,ell=ell,pairs=pairs,assembly=assembly_control(A,D,T)))
        descriptors.append((T,p,ell))
    assert sorted(descriptors)==[(3,4,5),(4,5,7),(5,6,9)]
    return descriptors,records


def compatible_complement(d,m,descriptors):
    best=geometric_complement(d,m);witness=None
    if m<3 or m>=2*d-1:return best,witness
    for T,p,ell in descriptors:
        for a in range(m-2):
            b=m-3-a
            for s in range(1,T*(d-1)//(T*b+ell)+1):
                q=min(T*(d-s)//(T*a+p),d-b*s-(ell*s+T-1)//T)
                if q>0 and s*q>best:
                    best=s*q;witness=dict(T=T,p=p,ell=ell,a=a,b=b,s=s,q=q)
    return best,witness


def complement_table(d,descriptors):
    F=[d*d,d*(d-1)//2,d*(d-1)//2];improvements=[]
    for m in range(3,2*d):
        value,witness=compatible_complement(d,m,descriptors);F.append(value)
        if witness:improvements.append(dict(m=m,old=geometric_complement(d,m),new=value,witness=witness))
    for m in range(2*d-2,-1,-1):F[m]=max(F[m],F[m+1])
    assert F[-1]==0
    return F,improvements


def refined_bin_bound(d,k,b,F,MG,H,left,right,scale):
    a0=(left*left*b).__floor__()+1;a1=min((right*right*b).__floor__(),d*d)
    if a0>a1:return None
    fb=F[b] if b<len(F) else 0
    R,C=max(d*d-a1,fb),d*d-b;e=max(0,a0+fb-d*d)
    if b==1:
        R=max(R,d*d-k*(k+1)//2);e=max(e,(k*k-3*k+4)//2)
    end=(right.numerator*H+right.denominator-1)//right.denominator
    nodes=[Q(j,H) for j in range(H,end)]+[right]
    base=1/(nodes[-2]+nodes[-1])
    terms=[base*(2*d*d+e)]
    for j,s in enumerate(nodes[:-1]):
        w=Q(1,2)-1/(1+nodes[1]) if j==0 else 1/(nodes[j-1]+s)-1/(s+nodes[j+1])
        assert w>0
        ss=s*s;mu=(ss*right*right-1)/(right**4-1);nu=(right*right-ss)/(right**4-1)
        g0=Q(MG[j][0],H*H)
        gb=(Q(MG[j][b],H*H) if b<len(F) else ss*b) if a0>ss*b else g0
        U=max(gb,min(a0,ss*b)+mu*R+nu*C);V=max(g0,b+nu*R+mu*C)
        terms.append(w*(U+V))
    lower=sum((term.numerator*scale)//term.denominator for term in terms)
    return dict(a_lower=a0,a_upper=a1,scaled_lower_bound=lower,terms=len(terms))


def dimension_certificate(d,F):
    floor=Q(16131,25);scale=10**9;H=16;bin_count=100
    assert 3363**2-2*2378**2==1
    assert floor.numerator**2378>d**3363*floor.denominator**2378
    G=weighted_budgets(d,F);MG4=mesh_budgets(d,F,4);MG16=mesh_budgets(d,F,H)
    total,base_count,integer_count,quarter_count=1,1,0,0;refined=[]
    for k in range(2,d+1):
        for b in range(1,(d*d-1)//((k-1)**2)+1):
            total+=1;*_,e=class_parameters(d,k,b,F)
            if Q(2*d*d+e,2*k-1)>=floor:base_count+=1;continue
            if class_bound(d,k,b,F,G)>=floor:integer_count+=1;continue
            qbound=rational_mesh_class_bound(d,k,b,F,MG4,4)
            if qbound>=floor:quarter_count+=1;continue
            bins=[]
            for j in range(bin_count):
                left=Q((k-1)*bin_count+j,bin_count);right=left+Q(1,bin_count)
                certificate=refined_bin_bound(d,k,b,F,MG16,H,left,right,scale)
                if certificate is None:
                    bins.append(dict(index=j,empty_integer_rank_interval=True));continue
                assert Q(certificate['scaled_lower_bound'],scale)>=floor,(d,k,b,j,certificate)
                bins.append(dict(index=j,**certificate))
            assert any('scaled_lower_bound' in z for z in bins)
            refined.append(dict(k=k,b=b,quarter_bound=fraction_record(qbound),bins=bins))
    assert total==1+sum((d*d-1)//((k-1)**2) for k in range(2,d+1))
    assert total==base_count+integer_count+quarter_count+len(refined)
    minimum=min((z['scaled_lower_bound'],r['k'],r['b'],z['index']) for r in refined for z in r['bins'] if 'scaled_lower_bound' in z)
    return dict(d=d,scalar_classes=total,base_term_classes=base_count,integer_mesh_classes=integer_count,
                quarter_mesh_classes=quarter_count,refined_ratio_classes=len(refined),
                rational_ratio_subintervals=bin_count,cut_mesh_denominator=H,downward_rounding_scale=scale,
                certified_common_floor=fraction_record(floor),upper_Pell_exponent=[3363,2378],
                minimum_scaled_subinterval_bound=dict(numerator=minimum[0],denominator=scale,k=minimum[1],b=minimum[2],bin=minimum[3]),
                refined_classes=refined,complement_table=F)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    parser.add_argument('--cores',type=Path,default=Path(__file__).with_name('compatible-compression-cores.json'))
    args=parser.parse_args();started=monotonic()
    descriptors,cores=verify_cores(args.cores)
    F,improvements=complement_table(97,descriptors)
    proof=dimension_certificate(97,F)
    result=dict(status='PASS',scope='All-field exact finite crossing rank-sum exclusion through97; no factoring exponent improvement',
                all_field_cutoff=97,previous_all_field_cutoff=96,
                small_source_tensor_power_constructions_not_excluded=True,
                core_certificates=cores,all_field_complement_improvements=improvements,dimension97=proof,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                literal_core_data_sha256=hashlib.sha256(args.cores.read_bytes()).hexdigest(),
                seconds=monotonic()-started)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('core_certificates','dimension97','all_field_complement_improvements')},indent=2))
    print(json.dumps({k:v for k,v in proof.items() if k not in ('refined_classes','complement_table')},indent=2))

if __name__=='__main__':main()
