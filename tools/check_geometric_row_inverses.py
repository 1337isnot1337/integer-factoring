#!/usr/bin/env python3
"""Exact row-inverse and finite-root checks for geometric crossings.

Uses only the standard library and the two adjacent nilpotent/geometric
checkers.  Numerical diagnostics are reported separately from exact checks.
"""
import argparse,json,math,time
from collections import Counter
from fractions import Fraction
from pathlib import Path
from check_nilpotent_flag_crossings import (
    nilpotent_family,compose_columns,apply_columns,pmul,smith_series,
    reverse_at_degree,column_cut,rank,det_poly,irreducible,
)
from check_geometric_transpose_crossings import (
    geometric_channel,powers_of,coefficient_columns,invert_columns,
    independent_basis,literal_endpoint,
)


def xor_sum(values):
    result=0
    for v in values:result^=v
    return result


def add_basis(pivots,value):
    while value:
        bit=value.bit_length()-1
        if bit in pivots:value^=pivots[bit]
        else:pivots[bit]=value;return True
    return False


def left_binary(M,K):
    n=len(K)
    selected=[[j for j in range(n) if M[j]>>i&1] for i in range(n)]
    return [[xor_sum(K[k][j] for k in selected[i]) for j in range(n)] for i in range(n)]


def right_binary(K,M):
    n=len(K);selected=[[i for i in range(n) if c>>i&1] for c in M]
    return [[xor_sum(row[k] for k in selected[j]) for j in range(n)] for row in K]


def polynomial_product(A,B):
    n=len(A)
    return [[xor_sum(pmul(A[i][k],B[k][j]) for k in range(n)
                     if A[i][k] and B[k][j]) for j in range(n)] for i in range(n)]


def row_inverse(K):
    n=len(K);degree=max(v.bit_length()-1 for row in K for v in row)
    ci=invert_columns(coefficient_columns(K,0))
    As=[compose_columns(ci,coefficient_columns(K,j)) for j in range(1,degree+1)]
    W=[1<<i for i in range(n)];spaces=[W]
    while W:
        W=independent_basis(apply_columns(A,v) for A in As for v in W)
        assert len(W)<len(spaces[-1]);spaces.append(W)
    basis=[];pivots={}
    for W in reversed(spaces[:-1]):
        for v in W:
            if add_basis(pivots,v):basis.append(v)
    assert len(basis)==n
    bi=invert_columns(basis)
    T=right_binary(left_binary(compose_columns(bi,ci),K),basis)
    assert all(T[i][j]==int(i==j) for i in range(n) for j in range(i+1))
    R=[[int(i==j) for j in range(n)] for i in range(n)]
    for i in range(n-1,-1,-1):
        for j in range(i+1,n):
            R[i][j]=xor_sum(pmul(T[i][k],R[k][j]) for k in range(i+1,j+1)
                            if T[i][k] and R[k][j])
    inverse=right_binary(left_binary(basis,R),compose_columns(bi,ci))
    identity=[[int(i==j) for j in range(n)] for i in range(n)]
    assert polynomial_product(K,inverse)==identity
    assert polynomial_product(inverse,K)==identity
    return inverse,[len(W) for W in spaces]


def adaptive_smith(K):
    n=len(K);degree=max(v.bit_length()-1 for row in K for v in row)
    limit=max(2,n*max(1,degree)+1);precision=max(4,degree+2)
    while True:
        orders=smith_series(K,precision)
        if max(orders)<precision:return orders,precision
        assert precision<limit
        precision=min(limit,2*precision)


def diagnostic(degrees,orders,q):
    def value(t):
        first=max(degrees)*t
        return first+math.log(sum(math.exp(x*t-first) for x in degrees))+math.log(sum(math.exp(-x*t) for x in orders))
    lo,hi=0.,8.
    for _ in range(120):
        a=(2*lo+hi)/3;b=(lo+2*hi)/3
        if value(a)<value(b):hi=b
        else:lo=a
    t=(lo+hi)/2
    return {'rho':value(t)/(4*math.log(q)),'s':math.exp(t),'B':math.exp(value(t))}


def ceil_root(value,power):
    lo=0;hi=1<<((value.bit_length()+power-1)//power)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**power>=value:hi=mid
        else:lo=mid
    return hi


def column_budget_certificate(q,degree,row_degrees):
    """Prove R(s)[1+(n-1)s^(-nL/(n-1))] exceeds the target for s>=1."""
    n=q*q;budget=Fraction(n*degree,n-1);den=budget.denominator
    target=ceil_root(q**99,35);terms=Counter()
    for d,count in Counter(row_degrees).items():
        terms[den*d]+=count
        terms[den*d-budget.numerator]+=(n-1)*count
    terms[0]-=target;shift=max(0,-min(terms));polynomial_degree=max(terms)+shift
    coefficients=[0]*(polynomial_degree+1)
    for exponent,count in terms.items():
        D=exponent+shift;binomial=1
        for j in range(D+1):
            coefficients[j]+=count*binomial
            binomial=binomial*(D-j)//(j+1)
    extra=0
    while any(c<0 for c in coefficients) and extra<512:
        coefficients=[coefficients[0]]+[a+b for a,b in
                     zip(coefficients,coefficients[1:])]+[coefficients[-1]]
        extra+=1
    assert all(c>=0 for c in coefficients) and coefficients[0]>0
    assert target**35>=q**99
    return {'column_budget':str(budget),'integer_target':target,
            'certificate_degree':polynomial_degree+extra,
            'extra_positive_factor_power':extra,
            'constant_coefficient':coefficients[0],
            'all_shifted_coefficients_nonnegative':True}


def field_multiply(a,b,modulus):
    out=0;degree=modulus.bit_length()-1
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if a>>degree:a^=modulus
    return out


def field_power(a,k,modulus):
    out=1
    while k:
        if k&1:out=field_multiply(out,a,modulus)
        k>>=1;a=field_multiply(a,a,modulus)
    return out


def evaluate_at_root(polynomial,modulus):
    out=0
    for i in reversed(range(polynomial.bit_length())):
        out=field_multiply(out,2,modulus)^((polynomial>>i)&1)
    return out


def extension_rank(matrix,modulus):
    matrix=[row[:] for row in matrix];n=len(matrix);pivot_row=0
    for column in range(n):
        pivot=next((i for i in range(pivot_row,n) if matrix[i][column]),None)
        if pivot is None:continue
        matrix[pivot_row],matrix[pivot]=matrix[pivot],matrix[pivot_row]
        inverse=field_power(matrix[pivot_row][column],
                            (1<<(modulus.bit_length()-1))-2,modulus)
        matrix[pivot_row]=[field_multiply(x,inverse,modulus)
                           for x in matrix[pivot_row]]
        for i in range(pivot_row+1,n):
            scalar=matrix[i][column]
            if scalar:
                matrix[i]=[x^field_multiply(scalar,y,modulus)
                           for x,y in zip(matrix[i],matrix[pivot_row])]
        pivot_row+=1
    return pivot_row


def finite_roots_q4(polynomial,source):
    factor={'geometric':21779,'inverse_superoperator':5139}[source]
    assert irreducible(factor)
    column=column_cut(polynomial,4)
    expected=(sum(((factor>>i)&1)<<(16*i)
                  for i in range(factor.bit_length())))<<96
    assert det_poly(column)==expected
    row_rank=extension_rank([[evaluate_at_root(x,factor) for x in row]
                            for row in polynomial],factor)
    col_rank=extension_rank([[evaluate_at_root(x,factor) for x in row]
                            for row in column],factor)
    assert (row_rank,col_rank)==(16,14)
    return {'zero_multiplicity':96,'other_irreducible_factor':factor,
            'factor_degree':factor.bit_length()-1,'factor_multiplicity':16,
            'row_rank_at_other_root':row_rank,'column_rank_at_other_root':col_rank}


def affine_pencil_check(q):
    n=q*q;ps=powers_of(geometric_channel(q),2*q);inverse_u=[0]*n
    for power in ps[:-1]:inverse_u=[a^b for a,b in zip(inverse_u,power)]
    from check_geometric_transpose_crossings import phi_rows
    rows=phi_rows(inverse_u,q)
    columns=[sum(((rows[i]>>j)&1)<<i for i in range(n)) for j in range(n)]
    swap=[1<<((j%q)*q+j//q) for j in range(n)]
    W=compose_columns(swap,columns)
    assert rank(W)==2 and not any(compose_columns(W,W))
    return {'q':q,'update_rank':2,'square_zero':True}


def run(qs,output):
    records=[]
    for q in qs:
        if q not in (4,8):raise ValueError('this finite checkpoint supports q=4,8')
        m=2*q;n=q*q;S=geometric_channel(q);ps=powers_of(S,m)
        Si=[0]*n
        for p in ps[1:-1]:Si=[a^b for a,b in zip(Si,p)]
        for name,channel in [('geometric',S),('inverse_superoperator',Si)]:
            start=time.monotonic();K=nilpotent_family(channel,q,m)
            old_orders,oprec=adaptive_smith(reverse_at_degree(K,m-1))
            inverse,dims=row_inverse(K);literal_endpoint(inverse,q)
            degree=max(v.bit_length()-1 for row in inverse for v in row)
            assert degree==max(old_orders)-(m-1)
            orders,precision=adaptive_smith(column_cut(inverse,q))
            degrees=[v-(m-1) for v in old_orders]
            lead=coefficient_columns(inverse,degree)
            sublead=coefficient_columns(inverse,degree-1)
            leadrows=[sum((lead[j]>>i&1)<<j for j in range(n)) for i in range(n)]
            subrows=[sum((sublead[j]>>i&1)<<j for j in range(n)) for i in range(n)]
            endpoint_product=n*orders.count(0)
            assert endpoint_product**35>q**99
            expected_orders=([0]*4+[8]*12 if q==4 else
                [0]*14+[8]*12+[16]*26+[24]*4+[32]*4+[40]*4)
            assert orders==expected_orders
            rec={'q':q,'source':name,'inverse_degree':degree,
                 'source_infinity_orders':old_orders,'inverse_row_degrees':degrees,
                 'inverse_column_zero_orders':orders,'column_precision':precision,
                 'leading_row_rank':rank(lead),'subleading_row_rank':rank(sublead),
                 'subleading_rank_on_leading_kernel':rank(leadrows+subrows)-rank(leadrows),
                 'protected_constant_product':endpoint_product,
                 'constant_twist_budget_certificate':column_budget_certificate(q,degree,degrees),
                 'source_constant_twist_budget_certificate':column_budget_certificate(
                     q,m-1,[m-1-v for v in old_orders]),
                 'joint_image_dimensions':dims,**diagnostic(degrees,orders,q),
                 'seconds':time.monotonic()-start}
            if q==4:rec['finite_column_roots']=finite_roots_q4(inverse,name)
            records.append(rec)
            print(json.dumps(rec),flush=True)
    result={'row_inverse_instances':records,
            'affine_pencil_checks':[affine_pencil_check(q) for q in [4,8,16]],
            'claim_scope':'exact finite row inverses and q4 finite roots; no exponent improvement'}
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--q',default='4,8');p.add_argument('--output',required=True);a=p.parse_args()
    run(list(map(int,a.q.split(','))),a.output)
