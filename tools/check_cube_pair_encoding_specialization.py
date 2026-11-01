#!/usr/bin/env python3
"""Exact negative control for diagonal specialization of an encoded cube tensor.

Uses only the standard library. Every retained Smith profile resolves all
16 pivots inside the stated power-series precision. The candidate itself
is defined by homogeneous degree-six Gaussian kernels, then restricted
to the local vacuum/pair alphabet.
"""

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path


FIELD_MODULUS=0x11B
PRECISION=64

def scalar_multiply(a,b):
    result=0
    while b:
        if b&1:result^=a
        a<<=1
        if a&256:a^=FIELD_MODULUS
        b>>=1
    return result

MUL=tuple(tuple(scalar_multiply(a,b) for b in range(256)) for a in range(256))
INV=(0,)+tuple(next(b for b in range(1,256) if MUL[a][b]==1) for a in range(1,256))

def power(a,k):
    result=1
    for _ in range(k):result=MUL[result][a]
    return result

def trim(a):
    while a and not a[-1]:a.pop()
    return tuple(a)

def add(a,b):
    result=list(a)
    if len(result)<len(b):result.extend([0]*(len(b)-len(result)))
    for i,x in enumerate(b):result[i]^=x
    return trim(result)

def multiply(a,b,limit=None):
    if not a or not b:return ()
    size=len(a)+len(b)-1
    if limit is not None:size=min(size,limit)
    result=[0]*size
    for i,x in enumerate(a[:size]):
        if x:
            row=MUL[x]
            for j,y in enumerate(b[:size-i]):
                if y:result[i+j]^=row[y]
    return trim(result)

def monomial(exponents,other):
    coefficient=MUL[power(other[0],exponents[1])][power(other[1],exponents[2])]
    return (0,)*exponents[0]+(coefficient,)

def kernel_parts(other):
    A_aug=[() for _ in range(8)];B_aug=A_aug.copy()
    for axis in range(3):
        exponent=[0]*3;exponent[axis]=6
        A_aug[1<<axis]=B_aug[1<<axis]=monomial(exponent,other)
    for a in range(3):
        for b in range(a+1,3):
            subset=(1<<a)|(1<<b)
            exponent=[0]*3;exponent[a]=exponent[b]=3
            A_aug[subset]=B_aug[subset]=monomial(exponent,other)
            for e,f in ((2,4),(4,2)):
                exponent[a]=e;exponent[b]=f
                A_aug[subset]=add(A_aug[subset],monomial(exponent,other))
            for e,f in ((1,5),(5,1)):
                exponent[a]=e;exponent[b]=f
                B_aug[subset]=add(B_aug[subset],monomial(exponent,other))
    A_aug[7]=B_aug[7]=monomial((2,2,2),other)
    s=()
    for exponents in ((4,1,1),(1,4,1),(1,1,4),(2,2,2)):
        s=add(s,monomial(exponents,other))
    def steps(aug):
        result=[() for _ in range(8)]
        for delta in range(1,8):
            for subset in range(1,8):
                if delta&subset==delta:result[delta]=add(result[delta],aug[subset])
        return result
    return steps(A_aug),steps(B_aug),s

def pfaffians(steps):
    result=[() for _ in range(256)];result[0]=(1,)
    for mask in range(1,256):
        if mask.bit_count()%2:continue
        bit=mask&-mask;first=bit.bit_length()-1;left=mask^bit;scan=left
        while scan:
            bit=scan&-scan;second=bit.bit_length()-1;scan^=bit
            result[mask]=add(result[mask],multiply(steps[first^second],result[left^bit]))
    return result

def encoded_tensor(other):
    A,B,s=kernel_parts(other);PA=pfaffians(A);PB=pfaffians(B)
    H=[]
    for mask in range(256):
        if not mask.bit_count()%2:H.append(multiply(PA[mask],PB[mask]));continue
        u=v=();scan=mask
        while scan:
            bit=scan&-scan;scan^=bit
            u=add(u,PA[mask^bit]);v=add(v,PB[mask^bit])
        H.append(multiply(s,multiply(u,v)))
    T=H.copy()
    for vertex in range(8):
        for mask in range(256):
            if mask>>vertex&1:T[mask]=add(T[mask],multiply(s,T[mask^(1<<vertex)]))
    return T,H,(A,B,s)

def literal_tensor(parts):
    """Independent direct Pfaffian recurrence on all sixteen mode indices."""
    A,B,s=parts
    @lru_cache(None)
    def visit(mask):
        if mask==0:return (1,)
        first_bit=mask&-mask;first=first_bit.bit_length()-1;left=mask^first_bit;scan=left;result=()
        while scan:
            bit=scan&-scan;second=bit.bit_length()-1;scan^=bit
            v,i=divmod(first,2);w,j=divmod(second,2)
            edge=() if v==w else (A[v^w] if i==j==0 else B[v^w] if i==j==1 else s)
            if edge:result=add(result,multiply(edge,visit(left^bit)))
        return result
    return [visit(sum(3<<(2*v) for v in range(8) if physical>>v&1))
            for physical in range(256)]

def smith_profile(tensor,precision=PRECISION):
    rows=[sum(((mask>>j)&1)<<(2*j) for j in range(4)) for mask in range(16)]
    columns=[sum(((mask>>j)&1)<<(2*j+1) for j in range(4)) for mask in range(16)]
    matrix=[[tensor[r|c][:precision] for c in columns] for r in rows]
    row_indices=list(range(16));col_indices=list(range(16));orders=[];witnesses=[]
    def valuation(a):return next((i for i,x in enumerate(a) if x),precision)
    for k in range(16):
        value,i,j=min((valuation(matrix[i][j]),i,j) for i in range(k,16) for j in range(k,16))
        assert value<precision,('unresolved precision',k,precision)
        matrix[k],matrix[i]=matrix[i],matrix[k];row_indices[k],row_indices[i]=row_indices[i],row_indices[k]
        for row in matrix:row[k],row[j]=row[j],row[k]
        col_indices[k],col_indices[j]=col_indices[j],col_indices[k]
        orders.append(value);length=precision-value;unit=matrix[k][k][value:]
        inverse=[INV[unit[0]]]
        for n in range(1,length):
            coefficient=0
            for j in range(1,min(n+1,len(unit))):coefficient^=MUL[unit[j]][inverse[n-j]]
            inverse.append(MUL[inverse[0]][coefficient])
        for i in range(k+1,16):
            quotient=multiply(matrix[i][k][value:],inverse,limit=length)
            if not quotient:continue
            for j in range(k+1,16):matrix[i][j]=add(matrix[i][j],multiply(quotient,matrix[k][j],limit=precision))
        witnesses.append({'rows':row_indices[:k+1].copy(),'columns':col_indices[:k+1].copy(),
                          'minor_valuation':sum(orders)})
    assert orders==sorted(orders)
    return orders,witnesses

DIAGONAL_PROFILE=[0,7,7,14,17,17,24,24,24,24,31,31,34,41,41,48]
UNEQUAL_F4_PROFILE=[0,2,3,3,3,3,4,6,13,13,13,13,22,24,26,28]
F256_PROFILE=[0,2,3,3,3,3,4,6,11,11,11,11,18,20,22,24]
BASELINE_PROFILE=[0,4,8,12,16,20,20,24,24,28,28,32,36,40,44,48]


def subset_sums(weights):
    return sorted(sum(value for j,value in enumerate(weights) if mask>>j&1)
                  for mask in range(1<<len(weights)))


def integer_digest(value):
    return hashlib.sha256(value.to_bytes((value.bit_length()+7)//8,'big')).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'evidence'/'cube-pair-encoding-specialization.json')
    args=parser.parse_args()
    assert all(MUL[a][INV[a]]==1 for a in range(1,256))
    roots=[a for a in range(256) if MUL[a][a]^a^1==0]
    assert roots==[188,189]
    assert subset_sums([7,7,17,17])==DIAGONAL_PROFILE
    assert subset_sums([4,8,16,20])==BASELINE_PROFILE

    records=[];literal_checks=0
    cases=[((1,1),DIAGONAL_PROFILE),((1,188),UNEQUAL_F4_PROFILE),((2,3),F256_PROFILE)]
    for other,expected in cases:
        T,H,parts=encoded_tensor(other)
        assert literal_tensor(parts)==T
        literal_checks+=256
        assert all(len(entry)<=6*mask.bit_count()+1 for mask,entry in enumerate(T))
        profile,witnesses=smith_profile(H)
        assert profile==expected,(other,profile)
        assert smith_profile(T)[0]==profile
        if other==(1,1):
            for mask,entry in enumerate(T):
                antipode=sum(((mask>>v)&1)<<(v^7) for v in range(8))
                value=0
                for coefficient in entry:value^=coefficient
                assert value==int(mask==antipode)
        records.append({'other_coordinates':list(other),'orders':profile,
                        'odd_physical_words_with_nonzero_polynomial':sum(bool(T[m]) for m in range(256) if m.bit_count()%2),
                        'attaining_prefix_minors_of_sheared_tensor':witnesses})

    f4_cases=[]
    for a in [1,*roots]:
        for b in [1,*roots]:
            T,H,_=encoded_tensor((a,b))
            expected=DIAGONAL_PROFILE if a==b else UNEQUAL_F4_PROFILE
            assert smith_profile(H)[0]==expected
            f4_cases.append({'other_coordinates':[a,b],'profile_kind':'diagonal' if a==b else 'unequal'})
    for r in (2,17):
        assert smith_profile(encoded_tensor((r,r))[1])[0]==DIAGONAL_PROFILE
    for other in ((5,7),(17,91)):
        assert smith_profile(encoded_tensor(other)[1])[0]==F256_PROFILE

    T,_,_=encoded_tensor((2,3))
    antipodal_minor=T[(1<<0)|(1<<7)]
    assert antipodal_minor==(0,0,31,0,0,0,0,0,20)
    field_leading=MUL[MUL[power(2,2)][power(3,2)]][power(2,6)^power(3,6)]
    assert field_leading==31
    assert all(power(a,6)^power(b,6)==0 for a in [1,*roots] for b in [1,*roots])

    margins={}
    for label,profile in [('unequal_f4',UNEQUAL_F4_PROFILE),('f256',F256_PROFILE)]:
        difference=[sum(BASELINE_PROFILE[:k])-sum(profile[:k]) for k in range(1,17)]
        assert min(difference)>=0
        margins[label]=difference

    # At t=41/60 and p=3/50, the diagonal profile gives rho<5/9.
    # This is exactly Z(t)^450 * t^(-432) * (3/50)^60 * (47/50)^940 < 1.
    z_numerator=(60**7+41**7)**2*(60**17+41**17)**2
    left=z_numerator**450*3**60*47**940
    right=60**21168*41**432*50**1000
    assert left<right
    assert 3*5**2<9**2  # 5/9 < 1/sqrt(3).
    assert 3*3**2>5**2  # 3/5 > 1/sqrt(3).

    report={
        'status':'passed',
        'scope':'Exact false-positive control for diagonal specialization; no improved factoring bound',
        'coefficient_field_modulus':hex(FIELD_MODULUS),
        'nonzero_field_inverses_checked':255,
        'power_series_precision':PRECISION,
        'literal_sixteen_mode_pfaffians_checked':literal_checks,
        'endpoint_physical_words_checked':256,
        'profiles':records,
        'f4_nonzero_specializations':f4_cases,
        'additional_diagonal_specializations':[[2,2],[17,17]],
        'additional_f256_specializations':[[5,7],[17,91]],
        'explicit_antipodal_minor_at_2_3':{'rows':[0,1],'columns':[0,8],'polynomial_nonzero_coefficients':{'2':31,'8':20}},
        'baseline_profile':BASELINE_PROFILE,
        'baseline_prefix_margins':margins,
        'false_positive_rational_witness':{
            't':'41/60','p':'3/50','strict_upper_bound':'5/9',
            'left_integer_bits':left.bit_length(),'right_integer_bits':right.bit_length(),
            'left_integer_sha256':integer_digest(left),'right_integer_sha256':integer_digest(right),
            'exact_comparison_passed':True},
        'generic_rejection':'Specialization and weak majorization place the unweighted common-type pressure above the existing baseline floor 3/5.'
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','literal_pfaffians':literal_checks,
                      'full_profiles':14,'false_positive_upper_bound':'5/9',
                      'generic_pressure_lower_bound':'3/5','evidence':str(args.output)},sort_keys=True))


if __name__=='__main__':
    main()
