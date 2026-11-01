#!/usr/bin/env python3
"""Exact finite checks for ordered matrix-valued seed amplification."""

from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import random
import time


def identity(n):
    return tuple(int(i==j) for i,j in product(range(n),repeat=2))


def madd(a,b,p):
    return tuple((x+y)%p for x,y in zip(a,b))


def mmul(a,b,n,p):
    return tuple(sum(a[i*n+k]*b[k*n+j] for k in range(n))%p
                 for i,j in product(range(n),repeat=2))


def ordered(values,n,p):
    answer=identity(n)
    for value in values:answer=mmul(answer,value,n,p)
    return answer


def inverse(a,n,p):
    rows=[list(a[i*n:(i+1)*n])+list(identity(n)[i*n:(i+1)*n]) for i in range(n)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if rows[i][j]),None)
        if pivot is None:return None
        rows[j],rows[pivot]=rows[pivot],rows[j]
        scale=pow(rows[j][j],-1,p)
        rows[j]=[x*scale%p for x in rows[j]]
        for i in range(n):
            if i!=j:
                scale=rows[i][j]
                rows[i]=[(x-scale*y)%p for x,y in zip(rows[i],rows[j])]
    return tuple(x for row in rows for x in row[n:])


def kron(a,b):
    return tuple(a[i*2+j]*b[k*2+l] for i,k,j,l in product(range(2),repeat=4))


def bits_rank(rows):
    pivots={}
    for value in rows:
        while value:
            bit=value.bit_length()-1
            if bit in pivots:value^=pivots[bit]
            else:
                pivots[bit]=value
                break
    return len(pivots)


def bits_mv(rows,value):
    return sum(((row&value).bit_count()&1)<<i for i,row in enumerate(rows))


def bits_bra(rows,value):
    answer=0
    for i,row in enumerate(rows):
        if value>>i&1:answer^=row
    return answer


def bits_mm(a,b):
    return tuple(bits_bra(b,row) for row in a)


def bits_inverse(a):
    rows=[row|(1<<(i+4)) for i,row in enumerate(a)]
    for j in range(4):
        pivot=next(i for i in range(j,4) if rows[i]>>j&1)
        rows[j],rows[pivot]=rows[pivot],rows[j]
        for i in range(4):
            if i!=j and rows[i]>>j&1:rows[i]^=rows[j]
    return tuple(row>>4 for row in rows)


def bits_kron(a,b):
    return tuple(sum((((a>>(2*i+j))&1)*((b>>(2*k+l))&1))<<(2*j+l)
                     for j,l in product(range(2),repeat=2))
                 for i,k in product(range(2),repeat=2))


def bits_reshuffle(a):
    return tuple(sum(((a[2*i+k]>>(2*j+l))&1)<<(2*k+l)
                     for k,l in product(range(2),repeat=2))
                 for i,j in product(range(2),repeat=2))


def f4mul(a,b):
    answer=0
    while b:
        if b&1:answer^=a
        b>>=1
        a<<=1
        if a&4:a^=7
    return answer


def f4dot(a,b):
    answer=0
    for x,y in zip(a,b):answer^=f4mul(x,y)
    return answer


def f4_mv_binary(rows,vector):
    return tuple(f4dot(tuple((row>>j)&1 for j in range(4)),vector) for row in rows)


def f4_bra_binary(vector,rows):
    return tuple(f4dot(vector,tuple((row>>j)&1 for row in rows)) for j in range(4))


def f4_product(vector):
    return any(vector) and f4mul(vector[0],vector[3])==f4mul(vector[1],vector[2])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    started=time.monotonic()
    rng=random.Random(20261030)
    p,b,d=2,2,2
    unit=identity(b)
    zero=(0,)*4
    first={indices:tuple(rng.randrange(p) for _ in range(4)) for indices in product(range(d),repeat=4)}
    first[(0,0,0,0)]=(0,1,0,0)
    first[(0,0,0,1)]=(0,0,1,0)
    second={indices:tuple((int(indices[0]==indices[3] and indices[1]==indices[2])*u-x)%p
                          for u,x in zip(unit,value)) for indices,value in first.items()}
    pieces=[first,second]
    assert mmul(first[(0,0,0,0)],first[(0,0,0,1)],2,2)!=mmul(first[(0,0,0,1)],first[(0,0,0,0)],2,2)
    left_entries=list(first.values())+list(second.values())
    for left in left_entries:
        for right in [zero,unit]:assert mmul(left,right,2,2)==mmul(right,left,2,2)

    crossing_checks=type_checks=factor_checks=scalarization_checks=0
    for k in (2,3):
        size=1<<k
        digits=[tuple((number>>(k-1-j))&1 for j in range(k)) for number in range(size)]
        for indices in product(range(size),repeat=4):
            local=[tuple(digits[a][j] for a in indices) for j in range(k)]
            dynamic=[unit]+[zero]*k
            for position,entry in enumerate(local):
                following=[zero]*(k+1)
                for count in range(position+1):
                    for piece in (0,1):
                        term=mmul(dynamic[count],pieces[piece][entry],2,2)
                        following[count+piece]=madd(following[count+piece],term,2)
                dynamic=following
            total=zero
            for value in dynamic:total=madd(total,value,2)
            assert total==tuple(int(indices[0]==indices[3] and indices[1]==indices[2])*x for x in unit)
            crossing_checks+=1
            if indices[0]<2 and indices[1]<2:
                for count in range(k+1):
                    expected=zero
                    for word in product(range(2),repeat=k):
                        if sum(word)==count:
                            expected=madd(expected,ordered([pieces[word[j]][local[j]] for j in range(k)],2,2),2)
                    assert dynamic[count]==expected
                    type_checks+=1
        for _ in range(12):
            local=[tuple(rng.randrange(2) for _ in range(4)) for _ in range(k)]
            word=[rng.randrange(2) for _ in range(k)]
            actual=ordered([pieces[word[j]][local[j]] for j in range(k)],2,2)
            for cut in ("row","col"):
                gathered=zero
                for channels in product(range(4),repeat=k):
                    left=[];right=[]
                    for j,(x,y,u,v) in enumerate(local):
                        a,c=divmod(channels[j],2)
                        index=(x,y,a,c) if cut=="row" else (x,a,u,c)
                        match=(a,c)==((u,v) if cut=="row" else (y,v))
                        left.append(pieces[word[j]][index]);right.append(unit if match else zero)
                    gathered=madd(gathered,mmul(ordered(left,2,2),ordered(right,2,2),2,2),2)
                assert gathered==actual
                factor_checks+=1
    for a,bits in product(range(16),repeat=2):
        left=tuple((a>>i)&1 for i in range(4));right=tuple((bits>>i)&1 for i in range(4))
        assert mmul(left,right,2,2)[0]==sum(left[j]*right[2*j] for j in range(2))%2
        scalarization_checks+=1

    x=(0,1,0,0);y=(0,0,1,0)
    interleaved=ordered([unit,x,y,unit],2,2)
    incorrectly_gathered=mmul(ordered([unit,y],2,2),ordered([x,unit],2,2),2,2)
    assert interleaved!=incorrectly_gathered and interleaved[0]!=incorrectly_gathered[0]

    trace_checks=0
    for _ in range(30):
        a=tuple(rng.randrange(3) for _ in range(4));bb=tuple(rng.randrange(3) for _ in range(4))
        while True:
            w=tuple(rng.randrange(3) for _ in range(16));wi=inverse(w,4,3)
            if wi is not None:break
        left=kron(a,identity(2));right=kron(identity(2),bb)
        assert mmul(left,right,4,3)==mmul(right,left,4,3)
        value=mmul(left,right,4,3)
        conjugated=mmul(mmul(w,value,4,3),wi,4,3)
        expected=((a[0]+a[3])*pow(2,-1,3)*(bb[0]+bb[3])*pow(2,-1,3))%3
        assert sum(value[5*i] for i in range(4))*pow(4,-1,3)%3==expected
        assert sum(conjugated[5*i] for i in range(4))*pow(4,-1,3)%3==expected
        trace_checks+=2
    assert sum(identity(4)[5*i] for i in range(4))%2==0

    products=[v for v in range(1,16) if ((v&1)*((v>>3)&1))==(((v>>1)&1)*((v>>2)&1))]
    assert len(products)==9
    invertible_count=needs_extension=0
    extension_example=None
    for w in product(range(16),repeat=4):
        if bits_rank(w)!=4:continue
        invertible_count+=1
        kets={bits_mv(w,v) for v in products if bits_mv(w,v) in products}
        bras=[v for v in products if bits_bra(w,v) in products]
        if not any((ket&bra).bit_count()&1 for ket in kets for bra in bras):
            needs_extension+=1
            if extension_example is None:extension_example=w
    assert invertible_count==20160 and needs_extension==432
    assert extension_example==(1,2,6,11)
    wi=bits_inverse(extension_example)
    density_candidates=0
    for a,bb in product(range(16),repeat=2):
        if ((a&1)^((a>>3)&1))!=1 or ((bb&1)^((bb>>3)&1))!=1:continue
        density_candidates+=1
        density=bits_kron(a,bb)
        transformed=bits_mm(bits_mm(wi,density),extension_example)
        assert bits_rank(bits_reshuffle(transformed))!=1
    assert density_candidates==64

    lines=[(1,a) for a in range(4)]+[(0,1)]
    product_vectors=[tuple(f4mul(x,y) for x in a for y in bb) for a,bb in product(lines,repeat=2)]
    kets=[f4_mv_binary(extension_example,ket) for ket in product_vectors
          if f4_product(f4_mv_binary(extension_example,ket))]
    bras=[bra for bra in product_vectors if f4_product(f4_bra_binary(bra,extension_example))]
    witness=None
    for ket,bra in product(kets,bras):
        pairing=f4dot(bra,ket)
        if pairing:
            inv=next(x for x in range(1,4) if f4mul(x,pairing)==1)
            bra=tuple(f4mul(inv,x) for x in bra)
            assert f4dot(bra,ket)==1
            assert f4_product(ket) and f4_product(f4_mv_binary(wi,ket))
            assert f4_product(bra) and f4_product(f4_bra_binary(bra,extension_example))
            witness={"field_polynomial":[1,1,1],"matrix_row_bitmasks":extension_example,
                     "bra":bra,"ket":ket,"pairing":1}
            break
    assert witness is not None

    result={"status":"passed",
            "scope":"Finite ordered matrix identities and scalarizations; no improved catalytic seed.",
            "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "ordered_crossing_entries":crossing_checks,"type_word_comparisons":type_checks,
            "gathered_factor_entry_checks":factor_checks,"matrix_entry_scalarization_checks":scalarization_checks,
            "normalized_trace_checks":trace_checks,"negative_interleaving_control":True,
            "gl4_f2_matrices":invertible_count,"no_basefield_common_bra_ket":needs_extension,
            "example_trace_one_product_densities_rejected":density_candidates,
            "f4_common_scalarization_witness":witness,"seconds":time.monotonic()-started}
    text=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:args.output.write_text(text)
    print(text,end="")


if __name__=="__main__":main()
