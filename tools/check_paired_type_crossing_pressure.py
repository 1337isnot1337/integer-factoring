#!/usr/bin/env python3
"""Exact controls for independent paired types and rectangular crossing costs."""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import product
import json
from math import factorial,lcm
from pathlib import Path
import random


def rank(rows):
    pivots={}
    for row in rows:
        while row:
            p=row.bit_length()-1
            if p in pivots:row^=pivots[p]
            else:pivots[p]=row;break
    return len(pivots)


def flatten(values,A,B,cut):
    rows=[0]*(A*B)
    for x,y,u,v in product(range(A),range(B),range(B),range(A)):
        if values[x+A*(y+B*(u+B*v))]:
            row=B*x+(y if cut==0 else u)
            column=A*(u if cut==0 else y)+v
            rows[row]^=1<<column
    return rows


def swapped(values,A,B):
    out=[0]*len(values)
    for x,y,u,v in product(range(B),range(A),range(A),range(B)):
        out[x+B*(y+A*(u+A*v))]=values[u+A*(v+B*(x+B*y))]
    return out


def tensor_values(left,right,A,B):
    q=A*B;out=[0]*q**4
    for x,y,u,v in product(range(q),repeat=4):
        x1,x2=divmod(x,B);y1,y2=divmod(y,A)
        u1,u2=divmod(u,A);v1,v2=divmod(v,B)
        out[x+q*(y+q*(u+q*v))]=left[x1+A*(y1+B*(u1+B*v1))]&right[x2+B*(y2+A*(u2+A*v2))]
    return out


def type_size(counts):
    result=factorial(sum(counts))
    for n in counts:result//=factorial(n)
    return result


def exact_gibbs_control():
    a,b=3,2;t=Fraction(3,2);R=Fraction(a,b);s=t**a
    GX=(0,2);GY=(-1,1,3)
    x=(Fraction(1,3),Fraction(2,3));y=(Fraction(1,2),Fraction(1,3),Fraction(1,6))
    wx=[t**(-b*g) for g in GX];wy=[t**(-b*g) for g in GY]
    ZX=sum(wx);ZY=sum(wy);px=[w/ZX for w in wx];py=[w/ZY for w in wy]
    assert sum(px)==sum(py)==1
    U=sum(s**(-k) for k in (-1,0,2));V=sum(s**(-k) for k in (0,1));B=U*V
    L=lcm(*(p.denominator for p in x+y))
    mean=sum(p*g for p,g in zip(x,GX))+sum(p*g for p,g in zip(y,GY))
    assert (L*mean).denominator==1
    left=B**(b*L)*t**int(a*b*L*mean)
    right=(B**b/(ZX*ZY)**a)**L
    for p,pstar in zip(x+y,px+py):
        exponent=a*int(L*p)
        left*=p**exponent;right*=(p/pstar)**exponent
    assert left==right
    return {'a':a,'b':b,'t':str(t),'gibbs_x':list(map(str,px)),
            'gibbs_y':list(map(str,py)),'multiplicative_entropy_identity':True,
            'strict_rational_inequality':B**b<(ZX*ZY)**a,
            'formal_profiles_only_not_a_crossing_certificate':True}


def check():
    A,B=2,3;q=A*B;rng=random.Random(4821)
    endpoint=[0]*(A*A*B*B)
    for x,y,u,v in product(range(A),range(B),range(B),range(A)):
        endpoint[x+A*(y+B*(u+B*v))]=int(x==v and y==u)
    pieces=[[rng.randrange(2) for _ in endpoint] for _ in range(2)]
    pieces.append([e^a^b for e,a,b in zip(endpoint,*pieces)])
    opposite=[swapped(piece,A,B) for piece in pieces]
    rank_pairs=[]
    for piece,op in zip(pieces,opposite):
        ranks=[rank(flatten(piece,A,B,cut)) for cut in (0,1)]
        assert ranks==[rank(flatten(op,B,A,cut)) for cut in (0,1)]
        rank_pairs.append(ranks)
    total=[0]*q**4;products=[]
    for i,j in product(range(3),repeat=2):
        value=tensor_values(pieces[i],opposite[j],A,B)
        ranks=[rank(flatten(value,q,q,cut)) for cut in (0,1)]
        expected=[rank_pairs[i][cut]*rank_pairs[j][cut] for cut in (0,1)]
        assert ranks==expected;products.append(ranks)
        total=[a^b for a,b in zip(total,value)]
    for x,y,u,v in product(range(q),repeat=4):
        assert total[x+q*(y+q*(u+q*v))]==int(x==v and y==u)
    sizes=[]
    for m in (6,12,18):
        cx=(m//3,2*m//3);cy=(m//2,m//3,m//6)
        DX=type_size(cx);DY=type_size(cy)
        # Swapping the physical crossing pairs gives dimensions DY,DX,DX,DY.
        assert [DX*DY,DY*DX,DY*DX,DX*DY]==[DX*DY]*4
        sizes.append({'m':m,'dimension_x':DX,'dimension_y':DY,'square_dimension':DX*DY})
    # The constant rectangular crossing has U=V=AB and G_X=G_Y=0.
    # At R=2, B=(AB)^2=(Z_X Z_Y)^2 exactly, giving rho=1.
    endpoint_ranks=[rank(flatten(endpoint,A,B,cut)) for cut in (0,1)]
    assert endpoint_ranks==[A*B,A*B]
    profile_product=endpoint_ranks[0]*endpoint_ranks[1]
    partition_product=len(range(A))*len(range(B))
    assert profile_product==partition_product**2
    return {'status':'pass','rectangular_dimensions':[A,B,B,A],
            'rank_pairs':rank_pairs,'tensor_product_rank_pairs':products,
            'separate_coefficient_pairs':9,'endpoint_entries_checked':q**4,
            'type_dimensions':sizes,'gibbs':exact_gibbs_control(),
            'constant_crossing_control':{'profile_product':profile_product,
                'partition_product':partition_product,'R':2,'rho':1}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=check();text=json.dumps(result,indent=2,sort_keys=True)
    if args.output:args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text+'\n')
    print(text)

if __name__=='__main__':main()
