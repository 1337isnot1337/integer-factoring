#!/usr/bin/env python3
"""Exact fixed-curve quartic jets and a common-type profile floor; no sweep."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
EXPECTED=ROOT/'evidence/fixed-curve-quartic-floor.json'
import cube_homogeneous_gaussian as cube
from check_quartic_charge_fusion import ELL2, QUARTIC

PRECISION=19
MASK=(1<<PRECISION)-1

def eval_curve(poly):
    out=0
    for a,b,c in poly:
        out ^= 1<<(2*a+c)
    return out

@lru_cache(None)
def pmul(a,b):
    out=0
    while b:
        bit=b&-b
        out ^= a<<(bit.bit_length()-1)
        b ^= bit
    return out&MASK

def val(p):
    return (p&-p).bit_length()-1 if p else None

def charge(word):
    out=0
    for a in word: out^=a
    return out

kernel=[[cube.multiply(e,ELL2) for e in row] for row in cube.kernel()]
left,right=cube.axis_parts(0)
pair=[[eval_curve(kernel[i][j])&MASK for j in right] for i in left]
h=eval_curve(QUARTIC)&MASK
assert min(val(p) for row in pair for p in row if p)==6 and val(h)==12
assert all(sum(e)==8 for row in kernel for p in row for e in p)
assert all(sum(e)==16 for e in QUARTIC)
# Exact complete one-particle minors before truncation.
cross=[[kernel[i][j] for j in right] for i in left]
prefix=[]
for size in range(1,5):
    values=[val(eval_curve(cube.minor(cross,r,c)))
            for r in combinations(range(4),size) for c in combinations(range(4),size)]
    prefix.append(min(v for v in values if v is not None))
orders=[prefix[0]]+[prefix[i]-prefix[i-1] for i in range(1,4)]
assert orders==[6,6,18,18],orders
# Endpoint and fixed superletter incidence, not adaptive reindexing.
assert all((len(kernel[i][j])%2)==(i^j==7) for i in range(8) for j in range(8))
assert len(QUARTIC)%2==0
assert all(((v&3)^((v^7)&3))==3 and ((v>>2)^((v^7)>>2))==1 for v in range(8))


def entry(aa,bb,fusion):
    total=len(aa)+len(bb)
    if not total: return 1
    if not aa or not bb or total>6: return 0
    out=0
    if len(aa)==len(bb):
        for perm in permutations(bb):
            if any(c!=d for (_,c),(_,d) in zip(aa,perm)): continue
            p=1
            for (i,_),(j,_) in zip(aa,perm): p=pmul(p,pair[i][j])
            out^=p
    if fusion and total==4:
        colors=[c for _,c in aa+bb]
        if len(set(colors))==4 and charge(colors)==0: out^=h
    if fusion and total==6:
        for i,c in aa:
            for j,d in bb:
                if c!=d: continue
                colors=[v for k,v in aa if k!=i]+[v for k,v in bb if k!=j]
                if len(set(colors))==4 and charge(colors)==0:
                    out^=pmul(h,pair[i][j])
    return out


def matrices(target,fusion):
    words=[w for w in product(range(8),repeat=4) if charge(w)==target]
    assert len(words)==512
    active=[tuple((i,c) for i,c in enumerate(w) if c) for w in words]
    blocks=[[0]*len(words) for _ in range(PRECISION)]
    entries=0
    for i,aa in enumerate(active):
        for j,bb in enumerate(active):
            p=entry(aa,bb,fusion)
            entries+=bool(p)
            while p:
                bit=p&-p; degree=bit.bit_length()-1
                blocks[degree][i]|=1<<j
                p^=bit
    return blocks,dict(dimension=len(words),occupancy=dict(Counter(len(a) for a in active)),nonzero_jet_entries=entries)


def ranks(blocks):
    n=len(blocks[0]); pivots={}; result=[0]
    for r in range(PRECISION):
        for i in range(n):
            row=0
            for d in range(r+1): row |= blocks[d][i]<<(n*(r-d))
            while row:
                bit=row.bit_length()-1
                if bit in pivots: row^=pivots[bit]
                else:
                    pivots[bit]=row
                    break
        result.append(len(pivots))
    cumulative=[result[i+1]-result[i] for i in range(PRECISION)]
    profile={i:cumulative[i]-(cumulative[i-1] if i else 0) for i in range(PRECISION)}
    assert all(v>=0 for v in profile.values())
    return dict(toeplitz_ranks=result,smith_prefix={i:v for i,v in profile.items() if v})

def column_rank(blocks,precision):
    """Independent column construction and low-pivot elimination at one truncation."""
    n=len(blocks[0])
    transposed=[[0]*n for _ in range(precision)]
    for d in range(precision):
        for i,row in enumerate(blocks[d]):
            while row:
                bit=row&-row; j=bit.bit_length()-1
                transposed[d][j]|=1<<i
                row^=bit
    pivots={}
    for c in range(precision):
        for j in range(n):
            column=sum(transposed[d][j]<<(n*(c+d)) for d in range(precision-c))
            while column:
                bit=column&-column
                if bit in pivots: column^=pivots[bit]
                else:
                    pivots[bit]=column
                    break
    return len(pivots)


def literal_creation_controls():
    """Quartic-support enumeration versus the fast matching formula.

    All placements of one zero-XOR quartet, with every possible extra
    equal-color cross pair, plus the complete binary occupied-word panel.
    The seven q8 zero-XOR quartets are one GL(3,2) color orbit.
    """
    edge={(a,b):pair[i][j] for i,a in enumerate(left) for j,b in enumerate(right)}
    edge.update({(b,a):value for (a,b),value in list(edge.items())})
    @lru_cache(None)
    def pfaffian(active):
        if not active: return 1
        a,color=active[0]; out=0
        for i,(b,other) in enumerate(active[1:],1):
            if color==other and (a,b) in edge:
                out^=pmul(edge[a,b],pfaffian(active[1:i]+active[i+1:]))
        return out
    def literal(active,fusion):
        value=pfaffian(active)
        if fusion:
            for indices in combinations(range(len(active)),4):
                quartet=[active[i] for i in indices]
                colors=[c for _,c in quartet]
                if (len(set(colors))==4 and charge(colors)==0
                        and any(v in left for v,_ in quartet)
                        and any(v in right for v,_ in quartet)):
                    rest=tuple(v for i,v in enumerate(active) if i not in indices)
                    value^=pmul(h,pfaffian(rest))
        return value
    panel=set(product(range(2),repeat=8))
    for support in combinations(range(8),4):
        if not (set(support)&set(left) and set(support)&set(right)): continue
        for colors in permutations((1,2,4,7)):
            word=[0]*8
            for v,c in zip(support,colors): word[v]=c
            panel.add(tuple(word))
            for a in set(left)-set(support):
                for b in set(right)-set(support):
                    for c in range(1,8):
                        extended=word.copy(); extended[a]=extended[b]=c
                        panel.add(tuple(extended))
    for word in panel:
        aa=tuple((i,word[v]) for i,v in enumerate(left) if word[v])
        bb=tuple((j,word[v]) for j,v in enumerate(right) if word[v])
        active=tuple((v,c) for v,c in enumerate(word) if c)
        for fusion in (False,True):
            assert entry(aa,bb,fusion)==literal(active,fusion)
    return 2*len(panel)


def pressure_certificate(combined):
    n=8**4
    assert combined==Counter({0:1,6:14,12:91,18:254})
    degree=16*(4*7*8**3)
    known=sum(combined.values())
    weighted=sum(k*v for k,v in combined.items())
    tail=n-known
    mean=Fraction(degree-weighted,tail)
    assert degree==56*n and mean<60
    # x=u^3. Subtract (1+7*u^32)^3 and factor u^18.
    coefficients={0:14,14:-21,18:91,36:254,46:-147,78:-343,162:tail}
    degree=200
    bernstein=[sum(Fraction(a*comb(k,j),comb(degree,j))
                   for j,a in coefficients.items() if j<=k) for k in range(degree+1)]
    assert min(bernstein)>0
    for j in range(degree+1):
        reconstructed=comb(degree,j)*sum(
            (-1)**(j-k)*comb(j,k)*bernstein[k] for k in range(j+1))
        assert reconstructed==coefficients.get(j,0)
    encoded=[[x.numerator,x.denominator] for x in bernstein]
    assert Fraction(3,4)**2>Fraction(1,2)
    return dict(whole_dimension=n,known_orders=known,known_order_sum=weighted,
                determinant_degree=56*n,tail_count=tail,tail_mean_upper=str(mean),
                tail_order_used=60,bernstein_degree=degree,
                factored_difference=coefficients,positive_coefficients=len(bernstein),
                minimum_coefficient=str(min(bernstein)),
                minimum_index=bernstein.index(min(bernstein)),
                coefficient_sha256=hashlib.sha256(json.dumps(encoded).encode()).hexdigest(),
                profile_floor='3/4',target='1/sqrt(2)')


def compute():
    record=dict(scope='Fixed (x^2,y^2,xy) quartic cube at q8; exact generic binary jets and an all-common-paired-type profile-envelope floor. No complete profile, coefficient-rank lower bound or factoring improvement.',
                precision=PRECISION,one_particle_orders=orders,
                registered_base='50b01ff3cce5138532c476fe856826e8bfcf1ef1',cases=[])
    combined=Counter()
    for fusion in (False,True):
        for target in (0,1):
            blocks,meta=matrices(target,fusion)
            result=ranks(blocks)
            checked={p:column_rank(blocks,p) for p in (7,13,18,19)}
            assert all(result['toeplitz_ranks'][p]==rank for p,rank in checked.items())
            case=dict(fusion=fusion,charge=target,**meta,**result,
                      independent_column_ranks=checked)
            record['cases'].append(case)
            if fusion:
                for k,v in result['smith_prefix'].items():
                    combined[k]+=(1 if target==0 else 7)*v
    record['literal_creation_controls']=literal_creation_controls()
    record['pressure_certificate']=pressure_certificate(combined)
    record['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                             for name in ['tools/check_fixed_curve_quartic_floor.py',
                                          'tools/check_quartic_charge_fusion.py',
                                          'tools/cube_homogeneous_gaussian.py']}
    return json.loads(json.dumps(record))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=compute()
    if args.write: EXPECTED.write_text(json.dumps(result,indent=2)+'\n')
    else: assert result==json.loads(EXPECTED.read_text())
    print('Fixed-curve quartic: full charge jets, literal/column checks and Bernstein 3/4 floor PASS')
