#!/usr/bin/env python3
"""Self-dual Pfaffian physical codes and their uniform equal-copy blocks.

All arithmetic is over F2. Polynomials are integer bit masks. The two
closed representative formulas are justified in the accompanying proof.
"""
from functools import lru_cache
from itertools import combinations,product
from math import comb

from check_cancellation_module_automata import pmul


@lru_cache(None)
def edges(n):return tuple(combinations(range(n),2))


@lru_cache(None)
def matchings(n,subset):
    if not subset:return (0,)
    if subset.bit_count()%2:return ()
    low=subset & -subset;i=low.bit_length()-1;rest=subset ^ low;result=[]
    positions={e:k for k,e in enumerate(edges(n))}
    for j in range(i+1,n):
        if rest >> j & 1:
            result.extend(v | (1 << positions[i,j]) for v in matchings(n,rest ^ (1 << j)))
    return tuple(result)


def code(n):
    if not isinstance(n,int) or n<2 or n%2:raise ValueError('Use a positive even vertex count n>=2.')
    subsets=[s for s in range(1 << n) if s.bit_count()%2==0]
    vectors=[matchings(n,s) for s in subsets]
    return dict(n=n,modes=comb(n,2),q=1 << (n-1),subsets=subsets,
                embedding=vectors,retraction=vectors,heights=[s.bit_count()//2 for s in subsets])


@lru_cache(None)
def one_plus(exponent):
    return sum((comb(exponent,i)%2) << i for i in range(exponent+1))


def direct_block(a,delta):
    """Sum all original edge-matching minor contributions, no closed formula."""
    n=a+delta;intersection=(1 << a)-1;difference=((1 << delta)-1) << a
    parts=[s << a for s in range(1 << delta) if s.bit_count()%2==a%2]
    vectors=[list(product(matchings(n,intersection | s),
                          matchings(n,intersection | (difference ^ s)))) for s in parts]
    matrix=[]
    for row in vectors:
        line=[]
        for column in vectors:
            value=0
            for x,y in row:
                twice=x & y;once=x ^ y
                for u,v in column:
                    if (u & v)!=twice or (u ^ v)!=once:continue
                    # Columns of the compound have their two physical labels swapped.
                    same=once.bit_count()-((((x & u) | (y & v)) & once).bit_count())
                    value^=one_plus(same) << (once.bit_count()+4*twice.bit_count())
            line.append(value)
        matrix.append(line)
    return matrix


def closed_block(a,delta):
    if a<0 or delta<0 or delta%2:raise ValueError('Use a>=0 and nonnegative even delta.')
    if delta==0:return [[1 << (2*a)]] if a%2==0 else []
    parts=[s for s in range(1 << delta) if s.bit_count()%2==a%2];matrix=[]
    for x in parts:
        row=[]
        for y in parts:
            odd_intersection=(x & y).bit_count()%2
            distance=(x ^ y).bit_count()//2
            if a%2==0:
                value=0 if odd_intersection else one_plus(distance) << (delta//2+2*a)
            else:
                value=one_plus(distance+1-odd_intersection) << (delta//2+2*a-1)
            row.append(value)
        matrix.append(row)
    return matrix


def blocks(n):
    for a in range(n+1):
        for delta in range(0,n-a+1,2):
            if delta==0 and a%2:continue
            yield a,delta,comb(n,a)*comb(n-a,delta)


def seed():
    # k01=1+t, k10=t+t², k11=t after setting y=1, x=t.
    result=[0]*16;result[0]=1
    steps=[0,3,6,2]
    for i,j in combinations(range(4),2):result[(1 << i) | (1 << j)]=steps[i ^ j]
    result[15]=pmul(steps[1],steps[1]) ^ pmul(steps[2],steps[2]) ^ pmul(steps[3],steps[3])
    return result


def literal_tensor(n):
    """Full original Gaussian tensor after the four physical embeddings."""
    data=code(n);local=seed();values={}
    for labels in product(range(data['q']),repeat=4):
        value=0
        for occupations in product(*(data['embedding'][i] for i in labels)):
            term=1
            for edge in range(data['modes']):
                mask=sum(((v >> edge)&1) << leg for leg,v in enumerate(occupations))
                term=pmul(term,local[mask])
                if not term:break
            value^=term
        values[labels]=value
    return data,values


def subset_tensor(n):
    """Uniform 2^O(n) construction in a square-zero vertex algebra.

    This avoids enumerating factorially many matching monomials as n grows.
    Each state has 4n vertex bits. Every local edge factor has eight terms.
    """
    local=seed();states={0:1}
    for u,v in edges(n):
        factors=[]
        for mask,polynomial in enumerate(local):
            if polynomial:
                vertices=sum(((1 << (leg*n+u)) | (1 << (leg*n+v)))
                             for leg in range(4) if mask >> leg & 1)
                factors.append((vertices,polynomial))
        updated={}
        for state,value in states.items():
            for vertices,polynomial in factors:
                if state & vertices:continue
                target=state | vertices
                updated[target]=updated.get(target,0)^pmul(value,polynomial)
        states={state:value for state,value in updated.items() if value}
    subsets=[s for s in range(1 << n) if s.bit_count()%2==0]
    return {labels:states.get(sum(subsets[index] << (leg*n)
                    for leg,index in enumerate(labels)),0)
            for labels in product(range(1 << (n-1)),repeat=4)}


def reverse(polynomial,degree):
    assert polynomial.bit_length()-1<=degree
    return sum(((polynomial >> i)&1) << (degree-i) for i in range(degree+1))
