"""Exact finite ranks of the Appendix B pieces; no asymptotic claim.

Python standard library suffices. NumPy, when present, only speeds up the
integer matrix products; a checked bound prevents int64 overflow.

Default: weights 2, 3, 4 at prime 65521. --full reproduces all 13 recorded
cases, including two primes and the actual source's lexicographic prefixes.
The operator recurrence is Appendix B (O2), independently transcribed from
anthropics/formal-math@16a5ce41304fcdaea4e0be7c132e1c361dc26da5.
"""
from functools import lru_cache
from collections import Counter, defaultdict
from math import comb, factorial, isqrt
import json
try:
    import numpy as np
except ImportError:
    np = None


def zeros(rows, cols):
    return [[0] * cols for _ in range(rows)]


def matmul_mod(left, right, prime):
    shared = len(right)
    assert all(len(row) == shared for row in left)
    assert shared * (prime - 1) ** 2 < 2 ** 63
    if np is not None:
        return (np.asarray(left, dtype=np.int64) @
                np.asarray(right, dtype=np.int64) % prime).tolist()
    answer = zeros(len(left), len(right[0]))
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(right[k]):
                    answer[i][j] += value * other
        answer[i] = [value % prime for value in answer[i]]
    return answer

@lru_cache(None)
def partitions(n, maxpart=None):
    if n < 0:
        return ()
    if not n:
        return ((),)
    if maxpart is None:
        maxpart = n
    return tuple((k,) + p for k in range(min(n, maxpart), 0, -1)
                 for p in partitions(n-k, k))

def remove(p, x):
    p = list(p)
    p.remove(x)
    return tuple(p)

@lru_cache(None)
def opcoef(mu, out, inp):
    """[t^out] (t^mu)[wt(out)-wt(inp)] t^inp, using source (O2)."""
    if not mu:
        return int(out == inp)
    n, rest = mu[0], mu[1:]
    value = sum(comb(l-1,n-1) * opcoef(rest, remove(out,l), inp)
                for l in set(out) if l >= n)
    # binom(-k-1,n-1)=(-1)^(n-1)*binom(k+n-1,n-1)
    value += sum((-1)**(n-1) * comb(k+n-1,n-1) * k * mult
                 * opcoef(rest,out,remove(inp,k))
                 for k,mult in Counter(inp).items())
    return value

@lru_cache(None)
def afam(mu):
    terms = {(): 1}
    for n in mu:
        nxt = defaultdict(int)
        for p,c in terms.items():
            for j in range(1,n+1):
                nxt[tuple(sorted(p+(j,),reverse=True))] += c*(-1)**(n-1)*comb(n,j)
        terms = nxt
    return tuple(terms.items())

def bcoeff(h, p):
    z = 1
    for k,m in Counter(p).items():
        z *= k**m * factorial(m)
    assert factorial(h) % z == 0
    return factorial(h)//z

def modrank(a, prime):
    a = [[x % prime for x in row] for row in a]
    rows, cols = len(a), len(a[0])
    rank = 0
    for j in range(cols):
        piv = next((i for i in range(rank,rows) if a[i][j]),None)
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        inv = pow(a[rank][j], -1, prime)
        a[rank][j:] = [value * inv % prime for value in a[rank][j:]]
        for i in range(rank+1,rows):
            value = a[i][j]
            if value:
                a[i][j:] = [(x - value*y) % prime
                            for x,y in zip(a[i][j:],a[rank][j:])]
        rank += 1
        if rank == rows:
            break
    return rank

def run(h,prime=65521,take=None):
    assert h >= 1
    assert prime > h and all(prime % d for d in range(2, isqrt(prime)+1))
    ps = tuple(sorted(partitions(h)))
    if take is not None:
        assert 1 <= take <= len(ps)
        ps = ps[:take]
    dim = len(ps)
    dc = [bcoeff(h,p) for p in ps]
    # Pairing is independently checked from opcoef, before pieces.
    for x,px in enumerate(ps):
        for v,pv in enumerate(ps):
            pairing = dc[v]*sum(c*opcoef(mu,(),pv) for mu,c in afam(px))
            assert pairing == factorial(h)*int(x == v)
    total = zeros(dim*dim,dim*dim)
    ranks=[]
    for s in range(-2*h,2*h+1):
        channels=[(i,p) for i in range(2*h+1) for p in partitions(s+i)]
        left=zeros(dim*dim,len(channels))
        right=zeros(len(channels),dim*dim)
        for j,(i,sigma) in enumerate(channels):
            for x,px in enumerate(ps):
                ax=afam(px)
                for y,py in enumerate(ps):
                    left[x*dim+y][j] = ((-1)**i * sum(c*comb(sum(mu)+h,i)*opcoef(mu,sigma,py)
                                                    for mu,c in ax if sum(mu)+h>=i))%prime
            for u,pu in enumerate(ps):
                for v,pv in enumerate(ps):
                    right[j][u*dim+v] = dc[v]*opcoef(pv,pu,sigma)%prime
        a=matmul_mod(left,right,prime)
        total=[[(x+y)%prime for x,y in zip(trow,arow)] for trow,arow in zip(total,a)]
        col=[[a[x*dim+y][u*dim+v] for y in range(dim) for v in range(dim)]
             for x in range(dim) for u in range(dim)]
        if s == 0:
            # Independently form source (55), not merely a reshape of (54).
            cchannels=[(i,p) for i in range(2*h+1) for p in partitions(2*h-s-i)]
            cleft=zeros(dim*dim,len(cchannels))
            cright=zeros(len(cchannels),dim*dim)
            for j,(i,sigma) in enumerate(cchannels):
                for x,px in enumerate(ps):
                    for u,pu in enumerate(ps):
                        cleft[x*dim+u][j] = ((-1)**i*sum(c*comb(sum(mu)+h,i)*opcoef(mu,pu,sigma)
                                                       for mu,c in afam(px) if sum(mu)+h>=i))%prime
                for y,py in enumerate(ps):
                    for v,pv in enumerate(ps):
                        cright[j][y*dim+v]=dc[v]*opcoef(pv,sigma,py)%prime
            assert matmul_mod(cleft,cright,prime) == col
        ranks.append([s,modrank(a,prime),modrank(col,prime),len(channels)])
    for x in range(dim):
        for y in range(dim):
            for u in range(dim):
                for v in range(dim):
                    assert total[x*dim+y][u*dim+v] == factorial(h)*int(x==v and y==u)%prime
    return {'h':h,'dimension':dim,'prime':prime,'domain':'full_partitions' if take is None else 'lexicographic_prefix',
            'ranks':[{'s':s,'row':r,'col':c,'row_channel_bound':b} for s,r,c,b in ranks],
            'max_rank_product':max(r*c for _,r,c,_ in ranks),
            'identity_pass':True,'central_column_factor_pass':True}

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('weights', nargs='*', type=int)
    parser.add_argument('--prime', type=int, default=65521)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--no-numpy', action='store_true')
    args = parser.parse_args()
    if args.no_numpy:
        np = None
    assert opcoef((1,1),(1,1,1),(1,)) == 1
    assert opcoef((1,1),(3,),(1,)) == 2
    cases = [(h,args.prime,None) for h in (args.weights or [2,3,4])]
    if args.full:
        cases = [(h,p,None) for p in (65521,1000003) for h in range(2,7)]
        cases += [(h,65521,d) for h,d in ((2,2),(4,4),(6,8))]
    for h,prime,take in cases:
        print(json.dumps(run(h,prime,take),sort_keys=True),flush=True)
