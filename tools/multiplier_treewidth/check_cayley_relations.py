#!/usr/bin/env python3
"""Angle 1: smooth relations as cycles in the class-group Cayley graph.

Vertices: reduced forms of discriminant D (the class group H, |H| = h).
Edges: f -- f*[p] for each split prime p <= B (prime form [p]).
A closed walk gives an exponent vector e in Z^r with prod [p_i]^{e_i} = 1,
i.e. a relation.  Checks, exactly:
 (C1) BFS-tree fundamental cycles generate the full relation lattice
      Lambda = ker(Z^r -> H): HNF determinant == h  (Schreier, KNOWN).
 (C2) cycle-space rank |E|-|V|+1 vs lattice rank r: almost all graph
      cycles are commutation cycles with exponent vector 0.
 (C3) shortest nonzero relation found by BFS collision, vs. the
      pigeonhole (Moore-type) bound 2L, L = min{L : #{|e|_1<=L} > h}.
 (C4) smallest B for which the graph is connected (= the primes generate H).
Uses ../smooth_relations_scope/bqf.py (exact form arithmetic).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'smooth_relations_scope'))
from bqf import reduced_forms, prime_form, compose, inverse, identity, primes_upto
from math import comb

def ball(r, L):  # #{e in Z^r : |e|_1 <= L}
    return sum(comb(r, k) * comb(L, k) * 2**k for k in range(min(r, L) + 1))

def hnf_det(vectors, r):
    rows = []
    for v in vectors:
        v = list(v)
        for i in range(r):
            if v[i] == 0: continue
            piv = next((R for R in rows if R[0] == i), None)
            if piv is None:
                if v[i] < 0: v = [-a for a in v]
                rows.append((i, v)); rows.sort(); v = None; break
            # gcd-combine v with pivot row
            _, P = piv
            while v[i] != 0:
                q = P[i] // v[i]
                P2 = [a - q * b for a, b in zip(P, v)]
                P, v = v, P2
            if P[i] < 0: P = [-a for a in P]
            rows = [(k, (P if k == i else R)) for k, R in rows]
        # v reduced to zero or inserted
    if len(rows) < r: return 0
    d = 1
    for i, R in rows: d *= R[i]
    return abs(d)

def analyse(D, B):
    H = reduced_forms(D); h = len(H)
    gens = []
    for p in primes_upto(B):
        f = prime_form(D, p)
        if f is not None: gens.append((p, f))
    r = len(gens)
    e0 = identity(D)
    # BFS
    vec = {e0: (0,) * r}; order = [e0]; fund = []; edges = 0; best = None
    frontier = [e0]
    while frontier:
        nxt = []
        for g in frontier:
            for i, (p, f) in enumerate(gens):
                for s, ff in ((1, f), (-1, inverse(f))):
                    h2 = compose(g, ff)
                    v = list(vec[g]); v[i] += s; v = tuple(v)
                    if h2 not in vec:
                        vec[h2] = v; nxt.append(h2)
                    else:
                        rel = tuple(a - b for a, b in zip(v, vec[h2]))
                        if any(rel):
                            fund.append(rel)
                            n1 = sum(abs(a) for a in rel)
                            if best is None or n1 < best: best = n1
        frontier = nxt
    connected = len(vec) == h
    E = h * r  # undirected edges f--f[p] (counting multi-edges/loops)
    det = hnf_det(fund, r) if connected else None
    L = 0
    while r > 0 and ball(r, L) <= h: L += 1
    return dict(D=D, h=h, r=r, connected=connected, det=det, cyc_rank=E - h + 1,
                best=best, pig=(2 * L if r else None), depth=None)

def min_B(D, Bmax=400):
    h = len(reduced_forms(D))
    for B in primes_upto(Bmax):
        a = analyse(D, B)
        if a['connected']: return B
    return None

if __name__ == '__main__':
    ok = True
    Ds = [-4*N for N in (15, 21, 33, 35, 77, 143, 221, 323, 437, 667, 899, 1147, 2021, 3127, 4087)] \
         + [-23, -47, -71, -199, -383, -991, -1999, -4027]
    print(f"{'D':>7} {'h':>4} {'Bconn':>5} {'r':>3} {'det==h':>6} {'cycRank':>7} {'shortRel':>8} {'pigeon2L':>8}")
    for D in Ds:
        B = min_B(D)
        a = analyse(D, max(B, 30))
        good = a['connected'] and a['det'] == a['h']
        ok &= good
        print(f"{D:7d} {a['h']:4d} {B:5d} {a['r']:3d} {str(good):>6} {a['cyc_rank']:7d} {str(a['best']):>8} {a['pig']:8d}")
        if a['best'] is not None and a['pig'] is not None and a['best'] > a['pig']: ok = False
    print("PASS" if ok else "FAIL")
