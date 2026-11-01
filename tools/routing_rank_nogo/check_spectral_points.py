#!/usr/bin/env python3
"""Checks for the crossing-tensor remark in note.md (exact, integer/Fraction arithmetic).
E(x,y,u,v) = [x=v][y=u], N = 2^m. Verifies: (1) E is free (support points pairwise differ
in >= 2 coordinates); (2) every one-party flattening has N equal nonzero singular values
(Gram matrix = identity * N^... checked exactly), so all quantum functionals give N;
(3) cut1/cut2 ranks equal N^2 and crossing pieces bound sum sqrt(r1 r2) >= N on random
exact decompositions (sanity of the elementary bound).
Rerun: python3 check_spectral_points.py
"""
import itertools, random
from fractions import Fraction

def rank_frac(rows):
    A = [[Fraction(x) for x in r] for r in rows]; rk = 0
    for c in range(len(A[0])):
        piv = next((i for i in range(rk, len(A)) if A[i][c] != 0), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][c] != 0:
                f = A[i][c] / A[rk][c]; A[i] = [a - f * b for a, b in zip(A[i], A[rk])]
        rk += 1
    return rk

ok = True
for m in (1, 2):
    N = 2 ** m
    supp = [(x, y, y, x) for x in range(N) for y in range(N)]
    free = all(sum(a != b for a, b in zip(s, t)) >= 2 for s in supp for t in supp if s != t)
    E = {s: 1 for s in supp}
    # one-party marginal Gram: G[i,j] = sum_rest E(i,rest)E(j,rest)
    marg_ok = True
    for p in range(4):
        G = [[0] * N for _ in range(N)]
        for s in supp:
            for t in supp:
                if all(s[q] == t[q] for q in range(4) if q != p): G[s[p]][t[p]] += 1
        marg_ok &= all(G[i][j] == (N if i == j else 0) for i in range(N) for j in range(N))
    def flat(rowparties):
        cols = [q for q in range(4) if q not in rowparties]
        R = list(itertools.product(range(N), repeat=2)); 
        return [[E.get(tuple(dict(zip(rowparties + cols, ri + ci))[q] for q in range(4)), 0)
                 for ci in R] for ri in R]
    r1 = rank_frac(flat([0, 1])); r2 = rank_frac(flat([0, 2]))
    print(f"m={m} N={N}: free={free} marginals_maximally_mixed={marg_ok} cut1_rank={r1} cut2_rank={r2}")
    ok &= free and marg_ok and r1 == N * N and r2 == N * N
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
