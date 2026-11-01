#!/usr/bin/env python3
"""Exact check of Lemma B (note.md): for q = floor(r/2)/r, S,S' iid Bin(Q,q), Q>=64, 1<=R<=Q/4,
   Pr[S - S' >= R] >= (1/2) exp(-9R^2/(2Q) - 72),
and of the sharper per-factor window bound Pr[Bin(Q,p) >= ceil(Qp+R/2)] >= 2^{-1/2} exp(-9R^2/(4Q)-36)
for p in {q, 1-q}. Probabilities are exact Fractions; exp compared via mpmath-free float of a
lower bound (the right sides are tiny, so float comparison has huge margin; we also print log ratios).
Rerun: python3 check_tail_lemma.py
"""
from fractions import Fraction
from math import comb, exp, log, ceil
ok = True; worst = 1e9
for r in (2, 3, 4, 5, 16):
    q = Fraction(r // 2, r)
    for Q in (64, 100, 256):
        pm = [comb(Q, s) * q**s * (1 - q)**(Q - s) for s in range(Q + 1)]
        cdf_tail = {}
        for R in range(1, Q // 4 + 1):
            pr = sum(pm[s] * pm[t] for s in range(Q + 1) for t in range(Q + 1) if s - t >= R)
            lhs = log(pr.numerator) - log(pr.denominator)
            rhs = log(0.5) - 9 * R * R / (2 * Q) - 72
            for p in (q, 1 - q):
                A = ceil(Q * p + Fraction(R, 2))
                tp = sum(comb(Q, k) * p**k * (1 - p)**(Q - k) for k in range(A, Q + 1))
                l2 = log(tp.numerator) - log(tp.denominator)
                ok &= l2 >= -0.5 * log(2) - 9 * R * R / (4 * Q) - 36
            ok &= lhs >= rhs; worst = min(worst, lhs - rhs)
        print(f"r={r} Q={Q}: checked R=1..{Q//4}")
print(f"min log-margin = {worst:.2f}")
print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
