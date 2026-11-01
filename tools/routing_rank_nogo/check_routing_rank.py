#!/usr/bin/env python3
"""Checker for Theorem A (routing-rank lower bound) in note.md.
Exact arithmetic only: integers, Fractions, and ranks over GF(p).
Rerun: python3 check_routing_rank.py
"""
import itertools, sys
from fractions import Fraction
from math import comb

P = 1_000_003  # prime for modular ranks

def rank_mod(rows, p=P):
    A = [r[:] for r in rows]; n = len(A); m = len(A[0]) if n else 0
    rk = 0; col = 0
    for col in range(m):
        piv = next((i for i in range(rk, n) if A[i][col] % p), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][col], p - 2, p)
        A[rk] = [(x * inv) % p for x in A[rk]]
        for i in range(n):
            if i != rk and A[i][col] % p:
                f = A[i][col]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk

def ham(x, y): return sum(a != b for a, b in zip(x, y))

def coords(r, Q):
    # coordinate u = (d, a) with d, a in [r]^Q
    words = list(itertools.product(range(r), repeat=Q))
    return words, [(d, a) for d in words for a in words]

# ---------- free-set lower bound (Theorem A) ----------
def free_bound(r, Q, k, key):
    words = list(itertools.product(range(r), repeat=Q))
    order = sorted(words, key=key); pos = {w: i for i, w in enumerate(order)}
    cnt = 0
    for delta in words:
        mx = max(pos[g] for g in words if ham(g, delta) <= k)
        cnt += sum(1 for a in words if pos[a] > mx)
    return cnt

def score_key(r):
    H = set(range(r // 2))
    return lambda w: (sum(x in H for x in w), w)

def prob_diff_gt(r, Q, k):
    """Exact Pr[S - S' > k], S,S' iid Bin(Q, floor(r/2)/r)."""
    q = Fraction(r // 2, r)
    pmf = [comb(Q, s) * q**s * (1 - q)**(Q - s) for s in range(Q + 1)]
    return sum(pmf[s] * pmf[t] for s in range(Q + 1) for t in range(Q + 1) if s - t > k)

# ---------- the construction M = S g(W0) ----------
def interp_poly(I):
    """Coefficients (Fractions) of the degree <= |I|-1 polynomial g with g(j)=(-1)^j on I."""
    pts = list(I); coeffs = [Fraction(0)] * len(pts)
    for i, xi in enumerate(pts):
        num = [Fraction(1)]; den = Fraction(1)
        for j, xj in enumerate(pts):
            if j == i: continue
            num = [Fraction(0)] + num  # multiply by x
            for t in range(len(num) - 1): num[t] -= xj * num[t + 1]
            den *= (xi - xj)
        for t in range(len(num)): coeffs[t] += Fraction((-1) ** xi) * num[t] / den
    return coeffs

def construction(r, Q, I):
    words, U = coords(r, Q); D = len(U); idx = {u: i for i, u in enumerate(U)}
    def swap_pos(u, t):
        d, a = list(u[0]), list(u[1]); d[t], a[t] = a[t], d[t]; return (tuple(d), tuple(a))
    # W0 = sum_t (I - Pi_t)/2 as a sparse dict-of-dicts with Fractions
    W = [dict() for _ in range(D)]
    for i, u in enumerate(U):
        for t in range(Q):
            W[i][i] = W[i].get(i, 0) + Fraction(1, 2)
            j = idx[swap_pos(u, t)]
            W[i][j] = W[i].get(j, 0) - Fraction(1, 2)
    def matvec_rows(A, B):  # A,B sparse row dicts -> A@B
        out = [dict() for _ in range(D)]
        for i in range(D):
            for k_, x in A[i].items():
                for j, y in B[k_].items():
                    out[i][j] = out[i].get(j, 0) + x * y
            out[i] = {j: v for j, v in out[i].items() if v != 0}
        return out
    c = interp_poly(I)
    G = [dict() for _ in range(D)]; Pw = [{i: Fraction(1)} for i in range(D)]
    for deg, ci in enumerate(c):
        if deg > 0: Pw = matvec_rows(Pw, W)
        for i in range(D):
            for j, v in Pw[i].items(): G[i][j] = G[i].get(j, 0) + ci * v
    # M = S G, S = full swap (d,a)->(a,d)
    M = [dict() for _ in range(D)]
    for i, u in enumerate(U):
        si = idx[(u[1], u[0])]
        M[i] = {j: v for j, v in G[si].items() if v != 0}
    return U, M

def check_construction(r, Q, I):
    U, M = construction(r, Q, I); D = len(U); k = len(I) - 1
    # support: M[u,w] != 0 only if dist(d_u, a_w) <= k
    bad = sum(1 for i in range(D) for j in M[i] if ham(U[i][0], U[j][1]) > k)
    rows = []
    for i in range(D):
        row = [0] * D
        row[i] = 1
        for j, v in M[i].items():
            row[j] = (row[j] - v.numerator * pow(v.denominator, P - 2, P)) % P
        rows.append(row)
    rk = rank_mod(rows)
    a, s = r * (r - 1) // 2, r * (r + 1) // 2
    # kernel of I - S g(W0) = sectors j (j antisymmetric positions) with g(j) = (-1)^j
    c = interp_poly(I)
    fixed = [j for j in range(Q + 1) if sum(ci * j**e for e, ci in enumerate(c)) == (-1) ** j]
    formula = D - sum(comb(Q, j) * a**j * s**(Q - j) for j in fixed)
    lb = free_bound(r, Q, k, score_key(r))
    lb_lex = free_bound(r, Q, k, lambda w: w)
    return bad, rk, formula, lb, lb_lex

def exhaustive_min_rank_r2Q1(p):
    # r=2,Q=1,k=0: D=4, M[u,w] allowed iff d_u == a_w
    words, U = coords(2, 1)
    allowed = [(i, j) for i in range(4) for j in range(4) if U[i][0] == U[j][1]]
    best = 99
    for vals in itertools.product(range(p), repeat=len(allowed)):
        M = [[0] * 4 for _ in range(4)]
        for (i, j), v in zip(allowed, vals): M[i][j] = v
        rows = [[((1 if i == j else 0) - M[i][j]) % p for j in range(4)] for i in range(4)]
        best = min(best, rank_mod(rows, p))
    return best

def main():
    ok = True
    print("== k=0: lex free-set bound vs r^Q(r^Q-1)/2 ==")
    for r, Q in [(2, 1), (2, 2), (3, 1), (3, 2), (2, 3), (4, 1)]:
        lb = free_bound(r, Q, 0, lambda w: w); N = r**Q; exp = N * (N - 1) // 2
        print(f"r={r} Q={Q}: bound={lb} expected={exp}"); ok &= lb == exp
    print("== exhaustive min rank(I-M), r=2 Q=1 k=0 (expect 1) ==")
    for p in (2, 3):
        b = exhaustive_min_rank_r2Q1(p); print(f"GF({p}): min rank = {b}"); ok &= b == 1
    print("== construction M=S g(W0): support, rank (GF(p)), eigen formula, lower bounds ==")
    for r, Q, I in [(2, 1, [0]), (2, 2, [0]), (2, 2, [0, 1]), (2, 3, [0, 1]), (2, 3, [1, 2]),
                    (3, 2, [0, 1]), (2, 4, [1, 2]), (2, 4, [0, 1, 2]), (4, 2, [0, 1])]:
        bad, rk, f, lb, lbl = check_construction(r, Q, I)
        print(f"r={r} Q={Q} I={I} k={len(I)-1}: support_violations={bad} rank={rk} formula={f} "
              f"score_LB={lb} lex_LB={lbl} D*Pr[S-S'>k]={prob_diff_gt(r,Q,len(I)-1)*r**(2*Q)}")
        ok &= bad == 0 and rk == f and lb <= rk and lbl <= rk
        ok &= lb >= prob_diff_gt(r, Q, len(I) - 1) * r**(2 * Q)
    print("== asymptotic comparison: smallest k with rank/D <= eps ==")
    print("   LB: need Pr[S-S'>k] <= eps (any design);  construction: best R-window of Bin(Q,(r-1)/2r)")
    for r, Q in [(2, 64), (2, 256), (4, 256), (16, 256), (2, 1024)]:
        q = Fraction(r - 1, 2 * r)
        for eps in (Fraction(1, 10**3), Fraction(1, 10**6)):
            # lower bound: S-S' distribution
            qq = Fraction(r // 2, r)
            pm = [comb(Q, s) * qq**s * (1 - qq)**(Q - s) for s in range(Q + 1)]
            diff = {}
            for s in range(Q + 1):
                for t in range(Q + 1): diff[s - t] = diff.get(s - t, 0) + pm[s] * pm[t]
            tail = 0; kLB = None
            for kk in range(Q, -1, -1):
                tail += diff.get(kk + 1, 0)
                if tail > eps: kLB = kk + 1; break
            # construction: minimal R so some window of R consecutive j has mass >= 1-eps
            pj = [comb(Q, j) * q**j * (1 - q)**(Q - j) for j in range(Q + 1)]
            Rc = None
            for R in range(1, Q + 2):
                if max(sum(pj[s:s + R]) for s in range(Q + 2 - R)) >= 1 - eps: Rc = R; break
            print(f"r={r} Q={Q} eps={float(eps):.0e}: k_LB>={kLB}  k_constr={Rc-1}  ratio={(Rc-1)/kLB:.3f}")
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
