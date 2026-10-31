#!/usr/bin/env python3
"""Optional certificate search for research/bare-seven-to-twelve-letters.md.

This FLOATING-POINT search (requires numpy and scipy) found the certificates
stored in evidence/bare-seven-to-twelve-letters.json.  It plays no role in
any proof: tools/check_bare_seven_to_twelve_letters.py verifies the stored
certificates exactly with the standard library.  Usage:

    python3 tools/generate_bare_seven_to_twelve_certificates.py OUT.json [q ...]

then copy the "certificates" and "thirteen_target" fields into the record
and run the checker with --write.  The certificates use the box budgets of
Theorem 1 only.  The integer witnesses come from a second solve that also
imposes both square-box adjugate inequalities, so that they satisfy every
constraint proved in the note.

The search is environment dependent: SLSQP and HiGHS results vary with the
numpy/scipy build and thread count, so a rerun generally produces different
(but, after the exact repair below, equally valid) multipliers and
witnesses.  Pin OMP_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1 for the most
stable behaviour.  A full run over q = 7..13 takes from about one to about
five minutes; run single alphabets (extra arguments) under a timeout.
"""
import json
import math
import sys
import time

try:
    import numpy as np
    from scipy.optimize import linprog, minimize
except ImportError:  # pragma: no cover - optional tool
    sys.exit("this optional search needs numpy and scipy")
from fractions import Fraction as Fr

def constraints(q):
    n = q * q
    cons = []
    for al in range(1, q + 1):
        for be in range(al, q + 1):
            N = al * be; r = al * al
            c = [max(0, N - j * r) for j in range(n + 1)]
            w = [c[i - 1] - c[i] for i in range(1, n + 1)]
            v = [1 if i <= N else 0 for i in range(1, n + 1)]
            cons.append((('C', al, be), w, v, N))
            cons.append((('R', al, be), v, w, N))
    return cons


def tails(vec):
    out = []; s = 0
    for x in reversed(vec):
        s += x; out.append(s)
    return list(reversed(out))  # out[i] = sum_{j>=i}


def matrices(q):
    n = q * q
    cons = constraints(q)
    # delta coords k=2..n (index 1..n-1 of tails)
    A = np.zeros((len(cons), 2 * (n - 1)))
    hi = np.zeros(len(cons))
    for c, (name, wm, wn, N) in enumerate(cons):
        A[c, :n - 1] = tails(wm)[1:]
        A[c, n - 1:] = tails(wn)[1:]
        hi[c] = N
    return cons, A, hi


def adjugate_rows(q):
    """Square-box adjugate inequalities (note, Section 5), normalized M = 1.

    (N-1) a_2 + N sum_{i<N} b_i <= N (N-1)   and the mirror with a and b
    exchanged, for N = alpha^2, alpha = 2..q, in increment coordinates.
    """
    n = q * q
    rows, rhs = [], []
    for al in range(2, q + 1):
        N = al * al
        wm = [0] * n
        wm[1] = N - 1
        wn = [N if i < N - 1 else 0 for i in range(n)]
        for x, y in ((wm, wn), (wn, wm)):
            rows.append(np.array(tails(x)[1:] + tails(y)[1:], dtype=float))
            rhs.append(float(N * (N - 1)))
    return np.vstack(rows), np.array(rhs)


def full_matrices(q):
    """Box budgets of Theorem 1 together with the adjugate inequalities."""
    cons, A, hi = matrices(q)
    Aa, ha = adjugate_rows(q)
    return np.vstack([A, Aa]), np.concatenate([hi, ha])


def lse(x):
    m = x.max(); return m + math.log(np.exp(x - m).sum())


def solve_t(q, t, A, hi, x0=None):
    n = q * q
    def prof(x):
        a = np.concatenate([[0.0], np.cumsum(x[:n - 1])])
        b = np.concatenate([[0.0], np.cumsum(x[n - 1:])])
        return a, b
    def f(x):
        a, b = prof(x); return lse(-t * a) + lse(-t * b)
    def g(x):
        a, b = prof(x)
        pa = np.exp(-t * a - lse(-t * a)); pb = np.exp(-t * b - lse(-t * b))
        ta = np.array(tails(list(pa)))[1:]; tb = np.array(tails(list(pb)))[1:]
        return np.concatenate([-t * ta, -t * tb])
    if x0 is None:
        x0 = np.full(2 * (n - 1), 1.0 / n)
    cons = [{'type': 'ineq', 'fun': lambda x: hi - A @ x, 'jac': lambda x: -A}]
    res = minimize(f, x0, jac=g, bounds=[(0, None)] * (2 * (n - 1)), constraints=cons,
                   method='SLSQP', options={'maxiter': 3000, 'ftol': 1e-15})
    return t + res.fun, res.x


def golden(q, A, hi, lo=3.0, up=6.0, it=30):
    cache = {}; x = [None]
    def F(t):
        v, xx = solve_t(q, t, A, hi, x[0]); x[0] = xx; cache[t] = (v, xx); return v
    g = (math.sqrt(5) - 1) / 2; a, b = lo, up; c = b - g * (b - a); d = a + g * (b - a)
    fc = F(c); fd = F(d)
    for _ in range(it):
        if fc < fd:
            b, d, fd = d, c, fc; c = b - g * (b - a); fc = F(c)
        else:
            a, c, fc = c, d, fd; d = a + g * (b - a); fd = F(d)
    t = min(cache, key=lambda k: cache[k][0])
    return cache[t][0], t, cache[t][1]


def rational_prob(p, D):
    k = [int(round(x * D)) for x in p]
    k = [max(1, x) for x in k]
    diff = D - sum(k)
    j = int(np.argmax(k)); k[j] += diff
    assert min(k) > 0 and sum(k) == D
    return k


def make_cert(q, D=10**7, verbose=True):
    n = q * q
    cons, A, hi = matrices(q)
    t0 = time.time()
    val, t, x = golden(q, A, hi)
    a = np.concatenate([[0.0], np.cumsum(x[:n - 1])]); b = np.concatenate([[0.0], np.cumsum(x[n - 1:])])
    if verbose:
        print('q', q, 'primal rho %.7f t %.4f (%.1fs)' % (val / (4 * math.log(q)), t, time.time() - t0), flush=True)
        print('  a-b max diff', np.abs(a - b).max())
    s = 0.5 * (a + b)
    p = np.exp(-t * s - lse(-t * s))
    k = rational_prob(p, D)
    pr = np.array(k, dtype=float) / D
    T = np.array(tails(list(pr)))[1:]
    c = np.concatenate([T, T])
    res = linprog(-c, A_ub=A, b_ub=hi, bounds=[(0, None)] * len(c), method='highs')
    lam = -res.ineqlin.marginals
    if verbose:
        print('  LP max', -res.fun, ' #positive duals', int((lam > 1e-12).sum()))
    lamr = {}
    for ci, (name, wm, wn, N) in enumerate(cons):
        if lam[ci] > 1e-13:
            lamr[name] = Fr(float(lam[ci])).limit_denominator(10**12)
    return dict(q=q, D=D, k=k, lam=lamr, t=t, a=list(a), primal=val), cons


def repair_and_pack(q, cert, cons):
    """Make the rationalized dual exactly feasible by raising the two determinant multipliers."""
    n = q * q; D = cert['D']; k = cert['k']
    lam = dict(cert['lam'])
    Pi = [Fr(x, D) for x in tails(k)]  # Pi[i] = sum_{j>=i} pi_j (0-indexed)
    tm = {c[0]: tails(c[1]) for c in cons}; tn = {c[0]: tails(c[2]) for c in cons}
    for side, key in (('mu', ('R', q, q)), ('nu', ('C', q, q))):
        tt = tm if side == 'mu' else tn
        worst = Fr(0)
        for i in range(1, n):
            lhs = sum(v * tt[nm][i] for nm, v in lam.items())
            d = Pi[i] - lhs
            if d > 0:
                worst = max(worst, d / tt[key][i])
        if worst > 0:
            lam[key] = lam.get(key, Fr(0)) + worst
    return lam


def integer_witness(q, a, M):
    # floor of a slightly shrunk feasible float profile; all constraint coefficients are >= 0
    return [int(math.floor(x * (1 - 1e-7) * M)) for x in a]


def witness_profile(q, M=10**6, lo=3.0, up=6.0, it=30):
    """Integer profile satisfying every proved constraint (box + adjugate)."""
    n = q * q
    A, hi = full_matrices(q)
    val, t, x = golden(q, A, hi, lo=lo, up=up, it=it)
    a = np.concatenate([[0.0], np.cumsum(x[:n - 1])])
    b = np.concatenate([[0.0], np.cumsum(x[n - 1:])])
    s = 0.5 * (a + b)  # the constraints and objective are symmetric and convex
    return dict(M=M, mu=integer_witness(q, list(s), M), t=repr(t)), val


def main():
    qs = [int(x) for x in sys.argv[2:]] or list(range(7, 14))
    out = {"certificates": {}}
    for q in qs:
        if q == 13:
            continue
        t0 = time.time()
        cert, cons = make_cert(q)
        lam = repair_and_pack(q, cert, cons)
        wit, val = witness_profile(q)
        print('  q', q, 'box+adjugate primal rho %.7f (%.1fs total)'
              % (val / (4 * math.log(q)), time.time() - t0), flush=True)
        out["certificates"][str(q)] = dict(
            D=cert["D"], k=cert["k"],
            multipliers=[[nm[0], nm[1], nm[2], str(v)] for nm, v in sorted(lam.items())],
            witness=wit)
    if 13 in qs:
        q, n = 13, 169
        cons, A, hi = matrices(q)
        val, t, x = golden(q, A, hi, lo=4.2, up=5.4, it=22)
        a = np.concatenate([[0.0], np.cumsum(x[:n - 1])])
        b = np.concatenate([[0.0], np.cumsum(x[n - 1:])])
        s = 0.5 * (a + b)
        M = 2000
        out["thirteen_target"] = dict(M=M, mu=[int(math.floor(v * (1 - 1e-7) * M)) for v in s],
                                      s=str(Fr(math.exp(t / M)).limit_denominator(10**6)))
    with open(sys.argv[1], "w") as handle:
        json.dump(out, handle)


if __name__ == "__main__":
    main()
