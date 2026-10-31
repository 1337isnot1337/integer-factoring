#!/usr/bin/env python3
"""Auxiliary binary crossings: exact controls for the contraction-domination floor.

For random matrix-valued binary Laurent crossings K(w) in Mat_D(F_p[w,1/w])
with K(1)e = E e, random integer physical heights at both places, and
independent auxiliary conjugation gauges at the two places (absent, diagonal
monomial, framed, or non-diagonal "conjugated" gauges built so that the
constant boundary contraction fails), this checker

* computes the rectangular Smith orders of both cut assignments at both
  places, with the auxiliary pair on the row side or on the column side of
  each cut independently (local elimination over F_p[[w]], cross-checked
  against the determinantal-divisor definition by all minors on a subset);
* constructs a Laurent "twisted" boundary pair (f~, e~) with e~(1)=e,
  f~(1).e=1 whose weights are regular at both places (linear algebra over
  F_p), and verifies that S~ = f~^T K e~ is a scalar binary strict crossing;
* checks coordinatewise domination of the auxiliary orders by those of S~
  (hence prefix domination), the exact Tomic--Weyl partition comparison,
  every normalized inequality of the scalar floor directly on the auxiliary
  profile, and outward dyadic pressure bounds against the analytic floor;
* replays the stored 421-leaf compact certificate of the scalar floor;
* records where the constant contraction f^T K e fails to dominate, per place
  and at the pressure level (exact rational comparison at e^{-t} rational);
* recomputes the explicit seed of the note's Section 5, in which every
  constant contraction fails per-place domination (exhaustively over F_2 and
  F_3; the note proves it over every field) and a mixed boundary contraction
  has strictly larger pressure than the auxiliary seed.

All conclusions are computed; exit status is nonzero on any failure.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random
import sys

from charged_gaussian_floor import Dyadic
from check_global_homogeneous_boolean_filter_floor import audit_cover

F = Fraction
ROOT = Path(__file__).resolve().parent.parent
LABELS = [((i >> 3) & 1, (i >> 2) & 1, (i >> 1) & 1, i & 1) for i in range(16)]  # (x,y,u,v)


def crossing(lab):
    x, y, u, v = lab
    return int(x == v and y == u)


# ---------------------------------------------------------------- Laurent polynomials over F_p
def ladd(a, b, p):
    out = dict(a)
    for k, c in b.items():
        s = (out.get(k, 0) + c) % p
        if s:
            out[k] = s
        else:
            out.pop(k, None)
    return out


def lscale(a, c, p, shift=0):
    c %= p
    return {k + shift: (v * c) % p for k, v in a.items()} if c else {}


def lmul(a, b, p):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] = (out.get(i + j, 0) + x * y) % p
    return {k: v for k, v in out.items() if v}


def lval(a):
    return min(a) if a else None


def leval1(a, p):
    return sum(a.values()) % p


def mat_mul(A, B, p):
    n, m, r = len(A), len(B), len(B[0])
    out = [[{} for _ in range(r)] for _ in range(n)]
    for i in range(n):
        for k in range(m):
            if not A[i][k]:
                continue
            for j in range(r):
                if B[k][j]:
                    out[i][j] = ladd(out[i][j], lmul(A[i][k], B[k][j], p), p)
    return out


def identity(D):
    return [[{0: 1} if i == j else {} for j in range(D)] for i in range(D)]


def transpose(A):
    return [list(row) for row in zip(*A)]


# ---------------------------------------------------------------- exact local Smith orders
def series_inverse(u, N, p):
    inv0 = pow(u[0], p - 2, p)
    out = [0] * N
    out[0] = inv0
    for k in range(1, N):
        s = 0
        for i in range(1, min(k, len(u) - 1) + 1):
            s += u[i] * out[k - i]
        out[k] = (-s * inv0) % p
    return out


def series_mul(a, b, N, p):
    out = [0] * N
    for i, x in enumerate(a[:N]):
        if x:
            for j in range(min(len(b), N - i)):
                if b[j]:
                    out[i + j] = (out[i + j] + x * b[j]) % p
    return out


def smith_orders(M, p, rank):
    """Sorted Smith orders at w=0 of a Laurent matrix of generic rank `rank`."""
    entries = [e for row in M for e in row if e]
    base = min(min(e) for e in entries)
    N = 8
    while True:
        rows = []
        for row in M:
            out = []
            for e in row:
                s = [0] * N
                for k, c in e.items():
                    if k - base < N:
                        s[k - base] = c
                out.append(s)
            rows.append(out)
        active_rows = list(range(len(rows)))
        active_cols = list(range(len(M[0])))
        orders = []
        while active_cols:
            best = None
            for i in active_rows:
                for j in active_cols:
                    s = rows[i][j]
                    v = next((k for k, c in enumerate(s) if c), N)
                    if best is None or v < best[0]:
                        best = (v, i, j)
            v, i, j = best
            if v >= N:
                break
            unit_inv = series_inverse(rows[i][j][v:] + [0] * v, N, p)
            for r in active_rows:
                if r == i:
                    continue
                entry = rows[r][j]
                if not any(entry):
                    continue
                factor = series_mul(entry[v:], unit_inv, N - v, p)
                for c in active_cols:
                    piv = rows[i][c][v:]
                    prod = series_mul(factor, piv, N - v, p)
                    target = rows[r][c]
                    for k in range(N - v):
                        if prod[k]:
                            target[k + v] = (target[k + v] - prod[k]) % p
            active_rows.remove(i)
            active_cols.remove(j)
            orders.append(v)
        if len(orders) == rank and max(orders) < N:
            assert orders == sorted(orders)
            return [o + base for o in orders]
        N *= 2
        assert N <= 4096, 'precision runaway (rank deficient?)'


def ldet(A, p):
    k = len(A)
    total = {}
    for perm in permutations(range(k)):
        sign = 1
        for i in range(k):
            for j in range(i + 1, k):
                if perm[i] > perm[j]:
                    sign = -sign
        term = {0: 1}
        for i in range(k):
            term = lmul(term, A[i][perm[i]], p)
            if not term:
                break
        if term:
            total = ladd(total, lscale(term, sign, p), p)
    return total


def determinantal_profile(M, p, rank):
    """Smith orders from the definition: minimum valuations of all k-minors."""
    prefixes = [0]
    m, n = len(M), len(M[0])
    for k in range(1, rank + 1):
        best = None
        for rows in combinations(range(m), k):
            sub_rows = [M[r] for r in rows]
            if any(all(not e for e in row) for row in sub_rows):
                continue
            for cols in combinations(range(n), k):
                d = ldet([[row[c] for c in cols] for row in sub_rows], p)
                if d:
                    v = min(d)
                    if best is None or v < best:
                        best = v
        prefixes.append(best)
    return [prefixes[k] - prefixes[k - 1] for k in range(1, rank + 1)]


# ---------------------------------------------------------------- cut matrices with heights
def physical_weight(heights, lab):
    return sum(heights[l][lab[l]] for l in range(4))


def cut_matrix(entries, D, place, cut, heights, aux_side='row'):
    """entries[lab] is a DxD Laurent matrix (D=1 for scalar).

    The auxiliary pair joins the physical row pair (aux_side='row', a
    4D^2-by-4 matrix) or the physical column pair (aux_side='col', a
    4-by-4D^2 matrix, returned transposed; Smith orders are unchanged by
    transposition)."""
    if cut == 'row':  # (x,y)|(u,v)
        rkey = lambda lab: (lab[0], lab[1])
        ckey = lambda lab: (lab[2], lab[3])
    else:  # (x,u)|(y,v)
        rkey = lambda lab: (lab[0], lab[2])
        ckey = lambda lab: (lab[1], lab[3])
    if aux_side == 'col':
        rkey, ckey = ckey, rkey  # build the transpose: aux joins the column pair
    rows = {}
    for lab in LABELS:
        w = physical_weight(heights, lab)
        for a in range(D):
            for b in range(D):
                e = entries[lab][a][b]
                if place == 0:
                    g = {k + w: c for k, c in e.items()}
                else:  # tau = 1/w, tau^{weight} * e(1/tau)
                    g = {w - k: c for k, c in e.items()}
                rows.setdefault((rkey(lab), a, b), {})[ckey(lab)] = g
    cols = [(0, 0), (0, 1), (1, 0), (1, 1)]
    return [[rows[r][c] for c in cols] for r in sorted(rows)]


def gauged(K, G, Ginv, p):
    return {lab: mat_mul(mat_mul(G, K[lab], p), Ginv, p) for lab in LABELS}


# ---------------------------------------------------------------- linear algebra over F_p
def solve(equations, nvars, p):
    """equations: list of (dict var->coef, rhs). Return one solution or None."""
    rows = [(dict(eq), rhs % p) for eq, rhs in equations]
    pivots = []
    solution_rows = []
    for eq, rhs in rows:
        eq = {k: v % p for k, v in eq.items() if v % p}
        for (pv, peq, prhs) in solution_rows:
            if pv in eq:
                c = eq[pv]
                for k, v in peq.items():
                    eq[k] = (eq.get(k, 0) - c * v) % p
                    if not eq[k]:
                        del eq[k]
                rhs = (rhs - c * prhs) % p
        if not eq:
            if rhs:
                return None
            continue
        pv = min(eq)
        inv = pow(eq[pv], p - 2, p)
        eq = {k: (v * inv) % p for k, v in eq.items()}
        rhs = (rhs * inv) % p
        new_rows = []
        for (qv, qeq, qrhs) in solution_rows:
            if pv in qeq:
                c = qeq[pv]
                for k, v in eq.items():
                    qeq[k] = (qeq.get(k, 0) - c * v) % p
                    if not qeq[k]:
                        del qeq[k]
                qrhs = (qrhs - c * rhs) % p
            new_rows.append((qv, qeq, qrhs))
        solution_rows = new_rows + [(pv, eq, rhs)]
    x = [0] * nvars
    for pv, eq, rhs in solution_rows:
        x[pv] = rhs  # free variables are zero; reduced echelon form
    return x


def twisted_pair(G0, G0inv, Gi, Giinv, e, p, window=6, twist_range=12):
    """Laurent vectors e~, f~ with e~(1)=e, f~(1).e=1, G0 e~, G0^{-T} f~ regular at
    w=0 and Gi e~ (exponents<=c), Gi^{-T} f~ (exponents<=-c) regular at infinity."""
    D = len(e)
    exps = list(range(-window, window + 1))
    var = {(b, j): idx for idx, (b, j) in enumerate((b, j) for b in range(D) for j in exps)}

    def forms(G):
        out = {}
        for a in range(D):
            for (b, j), idx in var.items():
                for k, c in G[a][b].items():
                    out.setdefault((a, k + j), {})
                    out[(a, k + j)][idx] = (out[(a, k + j)].get(idx, 0) + c) % p
        return out

    f0, fi = forms(G0), forms(Gi)
    g0, gi = forms(transpose(G0inv)), forms(transpose(Giinv))
    for c in sorted(range(-twist_range, twist_range + 1), key=abs):
        eqs = [({var[(b, j)]: 1 for j in exps}, e[b]) for b in range(D)]
        eqs += [(form, 0) for (a, k), form in f0.items() if k < 0]
        eqs += [(form, 0) for (a, k), form in fi.items() if k > c]
        etil = solve(eqs, len(var), p)
        if etil is None:
            continue
        eqs = [({var[(b, j)]: e[b] for b in range(D) for j in exps}, 1)]
        eqs += [(form, 0) for (a, k), form in g0.items() if k < 0]
        eqs += [(form, 0) for (a, k), form in gi.items() if k > -c]
        ftil = solve(eqs, len(var), p)
        if ftil is None:
            continue
        to_vec = lambda x: [{j: x[var[(b, j)]] for j in exps if x[var[(b, j)]]} for b in range(D)]
        return to_vec(etil), to_vec(ftil), c
    return None


def apply(G, vec, p):
    out = []
    for a in range(len(G)):
        total = {}
        for b in range(len(vec)):
            total = ladd(total, lmul(G[a][b], vec[b], p), p)
        out.append(total)
    return out


def contraction(K, ftil, etil, p):
    """Scalar tensor S(lab) = sum_ab ftil_a etil_b K(lab)_ab, as 1x1 matrices."""
    D = len(etil)
    weights = {(a, b): lmul(ftil[a], etil[b], p) for a in range(D) for b in range(D)}
    out = {}
    for lab in LABELS:
        total = {}
        for (a, b), wt in weights.items():
            if wt and K[lab][a][b]:
                total = ladd(total, lmul(wt, K[lab][a][b], p), p)
        out[lab] = [[total]]
    return out


# ---------------------------------------------------------------- random data
def random_laurent(rng, p, lo, hi, density):
    return {k: c for k in range(lo, hi + 1) if rng.random() < density and (c := rng.randrange(1, p))}


def random_gauge(rng, p, D, mode):
    if mode == 'none':
        return identity(D), identity(D)
    g = [rng.randrange(-3, 4) for _ in range(D)]
    diag = [[{g[i]: 1} if i == j else {} for j in range(D)] for i in range(D)]
    diag_inv = [[{-g[i]: 1} if i == j else {} for j in range(D)] for i in range(D)]
    if mode == 'diagonal':
        return diag, diag_inv
    # framed: diag(w^g) (I + w^j N) P with N strictly upper triangular, P constant invertible
    while True:
        P = [[rng.randrange(p) for _ in range(D)] for _ in range(D)]
        Pinv = constant_inverse(P, p)
        if Pinv is not None:
            break
    j = rng.randrange(-2, 3)
    N = [[{j: rng.randrange(p)} if b > a and rng.random() < 0.8 else {} for b in range(D)] for a in range(D)]
    N = [[{k: c for k, c in e.items() if c} for e in row] for row in N]
    U = [[ladd(identity(D)[a][b], N[a][b], p) for b in range(D)] for a in range(D)]
    Uinv = identity(D)
    power = identity(D)
    for k in range(1, D):
        power = mat_mul(power, N, p)
        sign = -1 if k % 2 else 1
        Uinv = [[ladd(Uinv[a][b], lscale(power[a][b], sign, p), p) for b in range(D)] for a in range(D)]
    Pm = [[{0: P[a][b]} if P[a][b] else {} for b in range(D)] for a in range(D)]
    Pim = [[{0: Pinv[a][b]} if Pinv[a][b] else {} for b in range(D)] for a in range(D)]
    G = mat_mul(mat_mul(diag, U, p), Pm, p)
    Ginv = mat_mul(mat_mul(Pim, Uinv, p), diag_inv, p)
    assert mat_mul(G, Ginv, p) == identity(D)
    return G, Ginv


def constant_inverse(P, p):
    D = len(P)
    A = [list(row) + [int(i == j) for j in range(D)] for i, row in enumerate(P)]
    for col in range(D):
        piv = next((r for r in range(col, D) if A[r][col] % p), None)
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        inv = pow(A[col][col], p - 2, p)
        A[col] = [(v * inv) % p for v in A[col]]
        for r in range(D):
            if r != col and A[r][col]:
                c = A[r][col]
                A[r] = [(x - c * y) % p for x, y in zip(A[r], A[col])]
    return [row[D:] for row in A]


def random_kernel(rng, p, D):
    lo = rng.randrange(-2, 1)
    hi = lo + rng.randrange(0, 4)
    density = rng.choice((0.35, 0.6))
    e = [rng.randrange(p) for _ in range(D)]
    if not any(e):
        e[rng.randrange(D)] = 1
    f = [rng.randrange(p) for _ in range(D)]
    i = next(k for k in range(D) if e[k])
    rest = sum(f[k] * e[k] for k in range(D) if k != i)
    f[i] = ((1 - rest) * pow(e[i], p - 2, p)) % p
    assert sum(a * b for a, b in zip(f, e)) % p == 1
    K = {}
    for lab in LABELS:
        M = [[random_laurent(rng, p, lo, hi, density) for _ in range(D)] for _ in range(D)]
        Me = [sum(leval1(M[a][b], p) * e[b] for b in range(D)) % p for a in range(D)]
        corr = [(crossing(lab) * e[a] - Me[a]) % p for a in range(D)]
        for a in range(D):
            for b in range(D):
                M[a][b] = ladd(M[a][b], lscale({0: 1}, corr[a] * f[b], p), p)
        K[lab] = M
    return K, e, f


# ---------------------------------------------------------------- pressure arithmetic
def point(ar, value):
    value = F(value)
    return ar.rational(value.numerator, value.denominator)


def log_one_plus(ar, value):
    return (ar.log_one_plus((value[0], value[0]))[0], ar.log_one_plus((value[1], value[1]))[1])


def log_partition(ar, orders, t):
    minimum = min(orders)
    total = ar.zero
    for value in orders:
        total = ar.add(total, ar.exp_negative(ar.multiply(ar.rational(value - minimum), t)))
    return ar.subtract(log_one_plus(ar, ar.subtract(total, ar.one)), ar.multiply(ar.rational(minimum), t))


def general_floor(ar, A, X, Y):
    """Phi(A,X,Y) of equation (18) of the universal binary floor (outward)."""
    a, x, y = (point(ar, v) for v in (A, X, Y))
    minimum = point(ar, min(X, Y))
    first = ar.exp_negative(ar.add(a, minimum))
    last_exp = ar.divide_integer(ar.subtract(ar.add(ar.multiply(point(ar, 3), a),
                                                    ar.multiply(point(ar, 4), ar.add(x, y))), minimum), 2)
    last = ar.multiply(point(ar, 2), ar.exp_negative(last_exp))
    positive = ar.add(a, ar.multiply(point(ar, 2), log_one_plus(ar, ar.add(first, last))))
    physical = ar.multiply(ar.root_two, ar.add(ar.softplus(ar.multiply(ar.root_two, x)),
                                               ar.softplus(ar.multiply(ar.root_two, y))))
    return ar.subtract(positive, physical)


def canonical(heights0, heightsi, sigma0, sigmai):
    """Letter-0 normalization and paired flips; returns canonical data."""
    diff = [[h[l][1] - h[l][0] for l in range(4)] for h in (heights0, heightsi)]
    H = [diff[0][l] + diff[1][l] for l in range(4)]
    EX, EY = H[0] + H[3], H[1] + H[2]
    flip = [EX > 0, EY > 0, EY > 0, EX > 0]
    shifts = [sum(h[l][1 if flip[l] else 0] for l in range(4)) for h in (heights0, heightsi)]
    canon = [[(-d if flip[l] else d) for l, d in enumerate(row)] for row in diff]
    Hc = [canon[0][l] + canon[1][l] for l in range(4)]
    a, b = -F(Hc[0] + Hc[3], 2), -F(Hc[1] + Hc[2], 2)
    s0 = [v - shifts[0] for v in sigma0]
    si = [v - shifts[1] for v in sigmai]
    return Hc, a, b, s0, si


def critical_pressure(ar, s0, si, a, b, t):
    tt = point(ar, t)
    actual = ar.add(log_partition(ar, s0, tt), log_partition(ar, si, tt))
    physical = ar.add(point(ar, 2 * (a + b) * t), ar.multiply(ar.root_two, ar.add(
        ar.softplus(ar.multiply(ar.root_two, point(ar, a * t))),
        ar.softplus(ar.multiply(ar.root_two, point(ar, b * t))))))
    return ar.subtract(actual, physical)


def prefixes(orders):
    out, s = [], 0
    for v in orders:
        s += v
        out.append(s)
    return out


# ---------------------------------------------------------------- compact certificate replay
def floor_interval(ar, rectangle):
    """F(A,t) of equation (25), outward; transcribed from check_binary_fixed_ray_floor.py (PR #39)."""
    lo_a, hi_a, lo_t, hi_t = rectangle
    a = (ar.rational(lo_a.numerator, lo_a.denominator)[0], ar.rational(hi_a.numerator, hi_a.denominator)[1])
    t = (ar.rational(lo_t.numerator, lo_t.denominator)[0], ar.rational(hi_t.numerator, hi_t.denominator)[1])
    first = ar.exp_negative(ar.add(a, t))
    second = ar.exp_negative(ar.divide_integer(ar.add(ar.multiply(ar.rational(3), a),
                                                      ar.multiply(ar.rational(7), t)), 2))
    argument = ar.add(first, ar.multiply(ar.rational(2), second))
    positive = ar.add(a, ar.multiply(ar.rational(2), log_one_plus(ar, argument)))
    negative = ar.multiply(ar.multiply(ar.rational(2), ar.root_two),
                           ar.softplus(ar.multiply(ar.root_two, t)))
    return ar.subtract(positive, negative)


def split_box(rectangle):
    la, ha, lt, ht = rectangle
    if ha - la >= ht - lt:
        m = (la + ha) / 2
        return ('0', (la, m, lt, ht)), ('1', (m, ha, lt, ht))
    m = (lt + ht) / 2
    return ('2', (la, ha, lt, m)), ('3', (la, ha, m, ht))


CERTIFICATE = ROOT / 'evidence' / 'binary-fixed-ray-floor.json'


def replay_cover(ar):
    record = json.loads(CERTIFICATE.read_text())
    compact = record['compact_certificate']
    margin = F(record['compact_margin'])
    assert margin == F(1, 16)
    # the scalar note's (28) needs F(A,t) > 1/16 on exactly [0,3/2] x [0,3/2]
    assert tuple(F(v) for v in compact['domain']) == (0, F(3, 2), 0, F(3, 2)), 'wrong certificate domain'
    function = lambda r: ar.subtract(floor_interval(ar, r), ar.rational(margin.numerator, margin.denominator))
    nodes = audit_cover(function, compact, split_box)
    rejected = 0
    for mutation in (dict(compact, cover=compact['cover'][1:]),
                     dict(compact, cover=[[path, str(int(v) + (i == 0))] for i, (path, v) in enumerate(compact['cover'])]),
                     dict(compact, domain=['0', '3', '0', '3'])):
        try:
            audit_cover(function, mutation, split_box)
        except AssertionError:
            rejected += 1
    assert rejected == 3, 'malformed cover accepted'
    # tail inequality constants used in (27): (5/3)^5 < 16 and 7/8 > 6/7.
    assert F(5, 3) ** 5 < 16 and F(7, 8) - F(6, 7) == F(1, 56)
    return {'leaves': compact['positive_leaves'], 'audited_nodes': nodes, 'malformed_rejected': rejected,
            'domain': compact['domain'], 'margin': str(margin)}


# ---------------------------------------------------------------- one control
SIDES = (('row', 'row'), ('row', 'col'), ('col', 'row'), ('col', 'col'))  # aux side at zero, at infinity
ASSIGNMENTS = (('row', 'col'), ('col', 'row'))  # cut at zero, cut at infinity
XS = (F(1, 2), F(9, 10), F(99, 100))  # exact values of e^{-t}
COUNT = {'rectangular_all_minor_replays': 0, 'scalar_all_minor_replays': 0}


def partition(orders, x):
    return sum(x ** v for v in orders)


def coordinatewise(mine, theirs):
    return all(x <= y for x, y in zip(mine, theirs))


def prefix_dominated(mine, theirs):
    return all(x <= y for x, y in zip(prefixes(mine), prefixes(theirs)))


def pressure_exceeds(aux, other):
    """Exact: G_other(t) > G_aux(t) at some tested rational e^{-t}. With equal heights the
    physical terms and the letter-0 shifts cancel, so this is a comparison of B(t)."""
    (a0, ai), (o0, oi) = aux, other
    return any(partition(a0, x) * partition(ai, x) < partition(o0, x) * partition(oi, x) for x in XS)


def nonconstant(vectors):
    return any(k != 0 for v in vectors for k in v)


def floor_down(ar, value, denominator=10 ** 6):
    """A rational at most value/scale (outward rounding down), as a string."""
    return str(F(value * denominator // ar.scale, denominator))


def check_instance(ar, p, K, e, f, gauges, heights0, heightsi, mode, minors_check):
    """All checks of one auxiliary seed; returns its record and constant-contraction statistics."""
    G0, G0inv, Gi, Giinv = gauges
    D = len(e)
    for lab in LABELS:  # literal endpoint: K(1)e = E e
        Ke = [sum(leval1(K[lab][a][b], p) * e[b] for b in range(D)) % p for a in range(D)]
        assert Ke == [(crossing(lab) * x) % p for x in e]
    assert sum(a * b for a, b in zip(f, e)) % p == 1
    K0, Ki = gauged(K, G0, G0inv, p), gauged(K, Gi, Giinv, p)

    pair = None
    for window in (6, 12, 24):  # widen the Laurent window until a regular pair appears
        pair = twisted_pair(G0, G0inv, Gi, Giinv, e, p, window=window, twist_range=2 * window)
        if pair is not None:
            break
    assert pair is not None, 'no regular twisted boundary pair found'
    etil, ftil, twist = pair
    # value conditions at w=1
    assert [leval1(x, p) for x in etil] == [x % p for x in e]
    assert sum(leval1(ftil[b], p) * e[b] for b in range(D)) % p == 1
    # regularity of every contraction weight at each place, checked directly
    weights = {}
    for name, G, Ginv, ok in (('zero', G0, G0inv, lambda k: k >= 0), ('infinity', Gi, Giinv, lambda k: k <= 0)):
        ehat = apply(G, etil, p)
        fhat = apply(transpose(Ginv), ftil, p)
        for a in range(D):
            for b in range(D):
                assert all(ok(k) for k in lmul(fhat[a], ehat[b], p)), (name, a, b)
        weights[name] = (fhat, ehat)
    S = contraction(K, ftil, etil, p)
    for lab in LABELS:  # S~ is a scalar strict crossing
        assert leval1(S[lab][0][0], p) == crossing(lab)
    # S~ equals the regular-weight contraction of each gauged matrix (input of Lemma 1)
    for name, Kg in (('zero', K0), ('infinity', Ki)):
        fhat, ehat = weights[name]
        alt = contraction(Kg, fhat, ehat, p)
        assert all(alt[lab][0][0] == S[lab][0][0] for lab in LABELS)
    naive = contraction(K, [{0: x} if x else {} for x in f], [{0: x} if x else {} for x in e], p)

    out = {'p': p, 'D': D, 'gauge_mode': mode, 'twist': twist,
           'twisted_pair_nonconstant': nonconstant(etil) or nonconstant(ftil),
           'heights_zero': heights0, 'heights_infinity': heightsi, 'assignments': []}
    stats = {'checked': 0, 'profile_failures': 0, 'pressure_failures': 0}
    for cut0, cuti in ASSIGNMENTS:
        ssig0 = smith_orders(cut_matrix(S, 1, 0, cut0, heights0), p, 4)
        ssigi = smith_orders(cut_matrix(S, 1, 'inf', cuti, heightsi), p, 4)
        # independent determinantal-divisor replay of the scalar profiles (always)
        assert determinantal_profile(cut_matrix(S, 1, 0, cut0, heights0), p, 4) == ssig0
        assert determinantal_profile(cut_matrix(S, 1, 'inf', cuti, heightsi), p, 4) == ssigi
        COUNT['scalar_all_minor_replays'] += 2
        nsig0 = smith_orders(cut_matrix(naive, 1, 0, cut0, heights0), p, 4)
        nsigi = smith_orders(cut_matrix(naive, 1, 'inf', cuti, heightsi), p, 4)
        _, _, _, sc0, sci = canonical(heights0, heightsi, ssig0, ssigi)
        for side0, sidei in (SIDES if D > 1 else SIDES[:1]):  # D=1: the side is immaterial
            M0 = cut_matrix(K0, D, 0, cut0, heights0, side0)
            Mi = cut_matrix(Ki, D, 'inf', cuti, heightsi, sidei)
            sig0, sigi = smith_orders(M0, p, 4), smith_orders(Mi, p, 4)
            if minors_check and side0 == sidei:  # rectangular replay by all minors (subset)
                assert determinantal_profile(M0, p, 4) == sig0
                assert determinantal_profile(Mi, p, 4) == sigi
                COUNT['rectangular_all_minor_replays'] += 2
            # (A) Lemma 1: coordinatewise, hence prefix, domination at both places
            for mine, theirs in ((sig0, ssig0), (sigi, ssigi)):
                assert coordinatewise(mine, theirs) and prefix_dominated(mine, theirs)
            # (B) exact partition comparison at rational e^{-t}
            for xval in XS:
                for mine, theirs in ((sig0, ssig0), (sigi, ssigi)):
                    assert partition(mine, xval) >= partition(theirs, xval)
            # constant contraction f^T K e: per-place and pressure-level comparison
            dominated = prefix_dominated(sig0, nsig0) and prefix_dominated(sigi, nsigi)
            exceeds = pressure_exceeds((sig0, sigi), (nsig0, nsigi))
            if mode == 'none':
                assert dominated  # constant contraction suffices without auxiliary gauges
            stats['checked'] += 1
            stats['profile_failures'] += not dominated
            stats['pressure_failures'] += exceeds
            # (C) direct normalized audit of the scalar proof on the auxiliary profile
            Hc, a, b, c0, ci = canonical(heights0, heightsi, sig0, sigi)
            mu0, mui = c0[0], ci[0]
            d = sum(Hc) - mu0 - mui
            r = [v - mu0 for v in c0]
            s = [v - mui for v in ci]
            C = 4 * d + 4 * (a + b)
            assert a >= 0 and b >= 0 and d >= 0 and r[0] == s[0] == 0
            assert sum(r) <= C and sum(s) <= C and r[1] + s[1] <= 2 * d + 2 * min(a, b)
            # (D) outward pressure: G_aux >= floor(d) > 0 and G_S~ > 0
            lower = []
            for t in (F(1, 10), F(1, 2), F(1), F(2)):
                g_aux = critical_pressure(ar, c0, ci, a, b, t)
                g_s = critical_pressure(ar, sc0, sci, a, b, t)
                bound = general_floor(ar, d * t, a * t, b * t)
                assert g_aux[0] >= bound[1] and bound[0] > 0 and g_s[0] > 0
                lower.append(floor_down(ar, bound[0]))
            out['assignments'].append({'zero_cut': cut0, 'infinity_cut': cuti,
                                       'aux_side_zero': side0, 'aux_side_infinity': sidei,
                                       'aux_zero': sig0, 'aux_infinity': sigi,
                                       'contraction_zero': ssig0, 'contraction_infinity': ssigi,
                                       'constant_contraction_zero': nsig0,
                                       'constant_contraction_infinity': nsigi,
                                       'strict_domination': sig0 != ssig0 or sigi != ssigi,
                                       'constant_contraction_profile_dominates': dominated,
                                       'constant_contraction_pressure_exceeds_aux': exceeds,
                                       'canonical_a_b_d': [str(a), str(b), str(d)],
                                       'floor_rounded_down_t_0.1_0.5_1_2': lower})
    return out, stats


def run_case(rng, ar, p, D, mode, minors_check):
    K, e, f = random_kernel(rng, p, D)
    heights0 = [[rng.randrange(-3, 4) for _ in range(2)] for _ in range(4)]
    heightsi = [[rng.randrange(-3, 4) for _ in range(2)] for _ in range(4)]
    G0, G0inv = random_gauge(rng, p, D, mode)
    Gi, Giinv = random_gauge(rng, p, D, mode)
    return check_instance(ar, p, K, e, f, (G0, G0inv, Gi, Giinv), heights0, heightsi, mode, minors_check)


def mono(k, c=1):
    return {k: c} if c else {}


def two_by_two(a, b, c, d):
    return [[a, b], [c, d]]


def conjugated_case(rng, ar, p):
    """Non-diagonal independent gauges built so that f^T K e can fail.

    K = V [[A,0],[E,0]] V^{-1} with A = E + (w-1)R a random scalar strict crossing and
    V = (permutation)(I + c w^j N); e = V(1)(1,1). Gauges G_0 = U_0 diag(1,w^g) V^{-1}
    and G_inf = U_inf diag(1,w^h) V^{-1}, with U_0 unipotent and regular at zero,
    U_inf unipotent and regular at infinity, g >= 4 and h <= -4."""
    one, zero = mono(0), {}
    E = {lab: (mono(0) if crossing(lab) else {}) for lab in LABELS}
    A = {lab: ladd(E[lab], lmul({1: 1, 0: p - 1}, random_laurent(rng, p, -1, 1, 0.5), p), p) for lab in LABELS}
    Kb = {lab: two_by_two(A[lab], zero, dict(E[lab]), zero) for lab in LABELS}
    j, c = rng.randrange(-2, 3), rng.randrange(1, p)
    V = two_by_two(one, mono(j, c), zero, one)
    Vinv = two_by_two(one, mono(j, (-c) % p), zero, one)
    if rng.random() < 0.5:  # swap the auxiliary basis
        Pm = two_by_two(zero, one, one, zero)
        V, Vinv = mat_mul(Pm, V, p), mat_mul(Vinv, Pm, p)
    assert mat_mul(V, Vinv, p) == identity(2)
    K = {lab: mat_mul(mat_mul(V, Kb[lab], p), Vinv, p) for lab in LABELS}
    e = [sum(leval1(V[a][b], p) for b in range(2)) % p for a in range(2)]
    while True:
        f = [rng.randrange(p) for _ in range(2)]
        if (f[0] * e[0] + f[1] * e[1]) % p == 1:
            break
    g, h = rng.randrange(4, 11), -rng.randrange(4, 11)
    i0, c0 = rng.randrange(0, 3), rng.randrange(1, p)
    i1, c1 = rng.randrange(0, 3), rng.randrange(1, p)
    U0, U0inv = two_by_two(one, mono(i0, c0), zero, one), two_by_two(one, mono(i0, (-c0) % p), zero, one)
    Ui, Uiinv = two_by_two(one, zero, mono(-i1, c1), one), two_by_two(one, zero, mono(-i1, (-c1) % p), one)
    D0, D0inv = two_by_two(one, zero, zero, mono(g)), two_by_two(one, zero, zero, mono(-g))
    Di, Diinv = two_by_two(one, zero, zero, mono(h)), two_by_two(one, zero, zero, mono(-h))
    G0, G0inv = mat_mul(mat_mul(U0, D0, p), Vinv, p), mat_mul(mat_mul(V, D0inv, p), U0inv, p)
    Gi, Giinv = mat_mul(mat_mul(Ui, Di, p), Vinv, p), mat_mul(mat_mul(V, Diinv, p), Uiinv, p)
    assert mat_mul(G0, G0inv, p) == identity(2) and mat_mul(Gi, Giinv, p) == identity(2)
    diagonal = all(not G[a][b] for G in (G0, Gi) for a in range(2) for b in range(2) if a != b)
    heights0 = [[rng.randrange(-2, 3) for _ in range(2)] for _ in range(4)]
    heightsi = [[rng.randrange(-2, 3) for _ in range(2)] for _ in range(4)]
    record, stats = check_instance(ar, p, K, e, f, (G0, G0inv, Gi, Giinv), heights0, heightsi,
                                   'conjugated', False)
    record['gauges_non_diagonal'] = not diagonal
    return record, stats


def explicit_negative_control(ar):
    """Section 5 seed: K = [[0, w^5 E], [E, 0]], e = (1,1), G_0 = diag(1,w^3), G_inf = diag(1,w),
    zero heights. Every constant contraction phi(K) with phi(K(1)) = E (exhaustive over F_2 and
    F_3) fails per-place prefix domination; the mixed boundary contraction f = (2,2) over F_3
    has strictly larger pressure for every t > 0; the twisted contraction (4) dominates."""
    gamma, eta, m = (0, 3), (0, 1), 5
    E = {lab: (mono(0) if crossing(lab) else {}) for lab in LABELS}
    K = {lab: two_by_two({}, (mono(m) if crossing(lab) else {}), dict(E[lab]), {}) for lab in LABELS}
    e = [1, 1]
    D = 2
    G0, G0inv = [[mono(gamma[a]) if a == b else {} for b in range(D)] for a in range(D)], \
        [[mono(-gamma[a]) if a == b else {} for b in range(D)] for a in range(D)]
    Gi, Giinv = [[mono(eta[a]) if a == b else {} for b in range(D)] for a in range(D)], \
        [[mono(-eta[a]) if a == b else {} for b in range(D)] for a in range(D)]
    zero_h = [[0, 0] for _ in range(4)]
    out = {'seed': 'K=[[0,w^5 E],[E,0]], e=(1,1), G_0=diag(1,w^3), G_inf=diag(1,w), zero heights'}
    for p in (2, 3):
        for lab in LABELS:
            assert [sum(leval1(K[lab][a][b], p) * e[b] for b in range(D)) % p for a in range(D)] == [crossing(lab)] * 2
        K0, Ki = gauged(K, G0, G0inv, p), gauged(K, Gi, Giinv, p)
        aux = {}
        for cut0, cuti in ASSIGNMENTS:
            for side0, sidei in SIDES:
                M0 = cut_matrix(K0, D, 0, cut0, zero_h, side0)
                Mi = cut_matrix(Ki, D, 'inf', cuti, zero_h, sidei)
                aux[(cut0, cuti, side0, sidei)] = (smith_orders(M0, p, 4), smith_orders(Mi, p, 4))
                assert determinantal_profile(M0, p, 4) == aux[(cut0, cuti, side0, sidei)][0]
                assert determinantal_profile(Mi, p, 4) == aux[(cut0, cuti, side0, sidei)][1]
                COUNT['rectangular_all_minor_replays'] += 2
        profiles = sorted({(tuple(a0), tuple(ai)) for a0, ai in aux.values()})
        # every constant functional phi on the auxiliary pair index with phi(K(1)) = E
        n_phi = per_place = pressure = 0
        for phi in product(range(p), repeat=D * D):
            Sphi = {lab: [[{}]] for lab in LABELS}
            for lab in LABELS:
                total = {}
                for idx, coefficient in enumerate(phi):
                    total = ladd(total, lscale(K[lab][idx // D][idx % D], coefficient, p), p)
                Sphi[lab] = [[total]]
            if any(leval1(Sphi[lab][0][0], p) != crossing(lab) for lab in LABELS):
                continue
            n_phi += 1
            fails_everywhere = exceeds_in_every_combination = True
            for (cut0, cuti, side0, sidei), (a0, ai) in aux.items():
                s0 = smith_orders(cut_matrix(Sphi, 1, 0, cut0, zero_h), p, 4)
                si = smith_orders(cut_matrix(Sphi, 1, 'inf', cuti, zero_h), p, 4)
                fails_everywhere &= not (prefix_dominated(a0, s0) and prefix_dominated(ai, si))
                exceeds_in_every_combination &= pressure_exceeds((a0, ai), (s0, si))
            assert fails_everywhere, ('a constant contraction dominates', p, phi)
            per_place += 1
            pressure += exceeds_in_every_combination
        # the mixed boundary contraction f^T K e with f = (2,2) over F_3 (f.e = 1)
        mixed = None
        if p == 3:
            f = [2, 2]
            assert (f[0] * e[0] + f[1] * e[1]) % p == 1
            Sf = contraction(K, [mono(0, x) for x in f], [mono(0, x) for x in e], p)
            exponent = set()
            for (cut0, cuti, side0, sidei), (a0, ai) in aux.items():
                s0 = smith_orders(cut_matrix(Sf, 1, 0, cut0, zero_h), p, 4)
                si = smith_orders(cut_matrix(Sf, 1, 'inf', cuti, zero_h), p, 4)
                assert all(len(set(v)) == 1 for v in (a0, ai, s0, si))  # constant profiles
                # then B(t) = 16 e^{-t(z+i)} and G_f(t) - G_aux(t) = t((z_aux+i_aux) - (z_f+i_f))
                exponent.add((a0[0] + ai[0]) - (s0[0] + si[0]))
            assert len(exponent) == 1 and min(exponent) > 0
            mixed = {'f': f, 'constant_contraction_zero': s0, 'constant_contraction_infinity': si,
                     'G_f_minus_G_aux_over_t': min(exponent)}
        # the explicit twisted contraction (4) and the solver's twisted pair
        cstar = max((b for b in range(D) if e[b] % p), key=lambda b: (eta[b] - gamma[b], -b))
        inv = pow(e[cstar], p - 2, p)
        S4 = {lab: [[{}]] for lab in LABELS}
        for lab in LABELS:
            total = {}
            for b in range(D):
                total = ladd(total, lscale(K[lab][cstar][b], e[b] * inv, p, gamma[cstar] - gamma[b]), p)
            S4[lab] = [[total]]
        assert all(leval1(S4[lab][0][0], p) == crossing(lab) for lab in LABELS)
        etil, ftil, twist = twisted_pair(G0, G0inv, Gi, Giinv, e, p)
        St = contraction(K, ftil, etil, p)
        twisted = {}
        for name, Sx in (('formula_4', S4), ('solver', St)):
            seen = set()
            for (cut0, cuti, side0, sidei), (a0, ai) in aux.items():
                s0 = smith_orders(cut_matrix(Sx, 1, 0, cut0, zero_h), p, 4)
                si = smith_orders(cut_matrix(Sx, 1, 'inf', cuti, zero_h), p, 4)
                assert coordinatewise(a0, s0) and coordinatewise(ai, si)
                for t in (F(1, 2), F(1)):
                    assert critical_pressure(ar, a0, ai, 0, 0, t)[0] >= critical_pressure(ar, s0, si, 0, 0, t)[1]
                    assert critical_pressure(ar, s0, si, 0, 0, t)[0] > 0
                seen.add((tuple(s0), tuple(si)))
            twisted[name] = [list(map(list, pr)) for pr in sorted(seen)]
        out['F_%d' % p] = {'aux_profiles_zero_infinity': [list(map(list, pr)) for pr in profiles],
                           'constant_functionals_with_endpoint_E': n_phi,
                           'constant_functionals_failing_per_place_domination': per_place,
                           'constant_functionals_with_larger_pressure_in_every_combination': pressure,
                           'mixed_boundary_contraction': mixed,
                           'c_star': cstar + 1, 'twisted_formula_4_and_solver_profiles': twisted,
                           'solver_twisted_e': [sorted(x.items()) for x in etil],
                           'solver_twisted_f': [sorted(x.items()) for x in ftil]}
    assert out['F_3']['mixed_boundary_contraction']['G_f_minus_G_aux_over_t'] > 0
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', type=int, default=10, help='cases per (p, D, gauge mode)')
    args = parser.parse_args()
    ar = Dyadic(80)
    cover = replay_cover(ar)
    rng = random.Random(20261031)
    records, stats = [], []
    plan = [(p, D, mode) for p in (2, 3, 5) for D in (1, 2, 3)
            for mode in ('none', 'diagonal', 'framed') if not (D == 1 and mode != 'none')]
    for p, D, mode in plan:
        for trial in range(args.cases):
            minors = D == 2 and trial < 2 and p != 5
            record, st = run_case(rng, ar, p, D, mode, minors)
            records.append(record)
            stats.append((mode, st))
    for p in (2, 3):
        for trial in range(args.cases):
            record, st = conjugated_case(rng, ar, p)
            records.append(record)
            stats.append(('conjugated', st))
    conj = [r for r in records if r['gauge_mode'] == 'conjugated']
    conj_stats = [st for mode, st in stats if mode == 'conjugated']
    needed = [r for r, st in zip(conj, conj_stats)
              if st['pressure_failures'] and r['gauges_non_diagonal'] and r['twisted_pair_nonconstant']]
    assert needed, 'no conjugated control where f^T K e has larger pressure under a nonconstant twist'
    rand_stats = [st for mode, st in stats if mode in ('diagonal', 'framed')]
    strict = sum(a['strict_domination'] for r in records for a in r['assignments'])
    negative = explicit_negative_control(ar)
    assert strict > 0
    assert COUNT['rectangular_all_minor_replays'] > 0
    result = {
        'status': 'Every auxiliary control passed: twisted scalar contraction found, coordinatewise and prefix '
                  'domination, exact partition comparison, direct normalized audit, outward pressure floor, '
                  'cover replay; the Section 5 seed defeats every constant contraction per place.',
        'scope': 'Random matrix-valued binary Laurent crossings over F_2, F_3, F_5 with D<=3, random integer '
                 'physical heights at both places, independent Laurent auxiliary conjugation gauges with '
                 'Laurent inverses, auxiliary pair on either side of each cut, hypothesis K(1)e=Ee. '
                 'Finite controls of the transcription; the general statement (rational gauges, real '
                 'heights, every field and D) is the written proof. No improved factoring exponent.',
        'cases': len(records),
        'cut_assignments': 2 * len(records),
        'checked_assignment_side_combinations': sum(len(r['assignments']) for r in records),
        'rectangular_profiles_by_local_elimination': 2 * sum(len(r['assignments']) for r in records),
        'rectangular_profiles_also_by_all_minors': COUNT['rectangular_all_minor_replays'],
        'scalar_contraction_profiles_by_all_minors': COUNT['scalar_all_minor_replays'],
        'strict_domination_combinations': strict,
        'random_gauged_cases': len(rand_stats),
        'random_gauged_combinations_where_constant_contraction_fails_per_place':
            sum(st['profile_failures'] for st in rand_stats),
        'conjugated_cases': len(conj),
        'conjugated_cases_non_diagonal_with_nonconstant_twist': sum(
            r['gauges_non_diagonal'] and r['twisted_pair_nonconstant'] for r in conj),
        'conjugated_combinations_where_constant_contraction_fails_per_place':
            sum(st['profile_failures'] for st in conj_stats),
        'conjugated_combinations_where_constant_contraction_has_larger_pressure':
            sum(st['pressure_failures'] for st in conj_stats),
        'conjugated_cases_with_larger_constant_pressure_nondiagonal_nonconstant_twist': len(needed),
        'combinations_per_conjugated_case': conj_stats[0]['checked'],
        'explicit_section_5_seed': negative,
        'compact_certificate_replay': cover,
        'critical_lower_bound': 'G_aux(t) >= G_contraction(t) > exp(-max(a,b) t)/28',
        'controls': records,
        'source_sha256': dict({name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'charged_gaussian_floor.py',
                                            'check_global_homogeneous_boolean_filter_floor.py')},
                              **{'evidence/binary-fixed-ray-floor.json': sha256(CERTIFICATE.read_bytes()).hexdigest()}),
    }
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    summary = {k: v for k, v in result.items() if k not in ('controls', 'source_sha256')}
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    try:
        main()
    except AssertionError as error:
        print('FAILED:', error, file=sys.stderr)
        raise
