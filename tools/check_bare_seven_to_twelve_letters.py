#!/usr/bin/env python3
"""Exact checks for research/bare-seven-to-twelve-letters.md.

The note proves box-minor Kronecker budgets for the bare two-place Smith
profiles of a strict polynomial crossing.  For every 1 <= alpha <= beta <= q,
with N = alpha*beta and r = alpha^2,

    sum_{j>=0} (mu_{j+1}-mu_j) (N-j r)^+  +  sum_{i<=N} nu_i  <=  N M,   (C)
    sum_{j>=0} (nu_{j+1}-nu_j) (N-j r)^+  +  sum_{i<=N} mu_i  <=  N M,   (R)

with mu_0 = 0 and increasing profiles.  It also proves the square-box
adjugate inequalities (N = alpha^2, 2 <= alpha <= q)

    (N-1) mu_2 <= N((N-1)M - sum_{i<N} nu_i)   and the mirror (mu <-> nu).

This script verifies, in exact rational arithmetic:

1. stored entropy-dual certificates (probabilities and nonnegative
   multipliers on the box budgets only) for q = 7..12, which prove
   rho_bare > 1/sqrt(2) for every real profile satisfying the budgets, at
   every degree and scale;
2. that the formal seven-letter optimizer of the shared-minor relaxation
   violates box budgets and both adjugate inequalities;
3. an integer witness profile for each q that satisfies every box budget and
   every adjugate inequality, which brackets the relaxation value from above
   (its pressure is evaluated in floating point and labelled as such);
4. an exact thirteen-letter integer profile satisfying every box budget and
   every adjugate inequality with B(s) < 13^(2 sqrt 2) at an explicit
   rational s, so these constraints do not close q = 13;
5. literal controls over F_p (consistency checks, not part of any proof):
   scalar crossings with rank-one jets, fixed-auxiliary crossings with both
   placements of the auxiliary indices, and unimodular nilpotent products.
   Complete local Smith profiles of both cuts and of C- and R-boxes are
   computed and checked against Lemma 1, interlacing, the box determinant
   budgets, (C), (R), and both steps of the adjugate refinement on 4x4
   boxes.  A platform-independent generator makes the counts reproducible.

The certificate data live in evidence/bare-seven-to-twelve-letters.json.
The stored derived fields are recomputed and compared (floats, which are
display values only, with a tolerance); any mismatch or failed inequality
exits with status 1.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "evidence" / "bare-seven-to-twelve-letters.json"


class CheckError(Exception):
    pass


def fail(message: str) -> None:
    raise CheckError(message)


# ---------------------------------------------------------------------------
# The constraint family
# ---------------------------------------------------------------------------

def box_weights(n: int, N: int, r: int) -> list[int]:
    """Coefficient of the i-th order (i = 1..n) in sum_j (x_{j+1}-x_j)(N-jr)^+."""
    c = [max(0, N - j * r) for j in range(n + 1)]
    return [c[i - 1] - c[i] for i in range(1, n + 1)]


def constraints(q: int) -> list[tuple[tuple[str, int, int], list[int], list[int], int]]:
    """(name, coefficients on mu, coefficients on nu, N); inequality <= N*M."""
    n = q * q
    out = []
    for alpha in range(1, q + 1):
        for beta in range(alpha, q + 1):
            N = alpha * beta
            r = alpha * alpha
            w = box_weights(n, N, r)
            v = [1 if i <= N else 0 for i in range(1, n + 1)]
            out.append((("C", alpha, beta), w, v, N))
            out.append((("R", alpha, beta), v, w, N))
    for _, w, v, _ in out:
        if min(w) < 0 or min(v) < 0:
            fail("a budget has a negative coefficient")
    return out


def tails(vec: list) -> list:
    out = [0] * len(vec)
    s = 0
    for i in range(len(vec) - 1, -1, -1):
        s += vec[i]
        out[i] = s
    return out


def satisfies(q: int, mu: list[int], nu: list[int], M: int) -> list[tuple]:
    """Return the list of violated budgets (empty when all hold)."""
    bad = []
    for name, wm, wn, N in constraints(q):
        lhs = sum(a * b for a, b in zip(wm, mu)) + sum(a * b for a, b in zip(wn, nu))
        if lhs > N * M:
            bad.append((name, lhs, N * M))
    return bad


def profile_ok(q: int, mu: list[int]) -> bool:
    return (len(mu) == q * q and mu[0] == 0
            and all(mu[i] <= mu[i + 1] for i in range(len(mu) - 1)))


# ---------------------------------------------------------------------------
# Rigorous rational logarithm bounds
# ---------------------------------------------------------------------------

LN_TERMS = 40


def _atanh_series_bounds(z: Fr) -> tuple[Fr, Fr]:
    """Bounds for 2*atanh(z) = log((1+z)/(1-z)), 0 <= z <= 1/3."""
    if not (0 <= z <= Fr(1, 3)):
        raise ValueError("series argument out of range")
    s = Fr(0)
    p = z
    z2 = z * z
    for j in range(LN_TERMS):
        s += p / (2 * j + 1)
        p *= z2
    tail = p / ((2 * LN_TERMS + 1) * (1 - z2))
    return 2 * s, 2 * (s + tail)


LN2 = _atanh_series_bounds(Fr(1, 3))


def ln_bounds(x: Fr) -> tuple[Fr, Fr]:
    """Rational lo <= log(x) <= hi for rational x > 0."""
    x = Fr(x)
    if x <= 0:
        raise ValueError("log of a nonpositive number")
    e = x.numerator.bit_length() - x.denominator.bit_length()
    y = x / Fr(2) ** e
    while y >= 2:
        y /= 2
        e += 1
    while y < 1:
        y *= 2
        e -= 1
    lo_y, hi_y = _atanh_series_bounds((y - 1) / (y + 1))
    if e >= 0:
        return lo_y + e * LN2[0], hi_y + e * LN2[1]
    return lo_y + e * LN2[1], hi_y + e * LN2[0]


def entropy_lower(D: int, ks: list[int]) -> Fr:
    """Lower bound for -sum (k/D) log(k/D) with sum ks = D."""
    cache: dict[int, Fr] = {}
    total = Fr(0)
    for k in ks:
        if k not in cache:
            cache[k] = ln_bounds(Fr(k))[1]
        total += k * cache[k]
    return ln_bounds(Fr(D))[0] - total / D


# ---------------------------------------------------------------------------
# 1. Entropy-dual certificates
# ---------------------------------------------------------------------------

def verify_certificate(q: int, cert: dict) -> dict:
    n = q * q
    D = int(cert["D"])
    ks = [int(k) for k in cert["k"]]
    if len(ks) != n or min(ks) <= 0 or sum(ks) != D:
        fail(f"q={q}: probabilities are not a positive rational distribution")
    lam = {}
    for side, alpha, beta, value in cert["multipliers"]:
        v = Fr(value)
        if v < 0:
            fail(f"q={q}: negative multiplier")
        lam[(side, int(alpha), int(beta))] = v
    cons = {c[0]: c for c in constraints(q)}
    for key in lam:
        if key not in cons:
            fail(f"q={q}: multiplier on an unknown budget {key}")
    # Tail sums of the probability vector (the same vector is used for both
    # profiles: the budget family is symmetric under exchanging mu and nu).
    Pi = [Fr(t, D) for t in tails(ks)]
    tail_mu = {key: tails(cons[key][1]) for key in lam}
    tail_nu = {key: tails(cons[key][2]) for key in lam}
    for i in range(1, n):
        lhs_mu = sum(v * tail_mu[key][i] for key, v in lam.items())
        lhs_nu = sum(v * tail_nu[key][i] for key, v in lam.items())
        if lhs_mu < Pi[i] or lhs_nu < Pi[i]:
            fail(f"q={q}: dual feasibility fails at increment {i + 1}")
    h = sum(v * cons[key][3] for key, v in lam.items())
    hp = max(Fr(1), h)
    H = 2 * entropy_lower(D, ks)
    lnq_hi = ln_bounds(Fr(q))[1]
    if not (H > 0 and H * H > 8 * hp * hp * lnq_hi * lnq_hi):
        fail(f"q={q}: certificate does not exceed q^(2 sqrt 2)")
    floor = H / (4 * hp * lnq_hi)
    # A decimal floor that is still above 1/sqrt(2) exactly.
    dec = Fr(math.floor(floor * 10**6), 10**6)
    if not (dec <= floor and 2 * dec * dec > 1):
        fail(f"q={q}: decimal floor bookkeeping")
    active = sorted(f"{k[0]}({k[1]},{k[2]})" for k, v in lam.items() if v > 0)
    return {
        "q": q,
        "dual_value_h": str(h),
        "h_minus_one_float": float(h - 1),
        "entropy_lower_float": float(H),
        "certified_rho_floor": str(dec),
        "active_budgets": active,
    }


# ---------------------------------------------------------------------------
# 2. The formal seven-letter optimizer of the shared-minor relaxation
# ---------------------------------------------------------------------------

def formal_seven() -> dict:
    q, M = 7, 12
    mu = [0] + [7] * 6 + [13] * 42
    old_ok = (sum(mu) <= 49 * M and sum(mu[:7]) + sum(mu[:7]) <= 7 * M)
    bad = satisfies(q, mu, mu, M)
    adj = adjugate_violations(q, mu, mu, M)
    if not old_ok or not bad:
        fail("formal seven-letter optimizer bookkeeping")
    worst = max(bad, key=lambda b: Fr(b[1], b[2]))
    return {"old_budgets_hold": old_ok, "violated_box_budgets": len(bad),
            "worst": f"{worst[0][0]}({worst[0][1]},{worst[0][2]}): {worst[1]} > {worst[2]}",
            "violated_adjugate_inequalities": [f"{a}: {l} > {r}" for a, l, r in adj]}


# ---------------------------------------------------------------------------
# 3./4. Feasible integer profiles
# ---------------------------------------------------------------------------

def check_witness(q: int, wit: dict) -> dict:
    M = int(wit["M"])
    mu = [int(x) for x in wit["mu"]]
    if not profile_ok(q, mu):
        fail(f"q={q}: witness is not a normalized increasing profile")
    bad = satisfies(q, mu, mu, M) + adjugate_violations(q, mu, mu, M)
    if bad:
        fail(f"q={q}: witness violates {bad[0]}")
    t = float(wit["t"])
    # FLOAT evaluation (upper bracket only; plays no role in any proof).
    sa = math.fsum(math.exp(-t * m / M) for m in mu)
    rho = (t + 2 * math.log(sa)) / (4 * math.log(q))
    return {"q": q, "M": M, "satisfies_box_and_adjugate": True,
            "float_rho_of_feasible_profile": round(rho, 7)}


def check_thirteen(target: dict) -> dict:
    q = 13
    M = int(target["M"])
    mu = [int(x) for x in target["mu"]]
    s = Fr(target["s"])
    if not profile_ok(q, mu) or s < 1:
        fail("q=13 target is not normalized")
    bad = satisfies(q, mu, mu, M) + adjugate_violations(q, mu, mu, M)
    if bad:
        fail(f"q=13 target violates {bad[0]}")
    counts: dict[int, int] = {}
    for m in mu:
        counts[m] = counts.get(m, 0) + 1
    # Exact pressure, then an exact dyadic upper rounding before the logarithm.
    U = sum(c * s ** (M - m) for m, c in counts.items())
    V = sum(c / s ** m for m, c in counts.items())
    B = U * V
    scale = 2 ** 96
    B_up = Fr(-((-B.numerator * scale) // B.denominator), scale)
    if B_up < B:
        fail("dyadic rounding")
    _, hi = ln_bounds(B_up)
    lnq_lo, _ = ln_bounds(Fr(q))
    # B < q^(2 sqrt 2)  <=>  log B < 2 sqrt 2 log q  <=  hi^2 < 8 lnq_lo^2.
    if not (hi > 0 and hi * hi < 8 * lnq_lo * lnq_lo):
        fail("q=13 target does not beat q^(2 sqrt 2)")
    adj_slack = min(N * ((N - 1) * M - sum(mu[:N - 1])) - (N - 1) * mu[1]
                    for N in (a * a for a in range(2, q + 1)))
    return {"M": M, "s": str(s), "rho_upper_float": float(hi / (4 * lnq_lo)),
            "satisfies_box_and_adjugate": True, "minimum_adjugate_slack": adj_slack,
            "distinct_orders": sorted(set(mu))[:6] + ["..."]}


# ---------------------------------------------------------------------------
# 5. Literal controls over F_p
# ---------------------------------------------------------------------------

class Rng:
    """Platform-independent 64-bit linear congruential generator.

    Python's random module may change its derived methods between versions;
    this keeps the literal controls, and hence the stored statistics,
    reproducible everywhere.
    """

    def __init__(self, seed: int) -> None:
        self.state = seed % 2**64

    def below(self, n: int) -> int:
        self.state = (self.state * 6364136223846793005 + 1442695040888963407) % 2**64
        return (self.state >> 33) % n

    def choice(self, seq):
        return seq[self.below(len(seq))]

    def sample(self, seq, k: int) -> list:
        pool = list(seq)
        return [pool.pop(self.below(len(pool))) for _ in range(k)]


def local_smith_orders(mat: list[list[list[int]]], p: int, prec: int) -> list[int]:
    """Smith orders at z = 0 of a rectangular matrix of truncated series mod p.

    mat[i][j] is a coefficient list.  Returns min(rows, cols) orders; orders
    >= prec are reported as prec.
    """
    n_rows, n_cols = len(mat), len(mat[0])
    A = [[(list(e) + [0] * prec)[:prec] for e in row] for row in mat]

    def val(e):
        for k, c in enumerate(e):
            if c:
                return k
        return prec

    def mul(a, b):
        out = [0] * prec
        for i, x in enumerate(a):
            if x:
                for j in range(prec - i):
                    if b[j]:
                        out[i + j] = (out[i + j] + x * b[j]) % p
        return out

    def inv_unit(a):
        inv0 = pow(a[0], p - 2, p)
        out = [0] * prec
        out[0] = inv0
        for k in range(1, prec):
            acc = sum(a[i] * out[k - i] for i in range(1, k + 1)) % p
            out[k] = (-acc * inv0) % p
        return out

    orders: list[int] = []
    rows = list(range(n_rows))
    cols = list(range(n_cols))
    size = min(n_rows, n_cols)
    while len(orders) < size:
        best = None
        for i in rows:
            for j in cols:
                v = val(A[i][j])
                if best is None or v < best[0]:
                    best = (v, i, j)
        v, pi, pj = best
        if v >= prec:
            orders.extend([prec] * (size - len(orders)))
            break
        orders.append(v)
        u = inv_unit(A[pi][pj][v:] + [0] * v)
        for i in rows:
            if i == pi or val(A[i][pj]) >= prec:
                continue
            f = mul(A[i][pj][v:] + [0] * v, u)  # entry / pivot
            for j in cols:
                prod = mul(f, A[pi][j])
                A[i][j] = [(a - b) % p for a, b in zip(A[i][j], prod)]
        rows.remove(pi)
        cols.remove(pj)
    return sorted(orders)


# Exact polynomial arithmetic over F_p for the adjugate controls.

def _padd(a, b, p):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] = x
    for i, x in enumerate(b):
        out[i] = (out[i] + x) % p
    return out


def _pmul(a, b, p):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] = (out[i + j] + x * y) % p
    return out


def _pdet(m, p):
    """Determinant of a small square matrix of polynomials (Laplace)."""
    size = len(m)
    if size == 1:
        return list(m[0][0])
    total: list[int] = []
    for j in range(size):
        if not any(m[0][j]):
            continue
        minor = [row[:j] + row[j + 1:] for row in m[1:]]
        term = _pmul(m[0][j], _pdet(minor, p), p)
        if j % 2:
            term = [(-x) % p for x in term]
        total = _padd(total, term, p)
    return total


def _deg(a):
    for k in range(len(a) - 1, -1, -1):
        if a[k]:
            return k
    return -1


def _ord(a):
    for k, x in enumerate(a):
        if x:
            return k
    return None


def _adjugate(m, p):
    size = len(m)
    adj = [[None] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            minor = [row[:j] + row[j + 1:] for k, row in enumerate(m) if k != i]
            c = _pdet(minor, p)
            if (i + j) % 2:
                c = [(-x) % p for x in c]
            adj[j][i] = c
    return adj


def adjugate_violations(q: int, mu: list[int], nu: list[int], M: int) -> list[tuple]:
    """Square-box adjugate inequalities of Section 5 and their mirrors.

    (N-1) mu_2 <= N((N-1)M - sum_{i<N} nu_i)  and the same with mu, nu
    exchanged, for N = alpha^2 and alpha = 2..q.
    """
    bad = []
    for alpha in range(2, q + 1):
        N = alpha * alpha
        for side, a, b in (("C", mu, nu), ("R", nu, mu)):
            lhs = (N - 1) * a[1]
            rhs = N * ((N - 1) * M - sum(b[:N - 1]))
            if lhs > rhs:
                bad.append((f"adj-{side}({alpha})", lhs, rhs))
    return bad


def make_crossing(rng: Rng, q: int, p: int, aux: int, M: int,
                  k_top: int, k_bot: int, place_R: str, place_C: str):
    """Random scalar or fixed-auxiliary crossing with rank-one jets.

    Coefficients K_j[(x,y,u,v,l,m)], j = 0..M, with auxiliary indices (l,m)
    attached to one physical side of each cut (place_R in {"xy","uv"},
    place_C in {"xu","yv"}).  The top k_top coefficients form an R-cut
    rank-one jet H_R = a(w) b(w)^T mod w^k_top and the bottom k_bot
    coefficients a C-cut rank-one jet; the auxiliary indices travel with the
    chosen side.  One free coefficient restores the endpoint of the scalar
    boundary contraction, sum_{l,m} f_l e_m K(1)[..,l,m] = E.
    """
    L = range(q)
    A2 = [(l, m) for l in range(aux) for m in range(aux)]
    keys = [(x, y, u, v, l, m) for x in L for y in L for u in L for v in L for (l, m) in A2]
    coeffs = [{k: rng.below(p) for k in keys} for _ in range(M + 1)]

    def split(key, place):
        x, y, u, v, l, m = key
        sides = {"xy": ((x, y), (u, v)), "uv": ((u, v), (x, y)),
                 "xu": ((x, u), (y, v)), "yv": ((y, v), (x, u))}[place]
        return sides[0] + (l, m), sides[1]

    def jet(k, place, positions):
        firsts = [dict() for _ in range(k)]
        seconds = [dict() for _ in range(k)]
        for i, j in enumerate(positions):
            for key in keys:
                f, g = split(key, place)
                total = 0
                for s_ in range(i + 1):
                    if f not in firsts[s_]:
                        firsts[s_][f] = rng.below(p)
                    if g not in seconds[i - s_]:
                        seconds[i - s_][g] = rng.below(p)
                    total += firsts[s_][f] * seconds[i - s_][g]
                coeffs[j][key] = total % p

    jet(k_top, place_R, [M - i for i in range(k_top)])
    jet(k_bot, place_C, list(range(k_bot)))
    free = [j for j in range(k_bot, M + 1 - k_top)]
    jf = rng.choice(free)
    f_vec = [1] + [rng.below(p) for _ in range(aux - 1)]
    e_vec = [1] + [rng.below(p) for _ in range(aux - 1)]
    for x in L:
        for y in L:
            for u in L:
                for v in L:
                    tot = sum(f_vec[l] * e_vec[m] * coeffs[j][(x, y, u, v, l, m)]
                              for (l, m) in A2 for j in range(M + 1)
                              if (j, l, m) != (jf, 0, 0))
                    coeffs[jf][(x, y, u, v, 0, 0)] = (int(x == v and y == u) - tot) % p

    def contracted(x, y, u, v):
        return [sum(f_vec[l] * e_vec[m] * coeffs[j][(x, y, u, v, l, m)] for (l, m) in A2) % p
                for j in range(M + 1)]

    return coeffs, keys, A2, contracted


def make_unimodular(rng: Rng, q: int, p: int, M: int):
    """Scalar crossing K_C(z) = E_C prod_k (z I + (1-z) N_k), N_k nilpotent.

    Each factor has determinant z^n and value I at z = 1, so the column
    determinant is a monomial of degree nM: the determinant budget
    (C) with alpha = beta = q is an equality (the inverse-degree setting).
    """
    n = q * q
    L = range(q)
    pairs = [(i, j) for i in L for j in L]

    def matmul(a, b):
        return [[sum(a[i][k] * b[k][j] for k in range(n)) % p for j in range(n)]
                for i in range(n)]

    def random_invertible():
        while True:
            S = [[rng.below(p) for _ in range(n)] for _ in range(n)]
            # Gaussian elimination for the inverse over F_p.
            A = [row[:] + [int(i == j) for j in range(n)] for i, row in enumerate(S)]
            ok = True
            for c in range(n):
                piv = next((r for r in range(c, n) if A[r][c] % p), None)
                if piv is None:
                    ok = False
                    break
                A[c], A[piv] = A[piv], A[c]
                inv = pow(A[c][c], p - 2, p)
                A[c] = [(x * inv) % p for x in A[c]]
                for r in range(n):
                    if r != c and A[r][c]:
                        f = A[r][c]
                        A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
            if ok:
                return S, [row[n:] for row in A]

    ident = [[int(i == j) for j in range(n)] for i in range(n)]
    poly = [ident]  # coefficient matrices of prod_k (z I + (1-z) N_k)
    for _ in range(M):
        S, S_inv = random_invertible()
        U = [[rng.below(p) if j > i else 0 for j in range(n)] for i in range(n)]
        N_k = matmul(matmul(S, U), S_inv)
        low = N_k
        high = [[(ident[i][j] - N_k[i][j]) % p for j in range(n)] for i in range(n)]
        new = [[[0] * n for _ in range(n)] for _ in range(len(poly) + 1)]
        for d, C in enumerate(poly):
            for shift, F in ((0, low), (1, high)):
                prod = matmul(C, F)
                for i in range(n):
                    for j in range(n):
                        new[d + shift][i][j] = (new[d + shift][i][j] + prod[i][j]) % p
        poly = new
    E_C = [[int(r[0] == c[1] and c[0] == r[1]) for c in pairs] for r in pairs]
    KC = [matmul(E_C, C) for C in poly]
    keys = [(x, y, u, v, 0, 0) for x in L for y in L for u in L for v in L]
    coeffs = [{(x, y, u, v, 0, 0): KC[j][pairs.index((x, u))][pairs.index((y, v))]
               for (x, y, u, v, _, _) in keys} for j in range(M + 1)]

    def contracted(x, y, u, v):
        return [coeffs[j][(x, y, u, v, 0, 0)] for j in range(M + 1)]

    return coeffs, keys, [(0, 0)], contracted


def literal_family(rng: Rng, q: int, p: int, aux: int, trials: int,
                   degrees: list[int], max_jet: int, stats: dict) -> None:
    n = q * q
    L = range(q)
    pairs = [(i, j) for i in L for j in L]
    for _ in range(trials):
        M = rng.choice(degrees)
        if max_jet == 0:
            place_R, place_C = "xy", "xu"
            coeffs, keys, A2, contracted = make_unimodular(rng, q, p, M)
        else:
            k_top = 1 + rng.below(max_jet)
            k_bot = 1 + rng.below(max_jet)
            while k_top + k_bot > M:
                if k_top >= k_bot:
                    k_top -= 1
                else:
                    k_bot -= 1
            place_R = rng.choice(["xy", "uv"])
            place_C = rng.choice(["xu", "yv"])
            coeffs, keys, A2, contracted = make_crossing(
                rng, q, p, aux, M, k_top, k_bot, place_R, place_C)
        if all(coeffs[M][k] == 0 for k in keys) or all(coeffs[0][k] == 0 for k in keys):
            continue
        prec = n * M + 2

        def entry(key, reverse):
            e = [coeffs[j][key] for j in range(M + 1)]
            return e[::-1] if reverse else e

        def cut(place, reverse):
            mat = []
            for a in pairs:
                for (l, m) in A2:
                    row = []
                    for b in pairs:
                        x, y, u, v = {"xy": (a[0], a[1], b[0], b[1]),
                                      "uv": (b[0], b[1], a[0], a[1]),
                                      "xu": (a[0], b[0], a[1], b[1]),
                                      "yv": (b[0], a[0], b[1], a[1])}[place]
                        row.append(entry((x, y, u, v, l, m), reverse))
                    mat.append(row)
            return mat

        mu = local_smith_orders(cut(place_R, True), p, prec)
        nu = local_smith_orders(cut(place_C, False), p, prec)
        if max(mu) >= prec or max(nu) >= prec or mu[0] or nu[0]:
            continue  # singular cut or uncancelled monomial: outside the interface
        stats["crossings"] += 1
        if aux > 1:
            stats["auxiliary_crossings"] += 1
        if mu[1] > 0 and nu[1] > 0:
            stats["unique_minimum_in_both_cuts"] += 1
        bad = satisfies(q, mu, nu, M) + adjugate_violations(q, mu, nu, M)
        if bad:
            fail(f"literal control violates {bad[0]}: mu={mu} nu={nu} M={M}")
        if sum(nu) == n * M:
            # Monomial column determinant: the adjugate inequality at
            # alpha = q reduces to the sharpened inverse-degree bound.
            stats["monomial_column_determinant"] += 1
            if (n - 1) * mu[1] > n * (nu[-1] - M):
                fail("literal control: sharpened inverse-degree bound fails")
        elif aux == 1 and mu[1] > 2 * (nu[-1] - M):
            # The (scalar) inverse-degree conclusion needs its hypothesis.
            stats["inverse_degree_form_fails_without_hypothesis"] += 1
        for alpha in range(1, q + 1):
            for beta in range(alpha, q + 1):
                N, r = alpha * beta, alpha * alpha
                for side, wm, wn in (("C", box_weights(n, N, r), [1] * N + [0] * (n - N)),
                                     ("R", [1] * N + [0] * (n - N), box_weights(n, N, r))):
                    lhs = sum(a * b for a, b in zip(wm, mu)) + sum(a * b for a, b in zip(wn, nu))
                    if lhs == N * M:
                        stats["tight_budgets"] += 1
                    if 1 < alpha < beta:
                        slack = Fr(N * M - lhs, N * M)
                        if stats["min_new_budget_slack"] is None or slack < Fr(stats["min_new_budget_slack"]):
                            stats["min_new_budget_slack"] = str(slack)
                S1 = rng.sample(list(L), alpha)
                S2 = rng.sample(list(L), beta)

                def box(rowf, colf, reverse):
                    out = []
                    for rk in rowf:
                        row = []
                        for ck in colf:
                            e = contracted(*rk(ck))
                            row.append(e[::-1] if reverse else e)
                        out.append(row)
                    return out

                # C-box: rows (x,u) in S1xS2, columns (y,v) in S2xS1.
                c_rows = [(lambda c, x=x, u=u: (x, c[0], u, c[1])) for x in S1 for u in S2]
                c_cols = [(y, v) for y in S2 for v in S1]
                # R-box: rows (x,y) in S1xS2, columns (u,v) in S2xS1.
                r_rows = [(lambda c, x=x, y=y: (x, y, c[0], c[1])) for x in S1 for y in S2]
                r_cols = [(u, v) for u in S2 for v in S1]
                for kind, rws, cls, layer_prof, layer_rev, inter_prof in (
                        ("C", c_rows, c_cols, mu, True, nu),
                        ("R", r_rows, r_cols, nu, False, mu)):
                    layered = local_smith_orders(box(rws, cls, layer_rev), p, prec)
                    other = local_smith_orders(box(rws, cls, not layer_rev), p, prec)
                    for j in range(n):
                        count = sum(1 for k in layered if k < layer_prof[j])
                        if count > j * r:
                            fail(f"literal control: Lemma 1 fails on a {kind}-box")
                        if count == j * r and 0 < j * r < N:
                            stats["tight_lemma1_layers"] += 1
                            if 1 < alpha < beta:
                                stats["tight_lemma1_layers_1<alpha<beta"] += 1
                    if sum(other) < sum(inter_prof[:N]):
                        fail(f"literal control: {kind}-box determinantal interlacing fails")
                    if sum(layered) + sum(other) > N * M:
                        fail(f"literal control: {kind}-box determinant budget fails")
                    stats["boxes_checked"] += 1
                if alpha == beta == 2:
                    adjugate_literal(q, p, M, mu, nu, S1, contracted, stats)


def adjugate_literal(q, p, M, mu, nu, S, contracted, stats) -> None:
    """The two steps of the square-box adjugate refinement, on 4x4 boxes."""
    N = len(S) ** 2
    for kind in ("C", "R"):
        if kind == "C":
            # rows (x,u), columns (y,v); adjugate indexed [(y,v)][(x,u)];
            # realigned rows (y,x), columns (v,u).
            P = [[contracted(x, y, u, v) for y in S for v in S] for x in S for u in S]
            idx = [(a, b) for a in S for b in S]
            adj = _adjugate(P, p)
            Rl = [[adj[idx.index((y, v))][idx.index((x, u))] for v in S for u in S]
                  for y in S for x in S]
            layer, other = mu, nu
        else:
            # rows (x,y), columns (u,v); adjugate indexed [(u,v)][(x,y)];
            # realigned rows (u,x), columns (v,y).
            P = [[contracted(x, y, u, v) for u in S for v in S] for x in S for y in S]
            idx = [(a, b) for a in S for b in S]
            adj = _adjugate(P, p)
            Rl = [[adj[idx.index((u, v))][idx.index((x, y))] for v in S for y in S]
                  for u in S for x in S]
            layer, other = nu, mu
        D = _pdet(Rl, p)
        if sum(D) % p == 0:
            fail(f"literal control: realigned adjugate of the {kind}-box is singular at one")
        deg, order = _deg(D), _ord(D)
        S_other = sum(other[:N - 1])
        if kind == "C":
            # divisibility at zero by the (N-1)-minors, and the Kronecker
            # layer at infinity.
            if order < N * S_other or N * (N - 1) * M - deg < (N - 1) * layer[1]:
                fail("literal control: C-box adjugate step fails")
        else:
            if order < (N - 1) * layer[1] or deg > N * ((N - 1) * M - S_other):
                fail("literal control: R-box adjugate step fails")
        stats["adjugate_boxes_checked"] += 1


def literal_controls(seed: int) -> dict:
    """Literal crossings over F_p: scalar and fixed-auxiliary families."""
    rng = Rng(seed)
    out = {}
    families = (
        # name, q, p, auxiliary dimension, trials, degrees, maximum jet length
        ("scalar_q3_F5", 3, 5, 1, 24, [2, 3], 1),
        ("scalar_jets_q3_F3", 3, 3, 1, 24, [3, 4], 3),
        ("scalar_jets_q2_F2", 2, 2, 1, 40, [3, 4, 5], 4),
        ("aux2_jets_q2_F3", 2, 3, 2, 40, [2, 3, 4], 3),
        ("aux2_jets_q3_F3", 3, 3, 2, 4, [2, 3], 2),
        ("unimodular_q3_F3", 3, 3, 1, 8, [1, 2, 3], 0),
        ("unimodular_q2_F5", 2, 5, 1, 16, [1, 2, 3, 4], 0),
    )
    keys = ("crossings", "auxiliary_crossings", "unique_minimum_in_both_cuts",
            "boxes_checked", "adjugate_boxes_checked", "tight_lemma1_layers",
            "tight_lemma1_layers_1<alpha<beta", "tight_budgets",
            "monomial_column_determinant", "inverse_degree_form_fails_without_hypothesis")
    total = dict.fromkeys(keys, 0)
    for name, q, p, aux, trials, degrees, max_jet in families:
        stats = dict.fromkeys(keys, 0)
        stats["min_new_budget_slack"] = None
        literal_family(rng, q, p, aux, trials, degrees, max_jet, stats)
        if stats["crossings"] < trials // 2:
            fail(f"literal control {name}: too few full-rank samples")
        out[name] = stats
        for k in keys:
            total[k] += stats[k]
    if not (total["tight_lemma1_layers_1<alpha<beta"] and total["tight_budgets"]
            and total["auxiliary_crossings"] and total["adjugate_boxes_checked"]):
        fail("literal controls lack power (no tight layer, budget, auxiliary or adjugate case)")
    out["total"] = total
    return out


def mutation_controls(record: dict) -> dict:
    """Corrupted certificates must be rejected."""
    cert = record["certificates"]["7"]
    out = {}
    # (a) Keep only the budgets known before this note (determinants, shared
    #     q-minor and square prefixes): the dual becomes infeasible.
    old = dict(cert)
    old["multipliers"] = [m for m in cert["multipliers"]
                          if int(m[1]) in (1,) or int(m[1]) == int(m[2])]
    # (b) Halve the largest multiplier.
    half = dict(cert)
    big = max(range(len(cert["multipliers"])), key=lambda i: Fr(cert["multipliers"][i][3]))
    half["multipliers"] = [list(m) for m in cert["multipliers"]]
    half["multipliers"][big][3] = str(Fr(half["multipliers"][big][3]) / 2)
    # (c) Move probability mass from the zero-order channel to the top one.
    skew = dict(cert)
    ks = list(cert["k"])
    move = ks[0] // 10
    ks[0] -= move
    ks[-1] += move
    skew["k"] = ks
    for name, bad in (("previous_budgets_only", old), ("halved_multiplier", half),
                      ("shifted_probability", skew)):
        try:
            verify_certificate(7, bad)
        except CheckError:
            out[name] = "rejected"
        else:
            fail(f"mutation control {name} was accepted")
    return out


# ---------------------------------------------------------------------------

def compare(stored, fresh, path: str):
    """Exact comparison, except floats, which are display values only."""
    if isinstance(stored, float) or isinstance(fresh, float):
        if not (isinstance(stored, (int, float)) and isinstance(fresh, (int, float))):
            return path
        return None if math.isclose(stored, fresh, rel_tol=1e-6, abs_tol=1e-6) else path
    if isinstance(stored, dict) and isinstance(fresh, dict):
        if stored.keys() != fresh.keys():
            return path
        for key in stored:
            bad = compare(stored[key], fresh[key], f"{path}.{key}")
            if bad:
                return bad
        return None
    if isinstance(stored, list) and isinstance(fresh, list):
        if len(stored) != len(fresh):
            return path
        for i, (a, b) in enumerate(zip(stored, fresh)):
            bad = compare(a, b, f"{path}[{i}]")
            if bad:
                return bad
        return None
    return None if stored == fresh else path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true",
                        help="rewrite the derived fields of the evidence record")
    args = parser.parse_args()

    record = json.loads(RECORD.read_text())
    derived = {"certificates": [], "witnesses": []}
    for q in range(7, 13):
        res = verify_certificate(q, record["certificates"][str(q)])
        derived["certificates"].append(res)
        print(f"q={q}: rho_bare >= {res['certified_rho_floor']} > 1/sqrt(2) "
              f"(h-1 = {res['h_minus_one_float']:.2e}; active {len(res['active_budgets'])} budgets)")
    for q in range(7, 13):
        res = check_witness(q, record["certificates"][str(q)]["witness"])
        derived["witnesses"].append(res)
        print(f"q={q}: a feasible integer profile has float rho {res['float_rho_of_feasible_profile']}")
    derived["formal_seven"] = formal_seven()
    print("formal q=7 optimizer:", derived["formal_seven"]["worst"],
          "; adjugate:", derived["formal_seven"]["violated_adjugate_inequalities"])
    derived["thirteen"] = check_thirteen(record["thirteen_target"])
    print(f"q=13: feasible profile with exact B < 13^(2 sqrt 2); float rho <= "
          f"{derived['thirteen']['rho_upper_float']:.6f}")
    derived["literal_controls"] = literal_controls(seed=20261031)
    for name, stats in derived["literal_controls"].items():
        print(f"literal controls {name}:", stats)
    derived["mutation_controls"] = mutation_controls(record)
    print("mutation controls:", derived["mutation_controls"])

    if args.write:
        record["derived"] = derived
        RECORD.write_text(json.dumps(record, indent=1) + "\n")
        print("wrote", RECORD.relative_to(ROOT))
    else:
        stored = record.get("derived")
        if stored is None:
            fail("the record has no derived fields; run with --write")
        mismatch = compare(stored, json.loads(json.dumps(derived)), "derived")
        if mismatch:
            fail(f"stored derived field differs from the recomputation: {mismatch}")
    print("PASS: box-minor Kronecker budgets close the bare criterion for q = 7..12")


if __name__ == "__main__":
    try:
        main()
    except CheckError as exc:
        print("FAIL:", exc)
        sys.exit(1)
