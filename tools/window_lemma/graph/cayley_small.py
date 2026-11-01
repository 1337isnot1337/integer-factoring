#!/usr/bin/env python3
"""Exact small-|Delta| experiment for the graph view of WL(Y, W).

For seeded semiprimes N = 3 mod 4 (Delta = -N, fundamental) whose class group is
CYCLIC (so Cl = Z/h with a PARI generator), it lists every reduced form as an
index i in Z/h (form = red(gen^i)), marks S = {i : A(i) is B_Y-smooth} (repo base,
repo primorial-style test), and measures:
  * Fourier: max_{chi != 1} |S^(chi)| against a random set of equal size;
  * single generator g_p: cycle length, first hit from the repo start K_p,
    largest S-gap over ALL starts (uniform window), clustering statistics,
    and the cusp lemma (A*2p < sqrt|Delta| and s in S  =>  g_p s in S);
  * Cayley graph Cay(Z/h, {+-t_p : p in base}): exact eigenvalues
    lambda_k = sum_p 2 cos(2 pi k t_p / h), BFS ball sizes from each repo start
    against the expander-mixing bound (lambda/d)^2 h / delta, covering radius,
    exact E_pi[tau_S] of the simple random walk against (1-delta)/(delta*gap).
Needs cypari2, numpy, scipy.  Repo tools are used for the base,
the prime forms, reduce_form and least_power_above.
Usage: python cayley_small.py --bits 26 30 34 37 --per-size 4 --u 3 4 --out data/cayley.jsonl
"""
import argparse, json, math, os, random, sys, time
from math import isqrt, gcd
sys.path.insert(0, os.environ.get('IF_TOOLS', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
from binary_quadratic_forms import reduce_form                   # noqa: E402
from prime_form_factor_base import build_prime_form_factor_base  # noqa: E402
from weighted_norm_relations import least_power_above            # noqa: E402
import numpy as np                                               # noqa: E402
import cypari2                                                   # noqa: E402
import scipy.sparse as sp                                        # noqa: E402
import scipy.sparse.linalg as spl                                # noqa: E402
pari = cypari2.Pari()
pari.allocatemem(512 * 10**6, silent=True) if hasattr(pari, 'allocatemem') else None


def semiprimes(bits, count, seed):
    rng = random.Random(f"cayley-{bits}-{seed}")
    out = []
    while len(out) < count:
        a = bits // 2
        p = int(pari.nextprime(rng.getrandbits(a) | (1 << (a - 1))))
        q = int(pari.nextprime(rng.getrandbits(bits - a) | (1 << (bits - a - 1))))
        N = p * q
        if p == q or N % 4 != 3:
            continue
        D = -N
        h, cyc, gens = (lambda r: (int(r[0]), [int(c) for c in r[1]], r[2]))(pari.quadclassunit(D))
        if len(cyc) != 1:
            continue  # need cyclic Cl for the simple index map
        out.append((N, p, q, h, gens[0]))
    return out


def max_zero_run_cyclic(bits):
    """Longest cyclic run of False in a 1-D bool array (len n), n if all False."""
    if not bits.any():
        return len(bits)
    idx = np.flatnonzero(bits)
    gaps = np.diff(np.concatenate([idx, [idx[0] + len(bits)]])) - 1
    return int(gaps.max())


def analyse(N, h, gen, u_target):
    D = -N
    H = isqrt(-D // 3)
    # index map: reduced form -> i
    index, A = {}, np.zeros(h, dtype=np.int64)
    f = pari.qfbred(pari.Qfb(1, 1, (1 - D) // 4))
    for i in range(h):
        a, b, c = (int(x) for x in pari.Vec(f))
        index[(a, b, c)] = i
        A[i] = a
        f = pari.qfbcomp(f, gen)
    assert len(index) == h, "generator order != h"
    Y = max(11, int(round(H ** (1.0 / u_target))))
    base = build_prime_form_factor_base(D, Y)['records']
    primes = [r.prime for r in base]
    t = {}
    for r in base:
        form = tuple(reduce_form(r.positive))
        assert form in index, ("normal-form mismatch", form)
        t[r.prime] = index[form]
    # S: A smooth over the base primes (repo test: P^bits(H) == 0 mod A)
    P = 1
    for p in primes:
        P *= p
    e = H.bit_length()
    S = np.array([pow(P, e, int(a)) == 0 for a in A])
    nS = int(S.sum()); delta = nS / h
    rng = np.random.default_rng(12345 + N % 1000)
    R = np.zeros(h, dtype=bool); R[rng.choice(h, nS, replace=False)] = True
    FSall = np.abs(np.fft.fft(S.astype(float)))
    kstar = int(np.argmax(FSall[1:]) + 1)
    FS = FSall[kstar]
    top10 = np.sort(FSall[1:])[::-1][:10]
    nbig = int((FSall[1:] > 3 * math.sqrt(nS * math.log(h))).sum())
    FR = np.abs(np.fft.fft(R.astype(float)))[1:].max()
    rec = dict(N=N, D=D, h=h, H=H, Y=Y, u=math.log(H) / math.log(Y), m=len(primes),
               nS=nS, delta=delta, maxFourier_S=float(FS), maxFourier_rand=float(FR),
               sqrt_S_ln_h=math.sqrt(nS * math.log(h)), kstar=kstar, top10_fourier=top10.tolist(),
               n_chars_above_3x_random=nbig)
    sqrtD = math.isqrt(-D)
    gens = []
    for p in primes:
        tp = t[p]; g = gcd(tp, h); n = h // g
        K = least_power_above(p, H)
        start = (K * tp) % h
        # orbit from repo start
        orb = (start + tp * np.arange(n)) % h
        hits = np.flatnonzero(S[orb])
        first = int(hits[0]) if len(hits) else None
        # all cosets of <t_p>: max gap over all starts, clustering
        maxgap = 0; nxt = 0; cl_sizes = []; fwd_mean = []
        for c in range(g):
            seq = S[(c + tp * np.arange(n)) % h]
            maxgap = max(maxgap, max_zero_run_cyclic(seq) if seq.any() else math.inf)
            nxt += int((seq & np.roll(seq, -1)).sum())
            if seq.any() and not seq.all():
                # cluster sizes (maximal cyclic runs of True)
                rot = np.roll(seq, -int(np.flatnonzero(~seq)[0]))
                runs = np.diff(np.flatnonzero(np.diff(np.concatenate([[0], rot.astype(int), [0]]))))[::2]
                cl_sizes.extend(runs.tolist())
                # mean forward hitting time (0 if start in S) from uniform start
                idx = np.flatnonzero(seq)
                gaps = np.diff(np.concatenate([idx, [idx[0] + n]])) - 1
                fwd_mean.append(float((gaps * (gaps + 1) / 2).sum() / n))
        # cusp lemma check, both directions
        a_s = A[S]
        cusp = (a_s * 2 * p) < sqrtD
        img = S[(np.flatnonzero(S) + tp) % h]
        imgm = S[(np.flatnonzero(S) - tp) % h]
        viol = int((cusp & ~(img & imgm)).sum())
        rand_gap = math.log(max(2, n * delta)) / -math.log1p(-delta) if delta < 1 else 0
        gens.append(dict(p=p, t=tp, cycle=n, cosets=g, K=K, first_hit=first,
                         uniform_window=(maxgap + 1) if maxgap < math.inf else None,  # None: some coset of <g_p> has no S point
                          rand_model_maxgap=rand_gap,
                         P_next_given_S=nxt / max(1, nS), cusp_frac=float(cusp.mean()) if nS else 0,
                         cusp_violations=viol, mean_cluster=float(np.mean(cl_sizes)) if cl_sizes else None,
                         mean_fwd_hit=float(np.mean(fwd_mean)) if fwd_mean else None))
    if len(primes) < 2:
        return None
    rec['gens'] = gens
    # phases of the most biased character on the generators: ||k* t_p / h||
    rec['kstar_phases'] = {str(p): abs(((kstar * t[p]) % h) / h - round(((kstar * t[p]) % h) / h)) for p in primes}
    rec['min_over_p_uniform_window'] = min((x['uniform_window'] for x in gens if x['uniform_window'] is not None), default=None)
    # Restrict the Cayley-graph part to the component of the identity, the subgroup
    # <g0> with g0 = gcd(h, t_p...); every repo start K_p t_p lies in it.
    g0 = h
    for p in primes:
        g0 = gcd(g0, t[p])
    rec['component_index'] = g0
    h_full = h
    h = h // g0
    S = S[::g0].copy()
    t = {p: t[p] // g0 for p in primes}
    nS = int(S.sum()); delta = nS / h
    rec.update(h_comp=h, nS_comp=nS, delta_comp=delta)
    for gi in gens:
        gi['t_comp'] = t[gi['p']]
    # Cayley graph spectrum (exact, characters of Z/h)
    k = np.arange(h)
    lam = np.zeros(h)
    for p in primes:
        lam += 2 * np.cos(2 * np.pi * k * t[p] / h)
    d = 2 * len(primes)
    lam_star = float(np.abs(lam[1:]).max()); lam2 = float(lam[1:].max())
    gap = 1 - lam2 / d
    rec.update(d=d, lam_abs=lam_star, lam2=lam2, gap=gap,
               eml_ball_bound=(lam_star / d) ** 2 * h / delta if delta else None,
               trace_lower=(h - d) / (d * (h - 1)), lam_over_d_sq=(lam_star / d) ** 2)
    # adjacency as index shifts
    shifts = sorted({s % h for p in primes for s in (t[p], -t[p])})
    # multi-source BFS distance to S
    dist = np.full(h, -1); dist[S] = 0; frontier = np.flatnonzero(S); r = 0
    while len(frontier):
        r += 1
        nb = np.unique(((frontier[:, None] + np.array(shifts)[None, :]) % h).ravel())
        nb = nb[dist[nb] < 0]; dist[nb] = r; frontier = nb
    rec['covering_radius'] = int(dist.max())
    # BFS from each repo start until S: explored ball size before the hit layer
    balls = []
    for gi in gens:
        x0 = (gi['K'] * gi['t_comp']) % h
        seen = np.zeros(h, dtype=bool); seen[x0] = True; front = np.array([x0]); size = 1
        while not S[front].any():
            nb = np.unique(((front[:, None] + np.array(shifts)[None, :]) % h).ravel())
            nb = nb[~seen[nb]]; seen[nb] = True; front = nb; size += 0 if S[nb].any() else len(nb)
        balls.append(size)
    rec['bfs_ball_max'] = int(max(balls)); rec['bfs_ball_mean'] = float(np.mean(balls))
    # exact E_pi tau_S for the simple random walk (non-lazy, d-regular multigraph)
    if nS == 0:
        return rec
    if nS < h:
        Sc = np.flatnonzero(~S); pos = -np.ones(h, dtype=np.int64); pos[Sc] = np.arange(len(Sc))
        rows, cols, vals = [], [], []
        for p in primes:
            for s in (t[p], -t[p]):
                y = (Sc + s) % h; ok = pos[y] >= 0
                rows.append(np.flatnonzero(ok)); cols.append(pos[y[ok]]); vals.append(np.full(ok.sum(), 1.0 / d))
        Pm = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(Sc), len(Sc)))
        M = sp.identity(len(Sc), format='csr') - Pm  # symmetric positive definite (S nonempty)
        tau, info = spl.cg(M, np.ones(len(Sc)), rtol=1e-10, maxiter=200000)
        assert info == 0, info
        rec['tau_residual'] = float(np.abs(M @ tau - 1).max())
        rec['E_pi_tau'] = float(tau.sum() / h); rec['max_tau'] = float(tau.max())
        rec['AF_bound'] = (1 - delta) / (delta * gap) if gap > 0 else None
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bits', type=int, nargs='+', default=[26, 30, 34, 37])
    ap.add_argument('--per-size', type=int, default=4)
    ap.add_argument('--u', type=float, nargs='+', default=[3.0, 4.0])
    ap.add_argument('--out', default='data/cayley.jsonl')
    a = ap.parse_args()
    with open(a.out, 'w') as fh:
        for bits in a.bits:
            for (N, p, q, h, gen) in semiprimes(bits, a.per_size, 0):
                for u in a.u:
                    t0 = time.time(); rec = analyse(N, h, gen, u)
                    if rec is None:
                        continue
                    rec['bits'] = bits; rec['secs'] = time.time() - t0
                    fh.write(json.dumps(rec) + '\n'); fh.flush()
                    print(bits, N, h, 'Y', rec['Y'], 'delta %.4f' % rec['delta'], 'gap %.3f' % rec['gap'],
                          'cov', rec['covering_radius'], 'Wunif_min', rec['min_over_p_uniform_window'],
                          '%.1fs' % rec['secs'], flush=True)


if __name__ == '__main__':
    main()
