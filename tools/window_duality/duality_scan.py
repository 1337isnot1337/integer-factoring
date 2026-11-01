#!/usr/bin/env python3
"""Numerical test of the Kronecker / Erdos-Turan duality for WL(Y, W).

For each (Delta, Y, p): exact first-hit window W of x_j = g_p^{K_p+j} into S, all
characters chi of Cl(Delta) (exact, via clgroup.ClassGroup), and
  * F(chi)   = sum_q |1-chi(g_q)|^2 / (ln q)^2    (ball-scale weights; 'F_ball')
  * Fp(chi)  = sum_q (1-Re chi(g_q)) / q          (pretentious distance^2 to 1)
  * s(chi)   = S^(chi)/|S|                         (normalized Fourier coefficient)
  * theta_p  = arg chi(g_p)/2pi.
Decomposition at W0 = ceil(1/delta) (exact identity):
  count_{W0}/(delta W0) - 1 = sum_{chi != 1} s(chi) chi(x_0) D_{W0}(theta_p)/W0 =: E_all,
and the parts E_L from character sets L (low-F, p-major arc, both).
Controls: a random S' (uniform, same size, inside H = <g_q>) analysed identically.
Also checks the PROVEN statements of duality.md (Theorems 1-3) exactly.
Usage: python3 duality_scan.py --count 120 --jobs 8 --out data/scan.jsonl
"""
import argparse, json, math, os, sys, time, random
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clgroup import ClassGroup, repo_base, smooth_mask, least_power_above, dist1, pari, np  # noqa
from math import isqrt, ceil

LOW_FRAC = 0.01      # low-frequency set L_F = lowest 1% of F among chi nontrivial on H


def subgroup_mask(G, tvecs):
    m = np.zeros(G.shape, dtype=bool)
    m[(0,) * len(G.shape)] = True
    while True:
        new = m.copy()
        for t in tvecs:
            new |= np.roll(m, shift=tuple(int(x) for x in t), axis=tuple(range(len(G.shape))))
        if (new == m).all():
            return m.ravel()
        m = new


def ball_vectors(primes, radius, ramified=()):
    """All v in Z^m with sum |v_q| ln q < radius (weighted l1 ball, strict).
    A ramified base prime (the repo's added 2 when Delta=-4N) only takes v_q in {0,1}."""
    logs = [math.log(q) for q in primes]
    ram = [q in ramified for q in primes]
    out = []
    def rec(i, rem, v):
        if i == len(primes):
            out.append(tuple(v)); return
        e = 0
        while e * logs[i] < rem - 1e-12 or e == 0:
            if e == 0:
                v.append(0); rec(i + 1, rem, v); v.pop()
            else:
                for s in ((e,) if ram[i] else (e, -e)):
                    if ram[i] and e > 1:
                        continue
                    v.append(s); rec(i + 1, rem - e * logs[i], v); v.pop()
            e += 1
            if e * logs[i] >= rem - 1e-12:
                break
    rec(0, radius, [])
    return np.array(out, dtype=np.int64)


def orbit(G, t, K, n):
    cyc = np.array(G.cyc)
    j = np.arange(n)[:, None]
    X = ((K + j) * t[None, :]) % cyc
    return np.ravel_multi_index(tuple(X.T), G.shape)


def order_of(G, t):
    o = 1
    for ti, c in zip(t, G.cyc):
        o = math.lcm(o, c // math.gcd(int(ti), c))
    return o


def decomposition(G, Sflat, sh, x0chi, th, W0, sets, index_H):
    """E parts: (1/[Cl:H]) sum_{chi in set} s(chi) chi(x0) D_W0(theta)/W0 (real part)."""
    z = np.exp(2j * np.pi * th)
    with np.errstate(divide='ignore', invalid='ignore'):
        D = np.where(np.abs(1 - z) < 1e-12, W0, (z ** W0 - 1) / (z - 1))
    terms = sh * x0chi * D / (W0 * index_H)
    return {k: float(np.real(terms[msk].sum())) for k, msk in sets.items()}


def analyse(job):
    N, Ys, tag = job
    t0 = time.time()
    G = ClassGroup(N)
    D, h = G.D, G.h
    res = []
    Kgrid = G.char_grid()
    cyc = np.array(G.cyc)
    for Y in Ys:
        base = repo_base(N, D, Y)
        primes = list(base)
        if len(primes) < 2:
            continue
        tv = {q: G.dlog(base[q]) for q in primes}
        ram = tuple(q for q in primes if D % q == 0)
        S, Hh = smooth_mask(G, primes)
        Hmask = subgroup_mask(G, list(tv.values()))
        hH = int(Hmask.sum())
        assert not (S & ~Hmask).any(), 'S not inside H'   # PROVEN: S subset H
        nS = int(S.sum()); delta_H = nS / hH
        theta = {q: G.phases(tv[q]) for q in primes}
        F = sum(np.abs(1 - np.exp(2j * np.pi * theta[q])) ** 2 / math.log(q) ** 2 for q in primes)
        Fp = sum((1 - np.cos(2 * np.pi * theta[q])) / q for q in primes)
        trivH = F < 1e-18          # characters trivial on H (Path B / Bach stratum)
        n_trivH = int(trivH.sum())
        assert n_trivH * hH == h
        Shat = np.fft.fftn(S.reshape(G.shape)).ravel()
        sh = Shat / nS
        rng = np.random.default_rng(N % (2 ** 32) + Y)
        # canonical order (by reduced form) so the control does not depend on PARI's random generators
        Hidx = np.array([f for k, f in sorted(G.index.items()) if Hmask[f]])
        Sr = np.zeros(h, dtype=bool); Sr[rng.choice(Hidx, nS, replace=False)] = True
        shr = np.fft.fftn(Sr.reshape(G.shape)).ravel() / nS
        nontriv = ~trivH
        # real (order <= 2) characters: genus characters, the Path B / least-nonresidue stratum
        real = nontriv & np.all((2 * Kgrid) % cyc[None, :] == 0, axis=1)
        Fp_real_min = float(Fp[real].min()) if real.any() else None
        Fp_min = float(Fp[nontriv].min()) if nontriv.any() else None
        Fp_med = float(np.median(Fp[nontriv])) if nontriv.any() else None
        cut = np.quantile(F[nontriv], LOW_FRAC) if nontriv.any() else 0
        lowF = nontriv & (F <= cut)
        # ---- Theorem checks: ball minorant (a = b = R/2)
        R = math.log(math.sqrt(-D) / 2)
        a = R / 2
        Ba = ball_vectors(primes, a, ram)
        Tm = np.array([tv[q] for q in primes])          # m x r
        img = np.ravel_multi_index(tuple(((Ba @ Tm) % cyc).T), G.shape)
        assert len(set(img.tolist())) == len(Ba)        # Lemma 1(ii) injectivity
        B2 = ball_vectors(primes, 2 * a, ram)
        img2 = np.ravel_multi_index(tuple(((B2 @ Tm) % cyc).T), G.shape)
        assert S[img2].all()                             # Lemma 1(ii): phi(B_R) in S
        ind = np.zeros(h); ind[img] = 1
        ca = np.abs(np.fft.fftn(ind.reshape(G.shape)).ravel()) / len(Ba)
        # boundary fractions eps_q(a)
        eps = {}
        for q in primes:
            lq = math.log(q)
            nb = len(ball_vectors(primes, a - lq, ram)) if a > lq else 0
            eps[q] = 1 - nb / len(Ba)
        # Theorem 2 check: c_a(chi)|1-chi(g_q)| <= 2 eps_q(a)
        thm2 = max(float((ca * np.abs(1 - np.exp(2j * np.pi * theta[q]))).max() / (2 * eps[q])) if eps[q] > 0 else 0.0 for q in primes)
        # Theorem 3 (fourth moment): sum_{chi|H != 1} c_a^4 >= [Cl:H](hH/|B_2a| - 1)  [over H-characters]
        sum4 = float((ca[nontriv] ** 4).sum()) / n_trivH
        rhs4 = hH / len(B2) - 1   # |phi(B_a)+phi(B_a)| <= |B_2a|
        maxc = float(ca[nontriv].max()) if nontriv.any() else 0
        rec_y = dict(N=N, D=D, tag=tag, h=h, cyc=G.cyc, certified=G.certified, Y=Y, m=len(primes),
                     H=Hh, u=math.log(Hh) / math.log(max(primes)), hH=hH, index_H=n_trivH, nS=nS,
                     delta=delta_H, R=R, Ba=len(Ba), B2a=len(B2), thm2_max_ratio=thm2,
                     thm3_sum4=sum4, thm3_rhs=rhs4, maxc_a=maxc,
                     maxFourier_S=float(np.abs(sh[nontriv]).max()) if nontriv.any() else 0,
                     maxFourier_rand=float(np.abs(shr[nontriv]).max()) if nontriv.any() else 0,
                     F_cut=float(cut), Fp_real_min=Fp_real_min, Fp_min=Fp_min, Fp_med=Fp_med,
                     n_real_nontrivH=int(real.sum()), n_lowF=int(lowF.sum()),
                     lowF_meanabs_S=float(np.abs(sh[lowF]).mean()) if lowF.any() else 0,
                     lowF_meanabs_rand=float(np.abs(shr[lowF]).mean()) if lowF.any() else 0,
                     hiF_meanabs_S=float(np.abs(sh[nontriv & ~lowF]).mean()) if (nontriv & ~lowF).any() else 0,
                     targets=[])
        W0 = max(1, ceil(1 / delta_H))
        for p in primes:
            if p in ram:
                continue
            t = tv[p]; K = least_power_above(p, Hh); n = order_of(G, t)
            orb = orbit(G, t, K, n)
            x0 = orb[0]
            x0vec = np.array(np.unravel_index(x0, G.shape))
            x0chi = np.exp(2j * np.pi * np.mod((Kgrid * x0vec[None, :] % cyc) / cyc, 1).sum(1))
            hits = np.flatnonzero(S[orb]); W = int(hits[0]) + 1 if len(hits) else None
            hits_r = np.flatnonzero(Sr[orb]); Wr = int(hits_r[0]) + 1 if len(hits_r) else None
            # ball hits: first j with x_j in phi(B_R) (minorant support)
            inball = np.zeros(h, dtype=bool); inball[img2] = True
            hb = np.flatnonzero(inball[orb]); Wball = int(hb[0]) + 1 if len(hb) else None
            th = theta[p]; dp = dist1(th)
            # Theorem 1 check on the ball minorant f = 1_{phi B_a} * 1_{phi B_a} / |B_a|
            Wm = (Wball - 1) if Wball else n   # number of initial misses of supp f
            if Wm > 0:
                with np.errstate(divide='ignore'):
                    mn = np.minimum(Wm, np.where(dp > 0, 1 / (2 * dp), np.inf))
                lhs = float((ca[nontriv] ** 2 * mn[nontriv]).sum())
                thm1 = lhs / (n_trivH * Wm)
            else:
                thm1 = None
            # p-major arc of width 1/W0
            maj = nontriv & (dp < 1 / W0)
            sets = {'all': nontriv, 'lowF': lowF, 'major': maj, 'lowF_major': lowF & maj,
                    'major_not_lowF': maj & ~lowF}
            E = decomposition(G, S, sh, x0chi, th, W0, sets, n_trivH)
            Er = decomposition(G, Sr, shr, x0chi, th, W0, sets, n_trivH)
            cnt = int(S[orb[:W0]].sum()) if W0 <= n else int(S[np.resize(orb, W0)].sum())
            cntr = int(Sr[orb[:W0]].sum()) if W0 <= n else int(Sr[np.resize(orb, W0)].sum())
            # exactness of the identity (counts over the H-coset; trivial-on-H chars give the main term)
            Eex = cnt / (delta_H * W0) - 1
            # S-independent geometry statistic: min F over chi nontrivial on H with ||theta_p|| <= delta
            arc = nontriv & (dp <= delta_H)
            Gp = float(F[arc].min()) if arc.any() else None
            Gpp = float(Fp[arc].min()) if arc.any() else None
            # coset (theta_p = 0) term: exact count over the full cycle
            zero = nontriv & (dp < 1e-12)
            rec_y['targets'].append(dict(p=p, K=K, ord=n, W=W, W_rand=Wr, W_ball=Wball, thm1_ratio=thm1,
                                         E=E, E_rand=Er, E_exact=Eex, count_W0=cnt, count_W0_rand=cntr, W0=W0,
                                         count_cycle=int(S[orb].sum()), G_p=Gp, G_p_pret=Gpp,
                                         n_p_trivial=int(zero.sum())))
        res.append(rec_y)
    return res, time.time() - t0


def random_N(rng, lo, hi):
    while True:
        N = rng.randrange(lo, hi) | 1
        if pari.issquarefree(N):
            return N


def adversary_N(p, bits, Y, rng):
    """Small instances of the repo's P3 family (long forced misses), squarefree N."""
    K = max(3, round((bits - 2) / (2 * math.log2(p))))
    W = K - 2
    while W >= 1 and (3 * p ** (K - W)) // 4 <= Y:
        W -= 1
    lo = max(Y + 1, -(-3 * p ** (K - W - 2) // 4)); hi = (3 * p ** (K - W)) // 4 - 1
    A = lo | 1
    while not (pari.isprime(A) and A % p):
        A += 2
    if A > hi:
        return None
    M = p ** (K + W)
    for x in range(1, A, 2):
        if x % p == 0:
            continue
        n = 4 * M * A - x * x
        if pari.issquarefree(n) and n < 4 * 10 ** 9:
            return n
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--count', type=int, default=120)
    ap.add_argument('--u', type=float, nargs='+', default=[2.0, 2.5, 3.0, 4.0])
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--lo', type=float, default=4.0)
    ap.add_argument('--hi', type=float, default=9.0)
    ap.add_argument('--no-adversary', action='store_true')
    ap.add_argument('--out', default='data/scan.jsonl')
    a = ap.parse_args()
    rng = random.Random(a.seed)
    jobs = []
    for i in range(a.count):
        e = a.lo + (a.hi - a.lo) * i / max(1, a.count - 1)   # log10|Delta| from lo to hi
        N = random_N(rng, int(10 ** e), int(10 ** e * 1.25))
        Hh = isqrt(N // 3)
        Ys = sorted({max(7, int(round(Hh ** (1 / u)))) for u in a.u})
        jobs.append((N, Ys, 'random'))
    seen = set()
    for p in (() if a.no_adversary else (3, 5, 7)):
        for bits in (24, 26, 28, 30, 32):
            for Y in (13, 31):
                N = adversary_N(p, bits, Y, rng)
                if N and N not in seen:
                    seen.add(N)
                    jobs.append((N, [Y], f'adversary_p{p}'))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with Pool(a.jobs) as pool, open(a.out, 'w') as fh:
        for res, secs in pool.imap_unordered(analyse, jobs):
            for r in res:
                r['secs'] = secs
                fh.write(json.dumps(r) + '\n')
            if res:
                print(res[0]['N'], res[0]['h'], [r['Y'] for r in res], '%.1fs' % secs, flush=True)


if __name__ == '__main__':
    main()
