#!/usr/bin/env python3
"""Summarize data/scan.jsonl from duality_scan.py.  Usage: python3 analyze_scan.py data/scan.jsonl"""
import json, math, sys
import numpy as np


def avgrank(v):
    v = np.asarray(v, float); o = np.argsort(v, kind='stable'); r = np.empty(len(v)); r[o] = np.arange(len(v))
    _, inv = np.unique(v, return_inverse=True)
    return (np.bincount(inv, r) / np.bincount(inv))[inv]   # ties get their average rank


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rx = avgrank(x); ry = avgrank(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def pearson(x, y):
    return float(np.corrcoef(np.asarray(x, float), np.asarray(y, float))[0, 1])


def partial_spearman(x, y, covs):
    r = avgrank
    Z = np.column_stack([np.ones(len(x))] + [r(c) for c in covs])
    def res(v):
        v = r(v); beta, *_ = np.linalg.lstsq(Z, v, rcond=None); return v - Z @ beta
    return float(np.corrcoef(res(x), res(y))[0, 1])


def boot_ci(x, y, f, n=400, seed=0):
    rng = np.random.default_rng(seed); x = np.asarray(x); y = np.asarray(y); out = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x)); out.append(f(x[i], y[i]))
    return np.percentile(out, [2.5, 97.5])


def main(path):
    allrecs = sorted((json.loads(l) for l in open(path)), key=lambda r: (abs(r['D']), r['Y']))  # order-independent
    recs = [r for r in allrecs if r['delta'] < 0.5 and r['maxFourier_rand'] > 0]
    print('dropped %d degenerate instances with delta >= 0.5 (Y too large for H)' % (len(allrecs) - len(recs)))
    T = [(r, t) for r in recs for t in r['targets']]
    print(f'instances (Delta,Y): {len(recs)}   targets (Delta,Y,p): {len(T)}')
    print('distinct Delta:', len({r["D"] for r in recs}), '  |Delta| range: %.2e .. %.2e' % (min(-r['D'] for r in recs), max(-r['D'] for r in recs)))
    print('h range:', min(r['h'] for r in recs), '..', max(r['h'] for r in recs),
          '  non-cyclic Cl:', sum(len(r['cyc']) > 1 for r in recs), '  [Cl:H]>1:', sum(r['index_H'] > 1 for r in recs),
          '  certified by reduced-form count:', sum(bool(r['certified']) for r in recs))
    print('m range:', min(r['m'] for r in recs), '..', max(r['m'] for r in recs), '  delta range: %.4f .. %.3f' % (min(r['delta'] for r in recs), max(r['delta'] for r in recs)))
    print('tags:', {k: sum(r['tag'] == k for r in recs) for k in sorted({r['tag'] for r in recs})})
    print()
    print('== PROVEN statements, exact checks ==')
    print('Lemma 1(ii) injectivity on B_{R/2} and phi(B_R) in S: asserted in every instance (script aborts otherwise)')
    print('S inside H = <g_q>: asserted in every instance')
    t1 = [t['thm1_ratio'] for r, t in T if t['thm1_ratio'] is not None]
    i1 = min(range(len(t1)), key=lambda i: t1[i])
    rr, tt = [(r, t) for r, t in T if t['thm1_ratio'] is not None][i1]
    print('  (Thm 1 tightest at D=%d h=%d Y=%d p=%d W_ball=%s)' % (rr['D'], rr['h'], rr['Y'], tt['p'], tt['W_ball']))
    print('Thm 1 (failure => sum c_a^2 min(W,1/2||th||) >= [Cl:H] W): min ratio LHS/RHS = %.3f over %d targets (must be >= 1); median %.1f' % (min(t1), len(t1), np.median(t1)))
    t2 = [r['thm2_max_ratio'] for r in recs]
    print('Thm 2 (c_a(chi)|1-chi(g_q)| <= 2 eps_q(a)): max ratio = %.3f over all chi, q, instances (must be <= 1); median of per-instance max %.3f' % (max(t2), np.median(t2)))
    t3 = [(r['thm3_sum4'], r['thm3_rhs']) for r in recs]
    print('Thm 3 (sum_{chi|H!=1} c_a^4 >= h_H/|B_2a| - 1): violations = %d / %d; median LHS/RHS = %.2f' % (sum(a < b - 1e-9 for a, b in t3), len(t3), np.median([a / b for a, b in t3 if b > 0])))
    w0 = [(4 * r['hH'] / r['Ba'], max(t['W'] for t in r['targets'] if t['W'])) for r in recs if r['targets'] and any(t['W'] for t in r['targets'])]
    print('Cor 1 minor-arc threshold W_* = h_H/(eta |B_a|) at eta=1/4: median W_*/max_p W = %.1f (min %.2f); W_* exceeds every observed window in %d/%d instances; median W_* / sqrt(h_H) = %.2f' % (
        np.median([a / b for a, b in w0]), min(a / b for a, b in w0), sum(a > b for a, b in w0), len(w0), np.median([4 * r['hH'] / r['Ba'] / math.sqrt(r['hH']) for r in recs])))
    mc = [r['maxc_a'] for r in recs]
    print('max_{chi|H!=1} c_a(chi): median %.3f, min %.3f  (a low-frequency character with c_a >= 1/2 exists in %d/%d instances)' % (np.median(mc), min(mc), sum(c >= 0.5 for c in mc), len(mc)))
    ex = max(abs(t['E']['all'] - t['E_exact']) for r, t in T)
    print('Fourier identity exactness max |E_all - E_exact| = %.2e' % ex)
    print()
    print('== Fourier bias of S on low-frequency characters (lowest 1% of F among chi|H != 1) ==')
    ratio = [r['lowF_meanabs_S'] / r['lowF_meanabs_rand'] for r in recs if r['lowF_meanabs_rand'] > 0]
    ratio2 = [r['lowF_meanabs_S'] / r['hiF_meanabs_S'] for r in recs if r['hiF_meanabs_S'] > 0]
    print('mean|s(chi)| on low-F chars, S vs random S\': median ratio %.2f (IQR %.2f-%.2f)' % (np.median(ratio), *np.percentile(ratio, [25, 75])))
    print('mean|s(chi)| on low-F vs other chars (S): median ratio %.2f' % np.median(ratio2))
    print('max|s| S vs random: median ratio %.2f' % np.median([r['maxFourier_S'] / r['maxFourier_rand'] for r in recs]))
    print()
    print('== Windows ==')
    fin = [(r, t) for r, t in T if t['W'] is not None]
    cens = [(r, t) for r, t in T if t['W'] is None]
    print('targets with orbit missing S entirely (coset miss): %d; random-control coset misses: %d' % (len(cens), sum(t['W_rand'] is None for r, t in T)))
    z = np.array([t['W'] * r['delta'] for r, t in fin])
    zr = np.array([t['W_rand'] * r['delta'] for r, t in T if t['W_rand'] is not None])
    print('W*delta (real S): median %.2f, 90%% %.2f, max %.1f;  random S\': median %.2f, 90%% %.2f, max %.1f' % (np.median(z), np.percentile(z, 90), z.max(), np.median(zr), np.percentile(zr, 90), zr.max()))
    print('W == W_ball (first hit already in phi(B_R)): %d / %d' % (sum(t['W'] == t['W_ball'] for r, t in fin), len(fin)))
    print()
    print('== Conjecture test A: long windows vs existence of a low-frequency chi in the p-major arc ==')
    print('G_p = min F(chi) over chi|H != 1 with ||theta_p(chi)|| <= delta (S-independent); "exists" := G_p <= F_cut (1% quantile)')
    rows = [(r, t) for r, t in fin if t['G_p'] is not None and t['W_rand'] is not None]
    g = np.array([t['G_p'] / r['F_cut'] for r, t in rows])
    lz = np.log([t['W'] * r['delta'] for r, t in rows]); lzr = np.log([t['W_rand'] * r['delta'] for r, t in rows])
    print('n=%d  Spearman(log W delta, G_p/F_cut): real S %.3f, random S\' %.3f' % (len(rows), spearman(lz, g), spearman(lzr, g)))
    gp = np.array([t['G_p_pret'] for r, t in rows])
    print('      Spearman(log W delta, pretentious G_p): real S %.3f, random S\' %.3f' % (spearman(lz, gp), spearman(lzr, gp)))
    ex_ = g <= 1
    long_ = lz >= np.percentile(lz, 90)
    longr = lzr >= np.percentile(lzr, 90)
    print('P(exists | long, top10%%) = %.3f   P(exists | not long) = %.3f   (real S)' % (ex_[long_].mean(), ex_[~long_].mean()))
    print('P(exists | long, top10%%) = %.3f   P(exists | not long) = %.3f   (random S\')' % (ex_[longr].mean(), ex_[~longr].mean()))
    print('counterexamples to "long => exists" (real S top-10%% windows with no 1%%-low chi in the p-arc): %d / %d' % ((long_ & ~ex_).sum(), long_.sum()))
    # control for confounding by delta, m, p, h
    covs = [np.log([r['delta'] for r, t in rows]), [r['m'] for r, t in rows], np.log([t['p'] for r, t in rows]), np.log([r['h'] for r, t in rows])]
    for name, gg in (('G_p/F_cut (ball weights)', g), ('pretentious G_p', gp)):
        pr = partial_spearman(lz, gg, covs); prr = partial_spearman(lzr, gg, covs)
        ci = boot_ci(lz, gg, spearman)
        print('partial Spearman(log W delta, %s | log delta, m, log p, log h): real %.3f, random %.3f   [raw Spearman 95%% bootstrap CI real %.3f..%.3f]' % (name, pr, prr, *ci))
    print('partial check: Spearman(G_p/F_cut, delta) = %.3f, Spearman(G_p/F_cut, m) = %.3f' % (spearman(g, [r['delta'] for r, t in rows]), spearman(g, [r['m'] for r, t in rows])))
    print()
    print('== Conjecture test B: which characters carry the deficit at W0 = ceil(1/delta) ==')
    print('E_all = count_{W0}/(delta W0) - 1 exactly; E_L = contribution of character set L')
    for key in ('lowF', 'major', 'lowF_major', 'major_not_lowF'):
        a = [t['E']['all'] for r, t in T]; b = [t['E'][key] for r, t in T]
        ar = [t['E_rand']['all'] for r, t in T]; br = [t['E_rand'][key] for r, t in T]
        print('  %-15s corr(E_all,E_L): real %.3f  random %.3f | mean|E_L|: real %.3f random %.3f | R^2 of E_all ~ E_L real %.3f random %.3f' % (
            key, pearson(a, b), pearson(ar, br), np.mean(np.abs(b)), np.mean(np.abs(br)),
            1 - np.mean((np.array(a) - np.array(b)) ** 2) / np.var(a), 1 - np.mean((np.array(ar) - np.array(br)) ** 2) / np.var(ar)))
    misses = [(r, t) for r, t in T if t['count_W0'] == 0]
    print('targets with zero hits in the first W0 steps: %d / %d (random S\': %d)' % (len(misses), len(T), sum(t['count_W0_rand'] == 0 for r, t in T)))
    if misses:
        fl = [t['E']['lowF_major'] for r, t in misses]
        fm = [t['E']['major'] for r, t in misses]
        print('  in those: E_major median %.3f (deficit carried by p-major arc if ~ -1), E_lowF_major median %.3f, frac with E_lowF_major < 0: %.3f' % (np.median(fm), np.median(fl), np.mean(np.array(fl) < 0)))
    print()
    print('== Path B stratum: real (genus) characters ==')
    rr = [r for r in recs if r['Fp_real_min'] is not None and r['targets'] and all(t['W'] is not None and t['W_rand'] is not None for t in r['targets'])]
    gr = [r['Fp_real_min'] / r['Fp_med'] for r in rr]
    mx = [max(math.log(t['W'] * r['delta']) for t in r['targets']) for r in rr]
    mxr = [max(math.log(t['W_rand'] * r['delta']) for t in r['targets']) for r in rr]
    gall = [r['Fp_min'] / r['Fp_med'] for r in rr]
    print('instances with a real character nontrivial on H: %d; real chars trivial on all of B_Y but not on Cl (exact Path B stratum): %d instances with [Cl:H]>1' % (len(rr), sum(r['index_H'] > 1 for r in recs)))
    print('Spearman(max_p log W delta, min pretentious dist of REAL chars / median): real S %.3f, random S\' %.3f' % (spearman(mx, gr), spearman(mxr, gr)))
    print('Spearman(max_p log W delta, min pretentious dist of ALL chars / median): real S %.3f, random S\' %.3f' % (spearman(mx, gall), spearman(mxr, gall)))
    cv = [np.log([r['delta'] for r in rr]), [r['m'] for r in rr], np.log([r['h'] for r in rr])]
    print('  partial (| log delta, m, log h): REAL chars real S %.3f random %.3f; ALL chars real S %.3f random %.3f' % (
        partial_spearman(mx, gr, cv), partial_spearman(mxr, gr, cv), partial_spearman(mx, gall, cv), partial_spearman(mxr, gall, cv)))
    print('fraction of instances whose most pretentious character is real (order 2): %.3f' % np.mean([abs(r['Fp_real_min'] - r['Fp_min']) < 1e-12 for r in rr]))
    idx = [r['index_H'] > 1 for r in recs if r['targets'] and all(t['W'] is not None for t in r['targets'])]
    zz = [np.median([t['W'] * r['delta'] for t in r['targets']]) for r in recs if r['targets'] and all(t['W'] is not None for t in r['targets'])]
    print('median W delta when [Cl:H]>1: %.2f (n=%d); when H=Cl: %.2f (n=%d)' % (np.median([z for z, i in zip(zz, idx) if i]), sum(idx), np.median([z for z, i in zip(zz, idx) if not i]), len(idx) - sum(idx)))
    ok = [r for r in recs if r['targets'] and all(t['W'] is not None for t in r['targets'])]
    for mm in sorted({r['m'] for r in ok})[:6]:
        a1 = [np.median([t['W'] * r['delta'] for t in r['targets']]) for r in ok if r['m'] == mm and r['index_H'] > 1]
        a0 = [np.median([t['W'] * r['delta'] for t in r['targets']]) for r in ok if r['m'] == mm and r['index_H'] == 1]
        if a1 and a0:
            print('  m=%d: median W delta [Cl:H]>1 %.2f (n=%d) vs H=Cl %.2f (n=%d)' % (mm, np.median(a1), len(a1), np.median(a0), len(a0)))
    print()
    print('== By size ==')
    for lo, hi in ((1e4, 1e6), (1e6, 1e8), (1e8, 1e10), (1e10, 1e12)):
        sub = [(r, t) for r, t in rows if lo <= -r['D'] < hi]
        if len(sub) < 10:
            continue
        gs = [t['G_p'] / r['F_cut'] for r, t in sub]
        print('  |Delta| in [%.0e,%.0e): n=%d  Spearman(log W delta, G_p/F_cut) real %.3f random %.3f; median W delta %.2f' % (
            lo, hi, len(sub), spearman(np.log([t['W'] * r['delta'] for r, t in sub]), gs),
            spearman(np.log([t['W_rand'] * r['delta'] for r, t in sub]), gs), np.median([t['W'] * r['delta'] for r, t in sub])))
    print()
    print('== Ten longest normalized windows (real S) ==')
    srt = sorted(fin, key=lambda rt: -rt[1]['W'] * rt[0]['delta'])[:10]
    for r, t in srt:
        print('  D=%d h=%d Y=%d m=%d p=%d W=%d Wdelta=%.1f W_rand=%s G_p/F_cut=%s E_major=%.2f E_lowF_major=%.2f %s' % (
            r['D'], r['h'], r['Y'], r['m'], t['p'], t['W'], t['W'] * r['delta'], t['W_rand'],
            '%.2f' % (t['G_p'] / r['F_cut']) if t['G_p'] is not None else None, t['E']['major'], t['E']['lowF_major'], r['tag']))
    if cens:
        print('== Coset misses ==')
        for r, t in cens[:10]:
            print('  D=%d h=%d cyc=%s Y=%d p=%d ord=%d count_cycle=%d n_p_trivial_chars=%d' % (r['D'], r['h'], r['cyc'], r['Y'], t['p'], t['ord'], t['count_cycle'], t['n_p_trivial']))


if __name__ == '__main__':
    main(sys.argv[1])
