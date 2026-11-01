#!/usr/bin/env python3
"""Summarise data/scan_*.jsonl: per-N density MLE, geometric-model checks, growth fits.
Usage: python3 analyze.py  (writes data/summary.json and prints markdown tables)"""
import glob, json, math, statistics as st
from collections import defaultdict
import numpy as np

def rho(u):            # Dickman rho by integrating u rho'(u) = -rho(u-1)
    if u <= 1: return 1.0
    h = 1e-3; n = int(u / h) + 1
    r = np.ones(n + 1)
    k1 = int(1 / h)
    for i in range(k1 + 1, n + 1):
        x = i * h
        r[i] = r[i - 1] - h * r[i - k1] / x if False else r[i - 1] - h * (r[i - k1 - 1] + r[i - k1]) / (2 * (x - h / 2))
    return float(r[min(n, int(round(u / h)))])

def load(path):
    return [json.loads(l) for l in open(path)]

def per_n(rec):
    hits = sum(len(t['hits']) for t in rec['targets'])
    steps = sum(t['scanned'] for t in rec['targets'])
    dens = hits / steps                              # geometric MLE with censoring
    firsts = [t['hits'][0] for t in rec['targets'] if t['hits']]
    cens_first = sum(1 for t in rec['targets'] if not t['hits'])
    gaps = [b - a for t in rec['targets'] for a, b in zip(t['hits'], t['hits'][1:])]
    w = [f + 1 for f in firsts]                      # repo's w_p = queries used
    return dict(bits=rec['bits'], dens=dens, m=rec['base_size'], firsts=firsts, gaps=gaps,
                wmax=max(w) if w else None, wmean=st.mean(w) if w else None,
                cens_first=cens_first, H_bits=rec['H_bits'], Y=rec['cutoff'],
                gcd_hits=sum(t['gcd_hit'] is not None for t in rec['targets']))

def fit(xs, ys):
    A = np.vstack([xs, np.ones(len(xs))]).T
    coef, res, *_ = np.linalg.lstsq(A, ys, rcond=None)
    pred = A @ coef
    ss = float(np.sum((ys - pred) ** 2)); tot = float(np.sum((ys - ys.mean()) ** 2))
    return dict(slope=float(coef[0]), intercept=float(coef[1]), r2=1 - ss / tot,
                rmse=math.sqrt(ss / len(xs)))

out = {}
for path in sorted(glob.glob('data/scan_*.jsonl')):
    name = path.split('scan_')[1][:-6]
    rows = [per_n(r) for r in load(path)]
    by = defaultdict(list)
    for r in rows: by[r['bits']].append(r)
    table = []
    for b in sorted(by):
        rs = by[b]
        inv = [1 / r['dens'] for r in rs]
        # normalised first hits j*dens should be ~Exp(1) under the random model
        z = [f * r['dens'] for r in rs for f in r['firsts']]
        g = [x * r['dens'] for r in rs for x in r['gaps']]
        u = (rs[0]['H_bits'] * math.log(2)) / math.log(rs[0]['Y'])
        table.append(dict(bits=b, nN=len(rs), Y=rs[0]['Y'], m_med=st.median(r['m'] for r in rs),
                          u=u, log2_inv_dens_med=math.log2(st.median(inv)),
                          log2_inv_dens_min=math.log2(min(inv)), log2_inv_dens_max=math.log2(max(inv)),
                          log2_inv_rho=-math.log2(rho(u)),
                          wmax_over_N=max(r['wmax'] for r in rs if r['wmax']),
                          wmax_over_inv_dens_med=st.median(r['wmax'] * r['dens'] for r in rs if r['wmax']),
                          z_mean=st.mean(z), z_cv=st.pstdev(z) / st.mean(z), gap_mean=st.mean(g) if g else None,
                          gap_cv=(st.pstdev(g) / st.mean(g)) if g else None,
                          censored_first=sum(r['cens_first'] for r in rs),
                          targets=sum(len(r['firsts']) + r['cens_first'] for r in rs),
                          gcd_hits=sum(r['gcd_hits'] for r in rs)))
    # growth fits on per-N log2(1/dens)
    B = np.array([r['bits'] for r in rows], float)
    L = np.array([math.log2(1 / r['dens']) for r in rows])
    U = np.array([(r['H_bits'] * math.log(2)) / math.log(r['Y']) for r in rows])
    fits = {'log2(1/dens) ~ b': fit(B, L),
            'log2(1/dens) ~ b/log2(b)': fit(B / np.log2(B), L),
            'log2(1/dens) ~ u*log2(u)': fit(U * np.log2(U), L),
            'log2(1/dens) ~ -log2 rho(u)': fit(np.array([-math.log2(rho(x)) for x in U]), L),
            'log2(1/dens) ~ sqrt(b)': fit(np.sqrt(B), L)}
    out[name] = dict(table=table, fits=fits)
    print(f'\n### {name}\n')
    print('| b | #N | Y | m (med) | u | log2(1/δ̂) med [min,max] | -log2 ρ(u) | max w over N | median w_max·δ̂ | mean(jδ̂) | CV(jδ̂) | mean(gap·δ̂) | CV gap | censored/targets |')
    print('|' + '---|' * 14)
    for t in table:
        print(f"| {t['bits']} | {t['nN']} | {t['Y']} | {t['m_med']:.0f} | {t['u']:.2f} | {t['log2_inv_dens_med']:.1f} [{t['log2_inv_dens_min']:.1f}, {t['log2_inv_dens_max']:.1f}] | {t['log2_inv_rho']:.1f} | {t['wmax_over_N']} | {t['wmax_over_inv_dens_med']:.2f} | {t['z_mean']:.2f} | {t['z_cv']:.2f} | {t['gap_mean'] if t['gap_mean'] is None else round(t['gap_mean'],2)} | {t['gap_cv'] if t['gap_cv'] is None else round(t['gap_cv'],2)} | {t['censored_first']}/{t['targets']} |")
    print()
    for k, v in fits.items():
        print(f"- {k}: slope {v['slope']:.4f}, intercept {v['intercept']:.2f}, R² {v['r2']:.4f}, RMSE {v['rmse']:.3f}")
json.dump(out, open('data/summary.json', 'w'), indent=1)
