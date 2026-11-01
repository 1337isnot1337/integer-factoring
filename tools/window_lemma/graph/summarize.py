#!/usr/bin/env python3
"""Summarise data/cayley.jsonl into data/cayley_summary.md (stdlib only)."""
import json, math, statistics as st, sys
src = sys.argv[1] if len(sys.argv) > 1 else 'data/cayley.jsonl'
out = sys.argv[2] if len(sys.argv) > 2 else 'data/cayley_summary.md'
R = [json.loads(l) for l in open(src)]
L = []
w = L.append
w(f"# cayley_small summary ({len(R)} (N, Y) instances from {src})\n")
w("| bits | h_comp | Y | m | delta | maxF(S)/maxF(rand) | max_p phase(k*) | (lam/d)^2 | EML ball bound / h | BFS ball max | cov radius | E_pi tau | AF bound | min_p Wunif | 1/delta |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
viol = cusp_ok = af_ok = pig_ok = 0; tot = 0
ratios = []; pn = {}; clus = []; NOHIT = []
for r in R:
    tot += 1
    g = r['gens']
    viol += sum(x['cusp_violations'] for x in g)
    cusp_ok += all(x['P_next_given_S'] >= x['cusp_frac'] - 1e-12 for x in g)
    if r.get('AF_bound') is not None and r.get('E_pi_tau') is not None:
        af_ok += r['E_pi_tau'] <= r['AF_bound'] + 1e-9
    # pigeonhole: uniform window >= 1/delta_c for some coset; check >= 1/delta - 1 per generator
    pig_ok += all(x['uniform_window'] is None or x['uniform_window'] >= 1 / r['delta'] for x in g if r['delta'] > 0)
    nohit = sum(x['uniform_window'] is None for x in g); NOHIT.append(nohit)
    ratios.append(r['maxFourier_S'] / r['maxFourier_rand'])
    for x in g:
        key = x['p'] if x['p'] < 30 else '>=30'
        pn.setdefault(key, []).append((x['P_next_given_S'] / r['delta'], x['cusp_frac'], (x['uniform_window'] or float('inf')) / max(1e-9, x['rand_model_maxgap'])))
        if x['mean_cluster'] and x['mean_fwd_hit'] is not None:
            clus.append((x['p'], x['mean_fwd_hit'], x['mean_cluster'] / r['delta']))
    ph = max(r['kstar_phases'].values())
    eml = r['eml_ball_bound'] / r['h_comp'] if r.get('eml_ball_bound') else float('nan')
    w(f"| {r['bits']} | {r['h_comp']} | {r['Y']} | {r['m']} | {r['delta_comp']:.4f} | {ratios[-1]:.1f} | {ph:.3f} | {r['lam_over_d_sq']:.3f} | {eml:.2f} | {r['bfs_ball_max']} | {r['covering_radius']} | {r.get('E_pi_tau', float('nan')):.1f} | {r.get('AF_bound') or float('nan'):.1f} | {r['min_over_p_uniform_window']} | {1/r['delta']:.1f} |")
w("")
w(f"* cusp-lemma violations (s in S, 2pA(s) < sqrt|D|, but g_p^{{+-1}} s not in S): **{viol}** over all instances/generators")
w(f"* instances with P(g_p s in S | s in S) >= cusp fraction for every p: {cusp_ok}/{tot}")
w(f"* instances with E_pi tau_S <= (1-delta)/(delta*gap) (Aldous-Fill): {af_ok}/{tot}")
w(f"* instances with uniform window >= 1/delta for every generator (pigeonhole, PROVEN bound; None = infinite counts as satisfied): {pig_ok}/{tot}")
w(f"* generators with a <g_p>-coset containing no S point (uniform window infinite): {sum(NOHIT)} (instances: {sum(1 for z in NOHIT if z)})")
w(f"* max|S^(chi)| / max|R^(chi)| (R random, same size): median {st.median(ratios):.1f}, min {min(ratios):.1f}, max {max(ratios):.1f}")
w("\n| p | n | median P(next in S given S)/delta | median cusp fraction | median uniform window / random-model max gap |")
w("|---|---|---|---|---|")
for k in sorted(pn, key=lambda z: (isinstance(z, str), z)):
    v = pn[k]
    w(f"| {k} | {len(v)} | {st.median(a for a,_,_ in v):.2f} | {st.median(b for _,b,_ in v):.3f} | {st.median(c for _,_,c in v):.2f} |")
if clus:
    rat = [a / b for _, a, b in clus if b > 0]
    w(f"\n* mean forward hitting time along g_p orbits (uniform start) / (mean cluster size / delta): median {st.median(rat):.2f}, IQR {st.quantiles(rat, n=4)[0]:.2f}-{st.quantiles(rat, n=4)[2]:.2f} (n={len(rat)})")
open(out, 'w').write("\n".join(L) + "\n")
print("\n".join(L))
