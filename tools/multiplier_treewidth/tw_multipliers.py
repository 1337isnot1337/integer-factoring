#!/usr/bin/env python3
"""Treewidth brackets for Tseitin primal graphs of n x n bit multipliers.

Circuits (gate level: AND, XOR, OR, 2 inputs each):
  array   : row-by-row ripple-carry accumulation of partial products
  wallace : Wallace-tree column compression (FA/HA layers) + final ripple adder
Primal graph: one vertex per input bit / gate output; each gate's scope
{out,in1,in2} is a clique. Product bits z_k are extra vertices joined to the
gate that drives them (constraint z_k = out_k), so G - Z is the circuit graph.

Upper bounds: greedy elimination orderings (min-degree, min-fill).  Every
elimination ordering's max clique-1 is a valid treewidth upper bound.
Lower bound: MMD+ (min-d) contraction heuristic (Bodlaender-Koster 2011),
a valid treewidth lower bound.  Pure Python, no dependencies.
"""
import sys, time, heapq, itertools

class Circ:
    def __init__(s): s.adj = {}; s.nv = 0
    def var(s):
        v = s.nv; s.nv += 1; s.adj[v] = set(); return v
    def clique(s, vs):
        for a in vs:
            for b in vs:
                if a != b: s.adj[a].add(b)
    def gate(s, a, b):
        o = s.var(); s.clique((o, a, b)); return o
    def ha(s, a, b):
        return s.gate(a, b), s.gate(a, b)          # sum, carry
    def fa(s, a, b, c):
        t = s.gate(a, b); sm = s.gate(t, c)
        g1 = s.gate(a, b); g2 = s.gate(t, c); co = s.gate(g1, g2)
        return sm, co

def ripple(C, A, B):
    """add bit lists (LSB first, None = 0); returns list of sum bits."""
    L = max(len(A), len(B)); out = []; c = None
    for i in range(L):
        bits = [x for x in (A[i] if i < len(A) else None,
                            B[i] if i < len(B) else None, c) if x is not None]
        if len(bits) == 0: out.append(None); c = None
        elif len(bits) == 1: out.append(bits[0]); c = None
        elif len(bits) == 2: s_, c = C.ha(*bits); out.append(s_)
        else: s_, c = C.fa(*bits); out.append(s_)
    if c is not None: out.append(c)
    return out

def multiplier(n, kind):
    C = Circ()
    x = [C.var() for _ in range(n)]; y = [C.var() for _ in range(n)]
    pp = [[C.gate(x[i], y[j]) for i in range(n)] for j in range(n)]  # pp[j][i] weight i+j
    if kind == 'array':
        acc = list(pp[0])
        for j in range(1, n):
            acc = ripple(C, acc, [None]*j + pp[j])
        prod = acc
    else:  # wallace
        cols = [[] for _ in range(2*n)]
        for j in range(n):
            for i in range(n): cols[i+j].append(pp[j][i])
        while max(len(c) for c in cols) > 2:
            new = [[] for _ in range(len(cols) + 1)]
            for k, c in enumerate(cols):
                m = 0
                while len(c) - m >= 3:
                    s_, co = C.fa(c[m], c[m+1], c[m+2]); m += 3
                    new[k].append(s_); new[k+1].append(co)
                if len(c) - m == 2:
                    s_, co = C.ha(c[m], c[m+1]); m += 2
                    new[k].append(s_); new[k+1].append(co)
                new[k].extend(c[m:])
            cols = new if new[-1] else new[:-1]
        A = [c[0] if len(c) > 0 else None for c in cols]
        B = [c[1] if len(c) > 1 else None for c in cols]
        prod = ripple(C, A, B)
    prod = (prod + [None]*(2*n))[:2*n]
    circ_vertices = set(C.adj)
    z = []
    for b in prod:
        v = C.var(); z.append(v)
        if b is not None: C.clique((v, b))
    return C.adj, x, y, z, circ_vertices

def copy(adj): return {v: set(nb) for v, nb in adj.items()}

def ub_min_degree(adj):
    g = copy(adj); h = [(len(nb), v) for v, nb in g.items()]; heapq.heapify(h)
    width = 0
    while g:
        d, v = heapq.heappop(h)
        if v not in g or d != len(g[v]): continue
        nb = g.pop(v); width = max(width, len(nb))
        for a in nb: g[a].discard(v)
        for a in nb:
            for b in nb:
                if a != b: g[a].add(b)
        for a in nb: heapq.heappush(h, (len(g[a]), a))
    return width

def fill(g, v):
    nb = list(g[v]); f = 0
    for i in range(len(nb)):
        ga = g[nb[i]]
        for j in range(i+1, len(nb)):
            if nb[j] not in ga: f += 1
    return f

def ub_min_fill(adj):
    g = copy(adj); key = {v: fill(g, v) for v in g}
    h = [(key[v], len(g[v]), v) for v in g]; heapq.heapify(h); width = 0
    while g:
        f, d, v = heapq.heappop(h)
        if v not in g or f != key[v] or d != len(g[v]): continue
        nb = g.pop(v); del key[v]; width = max(width, len(nb))
        for a in nb: g[a].discard(v)
        for a in nb:
            for b in nb:
                if a != b: g[a].add(b)
        touched = set(nb)
        for a in nb: touched |= g[a]
        for a in touched:
            key[a] = fill(g, a); heapq.heappush(h, (key[a], len(g[a]), a))
    return width

def lb_mmd_plus(adj):
    """MMD+ with min-d contraction: repeatedly take a min-degree vertex v,
    record its degree, contract it into its min-degree neighbour."""
    g = copy(adj); best = 0
    h = [(len(nb), v) for v, nb in g.items()]; heapq.heapify(h)
    while len(g) > 1:
        d, v = heapq.heappop(h)
        if v not in g or d != len(g[v]): continue
        best = max(best, d)
        nb = g.pop(v)
        if not nb: continue
        u = min(nb, key=lambda a: len(g[a]))
        for a in nb: g[a].discard(v)
        for a in nb:
            if a != u: g[u].add(a); g[a].add(u)
        for a in nb | {u}: heapq.heappush(h, (len(g[a]), a))
    return best

def fit(ns, vals):
    k = len(ns); mx = sum(ns)/k; my = sum(vals)/k
    sxx = sum((a-mx)**2 for a in ns); sxy = sum((a-mx)*(b-my) for a, b in zip(ns, vals))
    sl = sxy/sxx; ic = my - sl*mx
    ss_res = sum((b-(sl*a+ic))**2 for a, b in zip(ns, vals)); ss_tot = sum((b-my)**2 for b in vals)
    return sl, ic, 1 - ss_res/ss_tot if ss_tot else 1.0

def main():
    ns = [int(a) for a in sys.argv[1:]] or [4, 6, 8, 10, 12, 16, 20, 24, 32]
    MF_MAX = 24
    print("kind     n   |V|    lbMMD+  ubMinDeg ubMinFill  ub(G-Z)  ub/n")
    res = {}
    for kind in ('array', 'wallace'):
        for n in ns:
            t = time.time()
            adj, x, y, z, cv = multiplier(n, kind)
            gz = {v: adj[v] & cv for v in cv}
            lb = lb_mmd_plus(adj); md = ub_min_degree(adj)
            mf = ub_min_fill(adj) if n <= MF_MAX else None
            ub = min(md, mf) if mf is not None else md
            ubz = ub_min_degree(gz)
            if n <= MF_MAX: ubz = min(ubz, ub_min_fill(gz))
            res.setdefault(kind, []).append((n, ub, lb, mf is not None))
            print(f"{kind:8s} {n:3d} {len(adj):6d} {lb:6d} {md:9d} {str(mf):>9s} {ubz:8d}  {ub/n:5.2f}"
                  f"   ({time.time()-t:.1f}s)", flush=True)
    for kind, rows in res.items():
        mf_rows = [r for r in rows if r[3]]
        sl, ic, r2 = fit([r[0] for r in mf_rows], [r[1] for r in mf_rows])
        print(f"linear fit {kind} (rows with min-fill, n<={MF_MAX}): tw_ub ~ {sl:.3f} n + {ic:.2f}  (R^2={r2:.4f})")
        sl, ic, r2 = fit([r[0] for r in rows], [r[2] for r in rows])
        print(f"linear fit {kind} (all rows): tw_lb(MMD+) ~ {sl:.3f} n + {ic:.2f}  (R^2={r2:.4f})")
        # power-law exponent from the two largest min-fill points
        import math
        (n1, u1), (n2, u2) = [(r[0], r[1]) for r in mf_rows[-2:]]
        print(f"  local log-log slope of ub between n={n1} and n={n2}: {math.log(u2/u1)/math.log(n2/n1):.3f}")
    assert all(r[2] <= r[1] for rows in res.values() for r in rows)
    print("PASS (bounds computed; lb <= ub checked)" )

if __name__ == '__main__':
    main()
