#!/usr/bin/env python3
"""Graph-side checks for path A (routing-rank no-go, PR #67).

(A) Exact MAIS of the routing digraph on the full label array [r]^Q x [r]^Q.
    MAIS <= minrank(I-M) <= construction rank (from pathA/check_routing_rank.out).
    Also the score bound and the best product-potential bound.
(B) Occurrence/variable-level potential bound (Theorem A'): random instances,
    random M0 on the allowed support, exact ranks over GF(p).
(C) Cancellation-free crossing pieces: cut ranks equal support size (rho = 1).
(D) Potential-obstruction width on several host graphs (illustration only).
(E) Random sparse codes with Q >> log n: fraction of occurrences with no
    variable schedule within distance k of the clause schedule (unroutable).
Exact integer / GF(p) arithmetic throughout, except (D) which is counting.
"""
import itertools, random, sys
import pulp, networkx as nx

P = 1_000_003

def ham(x, y): return sum(a != b for a, b in zip(x, y))

def rank_mod(rows, p=P):
    rows = [dict(r) for r in rows if r]
    piv = {}
    rk = 0
    for r in rows:
        r = {c: v % p for c, v in r.items() if v % p}
        while r:
            c = min(r)
            if c in piv:
                pr = piv[c]; f = r[c]
                for cc, vv in pr.items():
                    nv = (r.get(cc, 0) - f * vv) % p
                    if nv: r[cc] = nv
                    else: r.pop(cc, None)
            else:
                inv = pow(r[c], p - 2, p)
                piv[c] = {cc: vv * inv % p for cc, vv in r.items()}
                rk += 1
                break
    return rk

# ---------------- (A) exact MAIS ----------------
def mais(nodes, succ, selfloop):
    """max acyclic induced subgraph via ILP with lazy cycle cuts (HiGHS)."""
    prob = pulp.LpProblem("mais", pulp.LpMaximize)
    x = {u: pulp.LpVariable(f"x{u}", cat="Binary") for u in nodes}
    prob += pulp.lpSum(x.values())
    for u in nodes:
        if selfloop[u]: prob += x[u] == 0
    for u in nodes:
        for w in succ[u]:
            if w > u and u in succ[w]:
                prob += x[u] + x[w] <= 1
    it = 0
    while True:
        it += 1
        prob.solve(pulp.HiGHS(msg=False))
        S = [u for u in nodes if x[u].value() > 0.5]
        Sset = set(S)
        G = nx.DiGraph(); G.add_nodes_from(S)
        G.add_edges_from((u, w) for u in S for w in succ[u] if w in Sset and w != u)
        cyc = []
        try:
            for _ in range(200):
                c = nx.find_cycle(G)
                cyc.append([e[0] for e in c])
                G.remove_edge(*c[0])
        except nx.NetworkXNoCycle:
            pass
        if not cyc:
            return len(S), it
        for c in cyc:
            prob += pulp.lpSum(x[u] for u in c) <= len(c) - 1

def routing_digraph(r, Q, k):
    pts = list(itertools.product(range(r), repeat=Q))
    nodes = list(range(len(pts) ** 2))
    lab = [(a, d) for a in pts for d in pts]
    succ = {u: [w for w in nodes if ham(lab[u][1], lab[w][0]) <= k] for u in nodes}
    selfloop = {u: ham(lab[u][1], lab[u][0]) <= k for u in nodes}
    return pts, lab, nodes, succ, selfloop

def potential_count(lab, phi, k):
    return sum(1 for a, d in lab if phi(a) - phi(d) > k)

def best_product_potential(r, Q, k, lab):
    best = 0
    subsets = [frozenset(s) for m in range(1, r) for s in itertools.combinations(range(r), m)]
    for Hs in itertools.product(subsets, repeat=Q):
        phi = lambda x: sum(1 for t in range(Q) if x[t] in Hs[t])
        best = max(best, potential_count(lab, phi, k))
    return best

def score_bound(r, Q, k, lab):
    # threshold potential sum_t [x_t >= ceil(r/2)]: one product potential
    phi = lambda x: sum(1 for v in x if v >= (r + 1) // 2)
    return potential_count(lab, phi, k)

CONSTR = {(2,1,0):1,(2,2,0):6,(2,2,1):1,(2,3,1):10,(3,2,1):9,(2,4,2):13,(2,4,1):94,(4,2,1):36}

def part_A():
    print("== (A) exact MAIS of routing digraph (full array) ==")
    print("   MAIS <= minrank(I-M) over every field <= construction rank")
    ok = True
    for (r, Q, k) in [(2,1,0),(2,2,0),(2,2,1),(2,3,1),(3,2,1),(2,4,2),(4,2,1),(2,4,1)]:
        pts, lab, nodes, succ, selfloop = routing_digraph(r, Q, k)
        m, it = mais(nodes, succ, selfloop)
        pp = best_product_potential(r, Q, k, lab)
        sb = score_bound(r, Q, k, lab)
        c = CONSTR.get((r, Q, k))
        good = (sb <= pp <= m) and (c is None or m <= c)
        ok &= good
        print(f"r={r} Q={Q} k={k} D={len(nodes)}: score={sb} best_product_potential={pp} "
              f"MAIS={m} construction_rank={c} ilp_rounds={it} {'OK' if good else 'FAIL'}")
    return ok

# ---------------- (B) Theorem A' ----------------
def part_B(trials=60, seed=1):
    print("== (B) occurrence/variable-level potential bound, random instances ==")
    rnd = random.Random(seed)
    ok = True; tight = 0; tot = 0
    for tr in range(trials):
        r = rnd.choice([2, 3]); Q = rnd.choice([3, 4, 5]); k = rnd.choice([0, 1, 2])
        n = rnd.randint(6, 14); mcl = rnd.randint(3, 8)
        a = [tuple(rnd.randrange(r) for _ in range(Q)) for _ in range(n)]
        dsch = [tuple(rnd.randrange(r) for _ in range(Q)) for _ in range(mcl)]
        occ = [(e, v) for e in range(mcl) for v in rnd.sample(range(n), 3)]
        # random M0 on allowed support: M0[o,o'] free iff dist(d_e(o), a_v(o')) <= k
        dens = rnd.choice([0.3, 0.7, 1.0])
        rho = []
        for (e, v) in occ:
            row = {v: 1}
            for (e2, v2) in occ:
                if ham(dsch[e], a[v2]) <= k and rnd.random() < dens:
                    row[v2] = (row.get(v2, 0) - rnd.randrange(1, P)) % P
            rho.append(row)
        # minimal sum_h c_{t,h}
        ct = []
        for t in range(Q):
            s = 0
            for h in range(r):
                Sh = {v for v in range(n) if a[v][t] == h}
                rows = [{c: x for c, x in rho[i].items() if c in Sh}
                        for i, (e, v) in enumerate(occ) if dsch[e][t] != h]
                s += rank_mod(rows)
            ct.append(s)
        # all product potentials
        subsets = [frozenset(sb) for m_ in range(1, r) for sb in itertools.combinations(range(r), m_)]
        best_lb = 0
        for Hs in itertools.product(subsets, repeat=Q):
            phi = lambda x: sum(1 for t in range(Q) if x[t] in Hs[t])
            A = [(e, v) for (e, v) in occ if phi(a[v]) - phi(dsch[e]) > k]
            for t in range(Q):
                lb = len({v for (e, v) in A if a[v][t] != dsch[e][t]})
                best_lb = max(best_lb, lb)
                if ct[t] < lb:
                    ok = False
                    print("FAIL", tr, t, ct[t], lb)
        tot += 1; tight += (max(ct) == best_lb)
    print(f"trials={tot}: sum_h c_(t,h) >= #{{v: some A_phi occurrence of v cut at t}} for every t, phi: "
          f"{'OK' if ok else 'FAIL'}; max_t bound equals max_t c_t in {tight}/{tot}")
    return ok

# ---------------- (C) cancellation-free crossing pieces ----------------
def part_C(seed=2):
    print("== (C) support-contained crossing pieces have cut ranks = support size ==")
    rnd = random.Random(seed); ok = True
    for N in [2, 3, 4]:
        supp = [(x, y) for x in range(N) for y in range(N)]  # E[x,y,u,v]=[x=v][y=u]
        for _ in range(20):
            S = rnd.sample(supp, rnd.randint(1, N * N))
            W = {s: rnd.randrange(1, P) for s in S}
            # cut1 rows (x,y) cols (u,v); entry at (x,y),(y,x)
            c1 = [{(y * N + x): W[(x, y)]} for (x, y) in S]
            # cut2 rows (x,u) cols (y,v): entry (x,u=y),(y,v=x)
            c2 = [{(y * N + x): W[(x, y)]} for (x, y) in S]
            r1 = rank_mod(c1)
            r2 = rank_mod(c2)
            ok &= (r1 == len(S) == r2)
    print(f"N=2..4, 60 random weighted supports: r1 = r2 = |S|: {'OK' if ok else 'FAIL'}"
          " => sum_j sqrt(r1 r2) = sum_j |S_j| >= N^2, rho = 1 without cancellation")
    return ok

# ---------------- (D) potential-obstruction width on host graphs ----------------
def part_D(eps=0.05, seed=3):
    print(f"== (D) host-graph potential width at eps={eps} (illustration; distance potentials only) ==")
    rnd = random.Random(seed)
    hosts = {
        "Hamming H(10,2)": nx.hypercube_graph(10),
        "Johnson J(12,6)": None,
        "random 4-regular, 1024": nx.random_regular_graph(4, 1024, seed=seed),
        "torus 32x32": nx.grid_2d_graph(32, 32, periodic=True),
        "cycle C_1024": nx.cycle_graph(1024),
    }
    J = nx.Graph()
    sets = [frozenset(c) for c in itertools.combinations(range(12), 6)]
    idx = {s: i for i, s in enumerate(sets)}
    for s in sets:
        for i in s:
            for j in set(range(12)) - s:
                J.add_edge(idx[s], idx[(s - {i}) | {j}])
    hosts["Johnson J(12,6)"] = J
    for name, G in hosts.items():
        G = nx.convert_node_labels_to_integers(G)
        V = G.number_of_nodes()
        roots = rnd.sample(range(V), 6)
        diam = 0; width = 0
        for z in roots:
            dist = nx.single_source_shortest_path_length(G, z)
            vals = sorted(dist.values()); diam = max(diam, vals[-1])
            # Pr_{x,y}[phi(y)-phi(x) > k] for phi = dist(z,.)
            from collections import Counter
            cnt = Counter(vals)
            for k in range(vals[-1] + 1):
                pr = sum(cnt[i] * cnt[j] for i in cnt for j in cnt if j - i > k) / V**2
                if pr <= eps:
                    width = max(width, k); break
        print(f"{name:26s} V={V:5d} diam>={diam:4d} width(eps)~{width:3d} width/diam={width/diam:.2f}")
    print("   Theorem A' forces the schedule metric to be Hamming; these numbers only show"
          " what a free host metric would allow, they are not realizable in the preprint.")
    return True

# ---------------- (E) sparse codes, Q >> log n ----------------
def part_E(seed=4):
    print("== (E) random sparse codes: clause schedules with no variable schedule within k ==")
    rnd = random.Random(seed)
    for n, Q in [(256, 8), (256, 16), (256, 32), (256, 64)]:
        a = [tuple(rnd.randrange(2) for _ in range(Q)) for _ in range(n)]
        dsch = [tuple(rnd.randrange(2) for _ in range(Q)) for _ in range(200)]
        import math
        k = int(math.ceil(2 * math.sqrt(Q)))
        iso = sum(1 for d in dsch if all(ham(d, x) > k for x in a)) / len(dsch)
        print(f"n={n} Q={Q} k=ceil(2sqrt(Q))={k}: fraction of clause schedules with empty k-ball = {iso:.2f}")
    print("   empty k-ball => every occurrence of that clause has rho_o = e_v (unroutable),"
          " so it is counted directly; random codes need structure to use Q >> log n")
    return True

# ---------------- (F) anticoncentration of the potential walk Z_w ----------------
def part_F():
    """Z_w = sum of w iid increments in {-1,0,+1}, +-1 each w.p. p_r =
    floor(r/2)ceil(r/2)/(r(r-1)) (random half-set H_t, one cut coordinate).
    Check exactly: Pr[Z_w >= k+1] >= 1/2 exp(-9(k+1)^2/(2w) - 72) for k+1 <= w/4."""
    import math
    print("== (F) Lemma B' for Z_w: Pr[Z_w > k] >= 1/2 exp(-9(k+1)^2/(2w)-72), k+1 <= w/4 ==")
    ok = True; worst = float("inf")
    for r in [2, 3, 4, 5, 16]:
        num = (r // 2) * ((r + 1) // 2); den = r * (r - 1)   # p = num/den
        for w in [64, 100, 256]:
            # integer weights: +1 -> num, -1 -> num, 0 -> den-2num ; total den^w
            dist = {0: 1}
            for _ in range(w):
                nd = {}
                for z, c in dist.items():
                    for dz, wt in ((1, num), (-1, num), (0, den - 2 * num)):
                        if wt:
                            nd[z + dz] = nd.get(z + dz, 0) + c * wt
                dist = nd
            tot = den ** w
            kmax = w // 4 - 1
            tail = sum(c for z, c in dist.items() if z > kmax + 1)
            for k in range(kmax, -1, -1):  # tail = Pr[Z >= k+1]
                tail += dist.get(k + 1, 0)
                lhs = math.log(tail) - math.log(tot)
                rhs = math.log(0.5) - 9 * (k + 1) ** 2 / (2 * w) - 72
                worst = min(worst, lhs - rhs)
                if lhs < rhs: ok = False
    print(f"r in {{2,3,4,5,16}}, w in {{64,100,256}}, all k: {'OK' if ok else 'FAIL'}; min log-margin {worst:.2f}")
    return ok

if __name__ == "__main__":
    ok = part_A() & part_B() & part_C() & part_D() & part_E() & part_F()
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)
