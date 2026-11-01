# Path A as graph problems: what known graph theory says

Follow-up to PR #67 (routing-rank no-go). Repository base commit `4c3d135`. Preprint `src/openai.tex`, sha256 prefix `b2e6159d4fb7ea6a`.
Labels: **PROVEN** (proof here or in `pathA/note.md`), **VERIFIED** (exact computation in this folder), **KNOWN** (published, cited), **CONJECTURE** (open or heuristic). "Not found" never means "new".

## 0. Verdicts

| Angle | Verdict |
|---|---|
| 1. Host graphs and minrank | Theorem A **is** the maximum-acyclic-induced-subgraph (MAIS) lower bound on minrank (KNOWN, Bar-Yossef–Birk–Jayram–Kol). The graph lens yields a stronger theorem, A′ (PROVEN): it covers every label set, every placement or merging of variables, and every support-respecting M_0. Door (c) closes, and door (b) closes whenever the number of evaluation points is Q = O(log n). θ, Haemers and orthogonal rank are upper-bound tools on capacity, so they cannot strengthen a rank lower bound. No host graph inside the preprint beats the √Q width: the accounting forces the Hamming metric on schedules, and product potentials are universal there. Among these geometry variants, (b′) sparse label sets with Q ≫ log n remains open. Point-dependent measurements, door (a), are outside this argument. |
| 2. Crossing lane | Every monotone or support-based graph parameter gives ρ = 1. That covers rectangle covers, Boolean rank, Körner graph entropy, Hansel/Fredman–Komlós and fractional covers (PROVEN: cancellation-free pieces give ρ = 1). No known graph theorem gives 1/√2 optimality or breaks it. The graded window (γ > 1/2, L < 4) has graph realizations of the growth: Fibonacci-type differential posets and the repo's own Fibonacci and Catalan path families. Whether any of them reaches energy range L ≤ 2^{1/γ} is open. |

## 1. Dictionary: routing rank is minrank of a digraph

- Take labels u = (a_u, d_u) ∈ ([r]^Q)². Build the routing digraph Γ_k with an arc u → w (w ≠ u) iff dist_H(d_u, a_w) ≤ k, and a self-loop at u iff dist_H(d_u, a_u) ≤ k.
  - Step 3 of the preprint (l.1236–1262) allows M_0[u,w] ≠ 0 exactly on the arcs.
  - So min rank(I − M_0) is the minrank of Γ_k, with free diagonal at the looped vertices. This is the index-coding formulation (KNOWN: Bar-Yossef, Birk, Jayram, Kol, FOCS 2006 / IEEE-IT 2011).
- **KNOWN:** for every field, minrank ≥ MAIS. An acyclic induced set gives a principal submatrix that is unitriangular up to permutation (same paper, and the inequality chain MAIS ≤ β ≤ minrk in arXiv:1711.06487).
- **PROVEN (identification):** Theorem A / A1 of PR #67 is this bound. The score set {u : s(a_u) > s(d_u) + k} is induced-acyclic, ordered by score.
  - The PR note already says it "may be folklore". It should now cite MAIS explicitly as the known technique. What is specific to the PR is the choice of acyclic set (a potential) and the Gaussian evaluation.
- **KNOWN, why the "graph bounds" do not help a lower bound:**
  - Haemers' minrank and Lovász θ are both *upper* bounds on Shannon capacity. θ and minrank are incomparable (C5: θ = √5 < 3; complement of the Schläfli graph: minrk ≤ 7 < 9 = θ; Haemers 1978/79).
  - There are graphs with constant θ and minrank ≥ n^δ (Golovnev–Regev–Weinstein, arXiv:1802.03920). So θ cannot lower-bound minrank, and orthogonal-rank variants sit in the same capacity sandwich.
  - The universal *lower* bounds are α (undirected) and MAIS (directed).
  - The fooling-set bound rank ≥ √(size) (KNOWN, Dietzfelbinger–Hromkovič–Schnitger 1996) is the only other pattern-only bound found. It gives rank ≳ √D, which is useless against the design budget rank ≤ εD.
  - Field dependence of minrank (KNOWN, e.g. Lubetzky–Stav) is irrelevant, because MAIS is field-free.
- **Finite controls (part A): HiGHS-reported MAIS optima, checked acyclic witnesses, and comparisons with the source bounds/construction values** (full array):

| r,Q,k | D | best product potential | PR #67 score_LB | **MAIS (HiGHS)** | construction rank S·g(W_0) |
|---|---|---|---|---|---|
| 2,2,0 | 16 | 5 | 6 | 6 | 6 |
| 2,3,1 | 64 | 7 | 8 | 8 | 10 |
| 3,2,1 | 81 | 4 | 6 | 8 | 9 |
| 2,4,2 | 256 | 9 | 10 | 10 | 13 |
| 4,2,1 | 256 | 16 | 27 | 30 | 36 |
| 2,4,1 | 256 | 37 | 45 | 45 | 94 |

  - On these cases MAIS improves the potential bounds only by constants. The true minrank lies between MAIS and the construction rank; the gap is open at (2,4,1).
  - The asymptotic degree ratio stays ≤ √2·(1+o(1)) as in PR #67.
  - The "best product potential" column uses only the half-set potentials Σ_t [x_t ∈ H_t]. The PR's score_LB comes from its own `free_bound` routine in `check_routing_rank.py`.

## 2. Theorem A′: the variable-level potential bound

**Setting.** Every design built from the preprint's Steps 2–4 fits this setting, with any number of copies of the array, any label subset, and any merging of coordinates that hold the same variable.
- Variables v have schedules a_v ∈ [r]^Q. Clauses e have d_e ∈ [r]^Q.
- Each occurrence o = (e,v) sits at a coordinate whose a-part is a_v. By definition, the group of ζ_v at point t is a_v(t).
- M_0 is any matrix with M_0[o,o′] ≠ 0 only if dist(d_{e(o)}, a_{v(o′)}) ≤ k.
- Substituting coordinates by variables gives rows ρ_o = e_v − Σ_{o′} M_0[o,o′] e_{v(o′)} ∈ K^V.
- At point t, the measurement condition (preprint condition 3/4) needs ρ_o|_{S_{t,h}} ∈ rowspace(C_{t,h}) for every group h ≠ d_e(t), where S_{t,h} = {v : a_v(t) = h}.

**Theorem A′ (PROVEN).** Take any potential φ: [r]^Q → ℝ with φ(y) ≤ φ(x) + k whenever dist_H(x,y) ≤ k (e.g. any 1-Lipschitz φ). Let A_φ = {o = (e,v) : φ(a_v) > φ(d_e) + k}. Then, for every t,

  Σ_h c_{t,h} ≥ #{ v : some o = (e,v) ∈ A_φ has a_v(t) ≠ d_e(t) }.

*Proof.*
1. For such v, pick o and put h = a_v(t) ≠ d_e(t). Then ρ_o|_{S_h} is required.
2. Its coefficient at v is 1. The occurrences of v at near coordinates would need dist(d_e, a_v) ≤ k, which gives φ(a_v) ≤ φ(d_e) + k, a contradiction.
3. Every other nonzero coefficient sits at a v′ with dist(d_e, a_{v′}) ≤ k, so φ(a_{v′}) ≤ φ(d_e) + k < φ(a_v).
4. Within S_h, order the chosen rows by φ(a_v). They are unitriangular, hence independent. Sum over h. ∎

**Averaging (PROVEN).** Use random half-sets H_t, put φ = Σ_t [x_t ∈ H_t], and let each variable occur ≤ b times. Then

  max_t Σ_h c_{t,h} ≥ (1/(bQ)) Σ_o w(o) · Pr[Z_{w(o)} > k].

- Here w(o) = #{t : a_v(t) ≠ d_e(t)} is the number of cut points.
- Z_w is a sum of w iid steps: ±1 each with probability p_r = ⌊r/2⌋⌈r/2⌉/(r(r−1)) ≥ 1/4, and 0 otherwise.
- The bound holds for **every** code X = {a_v} ∪ {d_e}, whatever its intrinsic geometry. This universality is the graph-theoretic content: product potentials are 1-Lipschitz on any subset of the Hamming graph.

**Corollary C′ (ordinary consequence under the stated expansion and accounting inputs).**
- Inputs:
  - The instance family has an expanding primal graph. Then, at every t, balanced groups (parts ≤ n/2) cut ≥ c_0 n occurrences.
    - KNOWN for random 3-CNF at constant density, up to a negligible non-giant part. An ETH algorithm must handle such instances.
    - Constants are not computed.
  - Lemma B′: Pr[Z_w > k] ≥ ½ exp(−9(k+1)²/(2w) − 72) for w ≥ 64 and k+1 ≤ w/4.
    - Integer distribution counts with a floating log comparison for r ∈ {2,3,4,5,16}, w ∈ {64,100,256}, all k (part F, reported min log-margin 71.96).
    - The [uniform companion proof](lazy-potential-tail-proof.md) now establishes this for all w>=64 and every real p in [1/4,1/2], with the stronger constant -70. The finite grid is not its general justification.
- Argument: Σ_o w(o) ≥ c_0 nQ. So ≥ c_0 n/2 occurrences have w(o) ≥ W := c_0 Q/(2C_1), where C_1 n bounds the number of occurrences. With the design budget Σ_h c_{t,h} ≤ 2τn:

  (k+1)² ≥ (2W/9)(ln(1/τ) − O(1)),  hence δ ≥ R/Q ≥ Ω(√(ln(Q/δ)/Q)).

- With Q = O(log n), this is exactly Corollary C's wall, δ = Ω(√(log log n / log n)), i.e. α ≤ 1/2. It now holds for **any label set and any merging**, not just the full product array.
- **Under these hypotheses, door (c) is closed** up to constants.
- **Under these hypotheses, door (b) is closed for Q = O(log n).** Corollary C's branch (i) already forces Q = O(log n) for the full array.

**VERIFIED (part B):** 60 random instances (r ∈ {2,3}, Q ∈ {3,4,5}, k ∈ {0,1,2}, random formulas, random M_0 on the allowed support at three densities). For every t and every product potential, the exact GF(p) value of Σ_h c_{t,h} is at least the bound. The bound is attained in 3/60 instances, so it is valid but not tight.

## 3. Host graphs: what would beat √Q, and why the preprint cannot use it

- **Graph form of the obstruction (PROVEN from A′).** The routing-rank lower bound is a concentration profile of the host metric: sup over k-Lipschitz potentials φ of Pr_{o}[φ(a) − φ(d) > k].
- **KNOWN concentration facts:**
  - Hamming cube: Gaussian width ≍ √Q (McDiarmid / Talagrand). This is the CLT √Q of PR #67, and it is tight there.
  - Expanders with spectral gap λ: width O(log(1/ε)/√λ), independent of size (Alon–Milman 1985; Gromov–Milman). Relative to the diameter (≍ log N), a free expander host would be strictly better than the cube.
  - Johnson J(Q,s) and Kneser-type slices: Gaussian like the cube.
  - Bounded-generator Cayley graphs of abelian groups (tori, cycles): poor concentration (width ≍ diameter).
  - HDX: Gaussian-type concentration (cube-like), no better than expanders.
  - Small independence / large chromatic number: necessary for small minrank (minrank ≥ α, resp. MAIS) but not sufficient. It is not a route to constructions.
- **VERIFIED illustration (part D, distance potentials only, N ≈ 1024).** Width/diameter at ε = 0.05: H(10,2) 0.40, J(12,6) 0.33, random 4-regular 0.38, torus 0.47, C_1024 0.68. At this size the expander/cube difference is invisible. These numbers do not bear on the asymptotics.
- **What the accounting requires of G (PROVEN from the preprint lines):**
  1. The multiplier ∏_{s∈E_{u,w}}(1 − T/s) has degree |E_{u,w}| = #{points where the groups differ} (l.1236–1262). So the degree metric **is** Hamming distance on schedules. It is not a graph distance one may choose.
  2. Condition 3 is per evaluation point and per group (r groups at each of Q points). So the host is a subset X ⊆ H(Q,r) with the **induced Hamming metric**.
  3. Balance 12D² ≤ M and hD ≤ 2M ≤ 6n (l.1293, l.1340–1348) give at most √(n/2) labels. The full array then has Q ≤ log n/(4 log r).
- **Consequence.**
  - By (1)–(2), the product potentials of Theorem A′ apply to every admissible host. So no Cayley, Johnson/Kneser, HDX or "expander-like" code X ⊆ H(Q,r) beats the √Q width (PROVEN, via A′).
  - An expander metric cannot even be realized: by (1) it would have to embed isometrically into a Hamming metric, hence into ℓ1. Expanders need distortion Ω(log N) there (KNOWN, Linial–London–Rabinovich 1995).

## 4. The one remaining door, as a graph problem

**(b′) CONJECTURE/open.**
- Requirements:
  - a code X ⊆ H(Q,r) with ≤ √n schedules and Q ≥ (log n)^{1+Ω(1)};
  - clause schedules d_e;
  - k = R−1 ≤ δQ;
  - routing digraph Δ (jump arcs a_v → d_e for occurrences; closeness arcs d_e → a_{v′} when dist ≤ k) with minrank of the occurrence system ≤ 2τn.
- Theorem A′ alone gives only δ ≥ Ω(√(log Q / Q)), so it does not exclude α > 1/2 here.
- Necessary graph conditions (PROVEN):
  - An occurrence whose closeness neighbourhood is empty has ρ_o = e_v and is counted directly.
  - More generally, the occurrences whose jump arc lies on no directed cycle of Δ form an acyclic set, so they count toward the rank (MAIS).
  - Small rank therefore needs almost every jump arc to lie on a cycle of length ≥ 2 through near pairs.
- Heuristic against (VERIFIED for random codes, part E): n = 256, k = ⌈2√Q⌉.
  - 98% of random clause schedules have an empty k-ball at Q = 64, against 0% at Q ≤ 32.
  - Random codes therefore stop routing once Q ≫ log n, because ball volume 2^{H(δ)Q} ≪ 2^Q/|X|.
  - A structured (non-random) code is needed. Repetition codes make things worse, since they scale k and Q together.
  - Whether any code gives a cycle-rich Δ with low minrank is open. The best candidate statement to prove next: "for every X ⊆ H(Q,r) with |X| ≤ √n, the minrank of Δ is ≥ τn unless Q = O(log n)" (CONJECTURE).

Door (a), a t-dependent non-single-matrix measurement, is untouched by the graph view and stays open as in PR #67.

## 5. Crossing lane (#52) as graphs and tensor networks

- **TN dictionary (PROVEN, trivial).** E(x,y,u,v) = [x=v][y=u] is two bonds, x–v and y–u, on the square of parties.
  - Both bonds cross cut 1 {x,y}|{u,v} and cut 2 {x,u}|{y,v}, so each cut rank is N².
  - Across the third cut {x,v}|{y,u}, E has rank 1.
  - A piece's TN bond dimensions only *upper*-bound its cut ranks. TN/min-cut arguments therefore give no lower bound on Σ√(r1 r2).
- **Cancellation-free pieces (PROVEN; VERIFIED part C for N = 2..4).**
  - A piece supported inside supp(E) = {(x,y,y,x)} with weights W_j has cut-1 rank = cut-2 rank = |supp W_j|. Both flattenings are partial permutation patterns.
  - So Σ_j √(r1 r2) = Σ_j |S_j| ≥ N², i.e. ρ = 1.
  - Nonnegative pieces summing to the 0/1 tensor E are automatically of this kind.
  - Hence every monotone/union parameter gives ρ = 1 and can say nothing below it. This includes rectangle and biclique covers, Boolean rank, Körner graph entropy with Hansel/Fredman–Komlós-type cover bounds, and fractional covers.
  - The repo's 0.710 decomposition uses signed cancellation, which these tools cannot see.
- **Log-rank / fooling sets (KNOWN, not applicable).**
  - They relate the rank of one matrix to its support or communication. The quantity here is a sum over pieces whose supports leave supp(E).
  - The only general bound for the sum remains subadditivity, i.e. the elementary ρ ≥ 1/2 of PR #67.
  - I found no graph theorem implying 1/√2 optimality, and none implying a construction below it.
- **Graded window (γ > 1/2, L < 4).**
  - KNOWN: rank sizes of r-differential posets (graded graphs with DU − UD = rI) lie between n^a e^{2√(rn)} (Stanley–Zanello 2012; γ = 1/2, Young's lattice) and the Fibonacci numbers F_r(n) (Byrnes 2012, sharp; γ = 1, Young–Fibonacci lattice, Stanley 1988 / Fomin 1988).
  - PROVEN, trivially: any bounded-degree graded graph has at most Δ^n paths of length n, so γ ≤ 1 for graph/path realizations.
  - The repo already has γ = 1 path families: bounded-spin and paired-path Fibonacci crossings, Catalan cup codes, fusion seams.
  - Caveat: Appendix B's algebra is sl₂-type, [D,U] = 2H (`homogeneous-duals.md` eq. 3), not the Heisenberg relation. So Young–Fibonacci is not a drop-in replacement.
  - Kept at L = 4, a γ = 1 grading gives the *worse* ρ = 2^{γ−1} = 1 (repo eq. 7).
  - The window needs L ≤ 2^{1/γ}, i.e. L ≤ 2 for Fibonacci growth. At γ = 1, L = 2 the rank budget (repo eq. 3) is tight: the extreme pieces must use essentially all of V_{2h}.
  - CONJECTURE/open: measure the energy range L of the repo's Fibonacci/Catalan path crossings in the sense of eq. (4). Only an L < 2^{1/γ} family can enter the window.

## 6. Rerun

```
pip install "pulp<3" highspy networkx
python tools/routing_rank_nogo/check_graph_minrank.py   # ~40 s, ends "ALL CHECKS PASS"
```

- Checker sha256 prefix `0341cda7b755f52b`. Output in `check_graph_minrank.out`.
- Arithmetic: GF(1000003) for ranks (a GF(p) rank lower-bounds the ℚ-rank, and the bounds hold over every field), exact integers for the Lemma B′ tails, ILP (HiGHS) plus an explicit acyclicity check for MAIS.
- Part D is a float illustration only.

## 7. Literature searched (this turn)

- Index coding, minrank and MAIS: Bar-Yossef–Birk–Jayram–Kol (FOCS 2006; IEEE-IT 57, 2011); arXiv:1711.06487; arXiv:1309.4166.
- Haemers / θ: Haemers 1978/79; Golovnev–Regev–Weinstein arXiv:1802.03920; arXiv:2002.02743.
- Fooling sets vs rank: Dietzfelbinger–Hromkovič–Schnitger (TCS 1996); Friesen–Theis arXiv:1208.2920.
- Differential posets: Stanley 1988 (JAMS); Fomin 1988; Byrnes 2012 (thesis / sharp upper bound); Stanley–Zanello (EJC 19(2) P13, 2012); arXiv:1202.3006.
- Cited from standard knowledge, not re-searched this turn: Alon–Milman 1985, Linial–London–Rabinovich 1995, Lubetzky–Stav (field dependence of minrank), Körner graph entropy / Hansel, expansion of random 3-CNF.

## 8. Disclosure

Written by an AI research assistant (Grok Bot) for Alkindus Arifuddin. Nothing was posted. The proofs are short and were checked only by the computations listed; there is no independent human review.

## Integration precision and scope

The MAIS selections are checked for acyclicity, and HiGHS must report an
optimal status. Its floating MILP optimality is not exported as an exact
independently checked upper-bound certificate. Part F computes integer
probability numerators but compares their logarithms numerically. The
separate uniform lazy-walk proof supplies the all-parameter lemma.
The corollary still assumes expanding primal instances, bounded occurrence,
balanced groups and the preprint's measurement accounting. It does not
exclude factoring-specific inputs or algorithms outside that interface.
