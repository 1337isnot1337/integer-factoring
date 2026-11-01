# Path A: the √log rate of the generic counting engine. A routing-rank no-go and the crossing lane

Base commit: `1337isnot1337/integer-factoring@4c3d135`. `origin` has since moved to `5b20436`; those two commits touch box/word-adjugate budgets only and do not affect this note.
OpenAI source: `main_full_version.tex` ("The Exponential Time Hypothesis is False", revised 25 Oct 2026), source: https://raw.githubusercontent.com/openai/math/main/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/main_full_version.tex (line numbers refer to that file as fetched on 2026-11-01).

## Scope

The question is whether α > 1/2 on the Big O Speedrun scale (time 2^{b/(log b)^α}) can be reached without factoring-specific structure, by beating the √log rate of the generic counting engine. The note attacks that question at three levels:

1. **Complexity level (main result).** In the OpenAI engine, the √log comes from Step 2, the auxiliary matrix M_0. Theorem A proves a lower bound on rank(I−M_0) for *every* matrix with the Step-2 support pattern. Within the paper's own design accounting, this bound forces the exponent to be ≥ n·Ω(√(log log n / log n)). So no choice of M_0 on the paper's auxiliary array gives α > 1/2.
2. **Algebraic level (crossing lane, issue #52).** No known universal spectral point certifies that 1/√2 is optimal for the crossing tensor. Known functionals only give ρ ≥ 1/2.
3. **Meta level.** The note states exactly what the no-go does not cover, so the remaining doors are explicit.

Out of scope: lower bounds against arbitrary algorithms. Nothing here is a circuit or algorithm lower bound.

## Setup (from the preprint)

- Points are t ∈ [Q] (l.671). Schedules are a_v, d_e ∈ [r]^Q. The auxiliary array is B = ([r]×[r])^Q, of size D = r^{2Q}, with labels u = (d_u, a_u) (Step 4, l.1265ff).
- Step 3 (l.1236–1262) turns any rational D×D matrix M_0 into replacement rows of degree ≤ R−1, using C with rowspace(C) ⊇ rowspace(I−M_0). This works provided
  (S) M_0[u,w] ≠ 0 only if dist_H(d_u, a_w) ≤ R−1. The multiplier ∏_{s∈E_{u,w}}(1−T/s) has degree |E_{u,w}| = dist_H.
- Accounting (l.1340–1348): c = rank(I−M_0) rows per copy, h copies with hD ≤ 2M ≤ 6n, and rhc ≤ 2τn total measurement rows, where 2τn is the (r,τ,D)-design budget (Def. at l.690, condition 3).
- The cost has three terms: the Taylor/interpolation term ≈ n·R/Q·log(Q/R), the measurement enumeration (2D+1)^{2τn}, and the group sums 2^{2n/r}.
- The paper's own M_0 = S·g(W_0) (S is the full swap (d,a)↦(a,d), W_0 = Σ_t (I−Π_t)/2, and g interpolates (−1)^j on R consecutive integers) has c/D = Pr[J ∉ window], J ~ Bin(Q,(r−1)/(2r)). Requiring this to be ≤ τ/(3r) gives R ≍ √(Q log(r/τ)), which is the source of the √log.

## Claims

### Theorem A (routing-rank lower bound). PROVEN

Let 𝔽 be any field and r, Q, k ≥ 0. Let M be any D×D matrix over 𝔽 on B with M[u,w] ≠ 0 only if dist_H(d_u, a_w) ≤ k. Then for every total order ≺ on [r]^Q,

  rank(I−M) ≥ #{(α,δ) ∈ [r]^Q×[r]^Q : α ≻ max_≺ Ball_k(δ)}.

*Proof.* Let ξ ∈ ker(I−M) and write Ξ[α,δ] = ξ(d=δ, a=α). The equation ξ = Mξ says that Ξ[α,δ] is a linear combination of entries in the rows γ ∈ Ball_k(δ), since w ranges over a_w ∈ Ball_k(δ). Let F = {(α,δ): α ⪯ max Ball_k(δ)}.

Claim: if ξ|_F = 0 then ξ = 0. Induct upward on α. Suppose every row γ ≺ α vanishes. If (α,δ) ∉ F, every γ ∈ Ball_k(δ) satisfies γ ≺ α, so Ξ[α,δ] = 0. If (α,δ) ∈ F, it is zero by assumption.

Hence dim ker ≤ |F|, which gives the bound. ∎

**Corollary A1 (score order).** Take H ⊂ [r] with |H| = ⌊r/2⌋, score s(γ) = #{t: γ_t ∈ H}, and order by score with ties broken arbitrarily. Changing ≤ k positions raises the score by ≤ k, so

  rank(I−M) ≥ D·Pr[S − S' > k], with S, S' iid Bin(Q, ⌊r/2⌋/r). PROVEN

**Corollary A2 (k = 0).** The lex order gives rank(I−M) ≥ r^Q(r^Q−1)/2. This is attained by M = S. PROVEN (VERIFIED for small cases).

**Tightness. PROVEN up to constants, VERIFIED numerically.**
- The construction needs k = R−1 ≈ 2z_{ε/2}·σ(J) for rank/D ≤ ε. The bound needs k ≳ z_ε·σ(S−S').
- As r→∞, σ(J)² → Q/4 and σ(S−S')² → Q/2, so the construction's degree is optimal up to a factor → √2. At r = 2 the factor is √(3/2).
- Exact table (`check_routing_rank.out`): the ratio is 1.19–1.49 for Q ≤ 1024 and ε ∈ {10⁻³, 10⁻⁶}.

### Lemma B (binomial anticoncentration). PROVEN; VERIFIED on a grid

Let q = ⌊r/2⌋/r ∈ [1/3,1/2], Q ≥ 64 and 1 ≤ R ≤ Q/4. Then

  Pr[S − S' ≥ R] ≥ ½·exp(−9R²/(2Q) − 72).

*Proof.*
1. Write S − S' + Q = Bin(Q,q) + Bin(Q,1−q), with independent summands. Then lower-bound Pr[Bin(Q,p) ≥ ⌈Qp+R/2⌉] for p ∈ {q, 1−q} and multiply.
2. Sum Pr[=k] over ⌈√Q⌉ consecutive k starting at A = ⌈Qp+R/2⌉. Each term is ≥ e^{−Q·KL(k/Q‖p)}/√(2Q), using binom(n,k) ≥ 2^{nH(k/n)}/√(8n·x(1−x)) (Cover–Thomas, *Elements of Information Theory* 2nd ed., Lemma 17.5.1; KNOWN).
3. Bound KL ≤ χ² = (x−p)²/(p(1−p)) (KNOWN), with p(1−p) ≥ 2/9.
4. The deviation is ≤ (R/2+1+√Q)/Q, and (R/2+1+√Q)² ≤ R²/2 + 8Q.
5. The index range stays ≤ Q because Q ≥ 64 and R ≤ Q/4.

Each factor is therefore ≥ 2^{−1/2}e^{−9R²/(4Q)−36}. The constants are not optimized. ∎

### Corollary C (no α > 1/2 within the paper's Step-2/Step-5 accounting). PROVEN, under the stated accounting

Assumptions:
- Any design built as in Steps 2–4: the auxiliary array B, any M_0 satisfying (S) with degree R−1, and measurement rows counted as rhc with c ≥ rank(I−M_0), h ≥ 1, hD ≥ M ≥ n, rhc ≤ 2τn.
- A normalized exponent δ with R/Q ≤ δ, τ·log₂D ≤ δ, 1/r ≤ δ. All three terms are ≤ (paper's exponent)/n up to constants.

Conclusion: δ ≥ Ω(√(log log n / log n)). So the time is 2^{n·Ω(√(log log n/log n))}, and the bound is consistent with the paper's n·log log n/√log n.

*Proof.* If Q < 64 then δ ≥ 1/64, so assume Q ≥ 64 and δ ≤ 1/4 (hence R ≤ Q/4). Set ε := c/D ≥ Pr[S−S' ≥ R] by Theorem A/A1, which is ≥ ½e^{−9R²/(2Q)−72} by Lemma B. Write τ ≤ δ/(2Q log₂ r).
- (i) From h ≥ 1: ε ≤ 2τn/(rD) ≤ n·r^{−2Q}. Then 9R²/(2Q) ≥ 2Q ln r − ln(2n) − 72. With R ≤ δQ this gives Q(2 ln 2 − 4.5δ²) ≤ ln(2n)+72, so Q ≤ 0.91(ln n + 73).
- (ii) From hD ≥ n: ε ≤ 2τ/r ≤ δ/2. Then R² ≥ (2Q/9)(ln(1/δ) − 72), so δ ≥ R/Q ≥ √(2(ln(1/δ)−72)/(9Q)).

Combining, δ² ≥ 2(ln(1/δ)−72)/(8.2(ln n+73)). This forces δ = Ω(√(log log n/log n)). ∎

Branch (i) also shows that a sparse or huge auxiliary array (D ≫ n) cannot help within this accounting. The k = R−1 support then forces R ≳ Q.

**What Corollary C does NOT cover (open doors). CONJECTURE / open:**
- (a) A measurement component that is not a single matrix M_0. Step 3 uses the decomposition ζ = (I−M_0)ζ + M_0ζ at every t. A design whose rowspace part varies with t non-polynomially is not covered.
  - Farkas reformulation: failure for label u ⇔ there exist K-valued κ(0..Q) with vanishing moments of order < R, κ(t) supported on {w: a_w(t) ≠ d_u(t)}, and κ(0)_u ≠ 0.
  - Single-subset necessary conditions are vacuous once R > r.
- (b) Auxiliary arrays other than ([r]×[r])^Q, or labels that are not schedule pairs.
- (c) Substitution effects. Merging coordinates when original variables are substituted can lower the true measurement count below rank(I−M_0). The corollary bounds the paper's accounting, not the true count.
- (d) Entirely different engines.

### Crossing lane, algebraic level

- **PROVEN (VERIFIED m = 1, 2):** the crossing tensor E(x,y,u,v) = [x=v][y=u], N = 2^m, is free and all one-party marginals are maximally mixed (`check_spectral_points.py`). Hence every CVZ quantum functional gives F^θ(E) = N. Strassen support functionals are spectral points on oblique tensors. As used here they also give N, because the support's coordinate marginals are uniform. Note this support-functional statement is not separately checked.
- **Elementary bound, PROVEN:** Σ_j √(r1_j r2_j) ≥ Σ√max ≥ √(Σ max) ≥ N, i.e. ρ ≥ 1/2. So no known universal spectral point certifies the 1/√2 barrier. That barrier is proven in the repo only for partition channels (≥ 1/√2) and site-symmetric binary pieces (≥ 3/4). Whether 1/√2 is optimal for general decompositions is OPEN.
- **Observation (PROVEN from repo eq. (4) of `graded-space-barrier.md`):** the dimension-only graded bound is ρ ≥ 2^{−γ} with energy range L = 2^{1/γ}. The repo excludes γ < 1/2. For γ > 1/2 and L < 4 it permits ρ < 1/√2 (e.g. γ = 1, L = 2 gives 1/2). Whether any algebra realizes such a grading compatible with the crossing is OPEN. The repo note itself says it "does not cover arbitrary intermediate spaces".
- The site-symmetric expansion SWAP = ⊗(I − 2P_singlet) grouped by singlet count gives ρ ≈ 0.896 at κ = 1/2. This is consistent with the repo's 3/4 barrier and gives nothing new.

### Analogy (heuristic, not a claim)

The Anthropic engine's α = 1 − γ comes from partition-type √E growth. The OpenAI engine's 1/2 comes from CLT width √Q of a score/parity statistic. Theorem A shows the second is intrinsic to Hamming-support routing on product arrays.

## Literature searched (2026-11-01)

- Christandl–Vrana–Zuiddam, "Universal points in the asymptotic spectrum of tensors", JAMS 36 (2023), arXiv:1709.07851. KNOWN: quantum functionals; Strassen support functionals for oblique tensors.
- Alman–Williams (FOCS'15 probabilistic polynomials; STOC'17 rigidity, interpolating (−1)^j on an interval). This is the construction's antecedent, cited by the preprint (l.441–453). I found no routing-rank lower bound there.
- Web queries for √log barrier / optimality of the measurement rank in the ETH preprint. Found:
  - the preprint itself, which presents its bound as algorithmic, not optimal;
  - O'Donnell's exposition (CMU), Kintali's simplified proof, Kominers' Hexagon exposition 2610.00574v1, the Papex verification blog, and the OpenAI Lean scope doc 378.md.
  - None of these states a lower bound on rank(I−M_0).
- Sign-rank / k-Hamming-distance literature (arXiv:2506.12022): related support-pattern rank questions with a different measure. Not this result.
- Repo greps at 4c3d135 for `r^{2Q}`, `Bin(Q`, measurement tuple, replacement row, `rank(I`, `I-M_0`, routing, quantum functional, asymptotic spectrum: no prior analysis of the preprint's partial-swap routing.

"Not found" ≠ new. Theorem A is an elementary free-set argument and may be folklore. Approximate-degree lower bounds for symmetric functions (Paturi; Nisan–Szegedy) are the natural one-dimensional analogue and were not searched exhaustively.

## Checkers (exact arithmetic: integers, Fractions, GF(p), p = 1000003)

| command | checks | result |
|---|---|---|
| `python3 tools/routing_rank_nogo/check_routing_rank.py` (~4 min) | lex bound = r^Q(r^Q−1)/2 (6 cases); exhaustive min rank over GF(2), GF(3) for r=2, Q=1, k=0 equals 1; construction S·g(W_0) for 9 (r,Q,I): zero support violations, GF(p) rank = eigen formula, score/lex bounds ≤ rank, D·Pr bound ≤ score bound; degree-ratio table | `ALL CHECKS PASS` (`check_routing_rank.out`) |
| `python3 tools/routing_rank_nogo/check_tail_lemma.py` (~1 min) | Lemma B and its per-factor bound, exactly, for r ∈ {2,3,4,5,16}, Q ∈ {64,100,256}, all R ≤ Q/4 | `ALL CHECKS PASS`, min log-margin 71.95 (`check_tail_lemma.out`) |
| `python3 tools/routing_rank_nogo/check_spectral_points.py` | freeness, maximally mixed marginals, cut ranks N² (m = 1, 2) | `ALL CHECKS PASS` (`check_spectral_points.out`) |

The GF(p) ranks are lower bounds on the ℚ-rank. Equality with the ℚ-rank follows from the eigen formula (proven), and the checker confirms they match.

## Verdict

- The generic OpenAI engine is **walled at α = 1/2** for every Step-2 matrix on the schedule-pair array (PROVEN within the accounting). The paper's construction is optimal in degree up to a factor ≤ √2.
- The path stays alive only via doors (a)–(d): t-dependent measurement components, non-product auxiliary arrays, or substitution-merging.
- On the crossing lane there is no functional obstruction below 1/√2, and the γ > 1/2 graded window is open but unrealized.

## AI-assistance disclosure

This note, the proofs, and the checkers were produced by an AI research agent (Grok Bot) and have not been reviewed by a human mathematician. Every finite claim was machine-checked by the scripts above. The asymptotic proofs are by hand and need independent review.
