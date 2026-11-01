# Path C, graph view: does a known graph theorem settle WL(Y, W)?

Box-only research note, follow-up to the window-lemma note in PR #69
(https://github.com/1337isnot1337/integer-factoring/pull/69). Repo checkout
`ab5f29b01e3a3dba1bbc3e1ebb74f49091b195d5`. The repo tools used here
(`binary_quadratic_forms.reduce_form`, `prime_form_factor_base`,
`weighted_norm_relations.least_power_above`) are byte-identical to `4c3d135`, the
commit Path C was measured on. Written 2026-11-01 (America/Chicago). Drafted with AI
assistance. Nothing was posted.

Labels: **PROVEN** = proved here, with the proof given. **VERIFIED** = checked
exactly by the scripts in this directory. **KNOWN** = in the literature, with a
reference. **CONJECTURE** = believed but not proved.

## 0. Short answer

Rewriting WL as a graph problem does not give a new unconditional bound. Every angle
either (a) reduces to the existing barriers B1/B3 of PR #69, or (b) gives only
randomized statements, whose deterministic versions are blocked by a named open
problem (circuit lower bounds for E). The rewrite does produce three facts that are
small but useful, and that PR #69 does not yet contain:

1. **Ball sandwich (PROVEN).** S lies between the injective image of a weighted-ℓ¹
   ball of radius ln(√|Δ|/2) and the image of the ball of radius ln H. The two radii
   differ by only ln(2/√3) ≈ 0.144. So S is "a ball in the relation-lattice torus"
   Z^m/Λ, not a random-looking set. This one fact explains the small-p clustering and
   the Fourier bias below, and it shows why the flatness and covering theorems are
   vacuous.
2. **Fourier bias (VERIFIED).** At every size tested, max_{χ≠1}|Ŝ(χ)| is 2.4–24×
   (median 10.8×) the value for a random set of the same size. The biased character
   has small phases on the base generators. So the uniform character-sum bound
   M ≈ √|S| assumed in B1 is false for S at these sizes. Any completion argument
   needs a major/minor-arc split.
3. **Pigeonhole (PROVEN).** If the window has to work for every start on the g_p
   cycle, then W ≥ 1/δ. A density-only argument that is uniform over starts cannot
   beat 1/δ. Separately, uniform-over-starts WL fails outright (VERIFIED, 6
   generators): some cosets of ⟨g_p⟩ contain no point of S.

## 1. The graph and the lattice

Notation as in PR #69: Δ, H = isqrt(|Δ|/3), base B_Y (split p ≤ Y, m = |B_Y|), S =
the classes whose reduced leading coefficient A is B_Y-smooth, δ = |S|/h, and
K_p = min{k : p^k > H}.

- **Cayley graph.** G_Y = Cay(Cl(Δ), {g_p^{±1} : p ∈ B_Y}), degree d = 2m. WL(Y, W)
  for p asks that the walk x_j = g_p^{K_p+j} (a geodesic of the single-generator
  cycle graph) enters S within W steps.
- **Lattice (KNOWN viewpoint).** This is the relation lattice of Hafner–McCurley and
  Buchmann. Let φ: Z^m → Cl, v ↦ ∏ g_p^{v_p}. Then Λ = ker φ has index h, and
  G_Y = Cay(Z^m/Λ, ±e_p) is the quotient of the grid Z^m by Λ. Write
  ‖v‖ = Σ|v_p| ln p and R = ln(√|Δ|/2).
- **Lemma 1 (PROVEN).**
  - (i) Every nonzero λ ∈ Λ has ‖λ‖ ≥ ln(|Δ|/4) = 2R.
  - (ii) φ is injective on B_R = {‖v‖ < R}, and φ(B_R) ⊆ S ⊆ φ(B̄_{ln H}).

  *Proof.* Let 𝔞_v = ∏_{v_p>0} 𝔭^{v_p} ∏_{v_p<0} 𝔭̄^{|v_p|}. It is primitive with
  N𝔞_v = e^{‖v‖}.
  - (i) If φ(v) = φ(w), then 𝔞_v·𝔞̄_w = (α) with N(α) = e^{‖v‖+‖w‖}. If α ∉ Z, then
    N(α) = (x² + |Δ|y²)/4 ≥ |Δ|/4. If α = n ∈ Z, comparing the 𝔭- and 𝔭̄-adic
    valuations gives v_p⁺ + w_p⁻ = v_p⁻ + w_p⁺ for every p, i.e. v = w. Taking w = 0
    gives (i). Taking ‖v‖, ‖w‖ < R gives injectivity.
  - (ii) For ‖v‖ < R the form of 𝔞_v has a < √|Δ|/2. After reducing b to |b| ≤ a we
    get c ≥ |Δ|/(4a) > a, so the form is reduced with A = N𝔞_v, which is smooth.
    Conversely, a smooth reduced A ≤ H factors as a norm of some 𝔞_v with
    ‖v‖ = ln A ≤ ln H. ∎
- **Reformulation (PROVEN, up to the 0.144 slack in the radius).** WL(Y, W) holds for
  p iff the tube T_W = {(K_p + s)e_p : 0 ≤ s < W} + B_r contains a point of Λ, for
  some r ∈ [R, ln H]. In words: does a rational line in the torus Z^m/Λ come within
  the packing radius of the lattice within W steps.

## 2. Verdict per angle

| Angle | Known theorem tried | Verdict |
|---|---|---|
| Single-generator cycle | Cycle C_n has gap 2π²/n² (KNOWN); Szemerédi | **No expansion; Szemerédi has the wrong quantifier.** S contains long APs along e_p trivially (chords of the ball, Lemma 1). WL needs syndeticity of S along a fixed start, which no density theorem gives. Pigeonhole gives W ≥ 1/δ for any start-uniform argument (PROVEN). |
| Multi-generator expander mixing (EML) | EML (KNOWN, Alon–Chung); JMV eigenvalues λ_χ = Σ_p 2Re χ(g_p) (KNOWN, arXiv:0811.0647) | **Reduces to B1.** If the BFS ball B_r is at distance ≥ 2 from S, then \|B_r\| ≤ (λ/d)²h/δ (PROVEN from the EML). Trace/Parseval gives (λ/d)² ≥ (h−d)/(d(h−1)) (PROVEN; VERIFIED 47/47). So the certified ball is ≳ h/(dδ), and with d ≤ 2Y and Y optimized this is the √(h/δ) of B1. Under GRH (JMV, Y ≥ (ln\|Δ\|)^{2+ε}) the bound is ≈ h·polylog/(Yδ), which is still exponential. A fixed deterministic path can be dodged by a density-δ set whatever the spectrum, so a spectrum cannot certify a path. |
| Random walk / hitting time | Aldous–Fill: E_π T_S ≤ t_rel(1−δ)/δ (KNOWN, Ch. 3); AKS expander-walk hitting | **Randomized only.** It gives ℓ ≈ ln√h/(δ·gap) steps for a random walk. That is a randomized relation finder, and Lenstra–Pomerance (JAMS 1992, doi:10.1090/S0894-0347-1992-1137100-0) is already unconditional randomized L[1/2, 1+o(1)]. VERIFIED: the AF bound holds in 47/47 instances but is 2–25× loose. |
| Derandomizing the walk | Impagliazzo–Wigderson STOC'97 (KNOWN, doi:10.1145/258533.258590) | **Blocked by a named open problem.** *Conditional consequence, PROVEN as a standard composition (I found no paper stating it for factoring):* if E needs circuits of size 2^{Ω(n)}, the IW generator with seed O(log T) fools the size-poly(T) circuit "Lenstra–Pomerance on input N with coins r". Trying all seeds and checking each output gives deterministic factoring in time T^{O(1)} = L_N[1/2, O(1)] for every N, which beats the repo's target exp(b/(log b)^α). So the repo goal follows from a standard hardness hypothesis. Proving it unconditionally by this route needs a circuit lower bound. |
| Abelian Cayley graphs, generator count | Alon–Roichman / Lubotzky–Weiss: abelian groups need Ω(log h) generators to expand (KNOWN); Diaconis/Saloff-Coste mixing ≳ k·h^{2/k} (KNOWN) | **Consistent, nothing new.** At the sizes tested (Y < (ln\|Δ\|)²) the graphs are not good expanders: (λ/d)² is 0.35–1.0, and the EML bound is vacuous (≥ h) in 47/47 (VERIFIED). Yet a BFS from the repo start finds S after exploring ≤ 23 vertices (VERIFIED). The structure of S, not the spectrum, is doing the work. |
| Bohr sets / Bogolyubov | Sanders, Bogolyubov–Ruzsa (KNOWN, arXiv:1011.0107) | **Reduces to B3 and is quantitatively worse.** Bogolyubov puts a Bohr set inside 2S−2S, which gives g_p^k ∈ S⁴ within exp(polylog(1/δ)·log(1/ε)) steps. But membership in S⁴ cannot be tested from the leading coefficient, and the step count is worse than trivial for subexponential δ. The new point (VERIFIED, Lemma 1): S itself correlates with a Bohr set (the biased character has small phases on all g_p). This hurts B1 (M ≫ √\|S\|), but it suggests a major/minor-arc split along the orbit. That idea is CONJECTURE only, not carried out. |
| Freiman / doubling | Freiman–Ruzsa (KNOWN) | **No structure to use.** \|SS\| is close to h, so the doubling is ~1/δ, which is useless. Lemma 1 gives the true structure: a ball image, which is a "GAP-like" set of rank m. |
| Flatness / covering minima on the tube T_W | Khinchine flatness; Kannan–Lovász covering minima (KNOWN) | **Vacuous (PROVEN).** T_W has width 2r/ln q along the integral dual vectors e_q, q ≠ p, and 2r/ln q ≤ Flt(m) once m ≳ ln\|Δ\|. A lattice-free T_W is therefore always "flat" and the theorem gives nothing. Covering-minima bounds need the body to be long in m directions, but T_W is long in one. |
| Small-prime clustering | — | **Explained (PROVEN, VERIFIED).** Corollary of Lemma 1, the cusp lemma: if s ∈ S and 2pA(s) < √\|Δ\|, then g_p^{±1}s ∈ S. There were 0 violations, and P(g_p s ∈ S \| s ∈ S) ≥ the cusp fraction in 47/47. The hits are chords of a ball along e_p, of length ≈ 2(r − dist)/ln p, which are longer for small p. Measured: P(next \| hit)/δ ≈ 7 for p ≤ 3 and 1.5 for p ≥ 30. The uniform window is 3.8× the random-model max gap for p = 2 and 1.1× for p ≥ 30. The mean forward hit ≈ 0.76·(cluster size)/δ (VERIFIED). So the effective event rate is δ/c_p. Short cycles: by Lemma 1(i) the only cycles shorter than weighted length 2R are commutator 4-cycles. The graph is locally the grid Z^m (PROVEN). |
| Multi-generator gap needed vs known | JMV under GRH (KNOWN) | **Reduces to B1.** Even optimal Ramanujan-type λ ≈ 2√(d−1) gives a deterministic BFS ball ≈ 4h/(dδ). Beating B1 needs d ≫ h^{1/2}, i.e. Y beyond the base the method can afford. Deterministic progress needs control of S along specific paths (B3), not a better gap. |

## 3. Experiment (VERIFIED, exact)

`cayley_small.py` uses 16 seeded semiprimes N ≡ 3 mod 4 of 30, 34, 38 and 42 bits
with cyclic Cl(−N), h up to 9.4·10⁵, and three values of Y each (u ≈ 2.5, 3, 3.5).
That gives 47 (N, Y) instances. Classes are enumerated exactly as powers of a PARI
generator. The base, the prime forms and K_p come from the repo tools. Summary in
`data/cayley_summary.md`, raw data in `data/cayley.jsonl`.

- Cusp lemma: 0 violations. P(next ∈ S | S) ≥ the cusp fraction in 47/47.
- Pigeonhole bound W_unif ≥ 1/δ: holds in 47/47. 6 generators (4 instances) have a
  ⟨g_p⟩-coset with no S point, so WL that is uniform over starts is false.
- Fourier: max|Ŝ|/max|R̂| has median 10.8 (range 2.4–24.3).
- Spectrum: (λ/d)² ranges over 0.35–1.0. The EML ball bound is ≥ h in 47/47
  (vacuous). The actual BFS ball is ≤ 23 and the covering radius is 2–39.
- Aldous–Fill bound holds in 47/47 and is 2–25× loose.

Rerun (about 5 min on 1 core, needs cypari2, numpy and scipy):

    cd tools/window_lemma/graph
    python cayley_small.py --bits 30 34 38 42 --per-size 4 --u 2.5 3 3.5 --out data/cayley.jsonl
    python summarize.py data/cayley.jsonl data/cayley_summary.md
    (cd data && sha256sum -c SHA256SUMS)

Limits:
- Only cyclic class groups, and |Δ| ≤ 2^42, so Y is far below the JMV regime.
- The spectral statements at cryptographic size are KNOWN-under-GRH, not measured.

## 4. What this adds to PR #69

- Lemma 1 (the ball sandwich, PROVEN) and its corollaries: the cusp lemma, clustering
  as chords, and local grid structure.
- The Fourier-bias measurement, which weakens B1's hypothesis.
- The pigeonhole 1/δ floor for start-uniform arguments.
- One table recording that the EML, hitting-time, Bogolyubov and flatness routes
  reduce to B1/B3 or are randomized-only.
- The conditional remark: IW hardness ⇒ deterministic L[1/2, O(1)]. This is labelled
  as a standard composition, not a new result.

None of it changes the Path C verdict: no proof, conjecture still alive.
