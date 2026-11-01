# Kronecker / Erdős–Turán duality for the window lemma (Path C), and its relation to Path B

Date: 2026-11-01 (America/Chicago, CST). Status: research note. Every claim is labelled
**PROVEN** (proof given here, also checked numerically), **KNOWN** (in the literature, cited),
**NUMERICAL** (measured, with the command that reproduces it), or **CONJECTURE**.

## 0. Verdict

The dual statement can be made precise and proved (Theorem 1). It is an Erdős–Turán inequality
on the finite group H = ⟨g_q : q ∈ B_Y⟩, using a minorant built from a ball. But the conjecture
as posed ("a failing window forces a nontrivial χ with χ(g_p) ≈ 1 and χ(g_q) ≈ 1 on B_Y; the
converse is the Path B wall, and GRH rules it out") does not hold up in that form, for three reasons:

1. **The conclusion is vacuous.** Characters with χ(g_q) ≈ 1 for every q ∈ B_Y always exist,
   unconditionally. This follows from pigeonhole and, quantitatively, from the additive-energy bound
   (Theorem 3). GRH therefore cannot rule them out. In every one of the 1168 numerical
   instances there is a nontrivial χ|_H with |σ_{R/2}(χ)| ≥ 0.52.
2. **GRH (Bach) controls only the exact stratum.** That stratum is the set of characters trivial on
   all of B_Y. Those characters are invisible to the window lemma, because S ⊆ H and the walk stays
   in H (Proposition 4). Bach's theorem is the Path B wall: it says H = Cl once Y ≥ 6 log²|Δ|.
   It says nothing about Path C windows inside H.
3. **What Path C actually needs is quantitative.** The requirement is a bound on the *mass* of the
   Fourier transform of S on the dual Bohr set {χ : ‖θ_p(χ)‖ < η} (Corollary 1). The trivial
   (Parseval) bound proves windows of length at most W_* = h_H/(η|B_{R/2}|) ≈ 6–10·√h_H.
   Observed windows are 7.5× (median, |Δ| ≤ 4e9) to 19× (median, |Δ| up to 4e11) shorter than that.
   Closing the gap is a short-orbit equidistribution statement about smooth classes. That is beyond
   GRH-type input, and it lies in the "sparse orbit" regime that the known subconvexity methods do
   not reach.

The two walls *are* the same object seen in two strata of the phase map
Θ : Ĥ → T^m, χ ↦ (θ_q(χ))_{q∈B_Y}. Path B (least nonresidue, generation, Bach) is the kernel
and the real (genus) characters. Path C is a mean-value bound over a neighbourhood of the kernel.
The numerics, however, show that the real characters do not drive Path C windows. They are the
most "pretentious" character in only 1% of instances (0% at large |Δ|), and their pretentiousness
does not correlate with window length after controls. So the unification is a change of language.
It does not reduce one wall to the other.

## 1. Setup (notation of `graphs/pathC-graph.md`)

- Δ < 0, Cl = Cl(Δ), h = |Cl|, H = isqrt(|Δ|/3), B = B_Y = split primes ≤ Y (the repo also adds the
  ramified 2 when N ≡ 1 mod 4, Δ = −4N; for that coordinate only v_2 ∈ {0,1} is allowed, and
  ramified primes are never targets). m = |B|.
- g_q = prime form of q. The subgroup they generate is H_B = ⟨g_q⟩ ≤ Cl, with h_H = |H_B|.
- S = classes whose reduced leading coefficient A is B-smooth. δ = |S|/h_H (density *inside* H_B).
- φ : Z^m → H_B, v ↦ Π g_q^{v_q}. Λ = ker φ. ‖v‖ = Σ |v_q| ln q. B_a = {v : ‖v‖ < a}. R = ln(√|Δ|/2).
- Walk for target p ∈ B: x_j = g_p^{K_p + j}, with K_p = min{k : p^k > H}.
  W = first j with x_j ∈ S (W_ball: first j with x_j ∈ φ(B_R)).
- For χ ∈ Ĥ_B write χ(g_q) = e(θ_q(χ)), with ‖·‖ the distance to Z.
  σ_a(χ) = |B_a|^{-1} Σ_{v∈B_a} χ(φ v) and c_a(χ) = |σ_a(χ)|.
- Lemma 1 (pathC-graph, **PROVEN** there; re-asserted on every instance by the script):
  min_{λ∈Λ∖0} ‖λ‖ ≥ 2R, so φ is injective on B_R, and φ(B_R) ⊆ S.

## 2. Theorems

### Theorem 1 (dual form of a failing window) — **PROVEN**
Let a, b > 0 with a + b ≤ R and a, b ≤ R/2 (a = b = R/2 is the natural choice). If
x_j ∉ φ(B_{a+b}) for all 0 ≤ j < W (in particular, if the window fails, x_j ∉ S), then

    Σ_{χ ∈ Ĥ_B, χ ≠ 1}  c_a(χ) c_b(χ) · min( W, 1/(2‖θ_p(χ)‖) )  ≥  W.

The same inequality holds with Ĉl in place of Ĥ_B if the right side is multiplied by [Cl : H_B];
the script checks this version.

*Proof.* Put f(x) = |B_b|^{-1} #{(v,w) ∈ B_a × B_b : φ(v − w) = x}. Its support lies in
φ(B_{a+b}) ⊆ φ(B_R) ⊆ S. Fix w. Injectivity of φ on B_a (because 2a < 2R ≤ min‖Λ∖0‖) gives at most
one v for each x, so 0 ≤ f ≤ 1_S. With f̂(χ) = h_H^{-1} Σ_x f(x) χ̄(x) we get
f̂(χ) = (|B_a|/h_H) · conj(σ_a(χ)) · σ_b(χ), and σ_b is real because B_b = −B_b.
So f̂(1) = δ_f := |B_a|/h_H and |f̂(χ)| = δ_f c_a c_b. Fourier inversion along the walk gives
0 = Σ_{j<W} f(x_j) = W δ_f + Σ_{χ≠1} f̂(χ) χ(x_0) Σ_{j<W} e(jθ_p(χ)). Now use
|Σ_{j<W} e(jθ)| ≤ min(W, 1/(2‖θ‖)) and divide by δ_f. ∎

### Corollary 1 (major-arc mass) — **PROVEN**
Parseval and injectivity give Σ_{χ∈Ĥ_B} c_a² = h_H/|B_a|. Splitting Theorem 1 at ‖θ_p‖ = η and
applying Cauchy–Schwarz on the minor arcs shows: if W > W_* := h_H/(η √(|B_a||B_b|)) and the window fails,
then

    Σ_{χ≠1, ‖θ_p(χ)‖<η} c_a(χ) c_b(χ) ≥ 1/2.

So a failing window longer than W_* forces at least half a unit of "ball-Fourier mass" onto the
p-major arc of the dual Bohr set. Since |B_{R/2}| is of size about √(h·(stuff)), W_* sits at the
familiar √h-type B1 barrier. NUMERICAL (η = 1/4, a = b = R/2): median W_*/√h_H = 6.3 (|Δ| ≤ 4e9)
and 9.9 (|Δ| ≤ 4e11). W_* exceeds *every* observed window in 993/1032 and 133/136 instances, with
median W_*/max_p W of 7.5 and 18.8. The corollary is true but too weak to explain the observed windows.

### Theorem 2 (large Fourier coefficient ⇒ χ ≈ 1 on g_q) — **PROVEN**
For q ∈ B split and every χ: c_a(χ) · |1 − χ(g_q)| ≤ 2 ε_q(a), where ε_q(a) = 1 − |B_{a−ln q}|/|B_a|.

*Proof.* σ_a(χ)(1 − χ(g_q)) = |B_a|^{-1}(Σ_{B_a} − Σ_{B_a+e_q}) χ∘φ, so its absolute value is at most
|B_a △ (B_a+e_q)|/|B_a|. Since B_{a−ln q} + e_q ⊆ B_a ∩ (B_a + e_q), the symmetric difference has
size at most 2(|B_a| − |B_{a−ln q}|). Injectivity is not needed. ∎

NUMERICAL: the maximum ratio is 1.000 (the bound is attained). However, in the factoring regime
(u = log H / log Y ≥ 2) ε_q(R/2) is close to 1, so the bound is close to vacuous. "Large c_a ⇒ χ(g_q) ≈ 1"
holds only weakly.

### Theorem 3 (low-frequency characters always exist) — **PROVEN**
Σ_{χ≠1} c_a(χ)⁴ ≥ h_H/|B_{2a}| − 1, and therefore
max_{χ≠1} c_a(χ)² ≥ (h_H/|B_{2a}| − 1)/(h_H/|B_a| − 1).

*Proof.* Σ_χ |σ_a|⁴ = h_H |B_a|^{-4} E, where E = #{v_1+v_2 ≡ v_3+v_4 mod Λ} = Σ_x r(x)² and
r(x) = #{(v_1,v_2) ∈ B_a² : φ(v_1+v_2) = x}. The support of r lies in φ(B_{2a}), so Cauchy–Schwarz gives
E ≥ |B_a|⁴/|B_{2a}|. Subtract the χ = 1 term, then use max c² ≥ Σc⁴/Σc². ∎

NUMERICAL: 0 violations out of 1168 instances (median slack ×2.2–3.0). The maximum c_a over χ ≠ 1 has
median 0.93 and minimum 0.52. **Consequence:** "failure ⇒ ∃ nontrivial χ near 1 on all g_q" carries
no information, and no hypothesis on L-functions can refute it.

### Proposition 4 (Path B stratum is invisible to the window lemma) — **PROVEN**
S ⊆ H_B, because a reduced form represents A, so its class is a product of the prime forms of the
primes dividing A. Those primes all lie in B or are ramified. The walk x_j also lies in H_B. Hence
characters trivial on H_B (the "exact" stratum, Bach's object) contribute only the constant
[Cl:H_B] to every identity above. They matter only for Path B coverage (H_B = Cl or not).
*Rank-1 Kronecker identity* (PROVEN, elementary): over a full period,
Σ_{j<ord(g_p)} 1_S(x_j) = ord(g_p) · Σ_{χ∈Ĥ_B : χ(g_p)=1} Ŝ(χ) χ(x_0). A full orbit misses S
iff the coset mass vanishes. NUMERICAL: this happens for the real S in 0 of 4985 targets, against
16 for random S' of the same size.

### Remark 5 (GRH scope) — KNOWN + assessment
Bach (1990, Thm 4 / §7) shows, under ERH, that every nonprincipal class-group character has χ(𝔭) ≠ 1
for some N𝔭 ≤ 12 log²|Δ|, and that prime ideals of norm < 6 log²|Δ| generate Cl in the quadratic case.
This is the *exact* stratum (χ ≡ 1 on B_Y) and is precisely Path B / B1b. The same explicit-formula
method, applied to Σ_{N𝔞≤x} Λ(𝔞)(1 − Re χ(𝔞))(1 − N𝔞/x), should give a quantitative statement:
for x ≥ C log²|Δ|, χ is not uniformly close to 1 on split primes ≤ x. That saves only a polylog factor
per character. Theorem 1 needs a saving summed over about h characters, so GRH-type input alone
cannot cross the W_* barrier. (Assessment, not proved here. Constants not derived.)

### What would close it — CONJECTURE
A *mean-value* (large-sieve) bound over the dual Bohr set:

    Σ_{χ≠1, ‖θ_p(χ)‖<η} c_{R/2}(χ)² ≤ 1/4   for η ≈ 1/W,  W ≈ C/δ,

or, equivalently, equidistribution of the short orbit segment {g_p^{K_p+j}}_{j<C/δ} with respect to
the ball measure on S. In the language of Michel–Venkatesh this is equidistribution of a very sparse
orbit (length |Δ|^{o(1)}). Known methods (Duke; Harcos–Michel; Michel–Venkatesh) handle orbits of
subgroups whose size is a fixed power of |Δ|, and need subconvexity even there. Exceptional (Siegel)
characters are not the issue: they are real characters, and real characters are numerically irrelevant
to Path C windows (§3).

## 3. Numerical verification — NUMERICAL

Exact class groups come from PARI `quadclassunit` (GRH-conditional for the group structure). h is
certified by an exhaustive reduced-form count for |Δ| ≤ 1e7 (411 instances). Factor base and prime
forms come from the repo's `tools/` (`build_prime_form_factor_base`, `least_power_above`,
`reduce_form`). Instances are random squarefree N with Y = H^{1/u}, u ∈ {2, 2.5, 3, 4}, plus small P3
adversary instances (p ∈ {3,5,7}). Control: a random S' ⊆ H_B of the same size, analysed identically.

| | small (`output_scan.txt`) | large (`output_scan_large.txt`) |
|---|---|---|
| instances (Δ,Y) / targets | 1032 / 4093 | 136 / 892 |
| \|Δ\| | 1.06e4 – 4.40e9 | 1.13e9 – 3.96e11 |
| h | 44 – 59928 | 6108 – 396256 |
| Thm 1 min LHS/RHS (must be ≥ 1) | 1.881 | 5.251 |
| Thm 2 max ratio (must be ≤ 1) | 1.000 | 1.000 |
| Thm 3 violations | 0 | 0 |
| Fourier identity exactness | 3e-13 | 6e-13 |
| W·δ real S: median / 90% / max | 1.24 / 4.84 / 27.0 | 1.10 / 4.20 / 18.7 |
| W·δ random S': median / 90% / max | 0.74 / 2.13 / 6.6 | 0.72 / 2.32 / 6.6 |
| mean\|Ŝ\| on lowest-1%-F chars, S / S' | 6.65× | 6.20× |
| R² of E_all by low-F chars, S vs S' | 0.53 vs 0.05 | 0.80 vs 0.00 |
| R² of E_all by p-major arc, S vs S' | 0.58 vs 0.91 | −0.32 vs 0.89 |
| raw Spearman(log Wδ, G_p) S / S' | −0.25 / 0.01 | −0.27 / 0.05 |
| **partial** Spearman (ctrl log δ, m, log p, log h) | +0.08 / −0.01 | +0.03 / +0.08 |
| partial Spearman, pretentious G_p | −0.07 / −0.02 | +0.01 / +0.05 |
| most pretentious χ is real (genus) | 1.0% | 0% |
| partial Spearman(max log Wδ, real-χ pretentiousness) | +0.07 | +0.16 |

Score definitions: F(χ) = Σ_q |1 − χ(g_q)|²/ln² q (ball-adapted).
The pretentious distance is D² = Σ_q (1 − cos 2πθ_q)/q.
G_p = min F over χ ≠ 1 on H_B with ‖θ_p(χ)‖ ≤ δ.
E_all = count_{W0}/(δ W0) − 1 at W0 = ⌈1/δ⌉, and E_L is the exact contribution of character set L.

Findings:
1. **S really is Fourier-concentrated on low-frequency characters** (6–7× random). Those characters
   carry the deviations of S along the walk (R² 0.53 → 0.80 as |Δ| grows, against about 0 for random S').
   This is the correct, positive part of the duality picture: S is ball-shaped, so Ŝ lives on the
   dual Bohr set.
2. **The per-target existence test is vacuous.** "∃ low-F χ in the p-arc" holds for 99.9–100% of
   targets for both S and the random S'. The raw correlation between window length and G_p (−0.25 to
   −0.40) is entirely confounding through δ and m (Spearman(G_p, δ) = 0.65, (G_p, m) = 0.46–0.81).
   After controlling, |partial ρ| ≤ 0.08 for small instances and ≤ 0.03 for large ones.
3. **The coarse p-major arc does not carry the deficit.** In zero-hit windows E_major has median −0.65
   (small instances) but **+0.22** (large instances): the major arc predicts *more* than average,
   and the miss comes from cancellation spread over many characters. Counterexamples to the simple
   form: D = −553863727, h = 12416, Y = 11, p = 2, W = 3078 (Wδ = 27, random-S' window 143),
   E_major = +5.57; D = −3166093599, h = 34060, Y = 13, p = 2, W = 3378, E_major = +7.75.
   The longest windows all have m ≤ 3 and p ∈ {2,3}. That is a two- or three-dimensional Kronecker
   effect of the ball shape, not one bad character.
4. **The Path B stratum does not drive Path C.** Real characters are almost never the most pretentious.
   Their pretentiousness has partial correlation +0.07 / +0.16 with window length (the "wrong" sign,
   and small). [Cl:H_B] > 1 instances have longer windows (median Wδ 4.4 against 1.5), but this is
   mostly the small-m effect (at m = 4: 1.52 against 1.39).

## 4. Literature searched (2026-11-01)

The web searches run are listed with what each one found. **No source found states Theorem 1 for this
walk.** Its ingredients are standard, and the theorems should be presented as routine applications
rather than new mathematics.

- "Erdős–Turán inequality finite abelian group discrepancy characters Bohr set": Erdős–Turán (1948),
  Koksma's multidimensional form, and the Erdős–Turán–Koksma inequality over finite groups (Math.
  Balkanica vol. 19). Bohr sets B(K,ρ) and the bound |B| ≥ ρ^{|K|}|G| are standard (e.g. arXiv
  1805.06847). Theorems 1–3 are a direct instance of this toolkit. **Known technique.**
- "Granville Soundararajan pretentious characters smooth numbers …": Granville–Soundararajan, *Large
  character sums* (JAMS 2001) and *Large character sums: pretentious characters and the
  Pólya–Vinogradov theorem* (JAMS 2007). Large character sums come from χ(p) misbehaving at small
  p, i.e. from sums over smooth numbers. Halász / Montgomery–Tenenbaum on mean values; *Pretentious
  multiplicative functions* (distance D). This is the right language for the "near-1 on small
  primes" stratum. **Known framework, not applied to WL.**
- "Ellenberg Michel Venkatesh Linnik's ergodic method …": EMV arXiv 1001.0897 (expander random walks,
  large deviations); Michel–Venkatesh ICM 2006, which states that sparse subgroup orbits
  equidistribute only above a power of |Δ| and that arbitrarily small subgroups are out of reach;
  Duke 1988.
- "Heegner points short orbits … Harcos Michel subconvexity Duke 1988": Harcos–Michel, Invent. Math.
  2006 (short Galois orbits via Rankin–Selberg subconvexity). Confirms that the CONJECTURE in §2
  sits in the open sparse-equidistribution regime.
- "Thorner Zaman least prime ideal …": Thorner–Zaman, ANT 11 (2017), arXiv 1604.01750: explicit
  Linnik/Weiss bound for the least prime in a Chebotarev class. Gives the least prime represented by
  a form ≪ |Δ|^694. Lagarias–Montgomery–Odlyzko 1979. Unconditional least-prime-in-a-class bounds
  are polynomial in |Δ|, which is useless for smooth-scale Y.
- "Bach 1990 explicit bounds …": Bach, Math. Comp. 55 (1990) 355–380, Thm 4 and §7: 12 log² and,
  for quadratic fields, 6 log² generation. This is the exact stratum only (Remark 5).
- "smooth ideals class group … Hafner McCurley rigorous": Hafner–McCurley JAMS 1989. The ERH-based
  rigorous analysis of the class-group algorithm uses Bach-type generation plus smooth-ideal
  distribution; IACR ePrint 2025/1514 (rigorous sampling). None addresses single-generator orbits.
- "Harper character sums over smooth numbers …": Harper arXiv 1208.5992 (Bombieri–Vinogradov for
  smooth numbers; uses the large sieve plus zero density) and arXiv 1103.2106 (smooth numbers in APs).
  The closest analogue to the needed mean-value input, but over Dirichlet characters, not class-group
  orbits.
- "distribution of powers of a fixed ideal class reduced forms leading coefficient smooth …":
  nothing relevant found.
- Not re-searched this session: Jao–Miller–Venkatesan arXiv 0811.0647 (expander graphs on isogeny /
  class-group Cayley graphs under GRH; their mixing is for random walks with many generators, not
  for a single-generator progression). Recorded from the earlier Path C note.

## 5. Placement recommendation

Keep this in the integer-factoring research repo, as `research/` notes next to the window lemma.
Do not start a separate repo. Reasons:
- Every object (S, B_Y, K_p, the walk) is defined by the factoring algorithm. The geometric question is
  the window lemma rewritten in dual coordinates, and it depends on Lemma 1, the repo's factor-base
  conventions, and the P3 adversaries.
- The general mathematical question (sparse-orbit equidistribution of smooth classes, a mean value of
  ball-Fourier mass on dual Bohr sets) is open and hard. A standalone repo would make sense only if
  someone decides to pursue that question for its own sake, independent of factoring. If that happens,
  spin it out then, with this note as its seed.

## 6. Files and rerun

    duality/clgroup.py         exact class group via PARI; enumeration; dlogs; repo factor base
    duality/duality_scan.py    per-instance theorem checks, Fourier data, windows, controls → jsonl
    duality/analyze_scan.py    summaries, partial correlations, genus stratum → output_*.txt
    duality/run_all.sh         reproduces everything (~20 s; PY=/path/to/python3 with cypari2+numpy)
    duality/output_scan.txt, output_scan_large.txt, SHA256SUMS (outputs are byte-reproducible)
    duality/data/scan*.jsonl, scan*.log

    cd duality && PY=python3 ./run_all.sh     # uses ../../integer-factoring/tools or $IF_TOOLS

Determinism: PARI's random generators are fixed with `setrand(1)`, the random control is drawn in
canonical form order, and the analysis sorts records and uses tie-averaged ranks. Two consecutive runs
produce identical `output_*.txt`.

## AI-assistance note
The theorem statements, proofs, scripts, and this text were drafted by an AI research assistant
(Grok Bot) and checked numerically on the box. The proofs are short and should still be read by a
human before being cited. Literature claims are based on the web searches listed in §4, not on a
complete survey.
