# Kronecker / Erdős–Turán duality for the window lemma (Path C), and its relation to Path B

Date: 2026-11-01 (America/Chicago, CST). Status: research note. Every claim is labelled
**PROVEN** (proof given here, also checked numerically), **KNOWN** (in the literature, cited),
**NUMERICAL** (measured, with the command that reproduces it), or **CONJECTURE**.

## 0. Verdict

The dual statement can be made precise and proved (Theorem 1). It is an Erdős–Turán inequality
on the finite group H = ⟨g_q : q ∈ B_Y⟩, using a minorant built from a ball. But the conjecture
as posed ("a failing window forces a nontrivial χ with χ(g_p) ≈ 1 and χ(g_q) ≈ 1 on B_Y; the
converse is the Path B wall, and GRH rules it out") does not hold up in that form, for three reasons:

1. **Existence is not the needed quantitative statement.** The fourth-moment bound in
   Theorem 3 concerns large ball Fourier coefficients. It does not alone give simultaneous
   closeness to one on every base prime, especially when the boundary estimate in Theorem 2
   is weak. A precise pigeonhole condition is given below. In all 1168 retained numerical
   instances, a nontrivial character has |σ_{R/2}|≥0.52; this is an observation about the
   coefficient, not a universal simultaneous-phase theorem.
2. **GRH (Bach) controls only the exact stratum.** That stratum is the set of characters trivial on
   all of B_Y. Those characters are invisible to the window lemma, because S ⊆ H and the walk stays
   in H (Proposition 4). Bach's conditional small-prime-ideal generation theorem controls
   this exact stratum when the base contains those generators. It does not itself supply
   the short-window mass estimate inside H.
3. **What Path C actually needs is quantitative.** The requirement is a bound on the *mass* of the
   Fourier transform of S on the dual Bohr set {χ : ‖θ_p(χ)‖ < η} (Corollary 1). The
   Parseval minor-arc estimate becomes smaller than W/2 beyond
   W_* = h_H/(η|B_{R/2}|) ≈ 6–10·√h_H. This threshold alone does not prove a window upper bound.
   Observed windows are 7.5× (median, |Δ| ≤ 4e9) to 19× (median, |Δ| up to 4e11) shorter than that.
   Closing the gap requires a correctly normalized short-orbit estimate, with minor arcs
   controlled as well. The estimates stated here and Bach's generation theorem do not
   provide that input. This is not an impossibility theorem for every use of GRH.

The two questions can be expressed through the phase map
Θ : Ĥ → T^m, χ ↦ (θ_q(χ))_{q∈B_Y}, but require different estimates.
Exact generation excludes characters trivial on the whole base; it does not automatically
control approximate phases or the mass near the target generator's annihilator.
The numerics, however, show that the real characters do not drive Path C windows. They are the
most "pretentious" character in only 1% of instances (0% at large |Δ|), and their pretentiousness
does not correlate with window length after controls. So the unification is a change of language.
It does not reduce one wall to the other.

## 1. Setup (notation of [the graph and lattice note](leading-power-window-graphs.md))

- Δ < 0, Cl = Cl(Δ), h = |Cl|, H = isqrt(|Δ|/3), B = B_Y = split primes ≤ Y (the repo also adds the
  ramified 2 when N ≡ 1 mod 4, Δ = −4N; for that coordinate only v_2 ∈ {0,1} is allowed, and
  ramified primes are never targets). m = |B|.
- g_q = prime form of q. The subgroup they generate is H_B = ⟨g_q⟩ ≤ Cl, with h_H = |H_B|.
- S = classes whose reduced leading coefficient A is B-smooth. δ = |S|/h_H (density *inside* H_B).
- Use Γ=Z^s for the unramified split coordinates. If the ramified prime 2 is included,
  use Γ=Z^s⊕Z/2, with its canonical exponent e_2 in {0,1}. Let
  φ:Γ→H_B, v↦Π g_q^(v_q), Λ=ker φ, and
  ‖v‖=Σ_(split q)|v_q|ln q+e_2 ln2. The balls B_a={‖v‖<a}
  are symmetric in this group and their norm is subadditive. Put R=ln(√|Δ|/2).
  When Γ has torsion, Λ is a relation kernel in Γ, not an integer lattice.
- Walk for target p ∈ B: x_j = g_p^{K_p + j}, with K_p = min{k : p^k > H}.
  W = first j with x_j ∈ S (W_ball: first j with x_j ∈ φ(B_R)).
- For χ ∈ Ĥ_B write χ(g_q) = e(θ_q(χ)), with ‖·‖ the distance to Z.
  σ_a(χ) = |B_a|^{-1} Σ_{v∈B_a} χ(φ v) and c_a(χ) = |σ_a(χ)|.
- The predecessor's unramified Lemma 1 gives min_(λ∈Λ minus0)‖λ‖≥2R,
  injectivity on B_R and φ(B_R)⊆S. With ramified 2, this statement uses the
  quotient Γ above. Indeed a principal primitive ideal of norm below |Δ|/4
  must have a rational-integer generator. Split-prime valuations then cancel,
  and its ramified-prime valuation must be even; the canonical exponent is
  zero or one, so it must be zero too. This excludes a nonzero relation in Γ.
  The same triangle-inequality and reduced-form argument proves injectivity
  and containment. The finite scripts check injectivity/containment, not a
  minimum over an infinite lattice.

The unquotiented integer-lattice statement is false: at Δ=-20 the ramified
form (2,2,3) squares to the principal form (1,0,5). Thus 2e_2 is a nonzero
integer relation with norm 2ln2<2R=ln5. Also e_2 and -e_2 both lie in B_R
and have the same image. In Γ they are the same element and 2e_2=0.
The [exact scope controls](../tools/check_window_orbit_duality.py) retain this counterexample.

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

### Theorem 3 (fourth moment and a conditional coefficient bound) — **PROVEN**
Σ_{χ≠1} c_a(χ)⁴ ≥ h_H/|B_{2a}| − 1. If |B_a|<h_H, then
max_{χ≠1} c_a(χ)² ≥ (h_H/|B_{2a}| − 1)/(h_H/|B_a| − 1).
If the numerator is nonpositive this lower bound need not force a large coefficient.

*Proof.* Σ_χ |σ_a|⁴ = h_H |B_a|^{-4} E, where E = #{v_1+v_2 ≡ v_3+v_4 mod Λ} = Σ_x r(x)² and
r(x) = #{(v_1,v_2) ∈ B_a² : φ(v_1+v_2) = x}. The support of r lies in φ(B_{2a}), so Cauchy–Schwarz gives
E ≥ |B_a|⁴/|B_{2a}|. Subtract the χ = 1 term, then use max c² ≥ Σc⁴/Σc². ∎

NUMERICAL: 0 violations out of 1168 instances (median slack ×2.2–3.0). The maximum c_a over χ ≠ 1 has
median 0.93 and minimum 0.52. This does not imply a universal bound on every phase θ_q.

A separate **conditional pigeonhole lemma** is elementary: if H is generated by m
listed elements, L>=2 is an integer and |H|>L^m, then a nontrivial character has
‖θ_q‖<1/L on every listed generator. Partition the phase torus into L^m half-open
cubes and divide two distinct characters in the same cube. The size condition is
essential to this argument. For example H=Z/2 with its generator has only one
nontrivial character, at phase 1/2; it has none with phase distance below 1/2.
A large coefficient for A={0} still equals one in that example.

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

Under ERH, Bach's prime-ideal generation bound is 12log²|Δ|; the
[primary treatment by Grenié and Molteni](https://arxiv.org/pdf/1607.02430),
Introduction, recalls this bound and proves sharper versions. If the factor base
contains the generating prime ideals (with any small ramified factor handled),
this rules out the exact nonprincipal annihilator. It does not by itself give
Theorem 1's aggregate approximate-phase estimate. The previous source's proposed
quantitative explicit-formula extension was an assessment, without constants or
a proof; no impossibility claim for all GRH-based methods is made here.

### What would close it — an unresolved mass estimate

A proposed bound at one radius,

    Σ_(χ≠1, ‖θ_p(χ)‖<η) c_(R/2)(χ)² ≤ 1/4,   η≈1/W,

controls only one part of the weighted Fourier sum. It is not by itself equivalent
to short-orbit equidistribution: the remaining arcs and their normalization still
need estimates. The [orbit-relative continuation](window-orbit-normalization.md)
absorbs characters trivial on K=⟨g_p⟩ into the positive orbit mean and derives a
sufficient all-scale inequality. A quantitative bound for the arithmetic balls at
an affordable W remains open. The sampled small contribution of genus characters
in Section 3 is not a universal exclusion of those characters.

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
- "Bach 1990 explicit bounds …": Bach, Math. Comp. 55 (1990) 355–380; conditional
  small-prime-ideal generation. The conservative 12log² bound is recalled in the
  primary treatment linked in Remark 5. This concerns exact generation, not the
  complete target-cycle approximate-phase mass.
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

From the repository root, a fresh scan is

    PY=python3 bash tools/window_duality/run_all.sh
    (cd tools/window_duality && sha256sum -c SHA256SUMS)

It needs cypari2 and NumPy. For analysis of the saved records only:

    python3 tools/window_duality/analyze_scan.py tools/window_duality/data/scan.jsonl
    python3 tools/window_duality/analyze_scan.py tools/window_duality/data/scan_large.jsonl

The original outputs and data are retained under `tools/window_duality/`.
The integration replay reproduces both analysis files byte-for-byte, but does not
rerun PARI because cypari2 is unavailable. Exact class arithmetic, conditional
large-group completeness, floating Fourier summaries and byte-identical analysis
are distinct scopes of evidence.

Determinism: PARI's random generators are fixed with `setrand(1)`, the random control is drawn in
canonical form order, and the analysis sorts records and uses tie-averaged ranks. Two consecutive runs
produce identical `output_*.txt`.

## AI-assistance note
The theorem statements, proofs, scripts, and this text were drafted by an AI research assistant
(Grok Bot) and checked numerically on the box. The proofs are short and should still be read by a
human before being cited. Literature claims are based on the web searches listed in §4, not on a
complete survey.

## Integration scope repairs

Imported from [PR #75](https://github.com/1337isnot1337/integer-factoring/pull/75),
exact source0462759d46485c439fc751c88f7d6e134c9a2fb3, with its author credit and
complete data preserved. The integration quotients the known ramified relation,
distinguishes the fourth-moment bound from simultaneous phase closeness, and
replaces the single-radius equivalence claim by an explicit open all-scale target.
The source commit preserves the original wording. These repairs do not promote
finite spectral evidence to an unconditional arithmetic window or stronger score.
