# Path C: a uniform window bound for prime-form power search

- **Base commit:** `4c3d135aab2f454d09887875d86674e6ed00569b` of https://github.com/1337isnot1337/integer-factoring (pulled on 2026-11-01).
- **Scope:** the open window theorem in `research/smooth-leading-power-selection.md` §3, for the selector `tools/leading_power_smooth_relations.py`.
- **Status:** integrated research checkpoint from PR69; author proof labels below do not denote independent mathematical review.
- **Labels:**
  - **PROVEN:** proved here, with the proof given.
  - **VERIFIED:** checked by an exact script in this directory.
  - **KNOWN:** in the literature, with a citation.
  - **CONJECTURE:** neither proved nor checked.

## 1. The exact lemma

The setup is copied from the tool at the base commit:

- N is odd with b bits. Δ = −N if N ≡ 3 (mod 4), otherwise Δ = −4N.
- H = isqrt(|Δ| // 3).
- B_Y is the split-prime base from `build_prime_form_factor_base(Δ, Y)`, plus 2 when N ≡ 1 (mod 4). Write m = |B_Y| and P = ∏ B_Y.
- For p ∈ B_Y, g_p is the base's reduced prime form, and K_p = min{k : p^k > H}.
- Q_{p,j} = red(g_p^{K_p+j}) has leading coefficient A_{p,j}, which satisfies 1 ≤ A_{p,j} ≤ H.

**WL(Y, W).** For every such N and every p ∈ B_Y there is some j ∈ [0, W) where either A_{p,j} divides P^{bitlen(H)}, so A is B_Y-smooth, or 1 < gcd(A_{p,j}, N) < N.

The repo's argument (note §3, read but not re-derived) is as follows. Suppose WL holds with Y = 2^{O((log b)^c)} and W = exp(O(b/(log b)^α)), and B_Y generates Cl(Δ). For generation, the repo uses Cureton–Thorner, https://hexagonmath.org/2610.00200v2, via `research/polylog-prime-form-generators.md`. Then the square, nonsingular matrix of rows k·e_p − v, fed to the ambiguous-class extractor, factors N deterministically in time m²·W·poly(b). The whole path therefore rests on WL.

Equivalent views:

- **Hitting problem (PROVEN).** Let S ⊂ Cl(Δ) be the classes whose reduced form has a B_Y-smooth leading coefficient. Reduced forms correspond one-to-one to classes, up to the sign of b, which does not change a. WL then says that every progression {[g_p]^{K_p+j} : 0 ≤ j < W} meets S. So WL is a *short-orbit equidistribution* statement in Cl(Δ).
- **Trivial bound (PROVEN).** The identity class has a = 1, which is smooth. Every orbit therefore hits within ord([g_p]) ≤ h(Δ) steps, so WL(Y, h(Δ)) holds for every Y. KNOWN: h(Δ) ≪ |Δ|^{1/2} log|Δ| (class number formula with L(1, χ) ≪ log|Δ|), so W ≈ 2^{b/2+O(log b)}. This is the strongest unconditional uniform bound I found. It gives α = 0 and is useless.

## 2. Lower bound: an explicit family with long misses (PROVEN, VERIFIED)

**Claim P3.** Fix an odd prime p, a cutoff Y and integers K > W ≥ 1. Take a prime A > Y, A ≠ p, and an odd x with 0 < x < A and p ∤ x, such that:

- 4p^{K+W}A − x² ≥ 3p^{2K−2}, and
- A < (3/4)p^{K−W}.

Put Δ = x² − 4p^{K+W}A, so N = −Δ ≡ 3 (mod 4). Then:

- p is split.
- K_p = K.
- A_{p,j} = A·p^{W−j} for 0 ≤ j ≤ W.
- A_{p,W+i} = A·p^i for every i with 2A·p^i < |Δ|^{1/2}, which includes all i ≤ W.
- No A_{p,j} is smooth and none shares a factor with N. So the first hit for p is at offset j ≥ 2W+1. Choosing W as large as possible gives j ≥ b/log₂p − 2·log_p Y − O(1).

*Proof.*

1. M = p^{K+W} and F = (M, x, A) has discriminant Δ. F is primitive because p ∤ x.
2. F corresponds to a primitive ideal of norm M with p split. That ideal is 𝔭₀^{K+W} for one prime 𝔭₀ above p. Applying S gives F ~ (A, −x, M), which is a primitive ideal 𝔄 of norm A. So [𝔭₀]^{K+W} = [𝔄]. If the repo's orientation is 𝔭 = 𝔭̄₀, replace 𝔄 by 𝔄̄, which has the same norm.
3. For 0 ≤ j ≤ W, [𝔭]^{K+j} = [𝔄·𝔭̄^{W−j}]. This ideal is primitive (the norms are coprime) with norm A·p^{W−j}.
4. That norm is below |Δ|^{1/2}/2, because A²p^{2W} < (3/4)MA ≤ MA − x²/4 = |Δ|/4. A primitive ideal of norm a < |Δ|^{1/2}/2 gives a reduced form, since c ≥ |Δ|/4a > a. So A_{p,j} = A·p^{W−j}.
5. The same argument applied to 𝔄·𝔭^i gives A_{p,W+i} = A·p^i.
6. Bounding H: |Δ| < 4MA < 3p^{2K}, so H < p^K. The hypothesis gives |Δ| ≥ 3p^{2K−2}, so H ≥ p^{K−1}. Hence K_p = K.
7. A > Y is prime, so no A·p^i is Y-smooth. gcd(A, N) = gcd(A, x²) = 1 and gcd(p, N) = 1. ∎

**VERIFIED:** `adversary.py` builds the family for p ∈ {3, 5, 7}, b ≈ 64, 96, 128, 160, at Y = 4093. Results are in `data/adversary_cutoff4093.json`. Every predicted leading coefficient equals the repo's actual orbit (`predictions_match: true`), and K_p = K in all 12 cases.

- With `--semiprime`, x is chosen so that N is a product of two primes. PARI is used only to pick x. This gives another 12 cases, in `data/adversary_semiprime_cutoff4093.json`. All match.
- Example: b = 160, p = 3 gives guaranteed misses 85. N = 2501117468911849 · 514943113216362068753106498087539.

So no constant or o(b) window is uniform. The family only forces linear windows, while the random model below predicts exponential ones. Explicit algebraic families appear to stop at O(b) misses (CONJECTURE).

## 3. Why the natural proof methods fail

In this section, h = h(Δ), δ = |S|/h, and Ŝ(χ) = Σ_{c∈S} χ(c).

- **B1. Uniform character-sum bounds (PROVEN).**
  - Take any weights w ≥ 0 on the window and the orbit f = Σ_j w_j·1_{g^{K+j}}. Then #hits = h^{-1} Σ_χ f̂(χ)·conj(Ŝ(χ)).
  - Suppose the only input is |Ŝ(χ)| ≤ B for all χ ≠ 1. The error can only be bounded by h^{-1} B Σ_{χ≠1}|f̂(χ)|.
  - Fourier inversion gives Σ_χ |f̂(χ)| ≥ h·max f, so the error bound is ≥ B·max w − O(B·W/h). The main term is ≤ W·δ·max w, so the method needs W ≥ B/δ.
  - Parseval gives Σ_χ |Ŝ|² = h|S|, so max_{χ≠1} |Ŝ| ≥ (|S|(1−δ)·h/(h−1))^{1/2}.
  - Hence W ≥ (1 − o(1))·(h/δ)^{1/2} ≥ 2^{b/4 − o(b)}.
  - This holds **whatever L-function hypothesis** produces B: GRH, the 31/32 strip, or zero-density.
- **B2. Cauchy–Schwarz / large sieve (PROVEN).** Bounding Σ|f̂||Ŝ| by ‖f̂‖₂‖Ŝ‖₂ needs W > h/δ, which is worse than trivial. Mixing the window with a second enumerated smooth set S_M (`covering` variant) needs W·|S_M| > h/δ.
- **B3. Pointwise covering (PROVEN reduction).** A deterministic "smooth ideal in every class" statement, used as Cl = S_H·S_M, needs every class to contain a smooth ideal of norm ≤ |Δ|^{1/2}·M with M subexponential. This is the class-group analogue of "every residue class mod q has a y-smooth element ≤ q^{1+o(1)}" with y = q^{o(1)}.
  - KNOWN: the Dirichlet record is q^{9/4+ε} with y = q^{1/(4√e)+ε} (Harman 1999) and q^{2+o(1)} with y = q^ε (Munsch–Shparlinski). https://arxiv.org/abs/1408.4515 says "all current methods require kϑ>2". Granville–Soundararajan discuss the least-nonresidue barrier: https://arxiv.org/abs/0707.0299.
  - The repo's existence result (`research/powered-prime-smooth-norms.md`) sits at norm |Δ|^{1+o(1)}, exactly at that barrier.
  - WL is stronger still: it asks for smooth elements in a *prescribed short geometric progression*, not just in each class.
- **B4. Saddle-point / pretentious savings (heuristic, argued).** Even under GRH, the known non-pretentious savings for smooth-supported character sums are about u^{−cu} relative to Ψ (KNOWN for Dirichlet characters, Harper/Soundararajan line). Beating the B1 threshold would need savings ≈ 1/h, which forces Y = O(1). That contradicts δ > 0 being subexponential.
- **Which parts of GRH are used.**
  - GRH and the strip give *generation* (Bach-type, and Cureton–Thorner) and equidistribution of *prime ideals of norm ≳ h²·polylog*. With only the strip, the norm needed rises to ≳ h^{32}.
  - Neither says anything about orbits of length ≪ h^{1/2}. By B1, no bound on |Ŝ(χ)| alone can work, so zero-density estimates do not help either.
  - KNOWN: no deterministic subexponential factoring algorithm is known even under GRH. The best deterministic rigorous bound is Õ(N^{1/5}) (Harvey, https://arxiv.org/abs/2010.05450; Hittmeir). The best GRH-conditional class-group method is randomized L[1/2, √2] (Hafner–McCurley, https://doi.org/10.1090/s0894-0347-1989-1002631-0). So WL would be a new result even under GRH.
- **Averaged-over-N angle (argued).** "WL for almost all N" might be approachable. Lenstra–Pomerance remark that the class-group method likely works unconditionally for almost all d. But the repo's score is worst-case over N, so the exceptions need a fallback with the same complexity. No efficient test recognizes the exceptional N, and the §2 family shows the exceptions have no simple algebraic shape. This angle cannot produce the score.
- **Closest literature on short orbits.**
  - KNOWN: Harcos–Michel prove equidistribution of subgroup orbits only for index o(D^{1/23042}), https://doi.org/10.1007/s00222-005-0468-6.
  - Shparlinski (*Character Sums with Exponential Functions*) and Konyagin–Shparlinski 2012 cover the Dirichlet analogue {g^x mod q} in short ranges. Nontrivial results there need ranges q^{ε}-type or longer, far above subexponential-in-log q windows. Matomäki–Radziwiłł (https://arxiv.org/abs/1501.04585, https://doi.org/10.4007/annals.2016.183.3.6) give smooth numbers in short *intervals* [X, X + C√X]. Orbit points are not intervals, so this does not transfer.

## 4. Empirics (VERIFIED data, CONJECTURE fits)

`window_scan.py` walks the repo's exact orbit using the repo's own `power`, `compose_with_prime_form`, `build_prime_form_factor_base` and primorial smoothness test.

- N are seeded, deterministic, balanced semiprimes: SHA-256 seeds, then next prime.
- The script records every smooth offset until a fixed number of hits (`--hits`) or the cap of 4·10⁶ steps.
- δ̂_N = hits / steps, the geometric MLE with censoring.

The runs were:

- **`scan_fixed4093`:** Y = 4093 (m ≈ 280); b = 60, 70, …, 160; 8 N per size, 16 targets spread over the base, 4 hits each.
- **`scan_bsq`:** Y = b² (m from 254 to 1426); same sizes.
- **`scan_alltargets`:** Y = 4093; b = 80, 100, 120; 16 N per size; every base prime is a target.
- **`scan_tail`:** Y = 4093; b = 90, 100, 110; 64 N per size; all targets.

In total that is 70,465 (N, p) orbits and about 4.0·10⁸ composition steps. No gcd hit occurred. Full tables are in `data/summary.md`.

**VERIFIED observations:**

1. **Orbit density equals Dickman's ρ(u), with u = ln H / ln Y.** Fitting per-N log₂(1/δ̂) against −log₂ρ(u) gives:
   - slope 1.011, intercept −0.08, R² 0.979 (fixed Y, 88 N);
   - slope 1.005, intercept −0.05, R² 0.972 (Y = b², 88 N).
   - Example: at b = 160 and Y = 4093, the median 1/δ̂ is 2^{18.2} and 1/ρ(6.58) is 2^{18.2}.
2. **First hits look geometric for p ≥ 32.** j·δ̂ has mean 0.98–1.03 and CV ≈ 1.0. Gaps between hits match too.
   - In the all-target runs, the median over N of (max over the base of the window used)·δ̂ is 5.9–6.5. The Exp(1) model predicts ln m + γ ≈ 6.2.
   - Across 67,649 pairs, P(z > 4) is 1.90·10⁻² against e⁻⁴ = 1.83·10⁻².
3. **Small primes are worse than random, because hits clump.** For p < 16 the mean z is 1.5–1.8. For p = 2, 74% of second hits come at the very next step. The tail is heavier as a result: the maximum z was 21.2 (p = 2), where the model predicts about 11.7. See `data/small_p.md`. A smooth a for small p tends to stay smooth as a·p or a/p, so hits cluster and arrive later. This is a constant-factor effect in the data, not a change in growth rate.
4. **Growth at fixed Y is exponential.** log₂(1/δ̂) ≈ 0.159·b − 7.2 (R² 0.976). The largest window used at b = 160 was 3,394,412, and 3 of 128 targets were censored at 4·10⁶.

**CONJECTURE (fit):** W_needed(N, Y) ≈ C·(ln m + γ)/ρ(u), with u = ln H / ln Y and C ≤ 3 (C reflects the small-p clumping). If this holds uniformly, then:

- A polylog base Y = b^c gives u ≈ b/(2c·log₂b), so 1/ρ(u) = 2^{Θ(b)}, which is exponential.
- Y = 2^{(log₂ b)^{1+α}} gives W = exp(O(b/(log b)^{α})) with a quasipolynomial base. So the conjectured uniform version **would** deliver the repo's target for any α. Balancing m² against W gives L[1/2], the known randomized exponent.
- The data support the random model on average and for typical worst targets. They cannot test the uniform-over-all-N statement. §2 shows that structured N exist where the model fails for a single target, though only linearly so far.

## 5. Verdict

- **Statement:** WL(Y, W) is stated exactly in §1.
- **Unconditional bounds:** the upper bound is h(Δ) and the lower bound is about b/log₂p (both PROVEN).
- **Data:** the orbits behave like random walks with density ρ(u) (VERIFIED on 70k orbits, 60–160 bits).
- **Status of a proof: walled.** No unconditional or strip/GRH-conditional route is visible.
  - Completion using only a uniform absolute bound for the nontrivial smooth-set character sums stops at W ≳ (h/δ)^{1/2} (B1, author proof). Phase information, character-dependent bounds and other Fourier arguments are outside this statement.
  - The covering reformulation sits past the open Dirichlet barrier "kϑ > 2" (B3).
  - Not even GRH is known to give deterministic subexponential factoring.
- **Status of the conjecture: alive.** As a heuristic it is consistent with all data.

## Literature searched (2026-11-01)

The searches covered:

- smooth numbers in short intervals (Matomäki–Radziwiłł);
- smooth numbers in every residue class (Harman 1999, Munsch–Shparlinski, Harman–Shparlinski arXiv:1408.4515, Granville–Soundararajan arXiv:0707.0299);
- equidistribution of class-group subgroup orbits (Harcos–Michel 2006);
- exponential-function character sums and consecutive-power gaps (Shparlinski's book, Konyagin–Shparlinski 2012);
- the rigorous class-group factoring and deterministic factoring records (Hafner–McCurley 1989, Lenstra–Pomerance, Harvey arXiv:2010.05450, Hittmeir);
- the quasi-RH strip (Cureton–Thorner, https://hexagonmath.org/2610.00200v2; Lean repo https://github.com/AxiomMath/DedekindZetaZeroFree/tree/09afbb4ce6570c1615d58509dc197fc9f0956b7d);
- Ivanyos–Karpinski–Saxena, which is not relevant.

No result giving WL, or any subexponential uniform short-orbit hitting bound in class groups, was found.

## Rerun

```
cd tools/window_lemma
rm -f data/*.jsonl && ./run_all.sh          # ~45 min on 8 cores (scans ~30 min, adversary ~15 min)
python3 -m venv /tmp/venv && /tmp/venv/bin/pip install numpy cypari2   # analyze.py / --semiprime only
```

- `window_scan.py`, `tail_check.py` and `small_p.py` need only the standard library and the repo's `tools/`. Set `IF_TOOLS` to point elsewhere.
- `analyze.py` needs numpy.
- `adversary.py --semiprime` needs cypari2, which is used only to *choose* x.

## AI disclosure

Grok Bot (an AI agent) wrote this note, the scripts and the proofs, working for Alkindus Arifuddin. The proofs in §1–3 have not been checked by a human. The numbers come from the scripts above.
