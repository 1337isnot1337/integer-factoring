# Deterministic integer factoring bounds

[![Verification](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the machine outputs its smallest prime factor. Runtime is measured in finite-alphabet Turing-machine steps, with `n` equal to the input bit length.

## Main theorem: logarithmic saving exponent 1/2

**Lean proves the exact pinned `Factoring.SqrtLogLogConjecture`:** one deterministic finite-alphabet machine computes the smallest prime factor for every input in time

```text
T(n) <= C * 2^(C * n / (Nat.sqrt (Nat.log 2 n) + 1))
     = exp(O(n / sqrt(log n))).
```

The constant `C` does not depend on the input. This gives score **`alpha = 1/2`** under the [Big O Speedrun rule](https://bigospeedrun.com/integer-factoring). The bound is asymptotic. It is not a practical factoring speedup, and independent review remains pending.

The complete Lean endpoint is [`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](formal-appendix/HalfTheorem.lean). It proves the source statement `Factoring.SqrtLogLogConjecture` and depends only on `propext`, `Classical.choice`, and `Quot.sound`. The proof covers the factoring-table constructor, deterministic grouping, value preservation, padding, prefix recovery, and the finite-alphabet machine bound.

The composition uses the existing clique evaluator from [Anthropic's Appendix B](https://github.com/anthropics/formal-math) and the multiplication-network reduction by [Isaac Rudich and Louis-Martin Rousseau](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring). See the [written proof](research/half-rate-factoring.md), the [constructor](research/half-top-constructor.md), and the [evaluator audit](research/half-evaluator-audit.md).

## Results and evidence

| Claim | Status | Details |
| --- | --- | --- |
| Saving exponent `1/2` | **Lean-verified** | [Final theorem](formal-appendix/HalfTheorem.lean), [proof modules](formal-appendix/README.md) |
| Saving exponent `1/30` for the k6 machine | **Lean-verified** | [Proof](research/verified-k6.md), [machine theorem](variants/k6/FactoringK6/Certified.lean) |
| Sharper k6 saving `log(12/11)/log(12) > 1/29` | Ordinary proof, not Lean-verified at this rate | [Rate derivation](research/k6-sharp-rate.md) |
| Saving exponent above `1/2` | **Not established** | [Research index](research/INDEX.md) |
| Score transfer with arbitrary vanishing width error | Ordinary conditional theorem | [Proof and effective hypotheses](research/vanishing-slack-score-transfer.md) |
| Uniform endpoint from arbitrarily sparse finite seeds | Ordinary conditional theorem | [Bounded search and type amplification](research/gap-free-seed-score-transfer.md) |
| One uniform finite-seed endpoint over all coefficient fields | Ordinary conditional theorem | [Input-wise selection](research/universal-field-seed-selection.md), [exact rank specialization](research/finite-field-rank-specialization.md) |
| Auxiliary divisor-count exponent with leading coefficient one | Ordinary unconditional bound; overall score unchanged | [Exact integer parameter](research/optimized-divisor-count-bound.md) |
| All-field binary-word and middle-layer obstruction | Ordinary theorem; restricted family only | [Uniform rank floor `rho >= 3/4`](research/all-field-binary-word-barrier.md) |
| Matrix-entry-oracle rank and affine queries | Quasipolynomial time and polylogarithmic auxiliary space; local algorithm only | [Polynomial weighting and recomputed bases](research/polylog-space-rank-recovery.md), [scalar Berkowitz foundation](research/small-space-berkowitz-rank.md) |
| Minimum-channel matrix factors over truncated polynomial rings | Same local time-space form after coefficient linearization; includes nonunit factors | [Smith counts and virtual generator selection](research/polylog-space-jet-factorization.md) |
| Border and finite-precision seed extraction | Ordinary conditional theorem; no stronger seed supplied | [Coefficient extraction and nonunit budgets](research/border-and-nonunit-seed-transfer.md) |
| Ordered matrix-valued seed transfer | Ordinary conditional theorem with cross-half commutation; no improving seed supplied | [One scalarization after type amplification](research/commuting-halves-matrix-seed-transfer.md) |

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.

An ordinary [second-order audit](research/half-rate-second-order-audit.md) determines the logarithmic term in the existing central channel width: `sqrt(2) lambda+(2sqrt(2)-1)log2(lambda)+O(1)` per cut. It uses the classical partition asymptotic and exact minimal-energy selection. The correction is positive, so sharpening that partition prefactor alone does not improve the achieved rate.

## Verification

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-further       # New transfers, arithmetic bounds, and word certificates
make verify-energy-full   # Optional larger weighted-energy certificate
make verify-search-native # Optional exact finite-state engines; requires C++17
make verify-witt-full     # Optional 202-dimensional Witt Gram calculation
make verify-weyl-orbit    # Optional full energy-32 Weyl orbit regeneration
make verify-docs          # Local Markdown-link check
```

The Python checks use Python 3 and its standard library. Some optional checks use NumPy for speed.
The finite-state search-engine checks also require a C++17 compiler.

For the Lean builds, install the stated Lean release, Lake, Git, and the prerequisites in the corresponding scripts. The scripts fetch pinned upstream revisions.

```sh
make verify-half-appendix # Lean 4.33.1; complete 1/2 theorem
make verify-k6           # Lean 4.35.0; complete 1/30 theorem
make verify-formal       # Lean 4.35.0; original local lemmas
```

The [GitHub Actions workflow](.github/workflows/verify.yml) runs these four groups on PRs and on `main`: exact Python checks, original Lean lemmas, k6, and the half-rate theorem. The half-rate script audits the axioms of its final theorem and supporting declarations.

## Repository map

| Directory | Contents |
| --- | --- |
| [`formal-appendix/`](formal-appendix/README.md) | Lean proof of the `1/2` theorem and its dependencies |
| [`variants/k6/`](variants/k6/FactoringK6/Certified.lean) | Earlier complete `1/30` proof and source patch |
| [`formal/`](formal/) | Supporting Lean lemmas and local checks |
| [`research/`](research/INDEX.md) | Mathematical proofs, candidate constructions, counterexamples, and open questions |
| [`evidence/`](evidence/) | Finite certificates and recorded verification outputs |
| [`tools/`](tools/) and [`tests/`](tests/) | Exact checkers, reproduction scripts, and unit tests |

## Current research frontier

No uniform factoring exponent above `1/2` is proved here. One open approach seeks effective crossing pieces with total rank width `4*rho*lambda + o(lambda)` for `rho < 1/sqrt(2)`. The [score-transfer theorem](research/vanishing-slack-score-transfer.md) permits this arbitrary vanishing relative error, together with piece entropy `lambda/(log lambda)^(delta-o(1))`, and gives score `min(log2(1/rho), delta)`. Its uniform construction and factor-query hypotheses are part of the theorem. It supplies no new crossing family; finite rank improvements alone do not establish a stronger factoring bound.

The [finite-seed endpoint theorem](research/gap-free-seed-score-transfer.md) removes any gap or generation-speed requirement from an enumerable improving sequence. Its [input-wise field extension](research/universal-field-seed-selection.md) gives one machine attaining the best finite rank-sum score over all coefficient fields. The [specialization proof](research/finite-field-rank-specialization.md) preserves every finite rank profile over a finite field, and type amplification makes prime fields sufficient for the same infimum. Finding a better finite certificate remains the mathematical requirement.

Existence of one machine attaining the supremum of achievable scores also follows from [universal verified search](research/universal-verified-factoring-search.md), using polynomial-time primality testing. That is a standard search principle; the explicit tensor selector additionally describes the construction and its local resource costs.

The [all-field binary-word theorem](research/all-field-binary-word-barrier.md) rules out crossing pieces invariant under simultaneous permutation of all binary positions, and includes middle Hamming layers and products of such layers. One piece must have both ranks at least `2^(3t/2-O(log^2 t))`, uniformly in the field, so this family cannot reach the half-rate threshold. This does not restrict arbitrary crossing constructions or factoring algorithms.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

The same word certificate applies to [symmetric binary-energy shells](research/binary-energy-symmetry-barrier.md), giving `rho >= 3/4` for all fields. For the [specific reciprocal pencil](research/weighted-energy-reciprocal-pencil.md), a separate exact block argument proves the sharper optimum `rho = 1`. Other constructions on these low-growth spaces remain open. The [border-transfer theorem](research/border-and-nonunit-seed-transfer.md) makes polynomial or sufficiently precise truncated-parameter identities usable as exact finite seeds, with their coefficient-extraction costs included.

The [finite-state polynomial families](research/finite-state-pairing-crossings.md) give exact order-dependent crossings outside the position-invariance assumption. An exhaustive screen covers 88,006 finite models and retains both maximum rank product and total rank sum. The [run-compressed rank algorithm](research/two-state-run-rank-compression.md) combines a local Jordan decomposition, interval types, and an exact quadratic-pivot identity to reach longer profiles. These algorithms and endpoint identities are proved; the finite profiles do not establish an improved asymptotic contraction.

The larger [four-feature matrix-product family](research/four-feature-parity-crossings.md) adds a local mixed parity term and retains an exact crossing endpoint with only linearly many coefficients. It has explicit order-dependent examples at every length. A separate exhaustive six-position screen covers 131,072 choices; its integer rank-sum bounds exclude every screened model. Independent full-tensor checks are included in `make verify-search-native`, and the note gives complete screen replay commands.

The [eight-feature extension](research/eight-feature-parity-crossings.md) adds local bit dependence and delayed auxiliary transitions. Its structural screens audit 365,482 records representing 365,426 distinct choices, and eleven selected candidates have complete profiles at lengths four, six, and eight. Every selected eight-position rank-sum exponent is at least one. For one order-dependent branch, an exact multiplicative law proves orthogonal coefficient projectors with ranks `2^t binom(t,k)` in both cuts, so its rank sum equals `4^t` at every length. The larger family remains open beyond the stated finite and algebraic results.

[Matrix-unit monodromy](research/matrix-unit-monodromy.md) supplies a direct crossing outside the local parity span. Its [exterior-algebra extension](research/exterior-monodromy-and-occupation-graphs.md) has exact exclusion walks and a graph-based rank certificate for fixed occupation boundaries. The certificate excludes all binary occupations and all 126 ternary half-occupations from the improving asymptotic target. It strengthens the original one-particle boundary restriction to `rho >= 1`. Independent exterior operators, full small cut ranks, and explicit acyclic witnesses are checked in `make verify-search-native`; general coherent boundaries remain open.

[Bivariate site weights](research/site-weighted-pairing-crossings.md) give another exact position-dependent family. A projection and complementation proof shows that centrally symmetric weight multisets have optimal leading contraction `H_2(3/4) = 0.811278...` throughout their coefficient span. Asymmetric weights remain outside that theorem. All seven reported profiles were regenerated, including every one of 16,824 exact coefficient rank pairs; square weights through sixteen positions still fail the rank-sum target. The note explains why the earlier small maximum-product exponents alone do not imply a score gain. Its independent finite checks are included in `make verify-search-native`.

The [Witt state-space construction](research/witt-duality-and-vertex-jets.md) proves an exact odd Gram determinant and effective homogeneous duals at every energy. It restores the pairing on a space with binary-partition growth. A separate rank budget prevents a direct small-window substitution, and the tested scalar exchange interface fails; the duality result supplies no stronger factoring exponent.

The [weighted Weyl construction](research/weighted-weyl-and-multivariate-extraction.md) gives a different exact crossing with ordinary occupation carries and explicit factor channels. A [detecting row filtration and parity-sector restriction](research/weighted-weyl-row-filtration.md) prove the original coefficient-span floor `rho >= 3/4` along every growing energy sequence. The common filtration is proved on power-of-two shells; restriction transfers its rank floor to arbitrary energies. At energy 32, an [exact protected minor](research/weighted-weyl-protected-minor.md) also excludes arbitrary recombinations of all external symmetry images. The same construction note extends finite-jet extraction to several parameters, with its coefficient cost included. The three standard-library checks are part of `make verify-further`.

The [odd-characteristic spectrum analysis](research/odd-characteristic-weyl-spectra.md) derives exact characteristic-three projectors. At energy 81, every original-span crossing decomposition has a piece with row rank at least 33,933 and column rank at least 1,617, excluding its finite target even over extension fields. A separate mixed-algebra central projector at energy 27 is constructed and measured; its split increases the rank-sum cost. [Mismatched characteristic and mode weights](research/mismatched-mode-weyl-crossings.md) produce a different, noncommutative coefficient algebra. At energy 63 a nilpotent-kernel and protected-block argument excludes every polynomial recombination of its constant operator; words in distinct coefficients remain outside that result.

The [constant alternate-order Weyl algebra](research/weyl-alternate-order-algebra.md) permits arbitrary products of the two zero-parameter endpoint operators. An elementary all-block involution, a square-zero coproduct, and a protected diagonal row prove the asymptotic floor `rho >= 3/4` for this entire word algebra. Its exact spectral core is much smaller than the original row matrix; the energy-32 core has minimal-polynomial degree 352. The broader algebra involving nonconstant coefficients remains open. The compact algebra checker is part of `make verify-further`; optional native profile and core commands are in the note.

The [exact coefficient multiplication laws](research/weighted-weyl-additive-row-law.md) simplify both single-orientation families. At power-of-two energy the original curve satisfies `P(z)P(w)=P(z+w+1)`. At every energy the opposite curve has the additive normalization `G(t)=C t^h Q(1/t)` and satisfies `Q(z)^2=z^(2h)I`. Integral Gauss coefficient identities prove these statements and give an elementary replacement for the original row-filtration argument. Mixed-orientation products remain a separate construction problem.

At energy 32, the [full-symmetry orbit certificate](research/weighted-weyl-rank-one-ideal.md) reaches arbitrary words in all eight external-symmetry coefficient families. The top rank-one pairing generates a square-zero ideal; its transpose-closed orbit has dimension 2,051 inside a 4,154-dimensional block, leaving a 52-dimensional middle quotient. The full certificate was regenerated from the binomial entry formula. A separate parity argument excludes a natural affine family at every power-of-two energy. These statements do not exclude the rest of the generated algebra or establish an all-energy radical theorem.

The [matrix-valued seed interface](research/commuting-halves-matrix-seed-transfer.md) allows noncommutativity within each factor half while retaining exact ordered amplification. A fixed auxiliary matrix dimension is paid once. Full tensor-factor halves collapse to an ordinary scalar seed when normalized trace is available; two-by-two halves in matrix order four collapse after a finite field extension in every characteristic. The remaining larger characteristic-dividing cases are constructive search targets. Exact ordered-product, field-extension, and negative interleaving checks are included in `make verify-further`.

The [boundary-stable extension](research/boundary-stable-catalytic-seeds.md) only requires the local sum to have the correct common boundary action, or to pass a finite reachable-quotient endpoint test. It retains the same ordered transfer. An explicit four-dimensional endpoint family admits no normalized product functional preserving that endpoint, even over extension fields, so the earlier scalarization argument does not automatically apply. CAR left/right multiplication gives concrete alternative commuting halves and asymmetric terminal ranks. The supplied toy is exact but has rank sum thirteen, so it is not competitive.

The [diagonal rank-one word transfer](research/diagonal-rank-one-word-transfer.md) gives exact individual word ranks without a cross-half commutation assumption. Grouping by ordered transition counts yields a polynomial-piece spectral construction. An explicit [reduction-map Fierz identity](research/reduction-map-fierz-crossing.md) improves its supplied base from `2d sqrt(d(d-1)/2)` to `d(1+sqrt(d(d-1)/2))`; it remains above the score threshold. An [all-field binary two-piece inequality](research/binary-two-piece-rank-sum-and-f4-boundary.md) and an [exhaustive binary transition control](research/three-piece-binary-boundary-control.md) resolve two unequal-terminal three-piece families without extrapolating finite tests to extension fields. A [four-piece extension with two free ket functions](research/four-piece-f4-boundary-pencil.md) also has exact spectral minimum four.

The [channel-walk extension](research/twisted-channel-walk-transfer.md) permits input-independent scalar commutation phases and discards channel words containing certified zero products. Reciprocal piece weights give a spectral bound with exact rational certificates. A small noncommutative example reduces its supplied amplification base from 12 to 8; that example is still worse than the direct scalar crossing and establishes no score gain. The theorem gives a concrete additional condition to search for, with all-input identity and negative-phase checks in `make verify-further`.

[Integer flow dual certificates](research/channel-flow-dual-certificates.md) give exact limits on this spectral optimization method. They prove that eight is already the best weighted-supervector base for the toy channel graph, and recover the ordinary rank-sum optimum in a complete-graph example. A disconnected-graph control demonstrates that these certificates do not lower-bound true matching-type counts or tensor ranks. The exact rational and integer checks are included in `make verify-further`.

The [component-path refinement](research/component-channel-walk-transfer.md) covers walks by paths of strongly connected components and permits separate reciprocal weights for each row/column path pair. The actual amplified tensor pieces stay unchanged. An exact abstract graph example improves the optimal certified base from two to one; applying the refinement to a competitive crossing tensor remains open. Its checker includes rational certificates, matching integer lower bounds, and a control against improperly separating connected components.

The [exact matching-type growth theorem](research/exact-channel-type-growth.md) identifies the exponential count rate using connected integer base walks and stationary-flow entropy. It proves that independent reciprocal weights on admissible exact edge supports lose no exponential factor in that counting bound. The proof includes zero limiting bridge flows and exact parity incompatibilities; a limsup is necessary. This is an exact analysis of a rank upper bound, so a large count rate does not by itself imply large tensor ranks. A [bounded binary carry algorithm](research/finite-carry-channel-feasibility.md) decides the required positive integer support condition and supplies finite witnesses for feasibility or infeasibility. Together these results compute the count rate to any prescribed rational precision.

[Product-subspace automata](research/scalarized-channel-automata.md) track longer ordered matrix products and their actual scalarized bilinear ranks. They can discard entire components that cannot reach a nonzero terminal form. A zero value at a prefix alone does not justify discarding it. The finite toy verifies these distinctions but retains exponential base eight.

The [generic nonscalarization theorem](research/generic-full-factor-nonscalarization.md) proves a complementary existence result: two generic full factor structures in matrix orders twelve or sixteen over characteristic two need not admit any common normalized product functional, even after extending the field. This shows that the order-four scalarization argument does not settle those larger cases. It does not supply an improving crossing identity; actual factors may still span smaller subspaces that do admit scalarization.

## Contribution and attribution

Start with [CONTRIBUTING.md](CONTRIBUTING.md) and [autoresearch.md](autoresearch.md). Coordinate new research through [issue #8](https://github.com/1337isnot1337/integer-factoring/issues/8). Document source revisions, proof status, reproduction commands, and remaining gaps.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).
