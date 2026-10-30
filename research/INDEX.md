# Research index

This index covers the research notes in this checkout. Open PRs can contain work beyond `main`. Each note states its own assumptions, proof scope, and evidence. An obstruction to a specific construction is not a factoring lower bound.

[Session report: 30 October 2026](session-2026-10-30-further-bounds.md) summarizes the strongest new results, the unchanged score, and the remaining synthesis obligations.

## Proved machine bounds and their foundations

- [Half-rate factoring composition](half-rate-factoring.md): ordinary proof of the `1/2` machine bound; its final Lean theorem is in [`formal-appendix/HalfTheorem.lean`](../formal-appendix/HalfTheorem.lean).
- [Top constructor and codec](half-top-constructor.md): deterministic table construction and exact input translation.
- [Evaluator and encoding audit](half-evaluator-audit.md): recursive evaluator, computable pad, and charged bit costs.
- [Sharp logarithmic channel correction](half-rate-second-order-audit.md): the minimal partition energy gives central width `sqrt(2) lambda+(2sqrt(2)-1)log2(lambda)+O(1)` per cut; the positive correction does not hide a further saving.
- [Verified k6 machine](verified-k6.md): Lean proof of the earlier `1/30` bound.
- [Sharper k6 rate](k6-sharp-rate.md): ordinary rate derivation above `1/29` for the k6 program; this sharper exponent is not Lean-certified.
- [k6 source audit](k6-source-audit.md): hashes, modified source scope, and compiler preservation.
- [Baseline audit](baseline-audit.md): source recurrence and the unchanged reference's machine-cost analysis.

## Transfers and proposed stronger evaluators

- [Weighted Weyl and multivariate extraction](weighted-weyl-and-multivariate-extraction.md): an exact occupation crossing, both supplied channel factorizations, and charged multivariate coefficient extraction; the direct profiles do not improve the score.
- [Exact weighted Weyl multiplication laws](weighted-weyl-additive-row-law.md): integral Gauss identities prove coefficient-algebra closure at every energy, an additive original curve at power-of-two energies, and an additive reciprocal normalization of the opposite curve at every energy.
- [Witt duality and vertex jets](witt-duality-and-vertex-jets.md): an explicit odd Hall determinant and length-unitriangular vacuum maps give effective homogeneous duals at all weights; the crossing and rank-budget obligations remain distinct.

- [Binary trace conservation kernels](trace-bilinear-conservation-kernels.md): explicit binary tensors have both cut ranks `(2^t-1)2^k` and ordinary kernel rank at least `binom(t,k)`, disproving a proposed general rank obstruction already at `t=33,k=11`; a competitive crossing synthesis remains open.
- [Small trace-kernel certificate](trace-kernel-small-certificate.md): complete multiplicative coefficient support gives the exact ordinary rank 94,183 at `t=19,k=7`; nineteen is the first possible extension degree for a counterexample in this trace-product family.
- [Restricted trace kernels and composition](trace-kernel-restriction-and-composition.md): exact syndrome-image intersections compute both cut ranks on arbitrary label subsets; low-degree full-cube kernel compositions vanish, including controlled versions.
- [Complete two-patch trace crossing](two-patch-trace-crossing.md): two complementary trace patches give an exact crossing of degree `2t-1`; five complete finite profiles and independent dense checks record the raw-rank limits.
- [Exact small trace-span optima](two-patch-span-recombination.md): every F2 recombination and every number of pieces in the two-patch spans at widths three and four has total rank cost at least `d^2`, certified by all-edge integer potentials and independent dense ranks.
- [Tagged trace-coordinate controls](tagged-two-patch-span-controls.md): all fourteen choices at widths three and four retain exact binary-span optimum `d^2`, after explicit basis reduction; a larger bivariate span remains untested.
- [Finite-state polynomial crossings](finite-state-pairing-crossings.md): explicit order-dependent identities, exact XOR-circulant rank reduction, and the complete 88,006-model finite screen.
- [Four-feature matrix-product crossings](four-feature-parity-crossings.md): a larger local tensor space has constant auxiliary dimension and order witnesses at every length; a complete 131,072-model screen and independent full-tensor checks record the finite limits.
- [Eight-feature parity crossings](eight-feature-parity-crossings.md): exact local bit dependence and delayed transitions enlarge the family; 365,482 screen records and 33 selected profiles are audited, and an all-length multiplicative law explains the exact `4^t` rank sum of one selected branch.
- [Matrix-unit monodromy](matrix-unit-monodromy.md): an auxiliary walk gives a direct crossing outside the local even-parity span, with exact histogram blocks and a protected product-boundary restriction.
- [Exterior monodromy and occupation graphs](exterior-monodromy-and-occupation-graphs.md): acyclic restrictions give full coefficient-span rank bounds; they exclude every binary fixed occupation and all 126 ternary half-occupations, while coherent boundaries remain outside the theorem.
- [Site-weighted pairing crossings](site-weighted-pairing-crossings.md): bivariate position weights give exact complementary rank profiles; symmetric weight multisets have the sharp leading floor `H_2(3/4)`, while seven asymmetric and comparison profiles through length sixteen retain every exact coefficient rank.
- [Polynomial shear diagnostics](polynomial-shear-diagnostics.md): 27 complete polynomial-conjugation profiles have exact rank-sum lower bounds above `d^2`; an independently checked engine and an explicitly partial width-18 restart checkpoint are preserved.
- [Run-compressed crossing ranks](two-state-run-rank-compression.md): exact local Jordan decomposition, interval-type grouping, and quadratic-pivot rank reduction accelerate longer constructive searches.

- [Border and nonunit seed extraction](border-and-nonunit-seed-transfer.md): finite-jet and polynomial identities yield exact seeds with rank-sum loss at most `(kt+1)` after `k` copies; weighted budgets still hold over finite chain rings.
- [Matrix-valued seeds with commuting halves](commuting-halves-matrix-seed-transfer.md): ordered type amplification pays a fixed matrix scalarization cost once; normalized trace and a common product-vector argument identify two cases that already collapse to scalar seeds.
- [Boundary-stable catalytic seeds](boundary-stable-catalytic-seeds.md): a common boundary action or a finite observable quotient suffices for the same transfer; an explicit order-four endpoint evades full-factor scalarization, and CAR operators give two concrete factor structures with unequal terminal ranks.
- [Diagonal rank-one auxiliary words](diagonal-rank-one-word-transfer.md): exact adjacent transition ranks give a spectral transfer without cross-half commutation; ordered transition counts are essential, and common terminal projections have a precise scalar collapse.
- [Observable transition realization](observable-transition-realization.md): any finitely certified observable matrix table lifts explicitly to diagonal rank-one crossing pieces, including an exact example with no common left or right crossing eigenvector.
- [Linear detectors for observable tables](linear-rank-detectors-for-observable-tables.md): a two-by-two linear detector proves spectral base at least four for every characteristic-two binary observable table with entries in the three-pairing span, at any auxiliary size; adding a separable coordinate tensor escapes that detector.
- [Local alphabet and pairing spectral bounds](local-alphabet-and-pairing-spectral-barriers.md): a universal detector requires at least five local labels for an improving word-rank spectral transfer; the three-pairing scalar space is excluded in characteristic two at every alphabet and auxiliary size.
- [Finite functional detector certificates](finite-functional-rank-detector-certificates.md): exact projective inequalities certify an entire scalar entry space before auxiliary search, with separate F2, F3, and F4 certificates and an explicit extension-field negative control.
- [Overlapping subalphabet certificates](overlapping-subalphabet-spectral-certificates.md): all six orbits of two or three three-label corrections on five labels have exact F2 spectral lower bounds at least thirteen, uniformly in auxiliary dimension.
- [Twisted-pairing catalogue](twisted-pairing-observable-preflight.md): 7,500 indexed F2 choices, allowing repeated underlying spaces, each have crossing-coefficient detector weight eleven; 225,000 exact cut-rank evaluations exclude them before auxiliary search.
- [Reduction-map Fierz crossing](reduction-map-fierz-crossing.md): an explicit boundary identity improves its supplied amplification base from `2d sqrt(d(d-1)/2)` to `d(1+sqrt(d(d-1)/2))`, still above the score threshold.
- [All-field binary two-piece inequality](binary-two-piece-rank-sum-and-f4-boundary.md): every two-piece binary crossing has rank sum at least four; a three-piece F4 boundary pencil has exact spectral minimum four over every extension.
- [Four-piece F4 boundary pencil](four-piece-f4-boundary-pencil.md): two independent free ket functions still have exact spectral minimum four over every extension; the proof treats zero transitions and all admissible binary rank profiles.
- [Three-piece binary boundary control](three-piece-binary-boundary-control.md): the `D,G,D+G` terminal family has exact word-rank spectral minimum four in every characteristic-two field, with complete F2 transition frontiers.
- [Scalar twists and channel walks](twisted-channel-walk-transfer.md): input-independent scalar commutation and uniformly zero channel products give exact weighted-walk bounds, certified by rational positive supervectors.
- [Integer flow dual certificates](channel-flow-dual-certificates.md): stationary edge counts with equal piece totals bound every reciprocal-weight supervector certificate, proving exact graph optima while distinguishing them from actual tensor ranks.
- [Component-path channel refinement](component-channel-walk-transfer.md): independent rational certificates on pairs of condensation paths improve the counting bound without changing the tensor pieces; an exact abstract example reduces the best certified base from two to one.
- [Exact matching-type growth](exact-channel-type-growth.md): connected integer base walks and stationary flows give the exact limsup entropy; the infimum over independent reciprocal weights on admissible exact supports equals that count rate, including connectivity and congruence restrictions.
- [Finite carry feasibility](finite-carry-channel-feasibility.md): a bounded binary carry graph decides exact positive support feasibility, with finite certificates for both outcomes; combined with entropy and rational dual bounds, the matching-type count rate is uniformly computable.
- [Product-subspace automata and scalarization](scalarized-channel-automata.md): a finite automaton tracks ordered factor spans, charges their terminal bilinear ranks, and removes precisely the states with no nonzero continuation.
- [Generic full-factor nonscalarization](generic-full-factor-nonscalarization.md): incidence dimensions prove that generic characteristic-two auxiliary structures of orders twelve and sixteen admit no normalized functional factorizing across both full tensor-factor pairs; constructing low-channel crossings remains open.

- [Vanishing-slack score transfer](vanishing-slack-score-transfer.md): arbitrary `o(lambda)` width error preserves the leaderboard exponent under stated uniform algorithms and entropy bounds; no improved unconditional exponent.
- [Gap-free finite-seed score transfer](gap-free-seed-score-transfer.md): bounded search attains an enumerable rank-sum infimum, with no seed-gap or generation-speed hypothesis.
- [Input-wise field and seed selection](universal-field-seed-selection.md): one uniform machine attains the finite-certificate endpoint over all coefficient fields, retaining one selected field throughout each input's recursion.
- [Finite-field rank specialization](finite-field-rank-specialization.md): every finite rank profile over an arbitrary field specializes exactly to a finite field; prime-field seeds approximate the same infimum, and ring factor-channel bounds specialize without increasing.
- [Universal verified factoring search](universal-verified-factoring-search.md): the standard verified-search principle already gives one machine attaining the supremum score; direct constructions retain additional resource information.
- [Optimized divisor-count bound](optimized-divisor-count-bound.md): a computable `q_*(n)=(1+o(1))n/log2 n` replaces the earlier auxiliary leading coefficient nine by one.
- [Exact optimal divisor parameter](exact-divisor-count-parameter.md): an exact bit-length oracle has `2^O(sqrt(n))` time and a proved `O(sqrt(n) log n)`-space construction; at 128 bits it improves `q` from 65 to the optimal 29. The fast shipped oracle and the complete theoretical low-space construction are distinguished.
- [Polylogarithmic-space rank and affine recovery](polylog-space-rank-recovery.md): polynomial weighting and recomputed prefix bases give `m^O(log m)` time and `O(log m (log m+log q))` auxiliary bits over every prime field.
- [Minimum-channel jet factorizations](polylog-space-jet-factorization.md): Toeplitz coefficient ranks and full-image radical quotients recover minimum factorizations and affine solutions over `F_q[t]/(t^H)` with polylogarithmic auxiliary workspace.
- [Small-space Berkowitz foundation](small-space-berkowitz-rank.md): scalar characteristic-coefficient queries, division-free identities, and the earlier linear-space pivot-mask algorithm.
- [Rank-sensitive digit-box Fourier contraction](hankel-digit-box-contraction.md): one restricted digit family reduces to an affine box sum; explicit outer traversal still has exponentially many full-rank frequencies.
- [Base-prime carry lifts](padic-digit-carry-transfer.md): exact lifting, digit-pivot enumeration, dense Fourier support, and independent small-characteristic determinant accounting.
- [Pieces-to-factoring transfer](pieces-to-factoring.md): conditional crossing-parameter transfer, with full bit-cost accounting.
- [Scale-one modular transfer](scale-one-modular-transfer.md): fixed-field transfer and machine costs.
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

## Obstructions and work beyond 1/2

- [Weighted Weyl row filtration](weighted-weyl-row-filtration.md): an actual common module filtration handles power-of-two energies; a parity-sector restriction extends the original coefficient-span floor `rho >= 3/4` to every growing energy sequence.
- [Weighted Weyl protected minor](weighted-weyl-protected-minor.md): an exact 1,821-pair witness excludes the entire coefficient and external-symmetry span on the 202 labels of energy 32.
- [Constant alternate-order Weyl algebra](weyl-alternate-order-algebra.md): an elementary endpoint involution and square-zero coproduct give a central normal form, a smaller spectral core, and a `rho >= 3/4` obstruction for the entire constant-generator word algebra.
- [Full-symmetry Weyl rank-one ideal](weighted-weyl-rank-one-ideal.md): exact energy-32 orbit closure proves a square-zero ideal for arbitrary words in all external-symmetry coefficient families, with a 52-dimensional middle quotient; a separate all-size parity minor excludes the natural affine perturbation search.
- [Odd-characteristic Weyl spectra](odd-characteristic-weyl-spectra.md): characteristic-three projectors give an exact all-recombination obstruction for the original 239-label coefficient span; a nondegenerate mixed orbit yields a separate central projector with checked cut ranks.
- [Mismatched mode Weyl crossings](mismatched-mode-weyl-crossings.md): characteristic-two radix-three weights give a noncommutative coefficient family; a nilpotent-kernel and zero-row certificate excludes every polynomial in its energy-63 constant operator.
- [Mismatched commutator projector](mismatched-commutator-projector.md): a short word in distinct normalized coefficients gives an exact projector and paired-sector tensor at every energy, with cuts `(361,1444)` in the 132-label instance.
- [Partial-mode Weyl sweeps](partial-mode-weyl-sweeps.md): frozen-mode crossings give ordered endpoints with `O(h log h)` degree and a unique-flow entry formula; exact profiles and constant-power controls identify the current finite limits.
- [All-field binary-energy shell symmetry](binary-energy-symmetry-barrier.md): an explicit balanced-layer embedding forces `rho >= 3/4` for every equal-weight-class-invariant crossing decomposition.
- [Weighted-energy reciprocal pencil](weighted-energy-reciprocal-pencil.md): genuine projected tilting blocks prove the exact leading optimum `rho = 1` for the full characteristic-two coefficient/symmetry span.
- [Lacunary kernel boundary](lacunary-kernel-boundary.md): binary-partition state growth, the exact scalar Artin--Schreier fixed field, and the limits of finite scalar locality; other exchange relations remain open.

- [All-field binary-word and layer barrier](all-field-binary-word-barrier.md): a common eight-monomial functional forces `rho >= 3/4` for the entire position-invariant space, uniformly in the field, with middle-layer and product-layer versions.
- [Characteristic-two word proof](modular-binary-word-barrier.md): the modular simple-module certificate and a composition-series rank argument underlying the all-field extension.
- [Directions beyond half](post-half-directions.md): experiments, scope, and unresolved questions.
- [Network tradeoff](half-top-structure.md): restriction for the current local-enumeration and wire-channel constructor.
- [Restricted-label counterexample](half-restricted-labels.md): failure of a proposed low-mode dimension reduction.
- [Restricted intermediate spans](half-restricted-span.md): quantitative obstruction, with precise selected-label conditions.
- [Half-rate recurrence barrier](half-recurrence-barrier.md): limits of the existing operators and recurrence.
- [Graded-space barrier](graded-space-barrier.md): lower bound under the stated graded-space assumptions.
- [Rectangular-grid frontier](rectangular-grid-frontier.md): limits of unequal splits under the declared-profile recurrence.
- [Top-channel rank](top-channel-rank.md): exact factor rank on the full partition domain.
- [Coefficient basis](coefficient-basis.md): scoped minimax calculation for the chosen coefficient span.

## Earlier grouping and crossing research

- [Integration notes](integration.md): earlier combined grouping/crossing proposal and open proof obligations.
- [Recurrence](recurrence.md): the original recurrence and cost accounting.
- [Grouping improvement](grouping-improvement.md): concentration argument and deterministic construction.
- [Grouping validation](grouping_validation.md): finite grouping examples.
- [Crossing compression](crossing-compression.md): explicit integer factorization and channel counts.
- [Grid-level validation](grid-level-validation.md): exact small-grid tests.

## Verification and proof artifacts

- [`../formal-appendix/README.md`](../formal-appendix/README.md) maps the Lean `1/2` proof.
- [`../variants/k6/FactoringK6/Certified.lean`](../variants/k6/FactoringK6/Certified.lean) gives the Lean `1/30` theorem.
- [`../evidence/`](../evidence/) stores finite certificates, source audits, and recorded CI output.
- [`../tools/`](../tools/) contains the checkers. Run `make verify` from the repository root.

No finite certificate alone proves a new asymptotic factoring exponent. Distinguish an exact finite identity from a completed uniform machine theorem.
