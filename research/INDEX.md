# Research index

This index covers research notes on `main`. Open PRs can contain later work. Each note states its own assumptions, proof scope, and evidence. An obstruction to a specific construction is not a factoring lower bound.

## Proved machine bounds and their foundations

- [Half-rate factoring composition](half-rate-factoring.md): ordinary proof of the `1/2` machine bound; its final Lean theorem is in [`formal-appendix/HalfTheorem.lean`](../formal-appendix/HalfTheorem.lean).
- [Top constructor and codec](half-top-constructor.md): deterministic table construction and exact input translation.
- [Evaluator and encoding audit](half-evaluator-audit.md): recursive evaluator, computable pad, and charged bit costs.
- [Verified k6 machine](verified-k6.md): Lean proof of the earlier `1/30` bound.
- [Sharper k6 rate](k6-sharp-rate.md): ordinary rate derivation above `1/29` for the k6 program; this sharper exponent is not Lean-certified.
- [k6 source audit](k6-source-audit.md): hashes, modified source scope, and compiler preservation.
- [Baseline audit](baseline-audit.md): source recurrence and the unchanged reference's machine-cost analysis.

## Transfers and proposed stronger evaluators

- [Weighted Weyl and multivariate extraction](weighted-weyl-and-multivariate-extraction.md): an exact occupation crossing, both supplied channel factorizations, and charged multivariate coefficient extraction; the direct profiles do not improve the score.
- [Witt duality and vertex jets](witt-duality-and-vertex-jets.md): an explicit odd Hall determinant and length-unitriangular vacuum maps give effective homogeneous duals at all weights; the crossing and rank-budget obligations remain distinct.

- [Finite-state polynomial crossings](finite-state-pairing-crossings.md): explicit order-dependent identities, exact XOR-circulant rank reduction, and the complete 88,006-model finite screen.
- [Run-compressed crossing ranks](two-state-run-rank-compression.md): exact local Jordan decomposition, interval-type grouping, and quadratic-pivot rank reduction accelerate longer constructive searches.

- [Border and nonunit seed extraction](border-and-nonunit-seed-transfer.md): finite-jet and polynomial identities yield exact seeds with rank-sum loss at most `(kt+1)` after `k` copies; weighted budgets still hold over finite chain rings.

- [Vanishing-slack score transfer](vanishing-slack-score-transfer.md): arbitrary `o(lambda)` width error preserves the leaderboard exponent under stated uniform algorithms and entropy bounds; no improved unconditional exponent.
- [Gap-free finite-seed score transfer](gap-free-seed-score-transfer.md): bounded search attains an enumerable rank-sum infimum, with no seed-gap or generation-speed hypothesis.
- [Input-wise field and seed selection](universal-field-seed-selection.md): one uniform machine attains the finite-certificate endpoint over all coefficient fields, retaining one selected field throughout each input's recursion.
- [Finite-field rank specialization](finite-field-rank-specialization.md): every finite rank profile over an arbitrary field specializes exactly to a finite field; prime-field seeds approximate the same infimum, and ring factor-channel bounds specialize without increasing.
- [Universal verified factoring search](universal-verified-factoring-search.md): the standard verified-search principle already gives one machine attaining the supremum score; direct constructions retain additional resource information.
- [Optimized divisor-count bound](optimized-divisor-count-bound.md): a computable `q_*(n)=(1+o(1))n/log2 n` replaces the earlier auxiliary leading coefficient nine by one.
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

- [Weighted Weyl row filtration](weighted-weyl-row-filtration.md): an actual common module filtration proves a full-row-rank detecting summand at every power-of-two energy, excluding the original coefficient span.
- [Weighted Weyl protected minor](weighted-weyl-protected-minor.md): an exact 1,821-pair witness excludes the entire coefficient and external-symmetry span on the 202 labels of energy 32.
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
