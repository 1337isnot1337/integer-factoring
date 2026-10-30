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

- [Pieces-to-factoring transfer](pieces-to-factoring.md): conditional crossing-parameter transfer, with full bit-cost accounting.
- [Scale-one modular transfer](scale-one-modular-transfer.md): fixed-field transfer and machine costs.
- [Global-block field transfer](global-block-field-transfer.md): varying number fields, raw construction cost `2^O(t^h)`, and reciprocal-logarithmic width error preserve the exact simultaneous time and space exponent.
- [Recovery in varying characteristic](characteristic-independent-recovery.md): short verified divisor descriptors extend the global-block transfer to varying finite fields, including characteristic two.
- [Exact exterior crossing algebra](bi-exterior-crossing-algebra.md): an explicit characteristic-two representation makes both crossing cuts exact square-zero multiplication ranks, with arbitrary generator changes and a directly checked `(36,36)` direction.
- [Quantitative raw rank-sum transfer](rank-sum-global-transfer.md): exponentially many raw pieces are permitted under a total channel-rank cost bound, preserving the exact time and space endpoint through one global block length.
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

## Obstructions and work beyond 1/2

- [Weighted cut-rank budgets](weighted-cut-rank-budget.md): an all-field harmonic lower bound and an exact exclusion of finite rank-sum seeds through 46 labels.
- [Residual cut-rank budgets](weighted-cut-residual-budget.md): a separate exact refinement excludes the same criterion through 50 labels.
- [Positive cut duals](weighted-cut-dual-certificate.md): six exact rational witnesses extend the finite rank-sum exclusion through 52 labels.
- [Two-product complement](two-product-complement.md): over every field, approximating the swap by two product operators leaves rank at least `d(d-1)/2`, including singular pencils.
- [Two-product seed cutoff](two-product-seed-cutoff.md): the stronger complement bound extends the exact finite rank-sum exclusion through 54 labels.
- [Flagged cut budgets](flagged-cut-rank-budget.md): a stronger complement quotient for every number of product operators and class-conditioned weighted budgets extend the all-field exact finite rank-sum exclusion through 82 labels.
- [Universal geometric compression](geometric-subspace-compression.md): explicit unimodular Jacobian minors force common input/output subspaces for every matrix tuple, in every characteristic, including moving-subspace formal limits.
- [Geometric cut budgets](geometric-cut-rank-budget.md): the stronger complement bounds and 167,984 exact rational scalar classes extend the exact finite rank-sum exclusion through 95 labels.
- [Limit of scalar cut budgets](weighted-cut-relaxation-limit.md): abstract dyadic rank profiles satisfy the current scalar constraints with cost `O(d log d)`; a tensor realization is not supplied.
- [Limit of geometric cut budgets](geometric-cut-relaxation-limit.md): the same abstract dyadic profiles satisfy every new geometric subset constraint; a proved upper envelope controls the whole split family.
- [Formal-degeneration stability](weighted-cut-border-stability.md): saturated quotients extend the weighted budgets and the harmonic 41-label exclusion to Laurent-series pieces with poles.
- [Geometric formal cutoff](geometric-border-cutoff.md): moving geometric subspaces extend the formal rank-sum exclusion through 64 labels, without using exact-only small-product bounds.
- [Factor kernels for formal limits](rank-one-formal-cutoff.md): class-dependent rank-one kernel quotients and exact rational enumeration extend that formal exclusion through 71 labels.
- [Quadratic Boolean-code oracle](quadratic-code-oracle.md): polynomial coefficient queries, a rank-six annihilator theorem, and finite raw-profile evidence; no improved contraction is proved.
- [Arbitrary word restrictions](word-restriction-barrier.md): the full 23-term coefficient span on any injective binary-word labels of subexponential length has one piece with both widths `2*lambda-o(lambda)`, independently of the piece count.
- [Dual code-section normal forms](code-section-normal-forms.md): exact reduction of arbitrary ambient lengths and complete finite screens through logical width five for coordinate sections and width three for general linear sections.
- [Common frame filtration](frame-filtration-barrier.md): a simultaneous rank quotient onto the elementary coefficient span gives the sharp `rho >= H_2(1/4)` obstruction for every affine-selector frame and arbitrary coefficient recombinations across frames.
- [Fixed-target Feistel quotient](fixed-target-feistel-barrier.md): arbitrary truth-table controlled translations with a common target set retain the elementary two-cut rank quotient; an exact two-round different-target example breaks that quotient.
- [All diagonal tensor words](diagonal-word-rank-sum-barrier.md): over characteristic two, every piece in the full unsymmetrized diagonal word space has `r*c>=d^2*support`, forcing rank sum at least `d^(3/2)` for any number of pieces.
- [Boolean-incidence support bound](incidence-kernel-rank-sum-barrier.md): over every field, all conservation kernels supported on `w subset S` force `S>=d^(3/2)`; any nonzero unrestricted conservation kernel has cost at least `d`.
- [Exterior degree budgets](bi-exterior-degree-budget.md): arbitrary inhomogeneous pieces in the full exterior model satisfy `S>=d*(lambda+1)`; an additional exact quadratic-rank certificate excludes all its binary seeds through 256 labels, including supplied formal limits.
- [Nonlinear crossing layers](nonlinear-layer-probe.md): an effective Toffoli-based polynomial family escapes the common linear convolution algebra; its complete thirteen-dimensional width-three coefficient span has the same finite optimum as the elementary pencil.
- [Pauli and Walsh dictionaries](pauli-walsh-dictionary.md): an exact support-to-rank interface and positive rank-sum certificates exclude two specified sixty-four-label dictionaries over characteristic-three extensions; unrestricted recombinations remain open.
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
