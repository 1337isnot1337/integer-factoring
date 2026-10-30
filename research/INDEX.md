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
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

## Obstructions and work beyond 1/2

- [Directions beyond half](post-half-directions.md): experiments, scope, and unresolved questions.
- [Network tradeoff](half-top-structure.md): restriction for the current local-enumeration and wire-channel constructor.
- [Restricted-label counterexample](half-restricted-labels.md): failure of a proposed low-mode dimension reduction.
- [Restricted intermediate spans](half-restricted-span.md): quantitative obstruction, with precise selected-label conditions.
- [Half-rate recurrence barrier](half-recurrence-barrier.md): limits of the existing operators and recurrence.
- [Graded-space barrier](graded-space-barrier.md): lower bound under the stated graded-space assumptions.
- [Rectangular-grid frontier](rectangular-grid-frontier.md): limits of unequal splits under the declared-profile recurrence.
- [Top-channel rank](top-channel-rank.md): exact factor rank on the full partition domain.
- [Coefficient basis](coefficient-basis.md): scoped minimax calculation for the chosen coefficient span.

## Current post-half campaign

The [frontier ledger](frontier-search-ledger.md) records scope and remaining obligations. These ordinary mathematical results await independent review; no saving exponent above `1/2` is established.

- [Modular recombinations](modular-recombination.md): exact three-label span optimum four and finite structured-label searches.
- [Bounded-weight Pfaffian families](hamming-pfaffian.md), [local tensor powers](hamming-tensor-power.md), and [higher sectors](hamming-sectors.md): effective degree, symmetry, and complement-division arguments.
- [Sparse-width interpolation](sparse-width-interpolation.md) and [the modular half endpoint](modular-half-endpoint.md): effective parameter schedules and a matching half-rate estimate.
- [Rational-function transfer](rational-family-transfer.md): explicit specialization and extension-field coordinate descent.
- [Joint-sector optimum](hamming-joint-optimum.md): exact cutoff 25 and an independently replayed RUP lower-bound certificate.
- [Relabeling tensors](hamming-relabeling.md) and [their exact optimum](hamming-relabeling-optimum.md): a six-piece cutoff-16 seed, a twenty-piece cutoff-four seed, and constant block models.
- [SL2 coefficient ranks](sl2-coefficient-ranks.md): polynomial-size characteristic-two tilting blocks and a tight reciprocal-pencil obstruction.
- [All-field pairing obstruction](all-field-pairing-barrier.md), [characteristic-two small seeds](pairing-small-seeds.md), [full pairing words](pairing-word-barrier.md), and [position symmetry](word-symmetry-barrier.md): restricted-family lower bounds, including arbitrary collective recombinations.
- [Affine-code obstruction](affine-code-barrier.md) and the [finite-core leakage record](frontier-search-ledger.md): limitations of two domain-restriction approaches.
- [Central ranks](modular-central-rank.md) and [particle-hole symmetry](modular-particle-hole.md): exact finite matrix and full coefficient-span identities.

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
