# Research index

This index covers research notes on `main`. Open PRs can contain later work. Each note states its own assumptions, proof scope, and evidence. An obstruction to a specific construction is not a factoring lower bound.

## Proved machine bounds and their foundations

- [Column-prefix budgets](cnot-cut-budgets.md): exact coefficient cutoffs and an infinite generating identity exclude the zero-basepoint Smith criterion for CNOT row profiles at 16, 32, and 64 labels, including constant basis changes that mix all labels.
- [Exterior majority-state precision witness](exterior-majority-precision-witness.md): an explicit single Cauchy minor, realized by rank-one `0/1` tables, forces `Omega(MD/u)` local precision for the proposed exterior lift. The original Boolean deformation and divisor tensors remain separate questions.
- [Rational cut meshes](rational-cut-mesh-budget.md): complete exact certificates extend the prior finite-seed cutoff from 95 to 96 over all fields, and to 97 in characteristic two. Small polynomial sources and their tensor powers remain separate candidates.
- [Interior Jacobian rank budgets](interior-jacobian-budget.md): explicit full-row-rank binary Jacobians and a proper-projection argument give stronger complementary rank bounds; the initial 96-label characteristic-two certificate is strengthened by the rational-mesh successor.
- [Half-rate factoring composition](half-rate-factoring.md): ordinary proof of the `1/2` machine bound; its final Lean theorem is in [`formal-appendix/HalfTheorem.lean`](../formal-appendix/HalfTheorem.lean).
- [Top constructor and codec](half-top-constructor.md): deterministic table construction and exact input translation.
- [Evaluator and encoding audit](half-evaluator-audit.md): recursive evaluator, computable pad, and charged bit costs.
- [Verified k6 machine](verified-k6.md): Lean proof of the earlier `1/30` bound.
- [Sharper k6 rate](k6-sharp-rate.md): ordinary rate derivation above `1/29` for the k6 program; this sharper exponent is not Lean-certified.
- [k6 source audit](k6-source-audit.md): hashes, modified source scope, and compiler preservation.
- [Baseline audit](baseline-audit.md): source recurrence and the unchanged reference's machine-cost analysis.

## Transfers and proposed stronger evaluators

- [Carry-based swap circuits](carry-swap-circuits.md): all-size direct, inverse, and binary-compressed polynomial crossings, with exact low-jet and finite-basepoint exclusions for the recorded instances.
- [Complete binary cycle-source screens](low-choi-cycle-searches.md): 19,520 quadratic-feedback cycle partitions and 4,285,440 joined-LFSR linear gauges, with independent first-jet controls and optional native reruns.
- [Larger exact divisor checkpoints](large-divisor-checkpoints.md): independent prime-order replays through 32,768 bits, saved complete maxima and witnesses, and a focused large-decimal CLI fix.
- [Factoring from scalar counts in every fixed prime characteristic](fixed-prime-factor-transfer.md): power-root tuples, cyclic factor orbits, and explicit linear-mass arithmetic/comparison networks remove the coefficient-algebra requirement for any fixed prime.
- [Scalar parity factoring transfer](parity-factor-transfer.md): an odd-witness family for every composite, exact proper-factor recovery, and an `O(n)`-width repeated-squaring constructor reduce factoring to polynomially many scalar `F2` clique evaluations.
- [Partial transpose and crossing power ranks](partial-transpose-power-rank.md): a characteristic-free proof of the partial-transpose rank inequality, product bounds for all powers, and a Smith-profile exclusion for every regular rank-two source with at least five labels.
- [Rank-one fixed pairs](rank-one-fixed-pairs.md): an all-field restriction on invertible rank-two channels, plus regular rank-two controlled-cycle examples at every binary label size.
- [Full CNOT tensor coefficient space](cnot-full-tensor-span.md): exact small row blocks, all-field rank-one and column-rank-two classifications, and an all-field exclusion of the complete eight-label tensor span.
- [Three-CNOT row-involutive paths](cnot-crossing-modules.md): an all-size construction, exact local modules, and an all-power obstruction for ordinary univariate coefficient grouping.

- [Rough-input sparse transfer](rough-sparse-transfer.md): deterministic small-factor removal reduces the conditional transfer overhead to `exp(O(sqrt n))`, preserving every fixed logarithmic saving exponent supplied by an algebra-valid weighted evaluator.

- [Arbitrary-bound divisor maxima and logspace membership](certified-divisor-oracle.md): an explicit polynomial oracle for every binary bound, strict highly-composite membership in `L`, certified prefix pruning, and independent exact outputs through 8,192 bits.

- [Explicit adaptive divisor oracle](adaptive-benefit-divisor-preprocessing.md): computed benefit certificates remove unknown algorithmic constants, giving polynomial time and a logarithmic-workspace version; an independent exact search agrees through 2,048 bits.

- [Regular roots with a rank-two terminal power](modified-transpose-root-target.md): an all-size paired-chain construction with the ideal column spectrum, plus an all-function Smith exclusion for full-cycle Feistel maps.
- [Complete normalized rank-two channel form](rank-two-channel-normal-form.md): every invertible two-term matrix pencil is regular; after scalar extension, identity and trace normalization force a commuting-pair form.
- [Geometric transpose channels and row inverses](geometric-transpose-crossings.md): exact row-unimodular instances through sixteen labels, complete finite polynomial-change screens, protected coefficient ranks, and column-degree certificates for constant twists.
- [Polynomial exact divisor preprocessing](polynomial-divisor-preprocessing.md): an unconditional uniform `n^O(1)`-time, `O(log n)`-space corollary of Nicolas's exponent theorem, with exact threshold-center code; fixed global constants are not numerically instantiated.
- [Quantum-spin arithmetic crossing paths](quantum-crossing-initial-forms.md): a rational deformation family, reduction at every prime, independent highest-weight controls, and 422 finite initial-form cases with exact local exclusions.
- [Nilpotent-superoperator crossing families](nilpotent-flag-crossings.md): a basis-free all-field construction with exact Jordan-block column spectra, cyclic-vector invariance, explicit finite controls, and an exact single-transpose row theorem.
- [Exact two-adic Cauchy channels](cauchy-two-adic-smith.md): explicit integral LDU factors, exact Smith exponents, and the remaining square-root mode-enumeration scale.
- [Tapered Hermite precision witness](tapered-hermite-precision-witness.md): one admissible polynomial proves the `MD/u` scale necessary for the source interpolation class and grid; special contraction polynomials remain outside the claim.
- [Arithmetic initial forms of rational crossing paths](arithmetic-crossing-initial-forms.md): exact reductions in every prime characteristic with rational weights and fractional slopes; independent source and Smith-profile controls, with no stronger exponent.
- [Exact divisor maxima by next-prime exchanges](exchange-divisor-preprocessing.md): `n^O(log log n)` time and `O(log n log log n)` workspace with cited uniform arithmetic; a practical oracle and separately bounded elementary streamed implementation.
- [Primorial-weighted pattern enumeration](primorial-weighted-divisor-preprocessing.md): a tight `2^Theta(sqrt(n/log n))` count of all canonical patterns and a smaller-space enumeration; the exchange algorithm skips most patterns.
- [Simultaneous finite-seed endpoint](simultaneous-seed-endpoint.md): one machine attains an enumerable seed family's limiting time and workspace exponent, with no generation-speed or field-uniformity assumption.
- [Pieces-to-factoring transfer](pieces-to-factoring.md): conditional crossing-parameter transfer, with full bit-cost accounting.
- [Scale-one modular transfer](scale-one-modular-transfer.md): fixed-field transfer and machine costs.
- [Sparse divisor nonzeroness](sparse-divisor-nonzeroness.md): detection of divisors in weighted polynomial counts.
- [Modular crossing](modular-crossing.md): normalized finite-field pieces and quotient construction.
- [Modular angle factors](homogeneous-central-modular-symbols.md): Boolean-zeta symbols, a three-block polynomial-time Pfaffian, reciprocal coefficients, and exact parity support.
- [Eighth-order modular residual](homogeneous-central-eighth-jet.md): an all-size even-energy identity, explicit determinant-minor pairing, diagonal state charge, and a sharp weight-32 witness.
- [Complete homogeneous profiles](homogeneous-central-profiles.md): exact larger cut ranks, Taylor coefficients, and independently checked native polynomial Pfaffians.
- [Complement-shell crossings](complement-shell-crossings.md): explicit identity-Gram labels, exact sparse support, a universal four-state tensor, and an all-size paired reciprocal-pencil identity.
- [Complete shell rank blocks](complement-shell-rank-blocks.md): exact syndrome direct sums through 496 labels and protected central minors for fixed-length full shells.
- [Three-parameter Pfaffian crossings](pfaffian-ternary-edge-family.md): arbitrary fixed-length codes, a complete uniform edge space, exact label retracts, and binary/four-element-field seed exclusions.
- [Exact particle-hole reduction](pfaffian-particle-hole.md): a local exterior-power identity preserving all grades and energies, refined signed-charge blocks, and complete profiles through 252 full-grade labels.
- [Stable Pfaffian modules](stable-pfaffian-modules.md): an explicit graded 24-dimensional intertwiner, an all-size two-factor decomposition, and polynomial-size bounds for complete source coefficients and bivariate jets.
- [Complete stable coordinate searches](stable-pfaffian-coordinate-search.md): general pencil and affine-chart evaluators, monomial ancillary shifts, exact small controls, and complete grids in the stated coordinate spaces.
- [Pfaffian boundary sections](noncommutative-crossing-pfaffian-boundaries.md): all-size section and Gram formulas, the universal correction factor, and an all-field rank-cost certificate for four pairs.
- [Sparse isotropic Pfaffian sections](noncommutative-crossing-pfaffian-isotropic.md): explicit endpoint cancellation at every size from six pairs, four-term lifts from seven pairs, exact active correction ranks, and a six-piece crossing curve.
- [Optimal common Pfaffian retracts](noncommutative-crossing-optimal-retracts.md): the largest possible common label dimension for two canceled vertices, attaining triangular maps, and a sparsity tradeoff over every characteristic-two field.
- [Complete retract rank profiles](pfaffian-retract-rank-profiles.md): exact binary and four-element-field exclusions for sparse and triangular maps, invariant-orbit maps, and all 512 orientations of the ten-label example.
- [Sparse symplectic radical codes](noncommutative-crossing-symplectic-radicals.md): closed Gram ranks, polynomial-support coordinate bases, two cancelled vertices, and preservation of the full leading label entropy.
- [Ward factors and overlapping forms](noncommutative-crossing-radical-ward.md): an all-size `(q+r)^k` factor, exact homogeneous product-image dimensions, overlapping quadratic-charge conditions, and the scalar-padding limit of nested pivots.
- [Complete symplectic code profiles](symplectic-boundary-code-profiles.md): nine binary and six four-element-field exclusions, exact coefficient polynomials, and a six-direction family on 26 labels for tensor-power analysis.
- [The 26-label coefficient algebras](symplectic-radical-seed-algebra.md): complete multiplication tables, independently checked radical filtrations and adapted bases, and actual character-minor lower bounds for normalized polynomial coefficient expansions.
- [Exact small tensor modules](symplectic-seed-tensor-modules.md): constant direct sums with largest blocks 22 and 6, an explicit certificate, and complete binary recombination of the common-parameter tensor square.
- [Complete field screens for shared-parameter powers](symplectic-seed-tensor-field-screens.md): powers 2, 3, and 4 over both binary and four-element fields, weighted direct-summand rank bounds, and scalar/rational span preservation.
- [Full weighted seed powers and mixed curves](symplectic-seed-weighted-powers.md): exact modules for all six generators, complete binary square and cube exclusions, and exact recombination for 28 pairs of distinct local curves.
- [All-power common-character theorem](symplectic-seed-character-barrier.md): explicit row minors and a shared 43-dimensional column compression exclude the full 24- and 26-label coefficient spans, all mixed tensor powers, and every characteristic-two extension field.
- [Direct common-functional compressions](common-functional-rank-compressions.md): a reusable subspace criterion gives an all-field, all-power obstruction for the 80-label seed and mixed products of the 24-, 26-, and 80-label spaces, with literal tensor regeneration and explicit constant maps.
- [Two-even-site sections and subcodes](two-even-site-section-search.md): all 728 four-site graph/section choices, independent exact ranks for all 42,536 binary coefficient vectors, and five-site binary and four-element-field profiles.
- [Uniform nonbinary two-even-site section](nonbinary-two-even-site-sections.md): all output rows changed over `F4`, with a complete six-dimensional field scan and independent direct field ranks for every affordable vector.
- [All-size two-even-site core factorization](two-even-site-core-factorization.md): at most four exceptional sites per entry, at most `2P-2` polynomial directions, and an exact 240-label evaluator with a complete ten-dimensional binary screen.
- [All-power 240-label compression theorem](two-even-site-P6-character-barrier.md): explicit constant maps of dimensions 15,688 and 718 exclude the full six-site coefficient span and all mixed products with the fixed 24-, 26-, and 80-label spaces over every characteristic-two field.
- [All-power 672-label Jordan certificate](two-even-site-P7-jordan-barrier.md): bounded core restrictions, a triangular row map of size 108,315, and a column compression of size 1,150 give exponent `10/7` for the seven-site span and all mixed products with the four smaller certified spaces.
- [All-power 1,792-label extension](two-even-site-P8-character-barrier.md): an inherited triangular row map and new bounded column modules give dimensions 769,060 and 2,106, excluding the eight-site span and mixed products with exponent `99/70`.
- [Small row corners and the 4,608-label frontier](two-even-site-row-corner-barrier.md): a reusable joint-nilpotence compression criterion strengthens eight sites, closes nine sites with exponent `17/12`, and gives exact all-size certificate formulas with an explicit ten-site shortfall.
- [Graded Ward rank descent](noncommutative-crossing-ward-rank-descent.md): homogeneous lower-grade preimages, a row factorization without channel multiplicity, and a separate column channel bound.
- [Computed Ward bounds and overlapping codes](ward-code-rank-bounds.md): independently audited rectangular stable matrices, full and orientation-half code bounds, and complete binary searches for two actual overlapping-form constructions.
- [Euler refinement and the determinant interface](noncommutative-crossing-ward-euler.md): two-rectangle bounds for a global form, a small determinant for several binary forms, and a joint vanishing-order theorem for syndrome-preserving outputs.
- [Dual incidence and charge filtration](noncommutative-crossing-ward-charge-filtration.md): exact both-cut descent without channel multiplicities for arbitrary characteristic-two form weights; a separate binary gauge proof gives interface entries as single minors.
- [Computed refined Ward bounds](ward-charge-bound-computations.md): exact integer costs, a large-code improvement to `rho = 0.83464`, literal small tensors, and independent determinant controls.
- [Ward code label entropy](ward-code-charge-entropy.md): exact total dimensions, the quadratic-charge penalty, a characterization of disjoint supports, and positive-parameter bounds for each homogeneous grade.
- [Computed Ward-bound envelope](ward-rectangular-bound-envelope.md): protected rank-one modules cover all implemented minima and force limiting `rho >= 29/40` in the stated fixed-ratio regime; actual restricted ranks remain open.
- [Stable Ward operator interface](ward-stable-operator-interface.md): explicit incidence and creation transport, an eight-dimensional ancillary algebra, and a bond-nine recurrence for one global form at every size.
- [Several-form operator transfer](ward-form-operator-transfer.md): arbitrary characteristic-two weights use the same eight ancillary channels over a square-free form ring; the separate exponential form-selection memory is retained explicitly.
- [Disjoint-support operator transfer](ward-disjoint-operator-transfer.md): width at most 16 for any number of disjoint forms, preserved global grades, and an exact streamed entry evaluator with linear field work and constant field-element workspace.
- [Overlapping-support operator transfer](ward-interval-operator-transfer.md): width `8*2^w`, an exact entry algorithm, and a row-basis reduction attaining the optimum active width for every fixed site order.
- [A width-two overlapping fan family](ward-fan-family.md): explicit forms at every size, exact graph and homogeneous-dimension recurrences, operator width 32, and a scoped first finite upper profile.
- [Polynomial-time physical Ward entry queries](ward-pfaffian-entry-queries.md): one streamed Pfaffian handles arbitrary overlapping forms and zero local entries in `O(P*k^2+k^3)` field work and `O(k^2)` field-element workspace, with an effective formal coefficient-query corollary.
- [Combined Ward operator ranks](ward-combined-operator-ranks.md): exact finite type decompositions combine derivative and creation before taking ranks, improving four upper bounds while preserving full literal controls.
- [Cubic radical output-map theorem](cubic-radical-output-maps.md): exact block formulas over every characteristic-two field, a rank-cost charge for all syndrome-preserving outputs, and cyclic/circulant computation records.
- [Untouched-output gate and perturbations](cubic-output-perturbations.md): an all-field gate based on unchanged logical labels and a complete 1,960-case search outside the original syndrome-preserving family.
- [Arbitrary dual cubic outputs](cubic-arbitrary-output-identity.md): an all-odd-size four-direction identity under the stated dual-map constraints, a universal two-piece upper bound, 680 multirow experiments, and an individually affordable tensor that does not complete a strict seed.
- [Nonzero complementary cubic outputs](cubic-complementary-output-maps.md): a universal five-direction identity, a three-pairing cost bound, and the exact seven-label `F4` optimum 36, with complete rank tables for 1,251 recorded maps.
- [Non-diagonal cubic outputs over F8](cubic-complement-F8-probe.md): 428 explicit maps and 2,003,468 projective coefficient vectors, with independent field ranks and direct physical coefficient controls; every recorded space misses the strict seed gate.
- [Full mixed cubic spaces over F8](cubic-mixed-F8-spaces.md): five structured map pairs with complete combined coefficient-space searches, 5,692,269 projective vectors, and an exact scalar-parameter span reduction.
- [Cyclic cubic all-field and all-power theorem](cubic-cyclic-all-power-barrier.md): exact Fourier blocks, a 14-dimensional row compression, and a multiaffine grid bound give cost at least `336^(n/2)` in the full scalar cyclic-output tensor space.
- [Jordan cubic all-field and all-power theorem](cubic-jordan-all-power-barrier.md): a 14-dimensional row compression and 28-dimensional triangular column restriction give cost at least `392^(n/2)`; an exact certificate separately rules out a full-leg weight degeneration.
- [Syndrome sumset rank-cost theorem](syndrome-sumset-rank-cost.md): both-cut block support gives `S >= |Sigma+Sigma|`, all tensor powers and same-support mixtures, including an all-field exclusion of the seven-label diagonal cubic family.
- [Mixing affordable cubic output pieces](cubic-mixed-output-pool.md): a single tensor entry separates the crossing from all 659 affordable binary tensors in the recorded 680-map pool; arbitrary new cross-space combinations remain open.
- [Mixed-grade endpoint contacts and repairs](symplectic-mixed-grade-codes.md): an explicit failure of the raw crossing endpoint, two exact repaired 28-label constructions, and complete finite field screens.
- [Additional control and bipartite screens](additional-crossing-screens.md): all 40320 D8 control permutations at three weight choices, different operator/state occupation codes, and complete small mixed-radix F4 spaces, with independent winning-profile audits.
- [High-mode residue folding](residue-folded-channels.md): exact ambient dimensions, finite incoming and outgoing factors, both complete channel counts, and eventual energy periodicity.
- [Complete folded Gram channels](folded-gram-channels.md): all-angle cancellation inside one Gram form and a constructive polynomial-time exact Johnson-rank algorithm.
- [Weighted Gram aggregation](folded-gram-aggregation.md): exact energy-preserving multiplicities, Gaussian-binomial evaluation, and 542 complete finite energy-bound calculations.
- [Unequal-grade folded Gram channels](unequal-grade-folded-gram.md): exact creator-annihilator offsets, odd-weight Frobenius slices, energy-preserving aggregation, and 217 complete finite cases.
- [Homogeneous modular selection](homogeneous-central-selection.md): one operator weight at every label width, selected rectangular vacuum minors, and an effective uniform half-rate construction.
- [Noncommutative group-word crossings](noncommutative-crossing-orbits.md): exact all-size construction, compressed constant cut, optimal rank-28 central shears on 64 labels, and complete-profile limitations.
- [Coded central quotients](coded-central-crossings.md): binary-code group construction, exact constant-cut compression, matroid bounds, and an all-coefficient row filtration estimate.
- [Graded Hopf crossings](noncommutative-crossing-graded-hopf.md): exact noncommutative families, homogeneous row bounds, compressed translation cuts, and complete finite recombination and gauge audits.
- [Cumulative bosonic crossing profile](bosonic-hermite-profile.md): exact moment identities, a Virasoro rank formula, and complete pieces with independently collapsed prefix and suffix channels.
- [Homogeneous duals](homogeneous-duals.md): alternative algebraic pieces and their checked identities.
- [Near-ones label families](near-ones-labels.md): label construction and rank experiments.
- [Binary-prefix top channel](binary-prefix-top-channel.md): structural analysis of channel widths.
- [Primary central Gram matrix](primary-central-gram.md): source-linked rank calculations.
- [Half-piece ranks](half-piece-rank.md): finite crossing-piece rank evidence.

## Obstructions and work beyond 1/2

- [Elementary mixed-cut budget](strict-mixed-cut-gate.md): a short all-field strict budget gives an independent cutoff through twelve labels; the [prior PR #27 geometric theorem](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md) has the stronger exact cutoff through 95.
- [Sixteen-label CNOT companion](cnot-sixteen-label-span.md): the mixed budget and rank-one classification give a short all-field certificate for the full fourth coefficient tensor power.
- [Universal finite-size gates](noncommutative-crossing-finite-gates.md): all-field rank-sum bounds exclude every seed on at most five labels and every non-elementary group orbit family of order at most eight.
- [Abelian midpoint blocks](noncommutative-crossing-abelian-midpoints.md): exact doubling-map column multiplicities, a real rank-floor lemma, and an all-field exclusion of every cyclic group.
- [Involution-free orbit minors](noncommutative-crossing-involution-support.md): a diagonal-minor theorem for arbitrary group kernels, excluding odd-order groups, odd-order index-two extensions, and unique-involution groups over every field.
- [Uniform six-pair boundary barrier](pfaffian-p6-boundary-orbit-barrier.md): the common active quotient and an invariant rank shift give an exact lower cost 834 for all cancelling sections with binary coefficients.
- [Stable-bound computation limitation](stable-boundary-envelope.md): an all-six-ancilla rank-one block and elementary fusion paths force computed `rho >= 29/40` for every fixed dense grade ratio; this is not a lower bound on actual ranks.
- [Directions beyond half](post-half-directions.md): experiments, scope, and unresolved questions.
- [Quadratic charges and label entropy](binary-quadratic-charge-shells.md): exact residue-syndrome characters and dimension bounds for joint eigenspaces of all low quadratic charges.
- [PBW degeneration of Hopf crossings](noncommutative-crossing-pbw-degeneration.md): a common four-leg specialization proves a finite rank-sum obstruction for the full PBW-domain coefficient spaces, including field extensions and filtered local gauges.
- [Arbitrary-central-control degeneration](noncommutative-crossing-group-degeneration.md): a common heavier central weight excludes full-domain class-two word crossings with every central lookup table, including mixed families; it also explains the scope limit of the earlier weight-two proof.
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
