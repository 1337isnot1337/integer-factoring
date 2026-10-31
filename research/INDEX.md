# Research index

Latest: [31 October ordinary-proof research session](frontier-session-20261031.md) — complete score unchanged at 1/2; exact constructors, global charge obstructions, completed finite searches, and final validation.

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

- [Dyadic inverse sections](dyadic-inverse-sections.md): exact arithmetic interval recursion, affine-fibre floor sums, integer-factor recovery, and the unresolved aggregation cost.
- [Dyadic orbit aggregation](dyadic-orbit-aggregation.md): polynomial bit time for complete fractional-linear cyclic families with interval selectors and arbitrary affine changes.
- [Interval orbit gcd](interval-orbit-gcd.md): exact correlated-root interface, two boundary recurrences, output-sensitive algorithm target, and source audit; its stronger factoring consequence is conditional.
- [Interval unit jets](interval-unit-jets.md): polynomial-bit-time local Taylor data for exponentially long interval products using parity, power sums, and exact Newton identities; global elimination remains open.
- [Interval packet circuits](interval-packet-circuits.md): global modular polynomial circuits with `O(n+p)` nodes and polynomial generation cost; exact artificial-common-factor controls state the remaining elimination obligation.
- [Normalized interval jets](normalized-interval-jets.md): exact binary content exponents, zero orders, and polynomial-cost unit jets at scaled integer centers; global root search remains open.
- [Arithmetic framework shift](arithmetic-framework.md): primary-source randomness audit and a new generator/relation construction contract, outside the crossing evaluator.
- [Deterministic torsion extraction](arithmetic-torsion-extraction.md): a black-box two-group basis algorithm from an annihilator, with no saturated relation lattice or random group elements.
- [Positive relation kernel](positive-relation-kernel.md): a Fourier-positive smooth-pair hitting lemma and diagonally dominant relation construction; efficient arithmetic pair extraction is open.
- [Biased words without coordinates](arithmetic-biased-words.md): a powering seed and published derandomized squaring supply biased samples from generators and a group-size upper bound.
- [Actual class-group controls](arithmetic-class-group-controls.md): exact composition, smooth prime-form words, unsaturated relations and complete two-torsion recovery; exhaustive searches produce factors for 18 finite targets.
- [Density/index transfer](arithmetic-density-lift.md): dense smooth forms and relations on a possibly proper factor-base subgroup give a global class-number multiple; missing ambiguous directions remain an explicit obligation.
- [Bounded smooth-norm interface](bounded-smooth-norm.md): exact equivalence between smooth-form quotients and balanced smooth represented norms; supplied witnesses convert at `2^O(b/log b)` cost, while fast witness search remains open.
- [Deterministic norm sieve](arithmetic-norm-sieve.md): dominant relations and exact determinant/torsion extraction factor 60-bit and 80-bit inputs without class enumeration; two smaller-box exhaustion controls retain the finite search limits.
- [Quantitative smooth-norm localization](smooth-norm-localization.md): the exact biased generator-power family, an explicit unproved short-box hypothesis and a conditional transfer above half; splitting coverage remains a separate hypothesis.
- [Polynomial divisor-gap extraction](smooth-norm-gap.md): integer multiplicative trimming converts promised inner-radius balanced norms at a specified outer radius, without enumerating all divisors.
- [Relations from small smooth norms](smooth-norm-relations.md): primitive smooth norms yield bounded prime-form words directly, removing the balanced-factor condition and strictly enlarging the accepted target set.
- [Exact search for a factored norm](factored-norm-search.md): ideal-lattice composition and deterministic meet-in-the-middle decide one primitive representation and recover coordinates in `Y poly(b)+2^O(b/log b)` bit time; selecting candidate norms remains open.
- [Reflected duals and a neutral crossing compiler](neutral-reflected-crossings.md): Lucas duality in every prime characteristic and a characteristic-two Pfaffian formula give complete balanced and unbalanced fixed-sector crossings.
- [All homogeneous dual completions](homogeneous-dual-completion.md): the complete right-inverse family at weight 32 has 120 dual parameters and exactly 88 independent tensor directions, with exact source and endpoint controls.
- [Crossing synthesis from arbitrary F4 frames](f4-recentered-crossing-lattices.md): full-entry affine coefficient search can recenter on nonbinary tensors; finite polynomial unit frames certify complete orders, with independent binary-jet replays.
- [Odd-cyclic matrix crossings](odd-cyclic-matrix-crossings.md): explicit matrix/Fourier cut normal forms, an all-degree unweighted `rho >= 3/4` floor, and a completed five-symbol weighted search with recorded numerical failures.
