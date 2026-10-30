# Integer factoring bounds

[![Verify research certificates](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml/badge.svg)](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-59636e)](LICENSE)

An open research effort on deterministic integer factoring, targeting the
[Big O Speedrun problem](https://bigospeedrun.com/integer-factoring). For an
`n`-bit integer `N >= 2`, the task is to output its smallest prime factor. The time
bound counts finite-alphabet Turing machine steps for every input.

## Ordinary proof: exponent 1/2

The [written composition proof](research/half-rate-factoring.md) gives one
deterministic machine with

```math
T(n)\le C\exp\!\left(\frac{Cn}{\sqrt{\log n}}\right).
```

The [modular time–space improvement](research/rough-input-transfer.md)
now gives this same time bound with **polynomial space**, by removing
small factors first and evaluating the grid through sequential table
queries. The [local factor proof](research/modular-factor-space.md)
uses streamed permanents and polynomial-size coefficient arrays. Both
the repeated-query cost and the exact half-rate endpoint are included.
This is an ordinary mathematical proof with finite checks, without new
Lean formalization.

The [explicit interpolation constructor](research/polynomial-top-queries.md)
now computes every top-table entry in polynomial time and space. Each
prescribed boundary admits at most one internal setting. An explicit
grouping replaces the initial grouping search and allows the larger
top schedule proved in that note. The tables remain available through
queries; explicitly writing all their entries still takes exponential time.

It combines Claude (Anthropic)'s existing Appendix B clique evaluator with
Isaac Rudich and Louis-Martin Rousseau's multiplication-network reduction.
The new connection includes a deterministic written-table constructor, an
exact value-preserving codec, a computable runtime pad, and all input,
arithmetic, recursion, and prefix-search costs. It gives score `alpha=1/2`
under the published rule and the rate of the reference's square-root-logarithm
conjecture.

**Evidence level: complete ordinary mathematical argument, with source audit
and finite checks.** The composition has not been formalized in Lean or
reviewed independently outside this research session. No leaderboard
acceptance or practical factoring speedup is claimed. The prior source
theorems supply the evaluator and reduction; this repository does not claim
to have discovered those constructions.

[Complete proof](research/half-rate-factoring.md) ·
[Top constructor and codec](research/half-top-constructor.md) ·
[Evaluator and encoding audit](research/half-evaluator-audit.md) ·
[Exact arithmetic checks](evidence/half-rate-arithmetic.json)

Run `make verify` for the finite checks, including the new table adapter,
parameter bounds, computable pad, and prefix driver. These tests do not
materialize the enormous evaluator or replace the infinite-input proof.

## Research beyond exponent 1/2

The [directions note](research/post-half-directions.md) records the initial
search beyond `1/2`. It includes restrictions on the operator interfaces,
an exact counterexample to projecting away high modes, and finite matrix
ranks. The results below further specify the mathematical frontier.

The [network tradeoff](research/half-top-structure.md) rules out improving
the explicit constructor through regrouping and separate wire compression
alone. The [restricted-label counterexample](research/half-restricted-labels.md)
and [quantitative span theorem](research/half-restricted-span.md) show why
another proposed dimension reduction fails on the full restricted spaces.

The [Fourier gate construction](research/fourier-multiplication-gates.md)
gives polynomial-time point queries for every original multiplication
gate after `exp(O(sqrt(n)))` character-field preprocessing. It includes
prefix weights, composite-modulus bilinear sums, and all carry endpoints.
A characteristic-two involution also reduces unit Kloosterman sums to
quadratic roots. These are explicit local contraction primitives; a
global evaluation method with a better rate is still missing.

The [general transfer theorem](research/pieces-to-factoring.md) proves that
an effective crossing family with rational parameter `rho` gives factoring
exponent `log2(1/rho)`. It includes a computable pad and the complete bit-cost
connection. For example, a family with `rho=7/10` would give exponent
`0.514573...`. This is a conditional target; that family is not constructed.

The [fixed-characteristic transfer](research/scale-one-modular-transfer.md)
opens another route. An effective scale-one crossing family over one fixed
finite field gives the same factoring exponent `log2(1/rho)`, up to and
including `alpha=1`. Its original [sparse divisor-polynomial argument](research/sparse-divisor-nonzeroness.md)
detects permitted divisors using cyclic polynomial rings, with additional
cost `exp(O(n/log n))`. The proof includes weighted input construction,
circuit generation, ring arithmetic, and finite-alphabet machine costs.

The [small-factor preprocessing extension](research/rough-input-transfer.md)
reduces that extra time to `exp(O(sqrt(n)))`. If the weighted arithmetic
evaluator costs `A(n)`, the resulting factoring bound is
`exp(O(log(A(n)+2)+sqrt(n)))`; any fixed logarithmic saving exponent
is preserved. A parameter tradeoff also gives polynomial-space recovery
for every exponent at most one, provided the evaluator has sequential
queries and polynomial-space local routines. The existing modular
half-rate construction meets those requirements.

The stronger [bit-tagged recovery theorem](research/bit-tagged-witness-recovery.md)
now achieves `exp(O(sqrt(n)))` recovery overhead with polynomial space
over **every fixed characteristic**, preserving any fixed logarithmic
saving exponent of a streamed evaluator. A different control field
provides deterministic isolation and additive-character projections;
square-zero tags return the isolated witness's bits. The proof includes
candidate verification, complete factorization, and uniform field
construction. [Streamed local factors](research/modular-factor-space-odd.md)
are also available in every fixed characteristic.

The [extended transfer](research/subexponential-piece-transfer.md) permits
more than polynomially many crossing pieces. Piece entropy
`log |I_lambda| = O(lambda/(log lambda)^delta)` gives factoring exponent
`min(log2(1/rho),delta,1)`, without an extra depth factor at the endpoint.
In particular, `2^(O(lambda^beta))` pieces, for fixed `beta<1`, preserve
the rate-derived exponent. The width sum can have error as large as
`O(lambda/loglog lambda)`. An explicit volume cutoff, a geometric cost
sum, exact recurrence checks, and independent internal review support
this conditional theorem. It supplies no new family with a smaller rate.

The [rational-height extension](research/rational-height-transfer.md)
also removes the old polynomial coefficient-height requirement. A uniform
`2^O(lambda)` bit-time routine can output rational factors with exponentially
many bits. The evaluator's degree is at most `P^6`, so exact rational
arithmetic still fits the rate. As a consequence, uniformly computable
piece entries and a proved asymptotic rank profile suffice: exact Gaussian
elimination constructs the factors within the allowed cost. This result
does not replace the missing rank profile with finite rank experiments.

The [finite-order transfer](research/finite-order-crossing-transfer.md)
also admits algebraic approximations to the crossing. If the crossing is
the coefficient of order `kappa`, extracting that coefficient multiplies
each generic matrix rank by at most `kappa+1`. Approximation orders with
`log(kappa+1)=O(lambda/loglog lambda)` preserve the leading width rate.
The proof includes regular factorization over formal power series and
exact effectivity after coefficient extraction. It supplies a target for
degeneration constructions; no better-rate approximation is asserted.

A [normalized modular construction](research/modular-crossing.md) supplies
such a family for every fixed `1/sqrt(2)<rho<1`. It removes the vanishing
factorial multiplier and permits rings with zero divisors. The endpoint
estimates in the time-space notes give the exact half exponent with
sublinear width slack. It supplies no family
with `rho<1/sqrt(2)`. A local one-bit improvement also fails to tensor into
a better asymptotic profile, as the note's exact rank argument shows.

The [coordinate-support theorem](research/coordinate-support-rank.md)
extends that last restriction to arbitrary binary label subsets, including
nonlinear codes. A tensor with the stated two-pairing support satisfies
`rank_row(B) rank_col(B) >= d(B)^2/m`, where `d(B)` counts its nonzero
crossing-diagonal entries. Any `M` such pieces that sum to the crossing
therefore have maximum rank product at least `m^3/M^2`. Subexponentially
many pieces still require `rho>=3/4`. The proof covers different fixed
orientation masks in different pieces. Combining masks within a single
piece can escape the support condition and remains a separate route.

On the full binary cube, the [exact coefficient-rank formula](research/binary-coefficient-entropy.md)
gives the elementary family the limiting rate `H_2(1/4)=0.811278...`.
The [symmetric three-pairing theorem](research/symmetric-three-matching-barrier.md)
extends this sharp lower bound to every decomposition in the complete
permutation-symmetric span of the three local pairings. It permits
arbitrary recombinations, any number of pieces, and controlled
degenerations within that span. Coordinate-dependent constructions and
restricted label sets are outside this theorem's scope.

The [reflected monomial construction](research/fermionic-piece-polynomial.md)
gives exact unit duals by complementing binary mode indices. Its dual
operators may have larger weight than the labels; the note proves the
corresponding locality and finite support bounds. On the three labels
`(1,11),(3,9),(5,7)`, four nonzero pieces have rank products
`9,8,20,12`. The maximum `20` improves the previous finite certificate's
`27` on that domain. An independent polynomial matching formula and
both original operator orders verify every entry. This finite gain
does not establish a better asymptotic factoring exponent. Separate
recombination work in [PR #24](https://github.com/1337isnot1337/integer-factoring/pull/24)
has since produced a rank-product-nine certificate on three labels;
the reflected construction here is a different family, not the current
best finite certificate.

Two [operator restrictions](research/half-recurrence-barrier.md) rule out
direct shortcuts. The current operators permit only linear changes of the
mode weights. Deleting higher modes breaks their identities, as exact
negative controls demonstrate. The existing coefficient recurrence also
has exactly the stated square-root decay.

The [graded-space rank theorem](research/graded-space-barrier.md) gives a
broader restriction. A polynomial number of crossing pieces cannot use
slower-growing spaces within the same energy range, even with new operators.
The proof uses the full rank of the crossing matrix. It identifies which
assumptions a replacement construction must change.

The [rectangular-grid theorem](research/rectangular-grid-frontier.md)
retains the geometry of every finite grid. It proves that unequal, adaptive
splits cannot improve exponent `1/2` within the declared-profile recursion
model. An analytic tail bound and 3,844 exact interval cases support it.

The [highest-channel rank theorem](research/top-channel-rank.md) proves a
sharp asymptotic rank for one factor on the full partition domain.
Its [binary-label extension](research/binary-prefix-top-channel.md) proves
`log2 rank(U_lambda)=sqrt(2) lambda+O(log lambda)` on the actual label set.
Consequently, separate lossless compression of the existing factor vectors
cannot improve their leading width sum. Joint refactorization remains open.
The [fixed-defect extension](research/near-ones-labels.md) gives the same
restriction for the natural family of partitions padded with many ones.
Its exact support calculation also shows why annihilation terms prevent
a direct cutoff by partition defect.
An [independent finite checker](research/half-piece-rank.md) also checks
the two operator orders, crossing sums, and small binary-label restrictions.
The theorem does not settle the rank of the complete central piece.
These results give no lower bound on factoring or on an evaluator's runtime.

The [integral Schur construction](research/integral-schur-crossing.md)
retains the full `p(h)`-dimensional state space in every characteristic.
Translated Schur functions give multiplier-one crossings over the integers.
At prime-power energies, the homogeneous vacuum pairing is invertible
modulo that prime, and its angle sum has only two or three terms.
A constructive capacity theorem fills every binary label length with
`O(log lambda)` prime-power blocks. It retains every block index, giving
piece entropy and width error `O(log^2 lambda)` without an energy overshoot.
This supplies another exact half-rate family and makes a rank improvement
on prime-power blocks alone sufficient for an improved exponent.

The [homogeneous-dual construction](research/homogeneous-duals.md) gives a
new effective crossing family. An explicit polynomial in two sparse
derivations replaces a large matrix inverse. The proof bounds coefficient
sizes and generation costs, and verifies the complete crossing interface.
It retains `rho=1/sqrt(2)` and exponent `1/2`; smaller ranks for its complete
pieces remain an open route. Its [exact checks](evidence/homogeneous-duals.json)
cover the dual pairing, operator identities, and complete crossing sums.
On primitive inputs, a [central-piece rank identity](research/primary-central-gram.md)
goes further: an invertible row swap gives a positive Gram form, so
alternating terms cannot cancel its primitive fusion image. The size of
that image still needs an asymptotic bound. This restricted identity does
not establish optimality of the full central piece.

Run `make verify-frontier` for these exact checks. The strongest complete
ordinary factoring argument here remains exponent `1/2`.
The [research verification record](evidence/beyond-half-verification.json)
lists the proof scope, reviewed claims, checks, and artifact hashes.

## Lean theorem: exponent 1/30

**The k6 variant has a complete Lean proof of the `1/30` bound:**

```math
\exists M\;\exists C>0\;\forall n\ge2,\qquad
T_M(n)\le C\exp\!\left(\frac{Cn}{(\log n)^{1/30}}\right).
```

One deterministic machine computes `Nat.minFac` on every natural-number input
in the reference's binary encoding. It retains `Nat.minFac 0 = 2` and
`Nat.minFac 1 = 1`.
The proof covers the modified program, its complete machine cost, and the
real-logarithm rate. Its final declarations use only `propext`,
`Classical.choice`, and `Quot.sound`.

The exponent `1/30 = 0.033333...` exceeds the explicit reference exponent
`log(12/11)/log(16) = 0.0313827205...`. Lean verifies this strict inequality.
The comparison concerns proved upper bounds. It establishes no lower bound
on the reference machine's actual runtime.

[Verified result and proof roadmap](research/verified-k6.md) ·
[Final Lean theorems](variants/k6/FactoringK6/Certified.lean) ·
[Verification record](evidence/k6-verification.json) ·
[Source audit](research/k6-source-audit.md)

**Contributing research or asking an agent to help?** Start with
[autoresearch.md](autoresearch.md), [AGENTS.md](AGENTS.md), and
[CONTRIBUTING.md](CONTRIBUTING.md). Check the
[coordination board](https://github.com/1337isnot1337/integer-factoring/issues/8)
and [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls)
for overlap; post your scope and handoff on the board.

Search the live board and open work alongside this checkout's research notes:

```sh
python3 tools/find_overlap.py grouping crossing
```

This read-only helper requires GitHub CLI access. Results are leads, not a
reservation or a novelty check; it reports which sources it did not search.

**Keep this README current.** Changes to results, evidence, verification, or
open obligations must update the relevant status, links, and reproduction
instructions here in the same change. Check the underlying evidence before
promoting a claim.

## Verified result: smaller grouping schedule

The variant replaces the reference's schedule by

```math
K_0=48,\qquad G_d=5K_d^6,\qquad K_{d+1}=25K_d^{12}.
```

The reference uses `G_d = 3 K_d^8` and `K_(d+1) = 9 K_d^16`.
A stronger grouping lemma supports the new schedule with the same label-width
guarantees. The patch changes the actual uniform program and its correctness
and cost proofs. The compiler, simulator, machine definitions, and target
function retain their reference definitions.

The quantitative proof covers all recursion levels. It uses contraction
`23/25` at positive heights and handles the bottom transition separately.
The exact inequality `12 * 23^30 < 25^30` supplies the saving exponent `1/30`.
An integer-root theorem covers every input length. A second theorem converts
that bound to the real-logarithm formula above.

The complete artifact consists of the [source patch](variants/k6/reference.patch)
and eight [Lean modules](variants/k6/FactoringK6). The
[reproduction script](tools/check_k6.sh) applies them to a fresh archive of the
pinned reference. It rebuilds the modified `Solution` and every new module.
Before compilation, it verifies all patched source hashes against the audit.
It then verifies the final theorem's axioms in a fresh Lean invocation.

## Sharper ordinary analysis of the k6 machine

The [sharper analysis of the existing k6 machine](research/k6-sharp-rate.md)
proves the ordinary mathematical bound with exponent
`beta_12 = log(12/11)/log12 = 0.0350159540186562...`, which is strictly
larger than `1/29` and `1/30`. It retains the already certified program and
machine bridge, and sums the per-level contraction errors instead of replacing
them by a fixed slack. The new rate argument is handwritten; the Lean endpoint
remains `1/30`. Run `python3 tools/check_k6_sharp_rate.py` for its finite exact
arithmetic checks.

The separate `1/2` composition above uses a different evaluator. It does not
change the k6 program or extend its Lean theorem.

## Earlier candidate: grouping and crossing compression

The proposed combination of a stronger grouping construction and a 56-channel
crossing factorization targets

```math
T(n)=\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right),
\qquad \alpha=\frac{\log(1/\rho)}{\log 4}=0.0755104681\ldots,
\quad \rho=\frac{8+\log_2 7}{12}.
```

Larger `alpha` is better. This earlier route improves on k6 but is weaker
than the ordinary `1/2` theorem above, which bypasses this proposed evaluator.

**Status of this earlier target: research candidate.** The
notes describe a uniform counting construction and its table-operation costs.
The generalized evaluator and complete Lean proof of this stronger factoring
machine remain open. Exact finite certificates and
checked local lemmas establish only their stated scope. These open obligations
do not apply to the completed k6 theorem with exponent `1/30`.

[Integration and open obligations](research/integration.md) ·
[Recurrence and costs](research/recurrence.md) ·
[Grouping proof](research/grouping-improvement.md) ·
[Crossing construction](research/crossing-compression.md) ·
[Formal verification evidence](evidence/formal-verification.json)

## Baseline

The starting point is Isaac Rudich and Louis-Martin Rousseau's
[SubexponentialDeterministicFactoring](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring),
pinned to commit `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
Its published Lean theorem establishes one deterministic machine with time
`2^{o(n)} = N^{o(1)}`. Its separate square-root-logarithm conjecture is
unproved in that source. The ordinary composition above supplies that rate.

## Research results

The repository contains a written `1/2` composition theorem, the complete
k6 Lean theorem, exact finite checks, and earlier research. The table
identifies the evidence for each result.

| Result | Value | Evidence and scope |
| --- | --- | --- |
| Appendix B factoring composition | `alpha = 1/2` | [Complete ordinary proof](research/half-rate-factoring.md), source and interface audit, finite adapter and arithmetic checks; no new Lean formalization |
| **Complete k6 factoring theorem** | **`alpha = 1/30`** | **Lean verifies one modified machine, all runtime costs, the real-logarithm conversion, and the strict exponent comparison** |
| Sharper rate for the same k6 machine | `beta_12 = log(12/11)/log12 = 0.0350159540... > 1/29` | [Complete ordinary derivation](research/k6-sharp-rate.md), using the existing machine bridge; new rate not Lean formalized |
| Explicit analysis of the unchanged baseline | `alpha = log(12/11)/log(16) = 0.0313827205...` | Handwritten derivation from the pinned source's inequalities |
| Explicit time bound for the reference machine | `C * 2^(C * (explicitZu(n) + 1))` | Lean-checked bound for one finite-alphabet factoring machine; logarithmic-power conversion remains handwritten |
| Grouping with independent accuracy | `g = O(D^2 K^2 log K)` | Matching decomposition, Holder's inequality, and an explicit deterministic search |
| Formally checked grouping improvement | `g > K^2 + 4D^2K^4`; original accuracy works at `g = 5K^6` | Lean kernel checked, including network value preservation |
| Smaller crossing factorization | 56 integer channels instead of a trivial 64-channel factor | Explicit integer factors and exhaustive finite verification |
| Grouping-only rate target | `alpha = log(12/11)/log(4) = 0.0627654410...` | New parameter schedule and runtime analysis; modified machine not Lean verified |
| Combined rate target | `alpha = log(1/rho)/log(4) = 0.0755104681...`, `rho = (8+log_2(7))/12` | Combined mathematical construction; full formal integration remains open |

These bounds have very large constants; the finite checkers are not practical
factoring implementations or runtime benchmarks.

The combined target uses exact integer widths: `ceil(log_2(7^i))` is
`(7^i - 1).bit_length()`. A simpler rational width formula gives the slightly
weaker target `alpha = log(312/281)/log(4) = 0.0754879493...`. Neither formula
requires real-number approximation in the algorithm.

## Proofs and certificates

- [Half-rate factoring theorem](research/half-rate-factoring.md): one uniform
  machine, the complete written composition, all costs, and the score implication.
- [Top constructor](research/half-top-constructor.md) and
  [evaluator audit](research/half-evaluator-audit.md): source interfaces,
  entry bounds, exact codec, computable pad, and actual encodings.
- [Half-rate source manifest](evidence/half-rate-source-audit.json): pinned
  source files inspected for the new composition, with Git and SHA-256 hashes.
- [Half-rate arithmetic checks](evidence/half-rate-arithmetic.json) and
  [table-adapter tests](tests/test_half_clique_codec.py): finite checks and
  negative controls, separate from the infinite-input proof.
- [Directions beyond the half rate](research/post-half-directions.md):
  scoped interface restrictions, [finite rank evidence](evidence/half-piece-ranks.json),
  and an [exact low-mode projection counterexample](evidence/half-projection-counterexample.json).
- [Multiplication-network tradeoff](research/half-top-structure.md): a lower
  bound for the present local-enumeration, separate-wire, and explicit-pad
  construction, with its scope stated explicitly.
- [Restricted-label counterexample](research/half-restricted-labels.md):
  an exact [nonzero minor](evidence/half-restricted-labels.json) disproves the
  literal low-mode product-span bound while preserving the crossing identity.
- [Effective homogeneous duals](research/homogeneous-duals.md): an explicit
  integer inverse polynomial, controlled coefficient sizes, and a complete
  alternative crossing family with the existing `1/sqrt(2)` parameter.
- [Quantitative restricted-label spans](research/half-restricted-span.md):
  an almost full intermediate-dimension bound, with exact constants and
  explicit conditions on any selected binary subset. It does not bound the
  ranks of the final composed pieces.
- [Actual binary-label channel rank](research/binary-prefix-top-channel.md):
  a sharp leading rank with logarithmic error, its lossless-compression
  consequence, and [finite exact checks](evidence/binary-prefix-channel.json).
- [Near-all-ones labels](research/near-ones-labels.md): the fixed-defect
  embedding, its factor-rank restriction, and exact support and small-rank
  [checks](evidence/near-ones-channels.json).
- [Primitive central-piece rank](research/primary-central-gram.md): an
  all-weight Gram identity after a row swap, with
  [exact rational checks](evidence/primary-central-gram.json). A large-image
  asymptotic theorem remains open.
- [Fixed-characteristic factoring transfer](research/scale-one-modular-transfer.md)
  and [sparse divisor nonzeroness](research/sparse-divisor-nonzeroness.md):
  a complete conditional route from effective modular crossings to one
  factoring machine, including every ring and divisor-recovery cost.
- [Normalized modular crossings](research/modular-crossing.md): exact
  scale-one pieces over a fixed finite field and all its commutative
  coefficient algebras. [Identity checks](evidence/modular-crossing.json)
  and [sparse-recovery checks](evidence/sparse-divisor-nonzeroness.json)
  cover their finite scope separately from the asymptotic proofs.
- [Baseline audit](research/baseline-audit.md): exact depth recurrence,
  contraction, all-level cost bound, and machine-cost bridge.
- [Verified k6 theorem](research/verified-k6.md): complete machine proof,
  integer and real rate statements, source changes, and reproduction steps.
- [k6 source audit](research/k6-source-audit.md): hashes of the protected
  reference definitions, actual runtime schedule, and compiler preservation.
- [Grouping theorem](research/grouping-improvement.md): complete concentration
  proof, integer tests, deterministic construction, and integration obligations.
- [Crossing compression](research/crossing-compression.md): explicit factors,
  coefficient identities, changed multiplier, and a limitation of this coefficient
  family.
- [Exact crossing checker](tools/check_crossing.py): finite identities, channel
  encodings, and rank/width comparisons.
- [Combined integration](research/integration.md) and [recurrence proof](research/recurrence.md):
  explicit integer parameters, common widths, all-level costs, and remaining
  formalization obligations.
- [Finite grid step](research/grid-level-validation.md): actual row and column
  contractions, signed tables, composite moduli, and the required `8^t` multiplier.
- [Coefficient-basis obstruction](research/coefficient-basis.md): an exact minimax
  result ruling out further gains from linear recombinations in the same coefficient
  span, with an exact small-instance search.
- [GroupingThreshold.lean](formal/GroupingThreshold.lean) and
  [CrossingComplement.lean](formal/CrossingComplement.lean): checked finite lemmas.
- [ExplicitChain.lean](formal/ExplicitChain.lean): a checked bound on every
  chain exponent and a complete time bound for the unchanged reference machine.
  This formal envelope retains the source's `23/24` contraction; the sharper
  `11/12` analysis remains a separate handwritten derivation. The completed
  k6 machine proof uses the separate modules linked above.
  [Verification evidence](evidence/formal-verification.json) records source hashes
  and the allowed axiom list.

## Reproduce the checks

The Python checks require only Python 3 and the standard library:

```sh
git clone https://github.com/1337isnot1337/integer-factoring.git
cd integer-factoring
make verify
```

Without Make, run the commands in the `verify` target of the [Makefile](Makefile).

The larger, optional rank suite is `make verify-post-half`. It checks 13
finite cases over two primes, including the source's first three binary
label domains. NumPy accelerates this suite when installed; the standard
library suffices. These rank checks test attempted extensions and do not
prove the asymptotic factoring result.

With Lean 4.35.0 and Lake on the executable path, run:

```sh
make verify-k6
```

This target verifies the complete k6 result from a fresh source tree.
It requires Python 3.11 or newer, Git, Bash, `flock`, and network access for
the pinned reference and Mathlib cache. It records logs and hashes under
`build/k6/`. It rejects final axiom reports outside `propext`,
`Classical.choice`, and `Quot.sound`.

For the original local lemmas and unchanged reference envelope, run:

```sh
make verify-formal
```

The GitHub Actions workflow runs both formal targets and the Python checks.
It pins the Lean archive hash and checkout action. It runs for pull requests,
pushes to `main`, and manual dispatch. Branch pushes without a pull request
do not start CI. No leaderboard submission accompanies this result.

The [workflow badge](https://github.com/1337isnot1337/integer-factoring/actions/workflows/verify.yml)
shows current GitHub Actions results. For the recorded local reference build and
axiom scope, see the [build record](evidence/upstream-build.json) and
[formal verification evidence](evidence/formal-verification.json). Historical
CI observations are kept in [evidence/ci-run.json](evidence/ci-run.json).
The complete k6 result has a separate [verification record](evidence/k6-verification.json)
and [final axiom output](evidence/k6-certified.log).

## Attribution

The factoring reduction and counting engine are due to the pinned upstream
development. Its README credits Isaac Rudich and Louis-Martin Rousseau and
describes its proof/code generation with Claude. Adapted source is identified
in file headers. This repository uses the upstream-compatible Apache-2.0 license.
The `1/2` proof additionally uses Claude (Anthropic)'s Appendix B pieces,
grouping, recursive evaluator, and finite-symbol implementation in
`anthropics/formal-math`; it credits those existing results explicitly.
The k6 patch, proofs, and research notes used substantial OpenAI assistant
support. Lean verification establishes the formal claims. It does not replace
independent review of their mathematical significance.
