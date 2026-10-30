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

It combines Claude (Anthropic)'s existing Appendix B clique evaluator with
Isaac Rudich and Louis-Martin Rousseau's multiplication-network reduction.
The new connection includes a deterministic written-table constructor, an
exact value-preserving codec, a computable runtime pad, and all input,
arithmetic, recursion, and prefix-search costs. It gives score `alpha=1/2`
under the published rule and the rate of the reference's square-root-logarithm
conjecture.

**Evidence level: complete ordinary mathematical argument, with source audit
and finite checks.** The complete composition has not been formalized in Lean or
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

The [work beyond `1/2`](research/post-half-directions.md) records the
directions tested next. It proves restrictions for the existing grid and
operator interfaces, gives an exact counterexample to projecting away high
modes, and checks finite matrix ranks. A [network tradeoff](research/half-top-structure.md)
rules out improving this explicit constructor by regrouping and separate
wire compression alone. A [restricted-label counterexample](research/half-restricted-labels.md)
shows why another proposed dimension reduction fails. These results establish
no exponent above `1/2` and no general factoring lower bound.

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
| Appendix B factoring composition | `alpha = 1/2` | [Complete ordinary proof](research/half-rate-factoring.md), source and interface audit, finite adapter and arithmetic checks. Lean checks the [value-level composition through `Nat.minFac`](formal-appendix/PureFactor.lean), with explicit grouping, tables and pads; the costed machine composition remains open. |
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

To rebuild the pinned Appendix B evaluator and factoring network dependencies
with Lean 4.33.1 and check the explicit-pad and grouped top-network theorems,
run:

```sh
make verify-half-appendix
```

This compiles the unchanged factoring network and grouping source against a
different Mathlib revision from the reference's Lean 4.35.0 build. The checked
table codec preserves the value of any reference clique, and the reference's
search supplies a valid top grouping. For `K^40 ≤ n`, Lean also proves that
the selected width satisfies `K^2 * width ≤ 6000n`, and bounds the selected
group's state count by an exponent whose size is at most `8000n/K`. Lean proves
the grouped network's weights and factors are zero or one, obtains the entry
bound, and supplies explicit table and evaluator pads. For a level satisfying
the selection inequalities, Lean bounds the evaluator-pad exponent by a
constant times `n / sqrt(log n)`, including for the explicitly selected level.
The level selection, selected `K`, network-size parameter, table width,
entry-bit bound, and pad exponents are polynomial-time in unary when the input
width is supplied in unary. One combined pad covers the grouping scan, table
writer, and evaluator and has an exponent on the same scale.
Lean also bounds the evaluator fragment's steps on the already-written input
on that scale. The composed mathematical function returns `Nat.minFac` on all
natural inputs. The reference challenge file contains an unrelated placeholder
theorem; the new theorem's axiom report does not depend on it. The costed
constructor and complete 1/2 factoring machine remain open. A separate
machine bridge constructs the reference finite-alphabet factoring machine
from a correct polynomial-time padded algorithm, but leaves the factoring
algorithm and its polynomial-time certificates as hypotheses. The prefix-search
polynomial-time theorem likewise assumes a polynomial-time count query.
The exhaustive grouping scan is polynomial-time given a polynomial-time
predicate and sufficient unit padding; its padding exponent is bounded by
`3n/r` for the selected schedule. Its group-size check is polynomial-time and
agrees with the reference; the crossing-width check remains to be certified.

The GitHub Actions workflow runs the reference Lean checks, the k6 rate checks,
and the Python checks. The Appendix bridge target above is run separately.
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
