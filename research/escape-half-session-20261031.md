# Ordinary-proof frontier session, 31 October–1 November 2026

## Status and reading path

The established deterministic factoring score remains **alpha = 1/2**:

```math
T(n)=\exp\!\left(O\!\left(\frac{n}{\sqrt{\log n}}\right)\right).
```

This session improves the complete arithmetic constructor and supplies
stronger explicit crossing inputs. It does not prove an exponent above
one half. All new mathematics is written as ordinary proofs and exact
computational certificates; no Lean source was changed or built.

The most useful reading order is:

1. [Width-sensitive prime schedule](width-sensitive-prime-constructor.md)
   and [rough-input observability](rough-observable-factor-constructor.md)
   for the strongest complete constructor result.
2. [Stable physical codes](stable-hot-reflection-codes.md),
   [complete occupation profiles](hot-code-charge-profile-constructor.md),
   and the [96-letter exact interval certificate](exact-natural-pressure-intervals.md)
   for the strongest certified crossing bound in this checkpoint.
3. [Finite-block amplification](finite-block-crossing-amplification.md)
   and the [asymmetric machine transfer](asymmetric-modular-transfer.md)
   for the exact remaining route to a higher score.
4. [Complete paired affine space](complete-paired-trace-affine-space.md)
   and [entrywise interpolation](entrywise-reflection-jet-interpolation.md)
   for a different source of explicit crossing tensors.

The work is on branch `research/escape-half-422a2e73`, based on
`a3c5bf15fcc80ea8ab9af7b4853a1a70ae6a88b0`, with
[ready-for-review PR 46](https://github.com/1337isnot1337/integer-factoring/pull/46).
The session began at `2026-10-31 21:13:54 UTC`. Three parallel research
workers cover crossing algorithms, arithmetic bounds, and proof audit;
the main worker integrates, independently checks, and publishes results.

## 1. A smaller complete factoring constructor

For every fixed rational epsilon > 0 and sufficiently large input bit
length n, the new construction has

```math
P=\Theta_\epsilon\!\left(\sqrt{\frac n{\log n}}\right),
\qquad \lambda=\Theta_\epsilon(\log n),
\qquad P^2\lambda<(12+\epsilon)n,
\qquad B_{\rm tab}=1.
```

Here P is the number of groups in the initial clique instance, lambda
is its label width, and B_tab bounds the bits needed for each table
entry. Every entry is Boolean. The total deterministic cost of the
small-factor pretest, prime search, coefficient-field construction,
local entry computation, and writing all tables is

```math
\exp\!\bigl(O_\epsilon(\sqrt{n\log n})\bigr).
```

The improvements proceed through explicit constructions, rather than
an uncharged change of parameters:

| Construction | Written top volume | Local entry bound or method |
|---|---|---|
| Adaptive transform | `P^2 lambda <= 99n` | `B_tab <= 117n/P` |
| Operand-length classes | `P^2 lambda <= 49n` | `B_tab <= 57n/P` |
| Observable boundary grouping | `P^2 lambda < 17n`, leading constant 12 | Boolean entries |
| Width-sensitive prime schedule | `P^2 lambda < (12+epsilon)n` | Polynomial-time interpolation after the pretest |

The [adaptive transform](adaptive-transform-constructor.md) removes
square-capacity padding, sizes the convolution modulus to the actual
coefficient range, and gives carries their own smaller domain.
[Operand-length conditioning](operand-length-transform-constructor.md)
then uses the true operand supports. The length classes are disjoint,
so summing their exact prefix counts preserves the original count.

The [observable grouping](observable-boundary-constructor.md) leaves
enough boundary Fourier observations to determine the internal digit
polynomial. The final rough-input construction first tests small
factors through an explicitly charged threshold. If the input survives,
both operands in every retained proper-factor class have sufficiently
small unknown local height. Their boundary observations determine them
by a Vandermonde system. All range, padding, multiplication, and carry
constraints are then checked directly.

This yields the exact Boolean local table, or the corresponding product
of local weights over an allowed coefficient algebra. The divisor
polynomial, cyclic nonzeroness test, and prefix descent recover the
least prime factor, with a separate prime-input branch. The pretest also
reduces the divisor-support bound to `exp(O(P+log n))` for the cyclic
recovery overhead. The pretest principle was independently recorded in
[PR 31](https://github.com/1337isnot1337/integer-factoring/pull/31);
the new use here forces the local observability inequalities.

The recursive evaluator's polynomial prefactor is controlled by label
width rather than an avoidable extra power of clique size. At the new
schedule its additional overhead is
`O(n log log n / log n)`, which is lower order than
`n/(log n)^alpha` for every fixed alpha < 1. The effective prime-in-
progression argument and explicit root-field construction include their
search costs and finite exceptions. See the primary references and
exact parameter statements in the linked proofs.

This makes the reduction cheaper. The complete evaluator still needs
a stronger crossing bound to improve the score.

## 2. Stronger finite crossing inputs

### Stable small physical codes

Retain a bounded-energy physical code inside a sufficiently large
reflection interval. Once its maximum energy satisfies the stated
stabilization condition, its primitive polynomial tensor no longer
depends on the larger ambient interval. This gives explicit small
crossings that differ from simply taking the full alphabet of a small
interval.

[Natural-cut duality](reflection-natural-cut-duality.md) identifies
the two complete weighted profiles up to opposite scalar shifts. For
coordinate downsets, triangular factors reduce the profile computation
to selected compound minors, separated by occupation number. The
[charge-block engine](hot-code-charge-profile-constructor.md) retains
every coefficient and every channel; insufficient precision is an
error, not a partial rank bound.

The following numerical values are complete-profile diagnostics at
the recorded parameters. They are not outward-rounded minima.

| Family | Alphabet | Recorded natural pressure |
|---|---:|---:|
| Original reflection interval | 16 | 0.716653659728 |
| Stable coordinate code | 16 | 0.714736631324 |
| Stable coordinate code | 24 | 0.712964178785 |
| Stable coordinate code | 32 | 0.7121199592 |
| Stable coordinate code | 48 | 0.7111215197 |
| Stable coordinate code | 64 | 0.710590358354 |
| Stable coordinate code | 96 | 0.7100058748811902 |

The exact certificates are separate:

- The [24-letter certificate](stable-hot-code-rational-pressure.md)
  proves `rho < 5/7`.
- The [64-letter certificate](hot64-rational-pressure.md)
  proves `rho < 32/45`.
- The [96-letter interval certificate](exact-natural-pressure-intervals.md)
  proves the exact rational bound `rho < 0.710005874881191`.
- A [coherent eight-letter boundary code](coherent-hot-boundary-crossings.md)
  proves `rho < 18/25`.

For the 64-letter code, the exact natural profile contains 4,096
channels. With rational base `459/460` and supporting parameter `64/45`,
the proof reduces to a strict integer inequality. Both positive
integers have 14,266,197 bits. Positive rational physical types and the
strict gap give an effective finite block, and finite amplification
supplies the uniform family at every external label width.

The 96-letter calculation completed all 9,216 channels, with largest
block 2,738, in about 25 minutes. Separate large-integer controls prove
`rho < 49/69` and `rho < 27/38`. The finer certificate uses exact outward
power intervals and rational logarithm tails, with 192-bit endpoints.
Its arithmetic takes less than one second on the saved complete profile;
that timing excludes the earlier profile construction. Every bound above
remains larger than `1/sqrt(2)`. None establishes a post-half factoring
exponent.

### Coherent boundary maps and full-letter gauges

The coherent construction permits arbitrary paired maps K,J on one
equal-energy boundary with `KJ=I`. All lower-energy states are retained.
This preserves the crossing endpoint and the regular triangular
factorization. The two natural profiles may differ, and both are
computed. Mixing occupation numbers requires a complete matrix profile
instead of the old occupation split.

The explicit eight-letter example has degree 56 and its exact
`18/25` comparison uses 720,752-bit integers. A declared finite grid
of 396 rank-one boundary idempotents is fully accounted for. The
archive explicitly joins an incomplete prefix with a separately
computed missing case, rather than presenting the prefix as complete.

The [weighted compound cone](weighted-hot-compound-profiles.md) allows
nonadditive full-letter physical heights under explicit integrality
inequalities. The [weighted search oracle](weighted-hot-compound-search.md)
checks the mirrored common-type pressure and its derivatives. These
are new search coordinates; they do not imply that every useful gauge
lies in the cone.

## 3. A complete affine space and direct jet construction

The paired-trace invariant proves scalar polynomial content one and
the degree floor M^2 for its declared reflection-frame family. An
[explicit uniform truncation](low-jet-reflection-degree-truncations.md)
attains degree M^2 from the original degree `3M^2/2`, retaining all
lower coefficients and both paired traces. The degree-16 four-letter
witness has complete polynomial provenance.

For the four-letter independent-paired closure, the exact affine
direction space is

```math
(1+w)^2(\mathfrak{sl}_4\otimes\mathfrak{sl}_4)[w]_{\le22},
\qquad \dim=15^2\cdot23=5175.
```

An invariant upper bound is attained by exact finite-field saturation.
Its canonical sparse basis is explicit. The degree-at-most-16 slice
has dimension `225*15=3375`, and every member has actual degree 16.
The subsidiary [polynomial-shift closure](polynomial-shift-affine-closure.md)
has dimension 3,674 and is contained in the larger completed space.

The direct entrywise solver handles each polynomial's Hermite
constraints first, then imposes the paired traces. The coupled stage
has rank at most `31(d-1)` at degree cap d. It supplies exact opposite
Bell-jet tensors without a dense elimination on all 5,175 directions.

Two ordinary theorems sharpen the search scope. The
[Bell transversality theorem](bell-jet-transversality.md) gives a unique
simultaneous eighth-order carrier for every pair of invertible Bell
forms. The [ramified-minor theorem](ramified-carrier-minor-domination.md)
shows that updates below the stated rank cannot improve any weighted
determinantal prefix at the existing ramified places. A separate
argument treats all carriers for the unweighted alternating form.
Neither theorem rules out arbitrary full-rank weighted changes or
additional singular places.

The [finite Smith-frame fibers](finite-jet-smith-fibers.md),
[general F4 iteration](f4-iterated-smith-fibers.md),
[paired restrictions](finite-fiber-paired-restrictions.md), and
[closed-place oracles](closed-place-weighted-oracles.md) provide
additional complete polynomial inputs. Their frame inequalities,
actual degree after scalar division, unequal paired alphabets, and
all closed-place channels remain explicit.

## 4. Arithmetic contractions and their remaining interface

[Original-field monoid convolution](unsplit-monoid-contraction.md)
decomposes prime-power multiplication without extending the coefficient
field. It gives explicit circuits and a near-linear-rank whole-middle
bundle in the declared smoothness regime. The complete linear maps,
not only the bilinear product count, are charged.

[Frobenius rational stars](frobenius-rational-trace-contraction.md)
contract products of scalar Frobenius-affine factors through exact
pole traces and a binary-carry coefficient recurrence. The
[infinity correction](frobenius-trace-infinity-correction.md) removes
the proper-degree restriction. [Fermat-torus boxes](fermat-torus-binomial-contraction.md)
have division-free fixed circuits over the stated coefficient algebras.
[Relative-linearized factors](relative-frobenius-trace-contraction.md)
retain an explicit dimension-dependent cost, and
[semisimple cyclic recovery](semisimple-cyclic-box-contraction.md)
handles componentwise inverse-or-zero values, including zero divisors.

These results do not yet close the full multiplication network. Its
jointly varying bilinear multiplication phase and the representation
after successive eliminations remain unbounded by the desired compact
class. A local fast sum or a small middle-layer rank cannot be applied
without charging its neighboring interfaces.

## 5. The exact gate for a higher score

Suppose a finite crossing on paired alphabets of sizes d_X,d_Y has a
complete decomposition into pieces with two ranks r_i,c_i. Put

```math
S=\sum_i\sqrt{r_i c_i}.
```

The [finite amplification theorem](finite-block-crossing-amplification.md)
constructs an effective all-width family with

```math
\rho\le\frac{\log S}{\log(d_Xd_Y)}.
```

It includes unequal alphabets, paired swapping, polynomially many
composition classes, factor queries, and coefficient reconstruction.
An exact strict bound below `1/sqrt(2)` would supply the required
balanced improvement. The complete
[asymmetric transfer](asymmetric-modular-transfer.md) also permits
a weighted supporting-line improvement with unequal child sizes.

The [finite typed adaptive checker](finite-typed-adaptive-certificates.md)
and [rational majorants](rational-finite-adaptive-bounds.md) retain all
endpoint pieces and the first positive jet after a free prefix. They
test integer inequalities without treating floating logarithms as
certificates. No tested finite block passes the post-half gate.

## 6. Completed searches and reproducibility

The published archives distinguish a finished time budget from an
exhausted finite search space. They preserve unresolved calls and
timeouts; neither is counted as a valid profile or a negative theorem.

| Search | Recorded scope | Main outcome |
|---|---|---|
| [Fused physical bases](fused-physical-search-checkpoint.md) | Two-hour saved budget, 1,049 attempted and 1,048 completed trials | No threshold candidate; both family winners freshly replayed |
| [Four-letter Laurent series](laurent-series-search-checkpoint.md) | All 65,536 parameter words | No gain; retained full tensors replayed |
| [Sixteen-letter one-direction affine search](filtered-reflection16-affine-search.md) | 1,386 choices, 211 distinct complete tensors | Baseline retained; both finite nontrivial F4 places screened |
| [Cross-direction affine combinations](filtered-reflection16-affine-combinations.md) | 849 distinct complete tensors | No gain; all polynomial entries reconstructed |
| [Weighted closed-place refinements](restricted-closed-place-search.md) | 1,978 partitions, 402 full refinements | Some weak seeds improve, but no competitive bound |
| [Finite-fiber basis and iteration](finite-fiber-search-checkpoint.md) | Two thirty-minute budgets | Full winners and forty retained fiber witnesses replayed |
| [Resumed physical families and targets](root-worker-search-checkpoint.md) | 6,722 saved candidate outcomes plus seven target/search records | No gain; explicit interruption and continuation lineage |

The interrupted physical workers retain immutable prefixes and new
declared continuation streams. Lost in-memory random states are not
claimed to have been restored. Unsaved interrupted work is not counted.
The affine refinement worker exhausted its selected queue before its
time ceiling; the Hasse target worker completed its time budget with
an incomplete target queue. Their notes state both distinctions.

The archives include exact source and input hashes. Some preserve the
complete raw bytes; the root physical archive deliberately keeps a
compact ledger of every candidate score and digest plus complete
winners. It does not claim to reproduce every discarded optimizer
trajectory. The original local reports remain intact.

## 7. Validation and next questions

The final default `make verify` replay passed all 50 main unit tests,
all five crossing unit tests, and its exact baseline certificates.
The documentation gate and direct Python module-dependency checks pass
at the integrated checkpoints. No Lean gate was run for this session.

Fresh main-worker replays also passed the full 64-letter profile and
integer certificate, 94 independent complete determinant controls,
the weighted oracle's 102 gradient coordinates, the 5,175-dimensional
affine saturation, the ramified-minor certificate, and the coherent
eight-letter certificate. The new 96-letter pressure arithmetic passed
911 exact power controls, 125 independent alternating-log controls,
24 integer-power sign comparisons, and separate executable audit.
Its full polynomial profile was constructed once and frozen; the fast
interval check does not claim a second large-matrix reconstruction.
Separate research-worker replays cover the
constructor, local contraction, finite-block, and archived-witness
packages. Each note states what is exact and what remains numerical.

The most concrete open questions are now:

- Find a stable physical code, coherent boundary, or complete finite
  adaptive block that passes an exact post-half gate. Larger alphabets
  alone have not yet done so.
- Use the full affine and finite-frame coordinates to change the
  complete two-cut pressure, rather than improving one selected rank
  or jet while losing elsewhere.
- Prove a compact class that survives successive arithmetic
  eliminations with the jointly varying multiplication phase, charging
  every external interface.

These questions identify the missing mathematical improvement. The
constructor and local identities are available for use once that
improvement is found.
