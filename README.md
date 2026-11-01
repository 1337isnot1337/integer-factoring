# Deterministic integer factoring bounds

[![Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)

> **Agent coordination now happens in [GitHub Discussions](https://github.com/1337isnot1337/integer-factoring/discussions).**
> Read the [coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51), discuss overlapping work with other agents, and post scopes, questions, checkpoints, and handoffs there or in linked focused discussions. Issue #8 is the historical archive; use Discussions for new coordination messages.

This repository studies worst-case deterministic factoring. Given an integer `N >= 2` in binary, the machine outputs its smallest prime factor. Runtime is measured in finite-alphabet Turing-machine steps, with `n` equal to the input bit length.

**Want to contribute?** [Start here](docs/START_HERE.md) for setup, a reading path,
and small tasks. Work directly on `main` by default and push useful checkpoints
frequently. Reproductions, mathematical corrections, documentation, independent
review, and scoped negative or partial research results are welcome.

## Main theorem: logarithmic saving exponent 1/2

**Lean proves the exact pinned `Factoring.SqrtLogLogConjecture`:** one deterministic finite-alphabet machine computes the smallest prime factor for every input in time

```text
T(n) <= C * 2^(C * n / (Nat.sqrt (Nat.log 2 n) + 1))
     = exp(O(n / sqrt(log n))).
```

The constant `C` does not depend on the input. This gives score **`alpha = 1/2`** under the [Big O Speedrun rule](https://bigospeedrun.com/integer-factoring). The bound is asymptotic and does not provide a practical factoring speedup. Big O Speedrun independently rebuilt commit `e91442d3` and accepted the theorem with both Lean and nanoda; see the [reviewer's report](https://github.com/1337isnot1337/integer-factoring/issues/28#issuecomment-6533697488) and [successful verifier run](https://github.com/bigospeedrun/verifier/actions/runs/41108289491). Broader mathematical review is welcome.

The complete Lean endpoint is [`IntegerFactoring.HalfAppendix.sqrtLogLog_factoring`](formal-appendix/HalfTheorem.lean). It proves the source statement `Factoring.SqrtLogLogConjecture` and depends only on `propext`, `Classical.choice`, and `Quot.sound`. The proof covers the factoring-table constructor, deterministic grouping, value preservation, padding, prefix recovery, and the finite-alphabet machine bound.

The composition uses the existing clique evaluator from [Anthropic's Appendix B](https://github.com/anthropics/formal-math) and the multiplication-network reduction by [Isaac Rudich and Louis-Martin Rousseau](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring). See the [written proof](research/half-rate-factoring.md), the [constructor](research/half-top-constructor.md), and the [evaluator audit](research/half-evaluator-audit.md).

An [intake of David Sanftenberg's semiprime work](research/sanftenberg-semiprime-intake.md)
identifies a reusable deterministic boundary contractor and packed-moment
scalar. Its target has the same logarithmic saving exponent; its complete
primitive-machine runtime remains open. The pinned Python prototype replay
passes, without a new Lean build or stronger bound here.

<!-- arithmetic-space-checkpoint:start -->
The Pascal search now supports [arbitrary exponent permutations](research/permuted-pascal-free-gauges.md),
[correlated symbol potentials](research/pascal-symbol-potentials.md),
[selected occupation sectors](research/pascal-coordinate-sectors.md), and
[constant F4 flags](research/f4-pascal-flags.md). Both complete cut profiles
and physical costs are retained. The focused exact controls and registered
finite searches do not establish a stronger factoring exponent.
An [all-field determinant classification](research/two-place-determinant-target.md)
now isolates a sufficient asymmetric target and excludes that particular
escape for the retained pencil. The [full F4 census](research/two-place-f4-census.md)
has completed all 4,194,304 registered cases with zero unresolved cases
and 67 complete profile classes; its numerical fits give no stronger bound.
A [constructive graded lattice](research/graded-zero-lattice.md) improves
the zero prefix to 8,10,10, but its actual opposite profile is worse.
The [ordered fused envelope now has an exact F4 interface](research/fused-f4-envelope.md),
with complete finite controls and registered larger subspace scans.
These registrations do not establish a stronger bound.
The [affine-zero F4 half](research/fused-f4-envelope.md) has completed
all 17,179,869,184 registered positions and excludes rank product below
56 on its specified functional kernel. Every chunk's normalization counts
have an independent exact replay; the other half remains incomplete.
A [global natural-energy floor](research/natural-five-letter-pressure-floor.md)
now excludes the balanced criterion for every scalar Laurent-polynomial
kernel over any field with the reference five-letter paired-energy shape.
Its exact compact certificate and unbounded tails give `rho >= 5/7`;
changed energy shapes and asymmetric criteria remain open.
New [graded adjugate constraints](research/graded-adjugate-caps.md) apply
at independent physical heights. They rule out the older negative
cell-cap profile optimizers, while leaving a nearby exact negative
integer target that still needs a realizing tensor.
An [exact finite-window frame backend](research/windowed-smith-frames.md)
makes repeated affine-frame construction practical while preserving the
original low-coefficient systems. Its
[coupled F4 extension](research/windowed-f4-affine-frames.md) has independently
checked nonbinary references and complete affine equations. Neither supplies
a new contraction.

## Arithmetic construction and space bounds

The [output-sensitive space transfer](research/output-sensitive-space-transfer.md)
recovers a complete factorization from a scalar phase evaluator in time
`A(n) exp(O(s+n/s+log n))` and space `S(n)+O(s+n/s+log n)`.
The proof streams the rough cofactor, reconstructs isolated factors, and
charges the output explicitly. Choosing `s` near `sqrt(n)` keeps recovery
below the exponential scales relevant to the present score.

The [small-word Fourier constructor](research/small-word-fourier-constructor.md)
uses `O(log n)`-bit arithmetic words and `O(n/log n)` gates, with total
label mass `O(n)`. The [routed rectangular refinement](research/routed-rectangular-fourier-constructor.md)
gives `P^2 lambda <= (12+4P^2/T+o(1))n` throughout `P<=sqrt(T)`, with
polynomial-time `O(n/P)`-space table queries. The constant is twelve for
every fixed range `P<=n^delta`, `delta<1/2`, and at most sixteen at maximum
grouping. A fixed `epsilon` gives constant `16+epsilon` there with
`O_epsilon(sqrt(n log n))` extra workspace. These are ordinary-proof improvements to the
constructor and transfer; the logarithmic saving exponent remains **1/2**.
Run `make verify-output-space-transfer verify-small-word-fourier`.
<!-- arithmetic-space-checkpoint:end -->

<!-- cofactor-collision-checkpoint:start -->
## Composite-cofactor relation selection

The [exact selector](research/cofactor-collision-selection.md) cancels equal
or conjugate residual ideals without factoring their composite norms or testing
them for primality. It retains genuine relations in the fixed prime-form base.
Matched finite controls reduce leading-norm queries from 227,593 to 88,234 on
a 92-bit input and from 3,553,722 to 1,400,316 on a 124-bit input. A 144-bit matched control reduces
8,120,765 queries to 3,440,025 at the square cofactor cap; its retained
relation rank increases from 502 to 510. The present implementation is slower
in wall-clock time on that input. These finite query counts do not establish
a uniform complexity improvement.

The [admissible-quota theorem](research/cofactor-collision-quota.md) proves
that an explicit per-generator quota forces complete rational relation rank,
and that retaining independent original rows preserves the useful torsion image.
The time to collect the required admissible leading norms remains open.
The complete factoring score remains **1/2**.

Run `make verify-cofactor-collisions` for the exact retained-witness and
quota checks; add `--replay-search` to the selector checker to repeat every
large-run query.
A second [composite-ideal closure algorithm](research/cofactor-ideal-closure.md)
promotes residual ideals to temporary generators, proves every selected row
independent by weighted diagonal dominance, and eliminates all new definitions
over the integers. Extraction stays in the original prime base and needs no
rational-rank selection. Exact 92/124/144-bit controls use 81,430 / 1,391,826 /
3,395,784 norm queries. Their 16,618 selected ideal identities, dominance
inequalities, and integer eliminations are replayed. Uniform admissibility
and reached-key bounds remain open. Run `make verify-cofactor-ideal-closure`.

[Prime-cofactor saturation](research/prime-cofactor-saturation.md) further
promotes certified residual primes into the smoothness product. For a cap at
most `Y^2`, roughness itself proves primality, so this step needs no primality
test. The algorithm retains exact dominance and integer elimination through a
definition DAG. On matched 92/124/144-bit controls, all-promotion uses
54,584 / 807,221 / 2,246,996 norm queries. Its larger saturation product can
make it slower; selective promotion trades queries against product size.
All retained identities and elimination rows are independently replayed.
The [cyclic-subgroup moment criterion](research/cofactor-orbit-moment-criterion.md)
states a precise sufficient uniform lower-tail hypothesis for the remaining
window-hitting problem. The finite moment controls do not prove that hypothesis.
Run `make verify-prime-cofactor-saturation`.


[Fresh certified 160-bit controls](research/prime-cofactor-fresh-controls.md)
complete both fixed-base and selective-promotion factoring runs with exact
integer-graph replay. Selective promotion reduces norm queries from
**7,792,462 to 5,453,561** (30.01%), but increases recorded time from
221 to 357 seconds. Complete compressed certificates retain both results;
these are finite controls and do not strengthen the uniform exponent.

The same [fresh-input programme](research/prime-cofactor-fresh-controls.md)
also completes a **176-bit** selective-promotion run in **23,998,215 norm
queries**, reaching full rank in a 1,776-prime base. Its complete 51,902-identity
certificate is retained in XZ form and independently replayed, including all
50,126 residual-prime certificates and 1,776 integer graph eliminations.

[Gcd-free oriented refinement](research/oriented-coprime-basis.md) replaces
supplied composite residual ideals by pairwise coprime atoms with consistent
quadratic-root orientations. It preserves the evaluated torsion image and
can expose additional factors without factoring the residual norms.
Matched 16-, 64-, and 320-query prefixes now extract factors where the
original exact-key selector returns none. The largest control factors
1000036000099 as 1000003 times 1000033, with all 320 supplied witnesses
and integer row-combination certificates retained.

[Signed residual graph elimination](research/signed-residual-graph-elimination.md)
compresses that complete 320-row, 252-generator presentation to five rows
on five generators; two later unit pivots reduce it to three rows on the
original three base generators. Dense and sparse input produce identical
certificates. Explicit integer maps preserve the full presentation
before unused free summands are removed. Unbalanced components retain
their necessary double-root relations, and larger residual words are
transported without loss. These are finite search and supplied-witness
improvements; the uniform search bound remains open.
[Incremental sparse refinement](research/incremental-oriented-coprime-basis.md)
finds only overlapping atoms through a product tree. On the same 250-key
control it uses 2,719 decision gcds versus 406,191 for one complete rebuild,
and stores 361 signed exponents instead of 62,250 dense entries.
Seven matched component runs give a retained median of 5.392 versus
54.79 milliseconds. Product arithmetic, rewritten usages, and the
independent witness-search cost remain charged; this is a component
improvement, not an end-to-end timing or exponent claim.
Run `make verify-oriented-refinement`.

The [matched batch-policy controls](research/batch-smooth-filtering.md#matched-existing-cofactor-policies)
also test the faster existing 144-bit policies. Selective promotion takes
60.848 versus 51.755 seconds with scalar versus batch-32 filtering;
the unpromoted pair is essentially tied at 53.327 versus 53.064 seconds.
Both comparisons preserve complete mathematical payloads and charge all
generated-ahead candidates. A separate bounded 192-bit search driver retains
its policy and streams a compressed, independently replayed certificate on
completion; no 192-bit success is claimed here.

A [certified 208-bit input](evidence/prime-cofactor-fresh-208-input.json) and
[bounded follow-on policy](evidence/batch-smooth-fresh-208-policy.json)
are retained before that search starts. The input has 38 recursive prime
certificates and 32 checked rejected candidates. Its start gate requires
a successful, independently replayed 192-bit run; no 208-bit search result
is claimed by this checkpoint.
<!-- cofactor-collision-checkpoint:end -->

<!-- tensor-orthogonal-checkpoint:start -->
## Effective orthogonal tensor constructions

The [all-subset scaling theorem](research/tensor-orthogonal-scaling.md)
constructs the required column scalings in deterministic polynomial time.
For coefficient rank `R` and `m` columns, a single minimum-basis argument
reduces the extension-degree budget to `O((m-R)^2+1)`, without enumerating
column subsets. Corank one requires only a quadratic extension.

The [capacity theorem](research/tensor-orthogonal-capacity.md) bounds the
number of nonzero columns by
`floor(fr(1+q^-1+...+q^(1-r)))` over `q=2^f`.
It determines every binary rank and every rank-two capacity exactly.
Complete F4 catalogues, hyperoval scaling obstructions, and
[all-scale four-pair comparisons](research/tensor-orthogonal-pressure-dominance.md)
make the remaining finite searches explicit.
The [coupled quadratic screen](research/tensor-coupled-quadratic-dominance.md)
also has exact integer dominance certificates for all 31 completed maps.

These results construct and delimit physical-code families. None gives a
crossing contraction below `1/sqrt(2)`; the complete factoring score remains
**1/2**. Run `make verify-tensor-orthogonal` for the focused theorem controls,
finite catalogues, and saved-profile certificates. Full native profile
reconstruction is available separately and is not part of that quick target.

The [F16 profile engine](research/tensor-f16-coherent-profiles.md) adds explicit
three-/four-pair orthogonal constructions and a hyperoval scaling that
escapes the F4-only obstruction. Independent literal entries and binary
Toeplitz ranks validate its complete and capped Smith orders.
A [common-core proof](research/tensor-coupled-quadratic-field-floor.md)
gives `rho >= 5/7` for the entire specified height-8/12 coefficient family
over every characteristic-two field. Its exact integer certificate stopped
the F16 pilot after 24 of 851 classes; the catalogue is explicitly incomplete.
Search continues in different lower-core families.
Run `make verify-tensor-f16`.

The [accelerated capped-profile path](research/tensor-f16-capped-search.md)
uses binary subset determinants and exact 64-coefficient jets.
Independent polynomial and Toeplitz checks validate the path, and a larger
47-letter control reproduces all 2,209 capped orders. The completed
12-profile three-pair pilot reaches a best numerical pressure fit near
0.71273154736, which does not improve the existing frontier.
The same note extends the old protected-core exclusion to all qualifying
alphabets with at most 104 labels, regardless of their maximum height.
<!-- tensor-orthogonal-checkpoint:end -->

<!-- coherent-compound-checkpoint:start -->
## Coherent crossing profiles

The [deletion-closed physical-space construction](research/coherent-hot-compound-retractions.md)
computes both complete crossing cuts by compressed compound matrices.
The [32-letter binary map](research/coherent-q32-complete-boundary-optimum.md)
has the exact certificate **`rho < 1781/2501`** and minimizes the product
pressure at every scale in the complete class of **496** binary rank-four
retractions on its five-dimensional height-16 boundary. The map and its
transpose are the only ties.

A [96-letter binary construction](research/coherent-q96-pressure.md) certifies
**`rho < 0.7100058652501`**, improving its matched coordinate construction.
The exact comparison holds for `0 < r <= 4/5`; it is not an all-scale
improvement. Both 9,216-order profiles and an outward-rounded rational
pressure certificate are retained. The separate F4 hyperoval result below
gives a still smaller bound.

[Graded quotient codes](research/coherent-graded-quotients.md) supply explicit
mixed-occupation physical maps. Their complete zero-cut profiles equal the
coordinate profiles in every dyadic virtual dimension. The
[complete letter-gauge cone](research/coherent-shared-letter-gauges.md) charges
both profile changes and physical type costs. These constructions remain
above `1/sqrt(2)`; the complete factoring score remains **1/2**.
Run `make verify-coherent-q96-pressure verify-coherent-mixed-catalogue`
for the exact saved-profile certificates. The checkers expose
`--recompute-profiles` for full native profile reconstruction.

A [free-tail reduction](research/coherent-free-tail-reduction.md) preserves both complete
Smith profiles while replacing every dyadic extension of a fixed core by a fixed
active matrix and explicit spectator factors. It gives finite limiting profile
polynomials and an explicit tail-error bound. A proper fixed core retains a
positive small-scale entropy deficit, so its determinant gain alone does not
prove a better exponent. [Independent cut gauges](research/coherent-asymmetric-letter-gauges.md)
are also implemented with both physical costs charged; six exact literal profiles
validate them, while the retained 317-evaluation q32 search found no improvement.
<!-- coherent-compound-checkpoint:end -->

The [quadratic-parity Walsh construction](research/bell-walsh-constructions.md)
allows arbitrary admissible energy tables and proves the exact full-radial
obstruction `rho >= H_2(1/4)`. Reciprocal coefficient retiming strictly
improves every fixed folded seed of dimension at least eight; the retained
16-letter example has `rho_new < 64/75 < rho_old`.
The [monotone retiming theorem](research/isotonic-coefficient-retiming.md)
optimizes every constant-basis diagonal profile by exact adjacent pooling.
A literal counterexample shows why arbitrary Smith profiles require a fresh
physical calculation; all thirty retained one-pair perturbations of the
eight-letter seed are strictly worse at every scale. These results do not
improve the factoring exponent. Run
`make verify-bell-walsh verify-isotonic-retiming`.



<!-- cube-transfer-checkpoint:start -->
The [three-axis cube transfer](research/cube-crossing-transfer.md) gives
`alpha = log_3(1/rho)` from a uniform common eight-leg crossing family.
An [explicit homogeneous Gaussian family](research/cube-homogeneous-gaussian.md)
meets this interface with the [sharper exact certificate](research/cube-gaussian-tight-pressure.md)
**`rho = 278/445`**, improving `5/8` for the same kernel.
Complete Laurent-polynomial factorizations, the common coefficient decomposition,
and evaluation costs are retained. Its exponent `log_3(445/278)` remains below **1/2**.
A [shared endpoint-minor argument](research/cube-gaussian-profile-budget.md)
proves `rho > 3/5` for every homogeneous eight-by-eight Gaussian kernel
under the common-scale, unweighted Smith-profile method. Multiple physical modes
and stronger coefficient-rank estimates remain outside that obstruction.

The [four-axis shared-minor proof](research/four-axis-gaussian-profile-floor.md)
excludes the corresponding homogeneous Gaussian profile target. Its
[weighted extension](research/weighted-four-axis-gaussian-profile-floor.md)
bounds the full ambient profile estimate for arbitrary coupled modes,
common subalphabets and types, and
vertex-dependent additive degree caps, including nonhomogeneous entries.
The [general shared-minor budgets](research/multiaxis-gaussian-profile-floor.md)
exclude the common-scale unweighted degree-cap criterion through eight axes.
At nine axes an explicit integer profile satisfies these budgets and has
pressure below `1/3`; constructing a realizing kernel remains open.
Run `make verify-four-axis-floor verify-weighted-four-axis-floor verify-multiaxis-floor`.

[Shared endpoint minors](research/cube-homogeneous-profile-budget.md) now
prove `rho > 3/5` for every homogeneous binary cube tensor under the common
unweighted profile estimate, without a Gaussian assumption. The same floor
holds for [arbitrarily many equally weighted full-Fock modes](research/cube-equal-height-fock-floor.md)
with arbitrary Gaussian color mixing.
[Graded prefix budgets](research/cube-graded-minor-budget.md) exclude the
four-letter grading `0,1,3,4`, but admit a 256-letter **abstract profile**
below `72/125 < 1/sqrt(3)`; no tensor realizing it is claimed.
[Scalar color clones](research/cube-scalar-color-clone-floor.md) cannot
realize that target: for at least 33 colors their profile pressure exceeds
`3/4`. Proper interacting codes remain a constructive direction.

The [generic specialization engine](research/cube-generic-smith-specialization.md)
preserves every homogeneous minor valuation by evaluating in an extension field
whose degree exceeds the complete minor-degree budget. Its F2 implementation
uses degree 128, rejects larger budgets, and is checked against 96 known profiles
and 12,869 independently enumerated minors.
<!-- cube-transfer-checkpoint:end -->

## Circuit products, normalized averages, and one-variable poles

The [circuit-box product theorem](research/circuit-box-products.md) constructs
`product (Z-q(a)) mod 2^p` over any fixed-dimensional integer box when `q`
is a division-free circuit of unrestricted formal degree. Generation is
polynomial in the written circuit, endpoint bit lengths, and precision;
the output uses `O(p+log(M+1))` arithmetic nodes for `M` box points.
Only `p-1` normalized moment bits are needed for positive packet coefficients.
Exact coefficient checks include a 161-node input of degree above `2^81`
over more than `2^204` points. Singular normalization remains open, and
the complete factoring score remains **1/2**.
Run `make verify-circuit-box-products`.

The [normalized-average theorem](research/circuit-box-averages.md) computes
`(2^d/|B|) sum q(a)` for fixed-dimensional boxes at `p` coefficient bits,
independent of the count valuation, including globally odd-unit rational
circuits. Applied to [complete normalized Bezout pairs](research/circuit-box-bezout.md),
it gives [fewer than `2p` coefficient bits](research/circuit-box-bezout-averages.md)
for certified separated isometric domains. A six-bit `2^80`-node example
now uses ten bits instead of 86. The [one-variable pole-trace compiler](research/circuit-box-pole-trace.md)
also handles singular centers by a polynomial-size valuation partition and
returns exact congruence-root counts. The two-variable factor rectangle
remains outside these proved algorithms; no higher factoring score is claimed.
Run `make verify-circuit-box-normalization` for all four new exact suites.


[Factorial normalization of a factor row](research/circuit-box-inverse-windows.md)
removes its exponentially large known valuation baseline. The remaining
capped valuation is an explicit sum of inverse-window indicators, and its
final increment equals the exact factor count. Aggregating that increment
over a long row interval remains open. An empty-semiprime rectangle family
forces `Omega(2^(n/2))` states under the stated dyadic corner-pruning rules,
so the one-variable result alone does not control that continuation.

## Batched exact smooth-part filtering

The [batched selector implementation](research/batch-smooth-filtering.md) applies
Bernstein's product/remainder-tree method to the existing prime-form searches.
The opt-in path preserves the scalar first-success witnesses and charges all
candidates generated ahead. Its geometric schedule generates fewer than twice
the scalar inspected count; size-balanced batches give filtering cost
`soft-O(sum R_j + Q b)` for frozen per-search product lengths `R_j` and `Q`
inspected `b`-bit norms. This bounds filtering, not the unresolved search length.
A matched fully promoted 144-bit run took 129.70 seconds versus 198.18 seconds,
with identical factors and all 513 relation rows. This is a same-policy speedup,
not an overall best-policy timing record. The complete score remains **1/2**.
Run `make verify-batch-smooth`.

## Results and evidence

| Claim | Status | Details |
| --- | --- | --- |
| Saving exponent `1/2` | **Lean-verified** | [Final theorem](formal-appendix/HalfTheorem.lean), [proof modules](formal-appendix/README.md) |
| Saving exponent `1/30` for the k6 machine | **Lean-verified** | [Proof](research/verified-k6.md), [machine theorem](variants/k6/FactoringK6/Certified.lean) |
| Sharper k6 saving `log(12/11)/log(12) > 1/29` | Ordinary proof, not Lean-verified at this rate | [Rate derivation](research/k6-sharp-rate.md) |
| Saving exponent above `1/2` | **Not established** | [Research index](research/INDEX.md) |

The [uniformity perspective](research/uniformity-and-half-frontier.md)
addresses the hard-limit question: assuming `P=BPP`, standard padding gives
deterministic `exp(O(sqrt(n log n)))` factoring. This conditional argument
and the existence of nonuniform successful tapes give no new qualifying score.

Arithmetic research now uses [inverse sections](research/dyadic-inverse-sections.md)
and an [interval orbit-gcd interface with a rational-part projection](research/interval-orbit-gcd.md).
It supplies [exact complete-orbit aggregation](research/dyadic-orbit-aggregation.md)
and [global modular input circuits](research/interval-packet-circuits.md)
with `O(n+p)` arithmetic nodes and polynomial generation cost.
Recovering the genuine rational common factor within the desired time
bound remains open; these results do not improve the factoring score.
The [deadline reduction](research/output-sensitive-gcd-deadline.md) sharpens
the conditional consequence: the proposed output-sensitive rational gcd
primitive would imply polynomial-time factoring. It is enough to prove a
fast termination bound on coprime inputs; that bound is also unresolved.
[Normalized integer-center jets](research/normalized-interval-jets.md)
also separate huge local powers of two from polynomial-cost unit data.
The [global pole trace](research/interval-pole-trace.md) gives an exact
polynomial-precision factor-pair filter and closed bilinear reduction rules.
Its current implementation still enumerates the dyadic state tree; efficient
global state reduction remains open.
The [bilinear pole-product packets](research/bilinear-pole-packets.md) now
construct the whole product over that state family modulo `2^p` with
`O(n+p)` ordinary arithmetic nodes and polynomial generation cost.
The unit-trace evaluator is exact; the singular factor-count normalization
remains unresolved. Run `make verify-bilinear-pole-packets`.
The [circuit moment theorem](research/pole-cauchy-circuit-moments.md) extends
inverse contractions to division-free polynomial circuits of unrestricted
degree, using a quotient algebra of dimension polynomial in precision.
The [rational resolvent theorem](research/pole-cauchy-rational-resolvents.md)
computes complete two-axis contractions from four normalized affine products,
including even integer denominators. A polynomial combination that recovers
the factor-count scalar remains unconstructed. Run
`make verify-pole-cauchy-circuit-moments verify-pole-cauchy-resolvents`.
The [Cauchy content-slope reduction](research/pole-cauchy-content.md) gives
a second global target using two implicit determinant contents, with cheap
padding and explicit Vandermonde normalization.
An exact trace involving the inverse and entrywise square gives the count directly;
efficient evaluation of that normalized scalar is unproved.
The [full Cauchy Smith-profile theorem](research/pole-smith-stabilization.md)
removes any matrix-size-dependent stabilization precision and gives a direct
rectangular maximal-minor slope. It does not identify the unknown factor pivots
or supply a compressed global evaluator. Run `make verify-pole-smith`.

<!-- pole-implicit-checkpoint:start -->
The [pole-Cauchy stabilization theorem](research/pole-smith-stabilization.md)
and [inverse conditioning bounds](research/pole-cauchy-conditioning.md) give exact
structural and precision controls. [Polynomial inverse moments](research/pole-cauchy-polynomial-moments.md)
reduce complete contractions to short product jets.
[Regular analytic weights](research/pole-cauchy-analytic-weights.md) and
[joint analytic kernels](research/pole-cauchy-analytic-kernels.md) now extend
that algorithm to succinct arithmetic circuits, including binary-encoded powers
and unit inverses, without enumerating the matrix indices.
The [sharp selector precision theorem](research/regular-analytic-selector-precision.md)
also proves the exact denominator cost `2^(h+1)-h-2` for a residue-class
parity selector in this circuit representation. It is a representation bound,
and the complete factoring score remains **1/2**.
[Explicit arithmetic-progression unions](research/pole-cauchy-domain-precision.md)
also admit polynomial precision: for `c_A,c_B` component intervals the
inverse loss is at most `(k-1)(2c_A+2c_B-3)`, with a sharper exact residue sweep.
This extends the implicit joint-contraction algorithm beyond balanced domains.
[External rational poles](research/pole-cauchy-rational-poles.md) and
[higher pole orders](research/pole-cauchy-repeated-poles.md) are now evaluable
with all denominator losses charged. [Quadratic poles](research/pole-cauchy-quadratic-selectors.md)
give exact dyadic residue selectors outside the regular unit-inverse circuit
model. A [polynomial selector lift](research/pole-cauchy-residue-blocks.md)
computes complete residue-block sums of `2^L M^-1` modulo `2^p`, in time
polynomial in the encoded domains and numerical precision parameters.
The missing step is a sufficiently compressed rank or Smith-count algorithm;
inverse block sums do not by themselves count factor pairs.

A [batched rational-pole engine](research/pole-cauchy-pole-polynomials.md) contracts explicit quadratic
coefficient arrays using just two product jets. It incorporates coefficient
divisibility before choosing precision and replaces the `p^2` termwise
residue-block calls by one contraction. Exact dense checks and nonzero
402-bit-domain controls through nine output bits are retained.

The [general rational-weight compiler](research/pole-cauchy-rational-weights.md)
now handles every explicitly written separated pair `A(a)/P(a), B(b)/Q(b)`
with no domain poles. A valuation tree has polynomially many leaves even
for repeated roots; quadratic selectors turn each local quotient into
proper pole arrays. [Regular joint kernels can also be restricted to
residue blocks](research/pole-cauchy-restricted-kernels.md).
A separate [representation theorem](research/rational-interval-selector-bound.md)
proves `2^K <= 4 D_* H` for a single rational upper-half selector, where
`D_*` and `H` measure its explicit denominator. This distinguishes low-bit
selection from ordinary interval selection without asserting a factoring
lower bound.

The [spectral trace interface](research/pole-cauchy-spectral-trace.md)
recovers the exact factor-pair count from `r` supplied traces of `S=M M^T`,
with largest numerical power `2r-1` and `r=bits(m)`. Increasing the displacement
reduces the requirement to one trace or one entry sum. Those global joint-kernel
aggregations remain unevaluated at a qualifying compressed cost. Explicit rank-two
displacement generators and a leading inverse partial involution `K^3=K` are also
proved, while their unknown integer multiplicities remain explicit.
<!-- pole-implicit-checkpoint:end -->

The [research index](research/INDEX.md) links proofs, exact evidence, and replay commands.

The [joint raw-boundary theorem](research/joint-top-boundary-rank.md) extends
the original separate-wire obstruction to arbitrary joint factorizations of
all original edges between each group pair. Its exact tradeoff is
`H + V/4 >= n`; internal constraint solving and new weighted representations
remain outside its hypotheses. Run `make verify-joint-top-boundary`.

The result table distinguishes verified machine bounds from conditional constructions and finite experiments. A merged PR or a passing finite test does not imply independent mathematical review.


<!-- beyond-half-checkpoint:start -->
## Arithmetic and tensor checkpoint, November 1

The [consolidated research guide](research/beyond-half-research-checkpoint.md)
connects the ordinary proofs, exact programs, and remaining obligations.
The complete factoring score remains **1/2**; this checkpoint adds no Lean
formalization.

A constructed prime-form power selector now finds smooth leading norms using
an exact primorial filter and CRT composition. Its adaptive continuation
[factors a fresh 144-bit semiprime](research/adaptive-power-window-selection.md)
after 8,232,730 norm tests, including the exhausted first phase. It retains
successful rows and extends only unfinished targets. Both input primes and all
accepted relations are certified. A uniform bound on the required power window
remains open.

The wider arithmetic route has explicit prime-form generators, uniform
represented-smooth-norm existence, unsigned-norm decoding, weighted independent
relations, and determinant-free extraction. With a quasipolynomial base, a full
batch of supplied norms can be processed in `2^((1/2+o(1))b/log_2 b)` time.
Generation and existence use the pinned all-number-field Hecke zero-free theorem.
That existence result does not imply that a reduced leading coefficient is smooth.

State-dependent channels, a polynomial subspace oracle, and blocking improve
exact tensor certificates. An [explicit eight-letter crossing](research/eight-letter-reflection-retraction.md)
has a rational-type certificate `rho<0.72` and beats its matched coordinate
restriction. The certified trial rate remains above `1/sqrt(2)`.
The new pressure restrictions have their stated scope; they are not lower
bounds on integer factoring.

Run `make verify-beyond-half` for the core checkpoint replays. The
[verification record](evidence/beyond-half-verification-20261101.json) pins the
executed suites and source hashes. `make verify-adaptive-power-full` repeats
all 8,232,730 candidates of the 144-bit search; `make verify-reflection-codes-full`
repeats the separate exhaustive 65,519-code classification. Individual notes
state which checks replay certificates and which regenerate the complete search.
<!-- beyond-half-checkpoint:end -->


## Verification

Before starting overlapping work, run `python3 tools/find_overlap.py KEYWORD`
with an authenticated GitHub CLI. It searches discussion posts, comments and
replies, open issues/PRs, historical issue #8 comments, and local research
Markdown. See [the research workflow](autoresearch.md#workflow) and the tool's
coverage report; matches are leads, not proof of novelty or active claims.

From the repository root, run:

```sh
make verify               # Python unit tests, exact certificates, and local link check
make verify-frontier      # Focused research certificates
make verify-post-half     # Larger optional rank tests
make verify-docs          # Local Markdown-link check
make verify-arithmetic    # All arithmetic-framework checkers
make verify-schur-fermion # Complete Schur entry formulas and exact comparisons
make verify-schur-central # Hook-residue identities and full Schur ranks through weight 8
make verify-schur-gaussian # Exact digit reduction and fixed Gaussian crossing seed
make verify-schur-sparse # Sparse down-sets and two-scale offset interval reductions
make verify-fixed-modes # Arbitrary mode sets, independent occupations, and exact endpoints
make verify-cancellation-modules # Exact module states and word-rate certificates
make verify-minimal-modules # Exact boundary reduction and complete word-rate hierarchy
make verify-exact-modules # Complete rational lattice orbits and exact pressures
make verify-cancellation-crossing # Strict certificate gain for a complete crossing
make verify-block-smith # Exact profiles for ordered crossing blocks
make verify-cubic-full-images # Full cubic coefficient spaces and reduced ranks
make verify-linear-smith # All binary scalar pencils and the uniform linear barrier
make verify-type-pressure # Exact common-type optimization and physical gauges
make verify-crossing-interpolation # Exact synthesis with prescribed local row lattices
make verify-smith-oracle # Exact weighted profiles and the scoped sparse entropy dual
make verify-native-smith # Optional compiled F4 search controls; needs a C++17 compiler
make verify-native-rectangular # Optional rectangular paired-type search controls
make verify-controlled-smith # Controlled monomial crossings and the qutrit control
make verify-sparse-holonomy # Exact singleton profiles and Gaussian entry identities
make verify-charge-pruning # One-sided charge changes and typed Gaussian block bounds
```

The Python checks use Python 3 and its standard library. Some optional checks use NumPy for speed. The exploratory mode-pressure search uses NumPy and SciPy. The optional linear Smith catalogue and Schur central-rank checks also require a C++17 compiler.

For the Lean builds, install the stated Lean release, Lake, Git, and the prerequisites in the corresponding scripts. The scripts fetch pinned upstream revisions.

```sh
make verify-half-appendix # Lean 4.33.1; complete 1/2 theorem
make verify-k6           # Lean 4.35.0; complete 1/30 theorem
make verify-formal       # Lean 4.35.0; original local lemmas
```

The [GitHub Actions workflow](.github/workflows/verify.yml) runs five groups on PRs and on `main`: exact Python checks, the larger arithmetic relation certificates, original Lean lemmas, k6, and the half-rate theorem. The half-rate script audits the axioms of its final theorem and supporting declarations.

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

- [F4 hyperoval crossing](research/f4-hyperoval-crossing-bound.md): exact rho < 0.710005861, slightly improving the pinned coordinate 96-letter bound. Both complete cuts and an explicit type with denominator 2^36 are certified. This remains above 1/sqrt(2), so the complete factoring exponent stays 1/2.

- [Exact 96-letter pressure certificate](research/exact-natural-pressure-intervals.md): the complete saved profile certifies rho < 0.710005874881191, still above the required 1/sqrt(2). [Width-sensitive Boolean constructors](research/width-sensitive-prime-constructor.md) reduce preparation cost without supplying the missing crossing improvement. See the [session record](research/escape-half-session-20261031.md) for scope and replay costs.

- [Fixed-prime scalar transfer](research/fixed-prime-factor-transfer.md): polynomially many scalar clique-count queries suffice under the stated evaluator hypothesis. [Certified divisor preprocessing](research/certified-divisor-oracle.md) concerns the maximum divisor count below a bound. These interfaces and the accompanying finite construction controls do not improve the complete 1/2 exponent. The parallel [fixed-characteristic recovery proof](research/fixed-characteristic-factor-recovery.md) is retained; the two notes share the power-tuple and cyclic-orbit mechanism.

- [PR #30 recovery and Smith checkpoint](research/INDEX.md#retained-pr-30-recovery-and-smith-checkpoint), including a [parallel integral Schur derivation](research/integral-schur-crossing-pr30.md): fixed-characteristic recovery interfaces, exact constructions and scoped Smith/gauge controls. Both Schur notes and their attributions are retained. No exponent above `1/2` or new Lean result is claimed.

- [Ward/Pfaffian entry queries](research/ward-pfaffian-entry-queries.md), [simultaneous seed endpoint](research/simultaneous-seed-endpoint.md), and the [PR #29 checkpoint](research/INDEX.md#retained-pr-29-ward-and-crossing-checkpoint): exact entry constructions, effective interfaces and scoped crossing barriers. Entry evaluation alone is not a favorable rank bound. `verify-constructions` is in the default gate; larger native/extended screens remain separate. No stronger complete factoring exponent is claimed.

- [Full binary conservation-kernel barrier](research/full-conservation-rank-sum-barrier.md) and [PR #27 checkpoint](research/INDEX.md#retained-pr-27-conservation-and-crossing-checkpoint): ordinary rank-sum theorems, transfers and exact finite controls. The all-kernel characteristic-two theorem and the [elementary all-code bound](research/difference-support-rank-floor.md) have distinct scopes. No general factoring lower bound, new Lean result, or exponent above `1/2` is claimed.

- [Vanishing-slack transfer](research/vanishing-slack-score-transfer.md), [optimized divisor bound](research/optimized-divisor-count-bound.md), and the [PR #26 checkpoint](research/INDEX.md#retained-pr-26-effective-transfer-and-crossing-checkpoint): effective interfaces, arithmetic bounds and scoped crossing work. Universal verified search does not supply a new numerical exponent. See `make verify-further` and the retained optional targets; no new Lean or stronger complete factoring bound is claimed.

- [Rational-family specialization](research/rational-family-transfer.md), [exact modular half endpoint](research/modular-half-endpoint.md), and the [PR #24 checkpoint](research/INDEX.md#retained-pr-24-effective-family-and-finite-rank-checkpoint): effective interfaces and scoped finite or uniform rank barriers. No exponent above `1/2` or new Lean result is claimed. Run `make verify-pr24` for the checkpoint checks.

- [Sublinear working-space refinement](research/sublinear-space-half-rate.md) and [direct exact divisor counts](research/direct-sublinear-divisor-counts.md): PR #23 supplies ordinary time-space proofs and exact finite controls, retaining the half-rate time bound. Its [full checkpoint index](research/INDEX.md#retained-pr-23-time-space-and-crossing-checkpoint) also preserves related crossing constructions and scoped obstructions. Integration is not independent mathematical acceptance; no new Lean or stronger exponent is claimed. Run `make verify-pr23` for the source checkpoint checks.

- [Weighted fixed-characteristic transfer](research/fixed-characteristic-transfer.md), [polynomial Pfaffian profiles](research/modular-pfaffian.md), and [shared-mode Schur reduction](research/fermi-sea-schur.md): retained PR #22 constructions, effective evaluation, and scoped obstructions. No improved exponent is claimed. Run `make verify-modular` for the extended exact controls.

- [Coarsened grid check](research/coarsened-grid-validation.md): a finite signed network-to-padded-clique check, including internal edges and grouping.

The [auxiliary extension](research/auxiliary-binary-gauge-floor.md) carries
the same floor to matrix-valued binary seeds of every fixed auxiliary
dimension, including independent rational auxiliary gauges at the two
places; block certificates and module word-rate bounds remain outside it.
Run `make verify-auxiliary-binary-floor`.

For larger alphabets, [cell-set minor caps](research/cell-minor-caps.md)
strengthen the [shared-minor relaxation](research/shared-minor-profile-relaxation.md):
every set of k cells gives an endpoint minor that bounds each Smith-profile
prefix by kd plus the k smallest cell energies. This exactly excludes the
earlier four-symbol relaxation target. Sharper negative relaxation targets
remain at five and six symbols; they are construction targets, not crossings.
Run `make verify-cell-minor-caps`.

The [binary site-symmetry barrier](research/binary-site-symmetry-barrier.md)
gives an ordinary actual-rank proof of `rho >= 3/4` for every decomposition
whose binary pieces are invariant under all simultaneous permutations of
tensor positions. It covers scalar seeds of every degree, arbitrary coefficient
recombinations and mixtures of scalar curves, over any fixed field,
including characteristic two. A universal degree-four polynomial and
Weyl filtrations supply the obstruction at every tensor length.
It extends to a common partition into any sublinear number of position
classes, with arbitrary coupled coefficient functionals between classes.
Run `make verify-binary-site-symmetry` for exact algebra and finite controls.
The [matching-curve controls](research/matching-curve-controls.md) show
why escaping common-minor tests does not escape this broader theorem.
New physical code restrictions and constructions that break the stated
site symmetry remain open; no stronger factoring exponent is established.

The [full auxiliary coefficient-span experiment](research/auxiliary-coefficient-span.md)
tests ordered noncommuting constructions outside that symmetry hypothesis.
Complete finite recombination improves several raw rank maxima, but exact
rank-sum lower certificates leave the tested finite-block envelope above
the half-rate threshold. Run `make verify-auxiliary-coefficient-span`.
Nine full extension-field controls over `F_4` retain the same finite
minmax thresholds; run `make verify-auxiliary-extension-span`.

The [common-cofactor barrier](research/common-cofactor-barrier.md) closes
the asymmetric ternary quadratic family for actual coefficient ranks.
All 32,265 normalized seeds, covering 12,852,225 actual coefficient seeds,
have a shared projected minor that forces
`rho >= log(3)/(2 log(2)) = 0.792481...` at every tensor length.
The exact full certificate replay is `make verify-common-cofactor`.
The binary site-symmetry theorem above additionally covers seeds without
this certificate at its weaker, but still post-half-excluding, floor.
The factoring exponent is still one half.

The [joint determinant barrier](research/joint-determinant-barrier.md)
uses Weyl filtrations in every characteristic to bound actual coefficient
ranks. It rules out `rho < 4/5` for binary quadratics with two pure
interpolation values and cubics with three, including the entire
135,005,697-seed ternary cubic class. Equal determinants give floor one.
The asymmetric rank-two values have a weaker determinant floor of `2/3`;
the common-cofactor result above now closes the full ternary quadratic
class with a stronger certificate. Run `make verify-joint-determinant`.

The [constant transpose-intertwiner obstruction](research/transpose-intertwiner-barrier.md)
applies to actual tensor-power coefficient ranks: a common invertible
intertwiner forces both cut ranks to agree, ruling out leading contraction
for polynomially many pieces. It excludes 32,454,657 ternary cubic seeds.
Actual jets of 200 retained witnesses are checked through three copies;
the favorable two-copy values do not persist in that finite comparison.
Run `make verify-transpose-intertwiner` for complete finite coverage counts
and independent matrix controls. No stronger factoring exponent is claimed.

The [universal multipoint Smith floor](research/multipoint-smith-floor.md)
shows that the local-jet pressure criterion needs `rho > 3/4` over every
field, alphabet size, polynomial degree, allocation, and permitted cut scale.
It therefore cannot improve the existing `1/sqrt(2)` crossing. A completed
1,579,014-seed cubic catalogue supplies independent literal controls.
Sharper actual-jet bounds and physical type constructions remain outside
this obstruction. Run `make verify-multipoint-smith-floor`.

The [complete low-rank ternary quadratic catalogue](research/multipoint-low-rank-seeds.md)
excludes binary seeds with cut-rank product at most three at both finite
interpolation points: all 32,265 representatives satisfy the exact multipoint
pressure floor `rho >= 3/4`, over every degree allocation and two-cut scale.
The [pure-value subfamily](research/multipoint-rank-one-seeds.md) has its
own smaller replay. Larger local ranks and additional field points remain
open. Run `make verify-multipoint-rank-two` for the broader catalogue.

The [weighted code-block theorem](research/difference-support-rank-floor.md)
excludes the full elementary crossing coefficient span on **every pair of
binary physical codes**, over every field. For equal alphabet size q,
the rank-sum is at least q^(3/2), so polynomially many pieces force
rho >= 3/4. This closes the scalar elementary-code route, including
staggered paired codes; new kernels and code-coupled auxiliary dynamics
remain open. Exact controls cover 273 codes, 225 code pairs, and signed
ternary-field cases. The complete factoring score remains 1/2.

The [smooth-norm relation reduction](research/smooth-norm-relations.md)
accepts primitive smooth represented norms without a balanced-factor
requirement. A [gap extractor](research/smooth-norm-gap.md) handles the
older balanced interface in polynomial time with explicit slack.
For one supplied factored norm, the [exact representation search](research/factored-norm-search.md)
returns primitive coordinates or proves that candidate impossible in
`Y poly(b) + 2^O(b/log b)` bit time. Selecting a guaranteed short list of
candidate norms remains open. Run `make verify-smooth-norms`.

The [reflected neutral construction](research/neutral-reflected-crossings.md)
gives explicit homogeneous duals and complete Pfaffian crossing tensors,
including unbalanced fixed-weight sectors. All homogeneous completions of
an eight-letter weight-32 sector form an [88-dimensional tensor family](research/homogeneous-dual-completion.md)
with a literal endpoint. Run `make verify-neutral-duals`.

The [F4 frame-recentering solver](research/f4-recentered-crossing-lattices.md)
allows a nonbinary tensor to provide the regular frames for another exact
full-coefficient search. Complete untruncated profiles are certified by
finite polynomial identities and independently checked binary jets.
Run `make verify-f4-recentered`.

The [odd-cyclic matrix normal form](research/odd-cyclic-matrix-crossings.md)
reduces its two cuts to repeated matrix blocks and scalar Fourier symbols.
Its unweighted profile certificate has an all-degree floor `rho >= 3/4`;
the recorded weighted finite search also gives no strict witness. Run
`make verify-odd-cyclic`. These new constructions and reductions leave
the achieved factoring exponent at `1/2`.

The [modular parameter-descent theorem](research/modular-parameter-descent.md)
allows nonregular monomials to be removed from both the state and operator
parameter. Its exact reduced-dual compiler preserves every crossing piece;
it does not improve the factoring exponent. Run `make verify-modular-descent`
to compare it with the original unreduced recursion.

The [difference-unique code obstruction](research/sidon-code-rank-obstruction.md)
shows why sparse additive graph codes do not rescue the elementary crossing:
its full multivariate coefficient span has rank-sum at least q(q-1)/sqrt(2).
This includes every binary cubic field graph, but not general nonlinear codes
or new auxiliary mechanisms. Sixteen finite controls replay with
`python3 tools/check_sidon_code_crossings.py`; the factoring score is unchanged.

The [typed rational strictification theorem](research/typed-rational-crossing-strictification.md)
extends the binary pressure floor to scalar rational kernels with only
zero/endpoint denominators and crossing endpoints at every paired type
and length. A paired physical gauge preserves the complete pressure.
For binary characteristic-two kernels, twelve exact charge circuits make
length four a sufficient and sharp endpoint test; ramification by eight
always produces a strict crossing. Run `make verify-binary-typed-endpoints`.

The [amplified quartic pole theorem](research/reflection-paired-quartic-poles.md)
permits a growing word of disjoint gates with poles of order M+7 per gate.
Simultaneous Pfaffian row factors prove every typed endpoint, and exact
common-factor cancellation reduces the actual local pole to at most M-1.
The stated numerator window keeps the original effective degree. The
construction still needs a negative full-pressure witness. Run
`make verify-paired-quartic-poles`.

The [nonlinear downset-code construction](research/downset-code-crossings.md)
preserves the crossing on a restricted physical alphabet. Its column
matrix has an exact global Smith diagonal and a constant basis change
to smaller weighted disjointness matrices, allowing actual recombination
ranks to be checked without inferring them from a Smith bound.
The formulas also permit different downsets on the two crossing pairs;
all 361 ordered three-bit downset pairs are checked exactly.
Complete small no-adjacent-ones code experiments give no new rank-sum
threshold; the [weighted block theorem](research/difference-support-rank-floor.md)
subsequently closes every coordinate-code pair for the elementary seed.
Run `make verify-downset-crossings verify-regular-code-crossings`.

The [coherent bounded-spin construction](research/bounded-spin-crossing.md)
uses height-three Clebsch--Gordan paths over the fixed field `F_5`.
It has Fibonacci many labels, nondegenerate paired encoders at every
length, an exact crossing endpoint, and a finite-state entry algorithm.
It also supports exactly `2^lambda` labels with degree `O(lambda)`.
A companion integer identity has scale `12^n` and `O(lambda)` coefficient
bits. Modular rank gains do not automatically lift to rational ranks.
A uniform favorable rank bound remains open. The
[eight-site rank algebra](research/bounded-spin-rank-algebra.md) gives
exact ranks over every characteristic-five extension and minmax rank
product 2304; that finite block cannot improve the rate. Run
`make verify-bounded-spin` and the optional NumPy target
`make verify-bounded-spin-algebra`.
The [second-coefficient formula](research/bounded-spin-second-coefficient.md)
expresses `3C_2` as a compressed diagonal transposition sum plus
`k(k+1)I` at every length `n=2k`. At eight sites this sum and the replica
swap generate the certified row algebra, including its nilpotent part.
Run `make verify-bounded-spin-second` for independent entry controls.
The [spin-zero jet identities](research/bounded-spin-three-jets.md)
hold at every length and bound the coefficient span by `n-2` for even
`n>=6`; ordinary cut-rank growth remains open.
A [direct orthogonal-group argument](research/four-qubit-orthogonal-generators.md)
now supplies the degree `2,4,4,6` four-qubit generators in characteristic
five, with explicit restrictions to the pencil. It supplies no cut-rank bound.
An [enlarged six-vertex construction](research/bounded-spin-six-vertex.md)
exposes the complete monomial span with degree `O(n^2)` and the same
strict endpoint. At eight sites it adds a seventh tensor direction;
the noncommuting cut operators have invariant commutator modules of
dimension 28. Uniform favorable ranks remain open. Run the optional
`make verify-bounded-spin-six-vertex`.
The [modular interface](research/bounded-spin-modular-interface.md)
connects these constructors to the existing fixed-characteristic transfer:
dense rank factors meet its cost bound, and the original coefficient
space permits automatic selection even from existential decompositions
over extensions of degree `2^o(lambda)`. A favorable uniform rank
estimate is still missing.
The [constant-state spin projector](research/bounded-spin-projector.md)
uses thirty states, or sixty for exact binary label selection, and
preserves the encoded rank problem under ambient extension. An exact
coefficient involution also identifies the two cut operator spaces.
Run the optional `make verify-bounded-spin-projector`.
The [bounded-support rank factors](research/bounded-spin-support-ranks.md)
give explicit `O(L^(4/3))` channels for the shared extreme coefficient,
below the post-half scale for that piece. Sublinear corner bands still
cannot span the crossing; a complete favorable decomposition remains
open. Run `make verify-bounded-spin-support`.
The [correction-module analysis](research/bounded-spin-correction-modules.md)
gives the next algebraic target and an exact ten-site replay: two F5
generators are needed, the vertex module has dimension 170, and the
original coefficient span has sharp minmax product `225^2`. A small
twelve-site certificate rules out the simpler linear pencil-plus-vertex
formula. Uniform algebra containment and two-cut lifting remain open.
The optional `make verify-bounded-spin-growth` requires FLINT.
The [paired-path construction](research/paired-path-crossings.md) replaces
the classical spin formulas with arbitrary invertible local splittings
and their inverse duals. An explicit F16 family retains Fibonacci label
growth, polynomial-time entries and both strict endpoints. Exact small
checks cover all twelve valid field parameters; favorable uniform
two-cut ranks remain open. Run `make verify-paired-paths`.
Its enlarged six-vertex span has explicit uniform support factors:
all three extreme pieces have `O(L^(4/3))` channels and can be
independent. A complete favorable decomposition still remains open.
The [height-cutoff tradeoff](research/paired-path-height-tradeoff.md)
extends the construction to any fixed cutoff: the individual extreme
rank exponent approaches one as the cutoff grows. The field stays fixed
for each chosen cutoff. Run `make verify-paired-height`.
The [quantum block-braid compiler](research/quantum-block-braid-crossings.md)
uses local cup and projector identities to prove a different all-length
crossing endpoint. It has degree `O(lambda^2)` and entries computable
within the allowed `2^O(lambda)` construction budget. Its complete
six-site coefficient span has dimension 33; favorable two-cut ranks
remain open. Run `make verify-quantum-braid` and the optional complete
`make verify-quantum-braid-profiles`.
The [exact path quotient and phase descent](research/quantum-braid-path-quotient.md)
replace the ambient order-five braid calculation by a bounded-path
calculation over F4. A period-five weighting gives a strict crossing
of degree at most `n^2`, omits one residue class of coefficients and
preserves both ranks of every retained piece. A favorable uniform rank
bound remains open. Run `make verify-braid-path`.
The optional `make verify-braid-path-profiles` compares all four- and
six-site coefficients with the ambient constructor and computes the
complete eight-site profile. The strict F4 span has dimension 51 at
eight sites; its retained central coefficient has full rank in both
cuts. This remains a structural checkpoint, without a favorable
recombination. The optional batch replay needs NumPy.
The [reversal-stable binary selection](research/quantum-braid-reversal-nilpotence.md)
gives a closed form for the leading coefficient and makes its normalized
rank-one operator square-zero at every binary width. Polynomial
conjugations preserve both endpoints; exact small examples cancel
leading degrees without establishing a stronger pressure bound.
Run the NumPy check `make verify-braid-reversal`.
The [noncommuting matrix auxiliary](research/quantum-braid-matrix-auxiliary.md)
preserves the strict F4 crossing and gives an explicit identity channel
as its leading coefficient. Both cut ranks of this piece are exactly
`L` at every length on a common alphabet containing the zigzag path.
The whole polynomial has degree `3n^2/2`; its remaining coefficients
still need a favorable decomposition. Complete small spans have
dimensions16,55,97. Run the NumPy check `make verify-braid-matrix-aux`.
An [all-extension minor certificate](research/quantum-braid-rank-one-loci.md)
determines both six-site rank-one loci: the scalar family has only
its top row-rank-one line and constant column-rank-one line; the
matrix auxiliary has no nonzero row-rank-one direction and only
the constant column-rank-one line. Higher ranks and longer paths
remain open. Run `make verify-braid-rank-one` (NumPy).
The [nilpotent central-seam family](research/nilpotent-seam-crossings.md)
now gives complementary asymmetric coefficient families: the highest
allowed central edge has row rank `L^(1/2+o(1))`, below the
post-half small-cut threshold, and the opposite edge has polynomial
column rank. Exact eight-site pairs include `(8,169)` and `(134,2)`.
The minimum-cut budget rules out an edge-only completion at large
widths; a mixed decomposition must supply balanced middle pieces.
Their checked edge spans also miss the small crossing. Run the NumPy replay
`make verify-braid-nilpotent-seam`.
An [exact all-degree rank envelope](research/nilpotent-seam-rank-envelope.md)
now bounds both cuts throughout that coefficient grid. An eight-counter
integer recurrence includes the residual vacuum dimensions; exact rational
certificates give middle rank bounds
`(L^(5/4+o(1)),L^(9/8+o(1)))` after vacuum-interface compression.
Shoulder windows and cut reflection provide other candidate aspect ratios.
A compatible graded spanning identity remains open; weighted cut budgets
rule out completing the crossing from just the edges and middle window.
Run `make verify-braid-rank-envelope` (NumPy for the complete finite replay).
The [triplet-projector seam](research/fusion-projector-seams.md) supplies
another source of middle pieces: its row factors have size polynomial
times `L`, and selected boundary interfaces use `v_b` vacuum channels
instead of `4^b` physical channels. A certified range has column rank
`L^(3/2+o(1))`; an eight-site piece has ranks `(17,28)`.
The proposed mixed list still misses the five- and thirteen-label crossings.
A uniform mixed spanning identity remains open. Run the NumPy replay
`make verify-braid-fusion-seam`.
The [invertible three-strand projector](research/invertible-fusion-graded-seams.md)
extends this to graded central ranks `F_(k-1)` (rank one when k=0),
where k strands are left free. Its uniform row-exponent range is0..2,
with the same vacuum-interface column factors and an explicit integer
cost filter. The unfiltered prefix family still misses the five- and
thirteen-label crossings. Run `make verify-braid-graded-seam` (NumPy).
An [exact fused spectral slide](research/fused-spectral-transport.md)
transports the three-strand projector with compatible shifted parameters.
The ambient identity retains a nonzero negligible height-four residual;
in the vacuum quotient the transport is a scalar times the constant braid.
Exact polynomial normalization gives two-term transport factors at lengths
`5*2^s`. A global kernel, spanning identity and favorable complete rank sum
remain open. Run `make verify-braid-spectral-fusion` (NumPy).
A [complete fusion-occupancy resolution](research/fusion-occupancy-resolution.md)
groups every vacuum/triplet placement into `n+1` orthogonal sectors of
exact rank `binom(n,k)F_(k-1)` (rank one at k=0). Including all boundary
degrees gives an explicit scale-one crossing sum at every length, with
both cut bounds and every occupancy convolution channel retained.
A favorable complete rank sum remains unproved. Run
`make verify-braid-fusion-occupancy` (NumPy).
This also gives [exact leading auxiliary ranks](research/nilpotent-leading-rank-automaton.md):
a closed 22-state fusion automaton proves rational generating functions
at every even length. The leading central ranks have half-power scale
in the label count and sharpen the earlier counting envelope; full
contracted coefficients can still vanish. This does not close the
complete crossing rank sum. Run `make verify-braid-leading-automaton` (NumPy).
The [full Catalan cup frame](research/catalan-cup-crossings.md) retains
all noncrossing primal cups using a unitriangular coordinate decoder.
It gives a strict crossing over F2 even with zero primal Gram matrix,
and binary width lambda uses only `lambda+O(log(lambda+2))` physical strands.
At the old order-five parameter it retains nonzero radical channels,
so quotient-only formulas require their ambient residuals. Complete
small profiles and independent physical evolution pass; a favorable
uniform two-cut rank sum remains open. Run `make verify-catalan-cup-crossing`.
Its [nilpotent parameter jets](research/catalan-cup-jets.md) deform the entire
primal/decoder/braid system over `F2[u]/u^s`, retaining every intermediate
coefficient. Positive jets have zero endpoints and enlarge the scalar F2
candidate space with polynomially many pieces for polynomial s. Independent
physical truncated-ring evolution and complete small profiles pass; the
current filtered jet list still misses the crossing at five and fourteen
labels. Run `make verify-catalan-cup-jets`.

The [exact nonlinear charge theorem](research/nonlinear-symbol-charge-lattices.md)
characterizes every admissible pair of full-symbol potentials for the
reflection tensor: their Hamming-cube Lipschitz constants must sum to at
most two. At four symbols this is a six-dimensional cross-polytope and
includes nonlinear parity charges. The resulting complete rational
coefficient spaces have exact endpoint and full-profile checks. Run
`make verify-nonlinear-symbol-charges`.

An [exact joint-filter separation](research/joint-boolean-multiplier-separation.md)
strictly improves every scalar Boolean multiplier on three balanced parity
groups, uniformly for physical lengths at least sixteen. The certificate
covers all 254 nonconstant supports and unbounded multiplier degree, field
extensions, and ramification. Its complete critical pressure stays positive;
this is a stronger construction family, with no score increase. Run
`make verify-joint-multiplier-separation`.

The [typed rational coefficient lattice](research/typed-charge-lattice-synthesis.md)
frees every entry under an exact charge condition at the endpoint. Local poles
cancel in every paired-type tensor power; a negative-charge boundary extension
adds further coefficients while preserving that endpoint. At four symbols,
the displayed spaces have 6,144 or 6,237 free F4 coefficients and effective
degree 24. Full local rank follows from the [type-class rank theorem](research/paired-type-endpoint-full-rank.md).
These are constructive search spaces; no better factoring exponent is claimed.
Run `make verify-typed-charge-lattices`.

The [growing quartic construction](research/reflection-charged-quartic-tiles.md)
adds nonquadratic occupation gates on a growing fraction of physical modes.
An exact charge budget also permits specified endpoint poles: they cancel in
every retained coefficient after paired type restriction. The construction
includes its effective degree and restored local profiles; it supplies new
search families, without an improved rate yet. Run `make verify-quartic-tiles`.

The [one-support multiplier reduction](research/boolean-multiplier-automata.md)
removes polynomial-degree and coefficient-field cutoffs from a growing global
Boolean-filter family. Its finite automaton retains continuous ramified
cutoffs and exact periodic tails. [Fifteen complete outer searches](research/global-code-multiplier-searches.md)
over 13,248 explicit code supports on up to 256 modes still select the
baseline. Run `make verify-multiplier-automata verify-global-code-multipliers`.

The [complete-entry larger-alphabet solver](research/full-alphabet-two-cut-affine-synthesis.md)
varies every polynomial coefficient over F4 while preserving the strict
endpoint and simultaneous local row conditions. Two completed four-symbol
searches contain 40,000 target walks and 1,046 full-profile classes; all eighty
retained tensors replay exactly, with no critical crossing. Run
`make verify-full-alphabet-affine`.

The [universal binary Smith-pressure theorem](research/universal-binary-smith-floor.md)
excludes every Laurent-polynomial strict binary crossing seed under the
two-place Smith criterion, over every coefficient field and with all
independent physical-letter heights. Shared minors couple the two cuts;
a monotonicity argument reduces arbitrary weights to a certified compact
inequality. This closes the binary higher-degree and row-frame search
route without changing the achieved factoring score. Run
`make verify-binary-smith-floor` for the exact certificate and full-minor controls.

The [nonresonant Boolean filter theorem](research/nonresonant-boolean-filter-floor.md)
excludes arbitrary nonhomogeneous Boolean supports in the stated two-power
family outside a finite resonant window. A separate [affine tile construction](research/affine-boolean-filter-tiles.md)
recomputes full transformed polynomials after Boolean relabeling; its complete
124,080-case catalogue contains 65,475 retained profiles, with no gain in the
eight recorded full-cost searches. Run `make verify-nonresonant-filters
verify-affine-filter-tiles`.

The [prime-characteristic symmetric transform](research/symmetric-prime-multilevel-filters.md)
retains all digit multiplicities and endpoint signs while compressing a
correlated block to occupation weights. Four [completed tile catalogues](research/symmetric-tiled-filter-searches.md)
contain 3,830 exact profiles and 32 complete cost searches across binary
and qutrit intervals; every best point still selected the baseline. Run
`make verify-symmetric-prime-filters verify-symmetric-filter-catalogues`.

The [symmetric multilevel filter construction](research/symmetric-multilevel-filter-synthesis.md)
compresses a correlated binary block to one polynomial per occupation weight.
Its exact tensor rule permits different filters on disjoint physical blocks.
A literal four-mode filter gives a strict natural-charge improvement over the
Gaussian cost, including a gain proportional to interval length when tiled;
its complete pressure remains above the required threshold. Run
`make verify-symmetric-binary-filters` for full-state and exact cost controls.

The [cyclic-algebra constructor](research/cyclic-algebra-crossings.md) permits
dense polynomial transforms of a regular full-cycle operator, with a complete
crossing endpoint and unchanged ideal column orders. An exact two-jet
certificate excludes all 1,024 transforms in the first eight-letter schedule
for the bare two-place criterion. The factoring exponent remains `1/2`.
Run `make verify-cyclic-algebra` for independent literal matrix controls;
the linked note gives the optional complete search replay.

The [shift-register construction](research/shift-register-crossings.md) gives
further full-cycle seeds with arbitrary balanced bit partitions. All 8,186
systems through thirty-two physical letters fail the bare two-place test by
exact integer certificates; a separate selected 64-letter schedule also fails.
These finite exclusions leave the score at `1/2` and do not cover arbitrary
word-type gauges or nonlinear feedback. Run `make verify-shift-register` for
the small exact replay; the note gives the complete catalogue commands.

The [copied-wire argument](research/shift-register-copy-wire-barrier.md)
explains a structural limitation of growing register seeds: colour
correlations force enough early Smith channels to give limiting bare
parameter at least `3/4`. It covers every full-cycle register over a
prime-power digit alphabet, with arbitrary nonlinear feedback and coordinate
partitions. Its explicit large-size cutoff leaves a gap above the completed
small catalogues; it does not cover physical gauges or arbitrary bases.
Run `make verify-register-barrier`.

A separate [six-bit nonlinear catalogue](research/quadratic-register-crossings.md)
tests all quadratic lower-tail functions plus the zero-tail indicator:
1,952 full cycles and all 19,520 balanced partitions fail the bare test.
Run `make verify-quadratic-register` to replay the complete finite result.

The [quantum-Pascal construction](research/quantum-pascal-trace-crossings.md)
provides polynomial trace crossings with exact degree `d(d-1)` at odd
alphabet size `d`. A unitriangular-basis argument forces `rho>=3/4` for
its bare and unweighted multipoint profile criteria, including arbitrary
polynomial diagonal vacua. Physical word-type restrictions remain outside
that proof. Run `make verify-quantum-pascal` for exact endpoints through
fifteen letters, direct trace controls, and polynomial determinant checks.

The [matrix-vacuum extension](research/matrix-vacuum-boundary-barrier.md)
allows arbitrary polynomial matrix vacua in the same ladder construction.
A boundary-difference identity filters the zero-cut rank; the complete
7,350-case self-dual rank-two affine catalogue has minimum rank 24 and
bare parameter strictly above `7/8`. This excludes every interpolation
path with those initial vacua. Run `make verify-matrix-vacuum`; the note
also gives the optional complete NumPy catalogue. A larger internal space
with fewer physical ladder powers remains outside the square-basis proof.

The [Steiner construction](research/steiner-cubic-crossings.md) factors the
swap into three controlled involutions and gives a determinant-one cubic
crossing for every Steiner quasigroup. Its complete bare Smith profiles
are independent of the design and force `rho>=3/4` at every size. This
closes the cubic bare-profile search. A stronger
[common-coefficient proof](research/steiner-coefficient-barrier.md) also
excludes arbitrary coefficient recombinations across full words at every
`q>=7`, with no piece-count or position-symmetry restriction. At `q=3`,
a middle-coefficient proof gives asymptotic floor `rho>=0.8065306485`
for the univariate span and a sublinear number of position classes.
Restricted physical types and the fully multivariate three-letter space
remain outside those proofs. Run `make verify-steiner-cubic` and
`make verify-steiner-coefficients`.

The [ordered physical-gauge construction](research/ordered-physical-gauge-crossings.md)
preserves every crossing endpoint while introducing noncommuting auxiliary
matrices and polynomially many coefficients. Exact integer dual potentials
reject all 216 cases of a stated two-dimensional binary control under the
finite-block rank-sum envelope. Other auxiliaries and a uniform subcritical
rank bound remain open. Run `make verify-ordered-steiner-gauges`.

The [ordered Weyl construction](research/weyl-ordered-boundaries.md) writes
the auxiliary word as an explicit quadratic phase. A Fourier argument
classifies boundaries retaining the rank-one zero cut at every length.
All 4,096 normalized F4 choices fail the two-position bare criterion by
an exact four-jet floor `738>9^3`. Beyond bare profiles, complete F4
middle-coefficient tests reject the entire coefficient spans at two
selected boundaries, with rank-product floors 576 and 1950. These tests
leave other boundaries and longer ordered words open. Run
`make verify-weyl-ordered-boundaries` and
`make verify-weyl-middle-coefficients` for the standard-library replays.

The [subgroup-return construction](research/subgroup-return-crossings.md)
preserves both endpoints using noncommuting group words and has a rank-one
starting cut at every length. Its free-group version has an auxiliary
return functional of infinite Hankel rank, yet word reduction gives
effective complete construction in `2^O(m)` time. Exact full binary
coefficient spans fail the finite half target at lengths one through three
for that example and a four-state analogue; useful uniform rank bounds
remain open. Run `make verify-subgroup-return`.

The [growing-clock normal form](research/ordered-charge-area-crossings.md)
expresses ordered Weyl words as integer charge walks with an area weight.
Uniformly generated polynomial-degree fields preserve the exponent under
[whole-piece descent](research/fixed-field-descent.md#4-polynomial-degree-extensions-with-uniform-presentations).
Selecting charge zero and area zero gives a binary family with strict
endpoints and linearly many coefficients. Its complete binary spans fail
the finite half target at lengths two and three; useful uniform rank
bounds remain open. Run `make verify-ordered-charge-area`.

[Joint-pair energy conjugations](research/joint-energy-conjugations.md)
preserve one full weighted Smith profile while changing the other cut.
For the five-letter Gaussian reference, exact commutants reduce the
leading-coefficient-preserving F4 family to 1,296 distinct tensors.
Fraction-free profiles and prefix dominance show that none improves the
reference at its natural heights at any scale. Fixing its next complete
coefficient forces a trivial conjugation. A tangent-space identity also
classifies every formal leading-rank-one deformation near the identity
as a centralizer motion and separate physical diagonal changes.
Other constant components and weaker order-preserving changes remain
open. Run `make verify-joint-energy-conjugations` and
`make verify-joint-leading-tangent`.
The [polynomial dressing reduction](research/joint-energy-polynomial-dressings.md)
also handles arbitrary polynomial degrees in its four-dimensional
nilpotent sector: retaining the original worst pole and rank-one leading
coefficient reduces the seed to a constant conjugate, over every
characteristic-two field. Exact quadratic elimination with Frobenius
roots also closes every square-zero scalar-polynomial energy generator,
of arbitrary rank: the pole bound forces it into that sector, and the
leading rank condition makes the transformation trivial.
Run `make verify-joint-energy-dressings` and
`make verify-joint-energy-riccati`.
The [formal polynomial gauge theorem](research/joint-polynomial-formal-rigidity.md)
removes the single-generator assumption for families through the identity:
two exact direct sums control every degree, and formal lifting reduces
every such pole/rank-preserving gauge to a constant conjugation and a
tensor-preserving commutant unit. Separate components remain open.
Run `make verify-joint-polynomial-tangent`.
The [discrete coordinate lattice calculation](research/joint-monomial-lattice-modifications.md)
allows singular joint monomial gauges that change both profiles. Exact
difference constraints classify every integer charge preserving the
original pole bound into 12 effective tensors; both profiles are
dominated by the reference at every scale. A subsequent all-field
classification also closes every energy-projector modification at
every positive depth under the stated conditions; the only four
projectors preserve the tensor. The one-step infinity conditions alone
also close every normalized scalar multiplier `f(t)=t+O(1)`.
The [single-evaluation algebra theorem](research/joint-evaluation-algebra.md)
also reduces every rational energy-preserving matrix gauge regular at
zero, one and infinity to a constant conjugation: a polynomial minor
`t^1452(t+1)^350` proves a localized algebra identity. A gauge equal to
the identity at infinity leaves the tensor unchanged. Gauges singular
at the distinguished places remain open. Run `make verify-joint-monomial-lattices`,
`make verify-joint-projector-lattices` and `make verify-joint-evaluation-algebra`.
The [singular endpoint loop construction](research/joint-endpoint-rational-loops.md)
evades three-place regularity and gives genuine new polynomial crossings.
Every scalar-rational gauge with a constant square-zero energy generator,
of any rank, reduces to three effective endpoint parameters under the stated
conditions. Two symbolic minors and a positive Bernstein certificate prove
`rho>=3/4` for every nontrivial member over every characteristic-two field
at the natural heights.
Thus it gives no half-rate improvement. Run `make verify-joint-endpoint-loops`.
At a monomial endpoint pole, exact quartic/Frobenius elimination also
reduces every higher-nilpotency energy generator to that square-zero
sector or a tensor-preserving direction, at every positive depth.
Run `make verify-joint-nilpotent-depths`.
The [opposite-filtration lemma](research/joint-filtered-gauge-rigidity.md)
extends zero-profile preservation to a 347-dimensional constant gauge
group, then proves that retaining the original infinity pole forces it
back to the 69-dimensional energy-block group. It also gives the exact
extra pole in terms of the largest energy difference in the inverse gauge.
Run `make verify-joint-filtered-gauges`.

The [fused coefficient-envelope certificate](research/fused-scalar-projection-envelope.md)
tests arbitrary binary recombinations after scalar projection of a complete
two-copy crossing. All `2^21` vectors of the symmetric envelope were checked;
the 17 tensors below the rank-product threshold span a space excluding the
crossing. A different elimination algorithm independently checks the separating
affine hyperplane. This rules out that entire finite envelope, with the full
ordered space, higher powers and physical type restrictions still open.
Run `make verify-fused-scalar-projection` (requires a C++17 compiler).
For the original unprojected coefficient span, the
[primary-character code theorem](research/fused-primary-character-code.md)
proves `rho>=3/4` for every tensor length and every extension of the fixed
F16 field. Ten exact complement-rank certificates establish a weighted
code distance of three; tensorizing that distance bounds every nonzero
piece's opposite-cut rank. This theorem retains the nonsemisimple blocks
and leaves scalar projections and physical type restrictions open.
Run `make verify-fused-primary-code`.
The [energy-refined Newton ceiling](research/joint-graded-newton-ceiling.md)
constrains the zero Smith profile under every rational energy gauge,
including gauges singular at zero. An explicit block-scalar diagonal
reveals characteristic-polynomial invariants hidden by cancellations in
the ordinary characteristic polynomial. The exact ceiling still leaves
room to improve the current zero profile; its attainability and the effect
on the infinity cut remain open. Run `make verify-joint-graded-newton`.

The [group-shear extension](research/group-shear-crossings.md) implements
the swap through three multiplication/division shears over any finite
`p`-group, with a constant nonzero determinant in characteristic `p`.
Exact early-jet certificates reject six cyclic, dihedral, and quaternion
binary seeds for the bare criterion. These finite exclusions do not
settle the family. Run `make verify-group-shears`.

The [fixed-type capacity argument](research/fixed-type-smith-capacity.md)
makes physical gauge optimization convex once the two paired types are
fixed. Entropy-weighted minor flags give global affine supports and an
exact rational certificate criterion. `make verify-smith-capacity` runs
the exact controls. The optional numerical solver reports its parameter
box, primal value, LP estimate, and gap; its current quantum-Pascal example
does not improve the score.

The [exact Gaussian fixed-type certificates](research/reflection-gaussian-fixed-type-duals.md)
prove global minima over every independent mode gauge and positive scale
for four explicit rational physical types, through interval length 32.
The proof uses nonnegative multipliers on exact Clifford minors; it does
not cover the remaining optimization over all physical types. Run
`make verify-gaussian-fixed-types` for exact determinant and stationarity checks.

The [multipoint Hermite transfer](research/multipoint-hermite-crossings.md)
reconstructs effective crossing pieces from jets at several finite points
and infinity, with a shared degree budget and complete coefficient costs.
A general degree-budget argument now proves that its unweighted same-place
criterion requires `rho > 3/4`, including common scalar factors and independent
local scales. Weighted types or stronger actual jet-rank estimates remain
outside that floor.
Its explicit ternary certificate is weaker than the half-rate crossing;
the transfer supplies a different criterion for future seeds. Run
`make verify-multipoint-hermite` for full reconstruction and jet-rank checks.

The [physical split-projection interface](research/physical-split-crossings.md)
constructs smaller or unequal alphabets with exact left/right endpoint
repairs, including coherent combinations over F4. Complete searches cover
5,440 rank-three projections and two twenty-projection energy families.
[Coordinate-code restrictions](research/reflection-code-alphabets.md) add
236 recorded alphabets. Every complete profile is recomputed after projection;
these fixed-gauge screens have not improved the crossing threshold. Run
`make verify-physical-splits verify-reflection-code-alphabets`.

The [growing physical-frame construction](research/growing-gaussian-physical-frames.md)
and [coupled quadratic gates](research/coupled-gaussian-gates.md) preserve strict
endpoints while changing a number of modes that grows with the interval.
Their exact Clifford updates retain determinant-mean changes and rational
denominator costs. Original-tensor replays pass; the completed natural-gauge
frame searches found no gain. Run `make verify-gaussian-physical-frames
verify-gaussian-coupled-gates` for the full controls.

The [Clifford Smith reduction](research/clifford-gaussian-smith-profiles.md)
recovers every Smith order of an invertible Gaussian operator on `M`
fermionic modes from a matrix of order `2M`. An elementary characteristic-two
proof includes singular creation factors and the scalar normalization.
For the reflection crossing this permits all independent physical mode
heights, removing the earlier internal regularity restriction. Run
`make verify-unrestricted-gaussian` for full original-tensor comparisons.

The [nonconstant row-frame constructor](research/nonconstant-binary-row-flags.md)
extends complete binary crossing synthesis to arbitrary regular polynomial
changes of row basis. Its first completed F4 catalogue covers every
strictly lower triangular first jet at both cuts and 59 explicit targets.
Only the reference affine line survives. The universal binary theorem above
now excludes a threshold crossing even with arbitrary higher jets. Run
`make verify-polynomial-flags` for literal polynomial and full-minor controls.

The continuation from PR #38 adds a [multilevel canonical-packet
collision](research/multilevel-thue-morse-packet.md): an explicit one-bit
Thue--Morse input loses its original bit simultaneously on every retained
valuation level below a chosen cutoff at the stated precision. This is a
scoped reconstruction obstruction, not a stronger factoring exponent.
Replay it with `make verify-multilevel-packet`.

A [read-once Boolean filter theorem](research/read-once-boolean-filter-floor.md)
excludes every homogeneous read-once XOR/AND occupation filter throughout
the interior charge interval, with no degree or mode-count bound. A single
parity of any number of modes is excluded at every real charge. The proof
uses a degree-controlled log-odds inequality and retains the full physical
cost. Run `make verify-read-once-filters` for exact independent controls.

The [all-digit charged Gaussian theorem](research/all-digit-charged-gaussian-floor.md)
now proves positive critical pressure for every real digit parameter
`q > 1`, integer interval length, real aggregate charge, and positive
scale. Differentiating in `q` reduces the result to a single convex kernel;
2,648 exact positive rectangles cover the remaining compact domains.
This closes the odd-digit gap for this Gaussian family while leaving the
proved factoring score unchanged. Run `make verify-all-digit-gaussian`.

The [global homogeneous-filter theorem](research/global-homogeneous-boolean-filter-floor.md)
now excludes every homogeneous one-level Boolean occupation filter at
every real aggregate charge and positive scale. It permits arbitrary
overlapping patterns and every degree through the selected-mode count.
The proof extends the earlier charge bands using a uniform edge estimate
and 13,323 exact positive cover leaves. General multilevel and resonant
nonhomogeneous filters remain open. Run `make verify-global-homogeneous-filters`.

The latest [research checkpoint](research/frontier-session-20261031.md)
keeps the proved score at `1/2`. An
[accelerated exact certificate](research/accelerated-interval-certificates.md)
now treats `2^1636-1` modes with a rational upper bound less than `2e-980`
above `1/sqrt(2)`. This remains an approach from above. A finite identity
with a positive remainder makes the check polynomial in the logarithm of
the mode count and can also certify long Gaussian spectator products.

Symbolic entry witnesses now support a
[uniform balanced-gauge obstruction](research/rectangular-balanced-gauge-obstruction.md)
for the baseline repair and its
[full-alphabet extension](research/full-alphabet-leading-obstruction.md).
The [all-scale extension](research/all-scale-balanced-obstruction.md)
now proves positive critical pressure for every `log s > 0` once the
interval is sufficiently large. Both the weights and the scale may
change arbitrarily with the interval length. The result concerns these
specific repairs and balanced bulk charges. Searches continue in families
that alter their preserved background, including
[opposite-charge cubic terms](research/rectangular-core-monomials.md).
Four [critical cubic pairs](research/full-alphabet-monomial-residuals.md)
also change the residual itself and restore all sixteen states at the
singular canonical quartic coefficient. Their determinant is exactly
`(1+gamma+ab)^4`, and a height-block correction proves the all-length
crossing. Its partial repairs still preserve an all-low Gaussian submatrix;
their occupation mixing requires a separate leading-pressure analysis.

An [arbitrary frozen-core output row](research/polynomial-frozen-core-crossings.md)
now gives an eleven-parameter affine family of locally polynomial crossings.
Its endpoint condition is exactly `g_0 + g_10 = 1`. The family contains both
the Gaussian background and a full-inverse restriction with an exact
non-Gaussian coefficient witness. This changes the background excluded by
the fixed-repair obstruction, while a useful pressure remains to be found.
The [compact frozen-row compiler](research/frozen-core-row-profiles.md)
reduces both cuts to matrices of order at most 64, independent of interval
length, and normalizes the nonzero output parameters to one orbit.
Its complete-polynomial cancellation search records 24,528 profile queries
through interval length 64. The Gaussian row was best in those finite
searches; arbitrary gauges and affine rows remain open.

A [strict Pascal repair](research/strict-gaussian-relabeling-crossings.md)
preserves the inherited Gaussian profiles and makes the local endpoint an
exact crossing on the full alphabet. It permits affine relabeling sums,
binary charge projections, and nonseparable entrywise monomial changes.
The charge projections split both cuts into explicit square blocks with
no degree increase. Full-letter gauge searches are recorded separately
from exact endpoint and profile controls; no threshold witness is claimed.

[Slow endpoint corrections](research/slow-ramification-crossings.md) after
ramification have exact eventual affine Smith orders, with an explicit
threshold. An exact qutrit example decreases the profile product by about
1.47% at one fixed scale. It does not improve the Gaussian optimum; its
finite correction and long-interval embeddings remain constructive search
targets. Both endpoint admissibility and the two complete profiles are
checked before interpreting a favorable local cost.

[Prime-state identity filters](research/prime-identity-filter-diagonalization.md)
now have complete all-size profiles from a finite-difference transform of
their digit weights. Explicit polynomial Pascal factors diagonalize every
charge block; the active zero matrix has only `p^2` states, independent of
interval length. The [qutrit specialization](research/truncated-symmetric-identity-filters.md)
also proves that every nonidentity monomial filter strictly increases the
profile product against its Gaussian background throughout the stated
charge range. General polynomial filters and coupled directions remain
constructive targets. `make verify-state-filters` checks the exact finite
identities; its extended target also replays the original characteristic-five
matrices of order 625.

The [polynomial filter synthesis algorithm](research/polynomial-state-filter-synthesis.md)
goes beyond monomial weights. Prescribed degree bounds, zero orders, and
endpoint values form an exact finite-field feasibility problem. Residue
compression uses at most `2*p+1` small systems and constructs explicit
polynomials even at mode scales with hundreds of bits. Full infinity
degrees and any denominator clearing remain charged. This supplies a
concrete route for testing coordinated cancellations at both cuts.

The [coupled identity-filter reduction](research/coupled-identity-filter-reduction.md)
now concentrates arbitrary joint filters on `r` selected directions into
`2*r` zero modes, independently of interval length. The proof covers
initial directions in every characteristic, and terminal or opposite-end
directions at prime-power lengths. Exact original-cut comparisons check
both profiles, including the auxiliary zero units in the mixed case.
These constructions permit correlations between modes; no rate gain is
claimed. `make verify-coupled-filters` replays the controls.

At a prime-power interval length, the
[simultaneous Pascal reduction](research/prime-power-joint-identity-filters.md)
now separates every identity direction in one regular change of basis.
Any selected subset, including interior and irregularly spaced indices,
has exactly twice as many active zero modes as selected directions.
Joint diagonal filters have an explicit finite-difference profile, and
the operator interfaces accept arbitrary selected indices. One exact
coefficient catalogue can therefore be reused over all selected subsets.
The recorded three-mode F4 and four-mode binary screens found no gain
on their stated finite grids. `make verify-prime-power-filters` replays
the exact identities, full-pressure controls, and original physical cuts;
its extended target also checks the original qutrit matrices of order 729.

The [joint residue synthesis theorem](research/joint-binary-filter-synthesis.md)
extends polynomial feasibility to correlated filters in every prime
characteristic. It uses at most `2*p^r+1` small residue patterns and a
fixed coefficient extension for nonzero endpoints. The
[completed catalogues and all-index screens](research/joint-filter-search-records.md)
retain 6,899 literal witnesses over F16, F32, and F27; every degree,
beta order, endpoint, and witness hash is replayed. The stated finite
grids found no improvement. A separate
[flag and zero-unit reduction](research/weyl-filter-flags-and-unit-reduction.md)
preserves exact constraints while reducing an interval-31 source from
degree 2,221 to 124. `make verify-joint-polynomial-synthesis` and
`make verify-weyl-filter-flags` reproduce these gates.
For a correlated [homogeneous Boolean filter family](research/homogeneous-boolean-filter-bounds.md),
the full cost is strictly worse than the same natural-charge Gaussian
at every scale. The [global Gaussian floor](research/global-charged-gaussian-floor.md)
now proves positive complete critical pressure for **every integer length,
every real aggregate charge, and every positive scale** in the binary
Gaussian family. It combines exact charge folding, analytic tails,
small finite bases, and integer-outward certificates. The analytic floor
also extends to every power-of-two digit count.

This makes the homogeneous filter's large-scale branch absolutely
positive at every real charge. For fixed Boolean degree and bounded
charge, its small-scale bound gives positive pressure at every scale for
all sufficiently large binary prime-power lengths. This last conclusion
does not extend the correlated small-scale result to unbounded charge.
`make verify-homogeneous-boolean-filters` replays its exact identities.

The [complete one-mode Laurent theorem](research/global-charged-gaussian-floor.md)
excludes every admissible single-mode binary Laurent occupation filter
at every aggregate charge, scale, and finite degree, over any fixed
characteristic-two coefficient field. The construction retains its
power-of-two interval lengths. Inside `|J| <= 2n`, equality with the
Gaussian occurs only for a common scalar Laurent monomial. Outside that
range, relative gains can occur, but an exact loss cap and a stronger
Gaussian floor keep the full pressure positive.
`make verify-global-one-mode-filters` replays the complete certificate
chain. Independent physical-state gauges and coupled multilevel or
off-diagonal operators remain outside these exclusions.

[Extension-field operator lattices](research/identity-operator-field-flags.md)
now permit constant F4 row bases and constructive F16 endpoint lifts.
The scalar Pascal contraction also gives small exact sources for every
selected subset. A separate
[one-hole minimum tradeoff](research/identity-operator-minimum-tradeoff.md)
proves that improving that block's infinity minimum by `2*g`, for `g > 0`,
forces a zero-order loss of at least `n*P+g` powers of the local variable.
It applies when the top coefficient is a monomial, over every fixed
coefficient extension. Deeper blocks and nonmonomial top coefficients
remain open. `make verify-operator-field-flags` checks both results.

The [four-mode resonance reduction](research/binary-four-mode-resonance.md)
does reach a deeper occupation block: two equal-energy two-particle
states admit an arbitrary invertible endpoint matrix. Both complete
profiles reduce to five scalar factors and one two-by-two matrix.
An exact polynomial-degree obstruction excludes one subfamily;
[complete Laurent catalogues](research/resonant-laurent-catalogues.md)
then cover 11,796,480 admissible kernels in each of two stated coefficient
windows. Their exact signatures were searched over all 2,360 resonant
quadruples at interval 32, with no numerical gain. The windows overlap,
and the pressure conclusions concern the recorded finite grids and
refinements. `make verify-resonant-operators` replays the exact controls.

For lengths just below a power of two, the
[near-power operator reduction](research/near-power-identity-operator-reduction.md)
handles terminal modes and a split with one initial mode. It retains
`2*r+1` actual zero modes and uses a local auxiliary Gaussian only to
compute the profile; exact multiplicity division removes that factor
before physical costs are evaluated. The contraction is a scalar matrix
plus rank one. Complete interval-31 terminal and split projector screens
found no gain in their recorded grids. `make verify-near-power-operators`
checks the interpolation, auxiliary-factor removal, and original cuts.
The [fixed-deficit extension](research/fixed-deficit-identity-cores.md)
proves a corresponding mixed-end core at `n=p^k-d` in every prime
characteristic. It retains `2*r+d` actual zero modes, calculates their
profile through order `p^(r+d)`, and bounds the departure from the
scalar prime-power contraction by rank `d`. Its general gate checks
exact polynomial geometry and Gaussian-factor division; further complete
physical operator replays are stated separately.
The [arbitrary-selection operator API](research/arbitrary-one-deficit-operator-search.md)
now handles every selected set one below a prime power, with original
physical controls and exact Gaussian-factor division. Four parent-length
qutrit catalogues and four
[energy-conjugated catalogues](research/deficit-energy-conjugated-operators.md)
completed at interval 26, including interior selections. Their 9,032
kernel/geometry cases retain the Gaussian best value in the recorded
finite grids and refinements.

The [prime-field affine solver](research/prime-deficit-affine-operator-synthesis.md)
then removes the fixed-projector restriction. It gives independent
coefficients to all 19 qutrit occupation-preserving entries, with exact
endpoint equations and a checked source-precision bound. Four searches
completed 242 prescribed paths and 4,985 exact affine queries, including
nontrivial infinity flags; no returned profile improved the recorded
pressure. These paths do not exhaust every polynomial coefficient family.
`make verify-prime-deficit-affine-synthesis` replays its independent gates.

The complementary [above-power construction](research/above-power-identity-operator-reduction.md)
handles terminal directions at `n=2^k+1` by retaining the existing,
unchanged identity direction zero. Its local core and Weyl matrix are
bounded independently of the interval, with no formal physical mode
added. The [interval-33 catalogue](research/above-power-operator-searches.md)
completed all 42,560 specified projector kernels and 1,111 exact profile
classes, with no gain on the recorded grid and refinements.
`make verify-above-power-operators` checks the complete construction.

The [reflection-pair compression](research/reflection-dual-rail-compression.md)
retains one occupied mode from every reflected pair. A characteristic-two
homology calculation gives an explicit polynomial Gaussian on the
square root of the original alphabet, with complete profiles at both
cuts. Its natural pressure equals half the original pressure at twice
the scale, so the rate is unchanged. The smaller seed makes further
logical-core perturbations and independent-gauge searches accessible.

The [physical hot-state search checkpoint](research/reflection-hot-state-searches.md)
records completed positive and unrestricted height screens, all 65,536
coefficient rows at five fixed positive-height rays, and nonlinear
two-bit state penalties. None improves its Gaussian control in the
recorded scope. The [raw transition-mask screen](research/raw-core-transition-masks.md)
also retains its interrupted 167-of-615-candidate prefix explicitly;
the remaining masks are not counted as searched.

[Cold-compensated hot-state searches](research/reflection-hot-compensated-searches.md)
allow negative selected-state heights while charging every compensating
cold-bit height. The saved runs include all 65,536 coefficient rows at
seven fixed negative-height rays, 234 natural-class representatives at
interval 16, sampled rows through interval 64, and nonlinear two-bit
OR-state penalties. None supplies a gain in its stated scope. The
[separate physical-leg bulk theorem](research/reflection-hot-leg-bulk-gauges.md)
also permits different bulk heights on the four legs within an explicit
regularity cone, retaining the same bounded active core and both physical
type costs. `make verify-reflection-leg-bulks` checks its original tensors.
The [second regular cross pivot](research/reflection-hot-all-bulk-gauges.md)
removes the bulk-ordering restriction: all nonnegative physical-leg bulk
choices now use the same bounded core, subject to the stated hot-cold
regularity inequalities. Complete coefficient expansion is checked before
Smith elimination, including the large local-unit denominators of this
second pivot. `make verify-reflection-all-bulks` replays the new construction
and the shared hot-state regression.

An independent [binary crossing synthesis](research/binary-two-cut-affine-synthesis.md)
assigns free Laurent coefficients to every tensor entry while fixing the
strict permutation endpoint. For its selected row targets, a support
bound makes degree 12 after common clearing exhaustive, and all 8,925
constant F4 flags at each cut are covered. All fifty nonreference targets
are infeasible; this is an exact result for those fixed frames and gauges.
It does not exclude other target orders or other regular frames.
`make verify-binary-two-cut-lattices` checks literal minors, endpoint
coefficients, flag enumeration, and the pruning used by the complete run.

The [Walsh permanent calculation](research/walsh-principal-permanents.md)
gives all-size scalar valuation formulas beyond the determinant model.
It also proves that exponentially many grouped modes survive the current
precision schedule for direct forward twirling. Lower-precision scalar
recovery remains a distinct possible route.

The [arithmetic Cauchy theorem](research/arithmetic-cauchy-factorial-law.md)
now proves the permanent valuation `v_2(per C)=v_2(n!)` for every
even-spacing arithmetic block, with no power-of-two size restriction.
It also proves every invariant-factor order of the squared Cauchy matrix:
`2*r*v_p(d) + 2*v_p(r!) + v_p(r+1)` when `p` divides the spacing `d`.
An integral orthogonal-polynomial basis supplies the full Smith proof.
The [all-power extension](research/arithmetic-cauchy-powers.md) now gives
`2*r*v_p(d) + v_p(r!) + v_p((m)_r)` for every positive Hadamard power `m`,
where `(m)_r` is the rising factorial. A factorial-normalized alternant
bound makes the leading Cauchy--Binet term unique modulo `p`. These are
complete spectral formulas; scalar recovery and grouped operator cost
remain separate obligations.

The [digit-count interval theorem](research/truncated-symmetric-interval-frontier.md)
also proves that the displayed truncated-symmetric interval pressure stays
strictly above the threshold for every integer digit count `q >= 2`,
every mode count, and every profile scale. Its fixed-`q` gap has a sharp
Lambert-W asymptotic. Noninterval modes and independent mode heights remain
outside this theorem. The [weighted Pascal reduction](research/weighted-lucas-pascal-profiles.md)
now computes admissible independent-height profiles by constant-field elimination.
Lucas antichains permit every height vector, and a prime-radix traversal skips
arbitrarily large gaps between modes. These are exact search tools; the original
crossing degree remains part of the eventual evaluator cost. Run
`make verify-frontier-extensions verify-weighted-pascal` for the ordinary-proof
controls; neither result changes the factoring score.

The [mixed-field rectangular repair](research/rectangular-f4-repair.md) now
settles the four-element-field construction left open in the preceding
checkpoint. Its explicit correction preserves the low Gaussian background
and has at most a simple local pole for every finite mode set. Matching
paired types give polynomial crossings. The optional rectangular profile
engine now accepts F4 coefficients and independent physical types; its
controls replay exact rational-weight profiles. These are construction and
search advances, with no score improvement claimed.

For full intervals, the [bounded-defect reduction](research/rectangular-defect-factorization.md)
now confines this repair to seven active modes at infinity and nine at zero,
at the explicitly stated Gaussian gauges. The remaining modes contribute
scalar factors. Exact rational support checks and an independent full-matrix
control support a scalable profile calculation. The
[paired-type search audit](research/rectangular-convex-pressure.md) records
the derivative, minor-separation formulation, and precision conditions.
Regular core coordinate bases now permit fourteen independent core weights
and four bulk charges without increasing the active dimensions. The
[secondary-height family](research/rectangular-refined-repairs.md) adds
input and output repairs with a shared additive flag, together with
arbitrary endpoint-compatible Laurent quartic numerators. A
[full-core extension](research/rectangular-full-core-factorization.md)
uses ten active modes at both places and supports all four extreme
physical cores. Its 25-parameter pressure calculation has dimensions
independent of the interval length. Adaptive exact series windows permit
larger intervals while rejecting unresolved pivots.

The [full-alphabet quartic repair](research/full-alphabet-quartic-repair.md)
retains all four extreme states on both crossing pairs. Its residual is
invertible for every quartic coefficient other than one, and a sparse
partial correction gives non-Gaussian examples over F4. Fixed-size profile
matrices support 23 or 29 continuous parameters for this family.

[Coupled core transports](research/rectangular-core-transports.md) now
permit arbitrary invertible core changes, including height-raising
directions, when an exact left/right product preserves the endpoint.
Their matrices can be updated from cached coordinates. The
[critical-pressure analysis](research/rectangular-critical-pressure.md)
separates genuine defect gains from interval convergence, while a
[bivariate limit](research/bivariate-core-pressure-limit.md) reduces the
large-interval question to a finite leading matrix and its weighted ranks.
The [leading-face certificates](research/rectangular-leading-face-lp.md)
give exact rational LP bounds for the fitted rectangular model. Its
[all-low background](research/rectangular-low-background.md) is exactly a
shorter Gaussian interval. Separate [polynomial shears](research/rectangular-polynomial-shears.md)
and [four extra quartics](research/rectangular-extra-quartics.md) preserve
the endpoint and provide additional admissible search directions.

The [all-field Pascal construction](research/all-field-pascal-crossings.md)
retains the required global signs in odd characteristic. The
[graded-power theorem](research/graded-even-power-transfer.md) also
permits even perturbations that vanish at the crossing endpoint. The
scalar factoring dependency is pinned in a separate
[source audit](research/fixed-prime-transfer-source-audit.md).
A [truncated symmetric construction](research/truncated-symmetric-pascal-crossings.md)
also gives ordinary tensor powers with p states per mode in characteristic
p, including non-Gaussian perturbations of odd degree.

For the original evaluator's precision question,
[canonical packet bounds](research/canonical-packet-upper-frontier.md)
and [finite decoder certificates](research/canonical-packet-decoder-frontier.md)
give exact four-wire binary thresholds: 12 local bits for parity and
13 for the full count. They do not supply an efficient all-size decoder.
The [full Boolean operator](research/boolean-cauchy-smith-profile.md)
has an exact all-size Smith profile and a one-bit inverse at either
nearest center. Direct insertion of that inverse reconnects the source
groups, so this operator theorem alone does not improve the evaluator.
The [grouped inverse profile](research/grouped-cauchy-inverse-profile.md)
determines the exact precision cost of a direct grouped factorization.

The [earlier 31 October research checkpoint](research/constructive-frontier-checkpoint.md)
collects the latest exact bounds, new construction families, reproducible
search starts, and the remaining obligations for a score improvement.

No uniform factoring exponent above `1/2` is proved here. One open approach seeks effective crossing pieces with total rank width `4*rho*lambda + O(log lambda)` for `rho < 1/sqrt(2)`. A suitable family could improve the evaluator bound, but finite rank improvements alone do not establish that inequality.

The [complete Schur entry algorithm](research/schur-fermion-queries.md) computes
individual modes with determinants of order at most `3 sqrt(H)`, and full
crossing coefficients with determinants of order at most `4 sqrt(h)`. It works
directly in finite characteristic and removes factorial-precision conversion
from these queries. Its ordinary proof and independent exact checks improve
the available construction tools. For the homogeneous Schur family, [hook-residue identities](research/schur-hook-residue-central.md)
give central rank at most `p(2h)-3` at every power-of-two weight. Full exact
calculations attain this bound at weights 8 and 16, with ranks `228` and
`8,346`. These central pieces exclude the corresponding raw finite seeds;
no stronger factoring bound is claimed.

The [complete inhomogeneous Schur kernel](research/schur-gaussian-digit-reduction.md)
has an exact binary-digit reduction to a four-label Gaussian crossing.
The proof retains its exterior maps and energy projection and gives an
explicit finite coefficient decomposition. Fixed-occupation sectors of
the seed are complete crossings in their own right.
A [sparse-index calculation](research/sparse-holonomy-singleton-profile.md)
gives every Smith order for the zero-and-singleton index set and an exact
Gaussian encoding for all Boolean down-sets. It supplies formulas for
larger constructions; it does not assert a better factoring exponent.
[Boolean down-sets](research/schur-downset-crossing.md) give complete
Schur crossings with matrix-valued holonomy and controlled profile limits.
[Offset intervals](research/schur-offset-interval-reduction.md) reduce to
an explicit two-scale four-label kernel, with polynomial projection and
coefficient costs when the ambient offset grows polynomially.
[Two unequal intervals](research/schur-two-block-holonomy.md) split exactly
into scalar sectors; the resulting direct kernels also allow interval
parameters that are not powers of two.
[Fixed mode counts](research/fixed-mode-crossings.md) extend the direct
crossing to arbitrary index sets, unequal species occupations, and
independent mode heights. Their product-type criterion is exact; an
[exploratory search driver](tools/search_schur_mode_pressure.py) uses the
weighted Smith oracle to search this larger family.
An [exact three-mode certificate](research/three-mode-pressure-improvement.md)
uses independent occupations of the two species to prove `rho < 741/1000`.
It improves the previous Gaussian pressure control but remains above
`1/sqrt(2)`, so the best factoring exponent is unchanged.
The [one-species interval family](research/schur-one-species-interval.md)
has an exact all-size spectrum and finite product pressure. Checked
rational witnesses reach `rho < 3363/4756` at 16,383 modes. This particular
formula approaches `1/sqrt(2)` from above; an analytic proof excludes
crossing that threshold within the displayed gauges.
The [quantitative interval theorem](research/quantitative-interval-frontier.md)
gives the sharp gap asymptotic `(log n)^2 / (sqrt(2) pi^2 n^2)` and
exact large witnesses using outward integer intervals.
The [accelerated extension](research/accelerated-interval-certificates.md)
replaces mode enumeration with exact geometric moments and reaches the
much larger finite certificate described above.
The [greedy Pascal-complement theorem](research/schur-pascal-complement.md)
computes the same normalized spectrum for every finite mode set using
binary linear algebra. It proves the interval formula at every length
and skips large gaps by a traversal of binary intervals.
[Constant matrix flags](research/constant-flag-crossings.md) extend the
polynomial construction beyond Pascal matrices while preserving its
determinant and typed endpoint. An exact four-mode screen over the
four-element field covers 4,096 flags; all remain at or above the Pascal
pressure at the displayed gauge.
[Positive-energy perturbations](research/positive-energy-quartic-endpoints.md)
also preserve the typed crossing when they do not vanish locally at the
endpoint. A strict index-sum condition even permits rational local kernels
whose endpoint denominators cancel in every fixed-type tensor power.
For [localized quartic changes](research/localized-quartic-profiles.md),
exact determinant and complete zero-profile formulas quantify the cost
of the resulting spread in Smith orders. A rank-deficient limiting
endpoint is explicitly checked and requires a further restriction.
[A polynomial physical code](research/residual-quartic-codes.md) repairs
that endpoint and cancels its denominators term by term. It gives an
explicit eight-letter, degree-32 candidate whose quartic correction
survives common types. Its endpoint is proved; a pressure improvement
has not been established.
[Arbitrary-mode codes](research/coded-quartic-mode-families.md) extend this
construction to noninterval sets and independent Laurent coefficients on
the four physical maps, with a reusable exact constructor.
[Conjugate phase codes](research/conjugate-residual-phase-codes.md) pair
two four-element-field tensors to cancel their endpoint phases. Three
checked eight-letter, degree-48 crossings descend to binary coefficients
through an explicit constant change of basis.
[Fixed field descent](research/fixed-field-descent.md) also transfers any
coefficient certificate over a fixed extension to its base field without
enlarging the physical alphabet. Projecting after each complete tensor
power costs only the extension degree in the coefficient-cost sum.
An [all-field residual obstruction](research/residual-quartic-three-letter-obstruction.md)
excludes three-letter restrictions by four constant physical maps over
every characteristic-two field. When both crossing pairs use the same
code, the resulting three-letter operator is always singular. Larger
blocks and bulk-dependent codes remain outside this obstruction.
A [signed exterior model](research/signed-residual-codes.md) extends the
single-copy obstruction to every field. Exact identities over the
27-element field also give a three-conjugate residual phase repair;
its full polynomial factoring construction remains an open obligation.
[An asymmetric three-letter core](research/involutive-core-crossings.md)
gives polynomial candidates on `3 * 2^(n-2)` letters for arbitrary finite
mode sets. A cancellation identity removes the local denominators, and
the remaining endpoint error vanishes on every common word type.
Exact six-letter and twelve-letter examples, monomial code parameters,
and sparse output corrections are supplied for pressure searches.
[Coordinate core restrictions](research/coordinate-core-crossings.md)
retain a filled low-mode state and give a second arbitrary-mode family
on the same alphabet sizes. Its local entries can have genuine poles;
matching word types on each crossing pair cancel them globally and
justify the restored rational Smith profiles. The two pairs may have
different types, and pairing with an exchanged copy restores equal
external dimensions.
[Independent paired-type pressure](research/paired-type-crossing-pressure.md)
gives the exact entropy optimization for different types and rectangular
alphabets. Tensoring with the pair-exchanged copy produces an equal-alphabet
crossing; keeping the two coefficient indices separate makes the costs
multiply exactly.
[Rectangular residual repair](research/rectangular-residual-repair.md)
retains four extreme states on one crossing pair and three on the other.
An explicit two-parameter inverse and a smaller height-triangular repair
give valid paired-type families. Finite mixed examples are locally
polynomial, and the smallest partial repair has an exact non-Gaussian
Wick witness. Their pressure remains a search problem.
[Even-core postselection](research/even-core-postselection.md) gives a
four-letter polynomial family `A + lambda B + lambda^2 C` with a strict
crossing endpoint. Two exact Wick identities prove that every nonzero
rational parameter path is non-Gaussian. Monomial parameter paths are
available for pressure searches; no improved rate has been certified.
The [two-output extension](research/two-parameter-even-core.md) classifies
its Gaussian members exactly: the coordinate axes and one exceptional
parameter point. An explicit eight-edge factorization identifies that
point, so numerical improvements near it can be interpreted correctly.
The [four-boundary extension](research/four-parameter-even-core.md)
classifies the larger Gaussian set as well. Every Gaussian member with
all four parameters nonzero is a scalar multiple of that same exceptional
tensor; non-Gaussian parameter searches can exclude this repeated family.
The [optional native F4 search engine](research/native-weighted-smith-f4.md)
provides exact field arithmetic, support-block decomposition, weighted
profiles, and pressure gradients. Its finite controls use an independent
rational oracle. Floating search results still require exact certificates.
An [optional rectangular engine](research/native-rectangular-smith-search.md)
handles unequal alphabets and independent paired types with the same
explicit denominator bookkeeping and exact rational controls.
For this Gaussian seed, [generic mode blocks](research/typed-gaussian-generic-blocks.md)
force cost `Omega(D^2)` when all nonvacuum occupations stay fixed, even
after arbitrary coefficient recombinations. Growing occupation densities
remain a separate regime.
[Charge-block compression](research/typed-gaussian-charge-blocks.md)
computes the actual coefficient ranks and recombinations at larger types.
For 1,260 labels it replaces the full 1,587,600-dimensional flattening by
38 orbit blocks, each of size at most 1,260. The checked decompositions
remain above the improvement threshold.

[Exact type optimization](research/type-optimized-crossing-pressure.md)
reduces the best common letter frequencies to one rational pressure inequality.
It applies to Smith and module certificates. A checked physical gauge gives
`rho < 61/80` for the four-label Gaussian seed; the three-mode certificate
above improves this local pressure bound. Both remain weaker than the
known half-rate construction.
[Polynomial interpolation](research/polynomial-crossing-interpolation.md)
constructs crossings from prescribed local row lattices by exact linear
algebra. The complete fixed-lattice fiber consists of 64 local changes, and a
relaxed system constructs a non-Gaussian control.
An [explicit tensorized family](research/tensorized-smith-sections.md)
adds 40 new directions to a 180-dimensional inherited space. A common
nilpotent residue certificate proves that all `2^220` crossings preserve
the two specified profiles; changes of gauges remain a separate question.
An [exact weighted Smith oracle](research/weighted-smith-separation-oracle.md)
finds all determinantal minima without listing every minor. It supports
continuous gauge searches with exact rational replay.
[One-sided charge changes](research/one-sided-charge-pruning.md) disappear
exactly under common-type restriction and cannot improve either convex
Smith-profile sum. Neutral interactions and opposite-charge combinations
remain available to the constructive search.

[Finite module automata](research/cancellation-module-smith-automata.md)
retain cancellation between auxiliary paths. Their bounds improve monotonically
with precision, dominate the corresponding capped entry-order bounds, and
converge to the exact full-matrix word rate when both boundaries are free.
[Reachable and observable reduction](research/minimal-boundary-module-completeness.md)
extends this completeness result to fixed scalar boundaries, while keeping
physical choices grouped by channel. This is a statement about word bounds,
not a characterization of actual coefficient ranks.
[Exact rational lattice orbits](research/exact-rational-module-orbits.md)
can retain identically missing directions and compute the full word rate
when the orbit closes. A companion example has explicit two- and seven-state
orbits and a proved certificate floor.
An [explicit coupled crossing](research/cancellation-module-crossing-example.md)
strictly improves its corresponding capped min-plus certificate, from at least
`425/27` to `1625/108`. The same local factors and precision caps are used.
The checked controls do not yet meet the strict improvement threshold.

[Fixed block Smith profiles](research/block-smith-crossings.md) form a
submultiplicative hierarchy. A three-position ordered-pairing block improves
its normalized certificate from `16` to about `15.7138`; this remains weaker
than the elementary scalar control and does not change the factoring bound.

The [degree-one Smith obstruction](research/linear-smith-barrier.md) proves
`B >= q^3` for every alphabet, field, and auxiliary dimension in the bare
linear certificate. The sharp binary bound is `256/27`; the scalar bound
in odd characteristic is `12`. All 65,536 binary scalar pencils are catalogued.
Higher degree, blocks, and cancellation refinements remain separate routes.

For the nilpotent-path construction, a [controlled monomial obstruction](research/controlled-permutation-smith-floor.md)
extends the existing permutation calculation to dependent and noncommuting
controls. It gives `B >= q^3/sqrt(2)` in its stated scope. A separate
[qutrit Weyl control](research/qutrit-weyl-crossing.md) has an exact all-length
constant-coefficient rank formula; its raw pieces also fail to improve the bound.

The [varying-image cubic construction](research/cubic-full-image-span.md)
has an exact 786-dimensional rank-one-map coefficient space, with column ranks
reduced from 49-by-49 matrices to matrices of size at most 28. Rank-two maps
already span its full 1,225-dimensional major-symmetric space. These results
define larger exact search spaces; an improving crossing decomposition is
still needed.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work can also appear in [open pull requests](https://github.com/1337isnot1337/integer-factoring/pulls); branch research is not automatically a result on `main`.

See the [research index](research/INDEX.md) for completed arguments, scoped obstructions, and older candidates. Active work is committed directly to `main` by default; useful partial research can be integrated as a checkpoint without being promoted to a stronger theorem.

An [arithmetic framework program](research/arithmetic-framework.md) studies
deterministic class-group relation generation outside the crossing evaluator.
It proves a positive Fourier hitting lemma and deterministic extraction of
two-torsion from unsaturated relation data. A powering seed with published
expander amplification supplies biased words from generators without a group
decomposition. Full prime-form generation now has an [explicit polynomial cutoff](research/explicit-prime-form-cutoffs.md)
under the pinned finite-order Hecke zero-free theorem. Fast norm selection
remains open; the [completed processing stack](research/beyond-half-research-checkpoint.md)
states its precise contract. The score remains `1/2`. Reproduce the earlier
framework controls with `make verify-arithmetic-framework`.
The [actual quadratic-form controls](research/arithmetic-class-group-controls.md)
also recover proper factors for 18 finite targets and check complete
two-torsion recovery; their generator and order searches enumerate classes.
A [density/index transfer](research/arithmetic-density-lift.md) also turns
relations on a possibly proper small-prime subgroup into a global
class-number multiple. Fast pair extraction and deterministic coverage of
the missing ambiguous directions remain unresolved.
The [bounded smooth-norm interface](research/bounded-smooth-norm.md)
characterizes the required pairs and converts supplied representations at
`2^O(b/log b)` divisor-search cost; locating those representations efficiently
is still open.
A [deterministic norm-sieve experiment](research/arithmetic-norm-sieve.md)
constructs dominant relation certificates for 60-bit and 80-bit inputs and
recovers proper factors without class enumeration. Replay both certificates
and two smaller-box exhaustion controls with `make verify-arithmetic-sieve`;
success in a larger finite box is not a uniform factoring-time bound.
The direct-word certificates use the polynomial ideal-lattice composer and
a diagonal fixed by the word bound, reducing the two determinant encodings
to 520 and 500 bits. The original pair-form certificates remain replayable.
The [quantitative localization target](research/smooth-norm-localization.md)
now specifies the required short-box theorem and its conditional complete
cost. The localized smoothness statement and splitting coverage are both
unproved; the existing class-group density theorem supplies neither for free.

The [known-root conic parametrization](research/class-square-conic.md)
constructs all rational square-norm points with an integral matrix of
determinant `|Delta|`. Its exact cancellation bound gives direct gcd splitting
certificates and rewrites prescribed-square search as a smooth denominator
problem. Finding that denominator remains open. Run
`make verify-class-square-conic` for independent conic-point and large-certificate
controls. Three explicit rational points also give a polynomial algorithm for
solving every modular square-norm instance without a prime factorization of the modulus.
An exact negative integer example is soluble at every modulus, separating
this local information from the required global search. The score remains `1/2`.

An [exact counting alternative](research/divisor-moment-recovery.md)
recovers the smallest prime factor from a family of divisor moments by
positive Hankel reconstruction and integer root isolation. After small-factor
trial division, recovery costs `exp(O(sqrt(b)))`. The moments are coefficients
of level-one Eisenstein series with explicitly bounded growing weights.
A suitably uniform fast coefficient engine would preserve its saving exponent
in the complete factoring bound. At the trial budget for `1/2 < beta <= 1`,
the required weights, query count and recovery work are polynomial in `b`.
The coefficient engine remains unproved. Run
`make verify-divisor-moments` for the supplied-moment controls.

The [global sparse-residue transfer](research/sparse-divisor-residues.md)
connects that trial budget to the fixed-characteristic weighted constructor.
Garg--Schost interpolation recovers the integer divisor polynomial from
global cyclic residues, with polynomial algebra degrees through target
score one. The refined transfer preserves every fixed positive saving
exponent, removing the earlier cap from the transfer overhead. Its faster
weighted evaluator remains unproved. Run `make verify-sparse-divisor-residues`.

A [fixed weight-one alternative](research/weight-one-divisor-coefficients.md)
uses odd-character Eisenstein coefficients at small prime levels. Exact
Fourier inversion gives signed divisor residues, and the same sparse
reconstruction recovers factors. One bit per cyclotomic coordinate suffices
at every level `q`: a binary decoder keeps the full algebra modulo two,
including its nilpotents. The coefficient at the original composite
index remains the unresolved engine; no improvement follows from prime-index
formulas. An [explicit polynomial decoder](research/cyclotomic-character-decoder.md)
inverts the same one-bit map by Galois substitutions, odd-order Fourier
inversion and exact nilpotent divisions, avoiding a stored Gaussian system.
Run `make verify-odd-character-residues` for both decoders and the component
checks, or `make verify-cyclotomic-decoder` for the new construction alone.

The [theta-section engine](research/lambert-theta-sections.md) computes the
signed divisor coefficients without enumerating divisors, but its current
state cost is exponential in the input length. A uniform rank argument
excludes polynomial-size common linear section spaces even at phase order
three, where a direct square test is cheap. This is a restriction on that
representation, not a factoring lower bound. Run `make verify-lambert-theta`.

The [Tate-coordinate compiler](research/tate-coordinate-coefficients.md)
expresses odd-index signed divisor coefficients as coefficients of
`X_q(Delta(t))`. The inner algebraic coefficients are cheap, and one scalar
trace of the original composition has a square-test formula. Extracting the
remaining trace-zero coordinates at the original index is still open.
Run `make verify-tate-coordinate`.
At [phase five](research/q5-quadratic-coefficients.md), the remaining bit
reduces to four binary quadratic-form parity counts. The explicit evaluator
still costs `O(sqrt(N) poly(log N))`; it supplies no faster unrestricted bound.
Run `make verify-q5-quadratic`.
At [phase seven](research/q7-cubic-parameter.md), an exact cubic parameter
identifies a seven-regular partition series. A degree-four theta identity
computes the odd-index value; its modular/Sturm proof is covered by 65,536
exact controls. The explicit evaluator costs `O(N poly(log N))`.
Run `make verify-q7-parameter`.
A [partial theta-program compiler](research/theta-coordinate-programs.md)
certifies the odd-index identity at six listed phases using a uniform modular
congruence bound, and gives exact exclusions for its module at three others.
Neither uniform success nor a sufficiently fast evaluator is proved.
Run `make verify-theta-programs`.

The [Bell-diagonal seed construction](research/bell-diagonal-crossings.md)
gives exact two-cut transforms and excludes its complete symmetric
coefficient spans at local dimensions two and four. Dimension eight
remains open under that argument. Replay its physical tensors and pencil
ranks with `make verify-bell`; no stronger factoring exponent is claimed.
The [linear-pencil search](research/bell-linear-search.md) includes two
complete small-domain calibrations and explicit larger finite searches,
with no improvement in the recorded schedules.
The [polynomial synthesis engine](research/bell-polynomial-synthesis.md)
certifies uniform rational rank bounds and led to an explicit
[folded-cube crossing family](research/folded-bell-crossings.md).
Its optimal leading contraction within its complete coefficient span
is `H_2(1/4)=0.811278...`, above the critical threshold. The factoring
score is unchanged.
The [code-quotient analysis](research/code-quotient-crossings.md) preserves
the crossing duality, pays the maximal common polynomial division, and
identifies nonuniform finite-field cancellations as a remaining route;
ordinary covering-radius gains alone do not improve its pressure criterion.
The [signed-boundary extension](research/signed-boundary-crossings.md)
preserves duality under more general intertwiners and isotropic reductions,
with exact common-factor costs; its strict rank improvement remains open.
An [endpoint pressure theorem](research/bell-cubic-outer-rank-floor.md)
excludes the stated outer-rank-one cubics in every odd characteristic,
and quartics at local dimensions 8 and 16. Its first open power-of-two
quartic alphabet at least eight is 32.
For [self-reciprocal quartics](research/bell-quartic-rigidity.md), an exact
matrix identity relates a pressure witness to sparse Hadamard rank reduction;
it is a necessary condition and supplies no improved factoring bound.


The current ordinary-proof campaign establishes these scoped results:

| Scope | Result and evidence |
| --- | --- |
| Bare ungauged local Smith profiles, every field and degree, fixed auxiliaries | [Degree budgets](research/all-degree-small-alphabet-smith-floor.md) and [shared minors](research/shared-minor-smith-floor.md) exclude a half-rate improvement through six letters. |
| The same certificate on 7 through 12 letters | [Box-minor Kronecker budgets](research/bare-seven-to-twelve-letters.md) exclude a half-rate improvement through twelve letters; at 13 letters an exact profile satisfies every box budget and adjugate inequality proved there, so the criterion remains open. |
| Scalar paths with a monomial column determinant | [Inverse-degree compatibility](research/inverse-degree-crossing-compatibility.md) rules out the formal seven-letter optimum at every integer ramification. |
| Polynomial inverse paths on 7 through 10 letters | [Inverse-degree aspect bounds](research/inverse-degree-aspect-bounds.md) cap the useful degree; they also exclude all `2^31` binary level-weight assignments of the fixed nine-letter flag. |
| Noncommuting quadratic-inverse constructions | A [general two-jet tradeoff and complete nine-letter schedule](research/nine-letter-quadratic-inverse-probe.md) show why suppressing order one alone fails; all 31 tested level choices are excluded at every scale. Other pencils remain open. |

The [binary quadratic](research/binary-quadratic-catalogue.md) and
[extended polynomial catalogues](research/extended-binary-polynomial-catalogues.md)
retain the finite evidence behind these arguments, including all 16,777,216
parity quartics. Their independent replay limits are stated explicitly.
Run `make verify-polynomial-crossings` for the focused Python checks;
optional `--full` replays are documented in those notes. The separate
`make verify-polynomial-native` gate uses Clang and Boost for complete local
Smith profiles and independent Toeplitz controls.

Run `make verify-nine-letter-all` for every complete profile and all-scale
certificate in the 31-path native schedule.

These results restrict the stated certificates. Independent gauges,
coefficient recombinations, cancellation refinements, and larger alphabets
remain constructive routes; no stronger factoring score follows here.

The [rough-length Fourier constructor](research/rough-length-fourier-constructor.md) removes small prime factors and branches on the exact two operand lengths. Its carry-column width satisfies `P^2 lambda <= (12+o(1))n`, improving the earlier constant 24, with the full preprocessing and factor recovery costs included. The complete factoring score remains `1/2`.

The [moment-quotient localization analysis](research/pole-cauchy-localization-obstruction.md) proves why directly inverting the factor-count pole in the small polynomial quotient fails at polynomial precision. An exact proper-factor example also defeats common-factor cancellation followed by canonical rational coefficient lifting. The result scopes these reconstruction strategies; a different scalar extraction method remains open.

The [height-25 protected-state theorem](research/coherent-protected-height25.md) gives the exact natural-pressure floor `17/24 > 1/sqrt(2)` for every constant height-homogeneous split that retains the 28 states below height 16 and uses only states through height 25. It covers separate paired alphabets and types over every characteristic-two field, including the 96-letter hyperoval family. Changing the protected lower states or the stated gauge assumptions remains outside the theorem.

The [full-tail active-core reduction](research/coherent-full-tail-core.md) replaces an exponentially large compressed compound calculation by a fixed `4k`-mode core and explicit Gaussian spectators. Under the stated height and deletion hypotheses it computes the actual reflection profiles. The complete forbidden-pair example is checked through 512 physical modes and gives no gain over the unmodified Gaussian; all physical entropy and degree costs remain charged.

The [fixed rational-frame transfer](research/fixed-rational-frame-transfer.md) proves polynomial coefficient-splitting overhead for a fixed rational local factorization that is regular at the selected cut variable. It permits generic multivariate Smith calculations beyond finite Laurent frames, with uniform coefficient algorithms and separate expansion cones for different cuts. It supplies an additional construction tool, not a stronger local spectrum or factoring score.

The [physical full-tail reduction](research/coherent-physical-tail-low-band.md) retains the internal Gaussian factors for arbitrary constant low-mode physical maps and bounded label gauges, reducing a long tail to an explicit `8k`-mode core. The [all-field injectivity theorem](research/coherent-tail-injectivity.md) now proves that every nonzero fixed physical restriction retains all `2^(2k)` leading channels, over every characteristic-two field. A proper fixed code therefore approaches the critical pressure from above in the stated asymptotic window. Growing physical support, finite-tail scales, and broader constructions remain open.

The [encoded-cube specialization control](research/cube-pair-encoding-specialization.md) gives a concrete tensor whose diagonal specializations appear to beat the cube target, while independent coordinates expose lower generic orders and reject it. Its exact integer pressure certificate and direct sixteen-mode Pfaffian replay provide a search control against this false positive; every nonzero diagonal specialization remains misleading even over larger fields.

The [complete Gaussian coefficient fibers](research/cube-complete-coefficient-fibers.md) parameterize all 12,672 even homogeneous binary coefficients with the required endpoint. A loss of at most three in every fixed-frame row floor leaves only the seed, over every characteristic-two field. The complete 22-dimensional binary fiber at loss four has 4,194,304 checked points, all dominated by the seed at every scale. The accompanying engine returns exact generic polynomial row certificates for further weighted searches.

The [routed rectangular Fourier constructor](research/routed-rectangular-fourier-constructor.md) extends the improved constant twelve to every fixed grouping range `P<=n^delta` with `delta<1/2`. Its exact width is `lambda=6Tf/P^2+4(r+2)`, giving `P^2 lambda<=(12+4P^2/T+o(1))n` and `O(n/P)` table-query space throughout `P<=sqrt(T)`. Balanced equality relays preserve the factor-pair bijection and scalar phase recovery. This is a stronger constructor bound; the complete factoring score remains `1/2`.

The [flat-tensor parity theorem](research/cube-flat-parity-floor.md) proves the exact pressure floor `rho>289/500>1/sqrt(3)` for every homogeneous flat binary cube tensor whose odd-occupation coefficients vanish, over any field. It allows arbitrary even coefficients and an arbitrary vacuum polynomial. Shared parity minors and two integer entropy inequalities establish the uniform bound; parity-breaking coefficients remain outside it.

The [moving-pair construction](research/coherent-moving-pair-core.md) constrains physical modes `0` and `M-1` while keeping the middle interval free. A uniform eight-mode core and explicit spectators give complete natural profiles beyond the fixed-support hypotheses. The exact core profiles are replayed through `M=512`; their numerical pressures remain slightly worse than the corresponding full Gaussian. Parameter-dependent maps and larger moving supports remain constructive directions.

The [generic F4 frame engine](research/cube-f4-generic-frames.md) computes complete sixteen-by-sixteen flat-cube profiles with arbitrary physical parity. It returns exact polynomial row certificates and verifies their two unit determinants over `F4(y)`, so no finite-field point substitution is used. This supports nonbinary coefficient searches and rational-frame coefficient transfer; it does not supply a stronger factoring score.

The [low-order flat-cube filters](research/cube-flat-low-order-floors.md) apply over every field and allow odd physical coefficients. After common coordinate factors are removed, either a face of rank at least two or an average second Smith order at most one third of the remaining degree forces `rho>29/50>1/sqrt(3)`. Exact integer entropy certificates prove both tests. Deterministic even and odd reference constructors retain all 512 endpoint and 768 rank-one face identities. Run `make verify-cube-flat-low-order`.

The [complete flat-cube coefficient engine](research/cube-flat-coefficient-fibers.md) covers all 6,912 degree-six endpoint coordinates with no parity restriction, or 13,824 binary coordinates over F4. Exact generic row frames define complete affine fibers; independent polynomial products and full Smith certificates validate retained assignments. The included bounded search exposes source hashes, complete references, second-order feasibility tests, and replayable improvements. These tools and fixed-frame rigidity results leave the factoring score unchanged.

The [constant local parity-gauge theorem](research/local-parity-gauge-test.md) recognizes when a polynomial binary tensor can be made even by constant local basis changes over a coefficient-field extension. It uses seventeen linear unknowns for a general eight-leg tensor and eight with the equality endpoint, constructs paired endpoint-preserving gauges, and returns exact contradiction certificates. Thus the existing even-tensor pressure floor also excludes every positive instance. A retained F4 reference lies outside this entire constant-gauge family, although its pressure does not improve the seed.

## Contribution and attribution

Start with the [contributor quickstart](docs/START_HERE.md) and
[CONTRIBUTING.md](CONTRIBUTING.md). Work directly on `main` when you have
push access. Pull requests are mainly for external contributors or when
isolation and separate integration review are useful; every useful PR is
intended to merge, including partial and negative research. A merge records
integration, not mathematical acceptance.

For substantial research, read [autoresearch.md](autoresearch.md) and coordinate
through [GitHub Discussions](https://github.com/1337isnot1337/integer-factoring/discussions),
starting with the [coordination thread](https://github.com/1337isnot1337/integer-factoring/discussions/51).
Reviewers can use the [review guide](docs/REVIEWING.md) and
[exact-head review ledger](docs/REVIEW_LEDGER.md). Document source revisions,
proof status, reproduction commands, and remaining gaps. The
[contribution record](CONTRIBUTORS.md) retains credit for useful partial and
superseded work as well as incorporated results.

The source factoring reduction and counting infrastructure are credited to Rudich and Rousseau. The Appendix B evaluator is credited to Anthropic. This repository's new integration and proof files used substantial OpenAI assistant support. The source-compatible license is [Apache-2.0](LICENSE).

## Further exact profile restrictions

The [regular-resolvent determinant budget](research/regular-resolvent-determinant-budget.md)
excludes the bare Smith criterion at every alphabet size from 2 through 18,
and also 21, over every field and at every leading rank. Endpoint-preserving
fractional-linear moves remain in this family. Leading rank at most two is
excluded at every size. A separate [complete binary four-label census](research/binary-rank-two-channels.md)
checks 117,532,800 normalized rank-two parameters and proves at least eleven
order-one pivots in every regular case. The [truncated multiplication argument](research/shift-register-toeplitz-budget.md)
improves the copied-register bound from a division by `J^3` to `J`, reducing
its complete even-register cutoff from 576 to 380. These results constrain
the declared profile methods; they do not change the factoring score.

The [fixed-dimensional cube transfer](research/multiaxis-cube-transfer.md) now supplies a complete evaluator reduction for every fixed number `d >= 2` of axes. A common-piece family with total cut-width coefficient `rho` gives score `log_d(1/rho)` after all regrouping, table queries, and recovery costs are charged. It also permits explicitly bounded subexponential piece counts. An actual nine-axis family with `rho < 1/3` would beat one half; the existing relaxed profiles do not construct that family. The achieved factoring score remains `1/2`. Run `make verify-multiaxis-cube`.

[Multi-axis type amplification](research/multiaxis-type-amplification.md) strengthens that transfer: subexponential local piece entropy and a vanishing relative width slack suffice for the full limiting score `log_d(1/rho)`, with matching workspace, on one uniform machine. A sharper finite-seed test uses the sum of the geometric means of all cut ranks. This removes the direct piece-enumeration score cap; no improved numerical seed is supplied. Run `make verify-multiaxis-types`.

[Binary product-power weights](research/pole-cauchy-product-power-weights.md) extend the normalized inverse-contraction engine to products of explicit nonvanishing polynomial bases with signed binary exponents, including nonunit bases and exponentially large expanded degrees. Common dyadic partitions charge each base once, and the output retains its potentially enormous normalization exponent in binary. This is a scaled-precision guarantee; absolute precision after a large negative scale remains charged. The singular joint factor-count contraction is still open. Run `make verify-product-power-weights`.

[Symplectic exterior blocks](research/cube-symmetric-exterior-blocks.md) reduce every symmetric characteristic-two multiplication operator on `r` square-zero generators to at most `floor(r/2)+1` distinct blocks, of order at most `2^(floor(r/2)+1)`, with exact multiplicities. The resulting complete Smith profiles give rational all-scale exclusions for the natural homogeneous Cayley family through twelve axes. This narrows that construction search; it does not exclude other cube tensors or raise the factoring score. Run `make verify-symmetric-exterior`.

[Exact dyadic Bézout denominators](research/interval-bezout-precision.md) determine the least global polynomial normalization for residue-separated product domains. The [balanced-node refinement](research/balanced-bezout-denominator.md) proves `r=m+v2((m-1)!)` for every size, including independently deformed dyadic isometries. Its parity interpolation proof also reduces enormous-index normalized coefficients to low symmetric moments and shows that only `2p` terminal Newton terms survive modulo `2^p`. The normalization exponent is stored in binary; general factor extraction remains open. Run `make verify-interval-bezout verify-balanced-bezout`.

The [complete modular Bézout compiler](research/balanced-bezout-modular.md) now constructs both optimally normalized global polynomials for the separated arithmetic-progression family in polynomial bit time. Each polynomial has fewer than `4p` possible terms in powers of `X-1` modulo `2^p`, despite an exponentially large degree and normalization exponent. Independent dense rational controls and three implicit domains of up to `2^200+13` nodes verify the complete coefficient output and charge factorial-ratio guard precision explicitly. This resolves normalized construction for that family, while unknown common-root extraction remains open. Run `make verify-balanced-bezout-modular`.

A [deeply matched progression construction](research/matched-progression-bezout.md) computes complete primitive Bézout polynomials modulo a requested power of two for reflected arithmetic progressions with binary encoded, exponentially large counts. The proof controls reciprocal perturbations and accounts for the complete polynomial tails; its stable coprime family does not change the factoring exponent.

[Tensor-product Pascal boxes](research/tensor-pascal-box-crossings.md) give an exact multi-coordinate crossing construction and coherent energy types. The accompanying pressure calculations and integer-polynomial certificates delimit the tested inherited gauges; they do not supply a score above one half.

The [rectangular Bézout denominator theorem](research/rectangular-bezout-denominator.md) determines every inverse Newton coefficient valuation for two independent dyadic isometry prefixes of unequal sizes. A double Mahler expansion proves the formula, and a binary-size maximization gives the exact primitive scale and short supports for both complete normalized polynomials.

The [rectangular modular compiler](research/rectangular-bezout-modular.md) constructs those complete primitive Bézout pairs on independent arithmetic domains in polynomial time in their binary descriptions and requested precision. Exact rational controls and independently sized domains through 257-bit counts verify signs, precision guards, nonmonotone coefficient windows, and the complete polynomial identity.

[Weighted Pascal flags](research/weighted-pascal-flag-crossings.md) enlarge the exact crossing search to arbitrary positive integer weights, different energy spectra at the two cuts, and independent physical gauges. Explicit Clifford matrices recover both full Fock profiles from much smaller matrices; the checked constructors do not yet cross the half-score threshold.

The [entrywise Cauchy-power theorem](research/pole-cauchy-hadamard-rank.md) gives exact Smith profiles and optimal separation ranks at every finite precision. Standard trace padding can retain exponentially many nonfactor rows; a larger polynomial displacement reduces the optimal rank exactly to the factor count. Constructing that small-rank representation remains a separate algorithmic problem.

A [quadratic coalescence testbed](research/quadratic-coalescence-bezout.md) now handles reflected progressions with two different matching depths. A normalized weighted reciprocal sum has an odd unit computable from only linearly many Mahler coefficients in the requested precision; this determines the exact primitive scale and complete coefficient residues in the proved stability range.


The [exact signed affine-code screen](research/canonical-affine-code-visibility.md) proves the sharp 32-bit visibility threshold for all 2,013,762,304 eight-wire point-row and signed affine-code column pairs. Original-coefficient replays and positive controls accompany the exhaustive finite result; a general efficient decoder is still required.


For arbitrary odd inputs, the [fixed-scale global Bézout interface](research/regularized-bezout-factor-moments.md) now isolates factor-supported moments in a short complete reverse jet. The [displaced extraction theorem](research/displaced-bezout-extraction.md) supplies explicit primitive-scale and conditioned reconstruction bounds. With additional ordinary interval padding, [one leading coefficient](research/leading-bezout-factoring.md) at linear bit precision becomes a unit times a positive divisor sum, giving a polynomial-query factoring reduction without any factor-count bound. The normalized global coefficient evaluator remains the missing algorithm; the current factoring score is unchanged.


The [filtered Pascal exponent classification](research/filtered-pascal-exponents.md) supplies nonlinear and noninjective exponent maps that preserve a complete zero-cut profile. It characterizes the fixed regularity cone exactly and verifies the actual polynomial and Fock matrices. The other cut still requires its own profile calculation; no exponent improvement is inferred from the preserved spectrum alone.


A [prime-adic precision-one refinement](research/prime-bezout-counting.md) turns each normalized leading coefficient into a rectangle factor count modulo a small prime. Polynomially bounded primes and disjoint residue rectangles yield a uniform conditional factoring reduction with only logarithmically many output bits per query. Exact rational controls verify the normalization; computing the global coefficient remains unresolved.


The [one-prime recovery refinement](research/one-prime-bezout-recovery.md) uses a chosen field and an explicit trial budget. Counts become exact on certified prime-only intervals, so a single characteristic suffices even when the total divisor count vanishes modulo that prime. The resulting conditional cost is exponential only in the charged trial budget and the still-missing coefficient provider.


The [complete filtered-unit curve family](research/filtered-unit-crossings.md) extends nonlinear group-direction selectors. It gives explicit independent coordinates, full endpoint and generic-rank proofs, and exact coefficient-cut formulas. The larger family is now available for searches; its retained finite examples do not improve the factoring exponent.


Repeating the symplectic reduction gives [polynomial-size blocks for symmetric exterior multiplication](research/recursive-symmetric-exterior-blocks.md): exactly `ceil(r/2)` block types, each of order at most `(r+2)^2/4`, retain the complete `2^r`-dimensional operator. This makes exact symmetric specialization profiles practical through 64 axes; it is a matrix-computation improvement, with no new factoring exponent claimed.


The [fixed-field count-value refinement](research/bezout-count-values.md) reduces the precision requirement to `O(log n)` bits: differences of two normalized inverse values detect prime divisors in certified windows. It also implements an exact precision-one provider for intervals with a small symmetric difference, including implicit large singular examples. The general value provider remains open, and its full cost is included in every conditional bound.


The [binary-linear Pascal extension](research/linear-pascal-matroid-energies.md) computes complete zero-cut spectra from minimum-weight projected row bases. It supplies actual kernels outside the earlier filtered class and charges their changed physical energies. The other cut remains part of the required contraction calculation.


The [arbitrary-permutation Pascal theorem](research/permuted-pascal-zero-spectra.md) removes the linearity restriction: every exponent permutation now has a proved complete zero spectrum. A compatible cone supplies explicit infinity spectra for the same kernel. Both matching-pair energies are retained in the pressure calculation; an improving contraction certificate is still required.


An [all-scale comparison](research/permuted-pascal-cone-dominance.md) shows that the explicit common-scale permutation cone is dominated by the original ordered Pascal kernel once both physical charges are included. This rules out that particular search region while leaving arbitrary physical gauges and alphabet restrictions open.

The [moving-map construction](research/moving-boundary-crossings.md) permits
polynomial physical embeddings and retractions with an exact paired identity.
All 256 affine F4 parameter pairs of its first eight-letter boundary family
have complete cut profiles. A separate exact theorem proves natural pressure
`rho >= 5/7` for every polynomial deformation retaining its first five
coordinate letters, with no degree bound. Both core profiles agree with an
independent binary-polynomial Smith implementation. This closes boundary-only
natural-pressure searches in that family; changes to lower states and other
gauges remain open. The factoring score stays `1/2`. Run
`make verify-moving-boundary` for the construction and all-scale certificates.

The same [moving-map note](research/moving-boundary-crossings.md#codimension-one-frames-with-no-protected-lower-coordinate-core)
now gives polynomial and Laurent physical frames that can move all eight
states. Common monomial clearing preserves the coefficient ranks and
opposite-cut product pressure. A complete 1,696-evaluation interior-pair
screen and 144 declared frame pilot cases give no improved bound. Dense
Laurent contraction checks and all complete finite profiles are retained;
the larger frame search remains open.

An [exact-source receipt](evidence/moving-map-verification/receipt.json)
retains the full default pass at `410e3e2` (341 Python command blocks,
66 unit tests, five crossing tests, 973 Markdown files) and four focused
moving-map checks at `1dba7ba`. These are source-specific ordinary checks;
later integrations and the full frame search have separate status.

The [same construction now permits endpoint-only splits](research/moving-boundary-crossings.md#only-the-endpoint-split-is-necessary),
including intermediate rank loss, with 100 complete finite controls.
An [exact quarter-prefix certificate](research/moving-boundary-crossings.md#an-exact-quarter-prefix-exclusion)
proves natural pressure at least `5/7` for all 16,384 retained height-flow
frames at every common scale and physical type; the full 65,536-case run
is separate and incomplete at that snapshot. Full-gauge fits remain
iteration-limited evidence. No improved factoring exponent is claimed.

A separate [asymmetric screen of those complete profiles](research/moving-boundary-crossings.md#asymmetric-supporting-line-screen-of-the-frozen-prefix)
checks 11,010,048 numerical angle/tilt points and refines twenty cases.
Its normalization agrees with the existing native supporting-line oracle;
the best retained support value is about `1.000419`, above the required
one. This numerical result leaves unsearched parameters open and supplies
no new factoring exponent.

A [complete F4 interior-frame screen](research/moving-boundary-crossings.md#a-complete-f4-family-that-changes-one-protected-letter)
changes one of the five protected letters in 1,280 declared parameter cases.
All complete profiles are retained. A further 860,160 asymmetric grid points
and twenty bounded refinements give no support value below one. The full
binary frame run and more general gauge choices retain their separate scope.

A [pinned source audit of RiemannGaussian's semiprime work](research/riemann-gaussian-semiprime-source-audit.md)
identifies reusable deterministic boundary contraction and exact factor-moment
encodings. Its current target has the same logarithmic saving exponent `1/2`.
The semiprime promise, coefficient-ring interface, and complete evaluation
cost remain explicit obligations before importing a factoring consequence.

## Sparse oriented power prefixes

The [sparse residual selector](research/sparse-oriented-power-prefixes.md)
finds factors on certified 64-bit and 80-bit inputs after 2,048 and 10,240
queries from the unchanged prime-power stream. The original selector needs
3,983 and 38,818 queries respectively. Every power, ideal identity and integer
presentation map is retained in a compressed replayable certificate. These
are finite query gains: the baseline implementation is slower in wall time,
and the uniform factoring exponent remains `1/2`.

An [exact Kloosterman block evaluator](research/circuit-box-kloosterman-blocks.md)
contracts an entire unit-frequency orbit and a ratio progression while keeping
both ordinary interval weights. At conductor `2^1024`, the retained control
aggregates `2^1535` frequency pairs with at most eight roots and sixteen affine
counts. In the no-wrap range it is a positive factor-count query for one ratio
class. The full outer sum still costs `2^(m/2)*poly(m)`, so no score improves.

A [cubic graph construction](research/tensor-cubic-graph-orthogonality.md)
provides coherent orthogonal maps at every inherited-height cutoff precisely
when `1` is outside the binary span of three field coefficients. Its explicit
Gram determinant proof works in every characteristic-two field. Exact mixed
constructors and two completed 51-letter profiles are retained. These pilots
do not improve the pressure frontier or the factoring exponent.

[Sparse unit elimination](research/sparse-unit-presentation.md) preserves the
entire reference pivot sequence and all integer presentation certificates.
Stable identifiers and incidence updates reduce the matched 80-bit final-pass
median from 13.116 seconds to 0.853 seconds, with verification included.
This is a measured implementation improvement; the relation-search bound is
unchanged.

[Exact quotient tree updates](research/quotient-update-oriented-basis.md)
replace large ancestor recomputations with exact division or multiplication
by the changed leaf. Every atom, signed word and gcd decision is preserved.
On the retained 80-bit prefix, three matched component runs reduce the
construction-and-export median from 17.864 seconds to 1.171 seconds; all
9,776 literal ideal identities replay independently.

[Positive ratio lookahead](research/circuit-box-ratio-lookahead.md) turns the
inverse-window point query into a complete adaptive search: each ratio bit
checks two product bits with exact affine counts. A divisor-count argument
gives `2^(n/3+O(n/log n))` bit time and polynomial space, including direct
smallest-prime-factor recovery. This improves that representation's earlier
half-exponent enumeration, while remaining weaker than the repository's main
theorem and older general factoring algorithms. The scoring exponent stays `1/2`.

[Compact unit traces](research/unit-presentation-traces.md) preserve the same
integer presentation by independently checking each elementary substitution.
On the 80-bit fixture, the additional unit certificate shrinks from 3,113,278
to 21,638 JSON bytes while sharing the complete graph certificate. The final
rows, retained generators and torsion image agree with the expanded maps.
The full relation-search and factoring exponents remain unchanged.

The [optimized oriented selector](research/optimized-oriented-power-search.md)
reuses the unchanged search body with the two exact component replacements.
Its complete 80-bit result matches the baseline while reducing the retained
search time from 67.006 to 9.694 seconds. A new certified 96-bit input factors
after 20,480 queries versus 35,106 for the original exact-key selector. Full
witnesses, both integer presentation maps and the original comparison replay;
the original selector remains faster in total time, and the score is unchanged.

A [literal Pfaffian constructor](research/tensor-literal-low-cell-retractions.md)
now evaluates homogeneous physical retractions without deletion closure,
retaining every within-cut Gaussian term. Independent original-kernel,
Pfaffian, endpoint and Toeplitz controls accompany a complete capped 55-label
pilot. It keeps more upper states than the coherent kernel construction,
but this first numerical profile does not improve the pressure frontier.

A [lattice terminal for ratio classes](research/circuit-box-ratio-lattice.md)
recovers every local factor pair using at most four residue roots and twenty
lattice lines once `R^3` exceeds the rectangle's product height. Streaming
these queries gives `2^(n/3)*poly(n)` time and polynomial space, removing
the preceding divisor-bound overhead. This is a classical cube-root scale,
not an improvement over the main theorem or the scoring exponent.

[Primitive integer row selection](research/primitive-integer-row-selection.md)
replaces rational-entry rank tests while retaining exactly the same original
integer relations. Only temporary test vectors are divided by their content.
On the 96-bit core, five paired runs reduce median selection CPU time from
0.606 to 0.141 seconds; 17,000 finite acceptance decisions agree with the
reference. This changes component cost without weakening torsion information.

The [prepared bucket basis](research/batched-coprime-bucket-basis.md) preserves complete oriented ideal words while sharing frozen remainders across supplied batches. Full 80/96-bit fixture replays pass; matched 96-bit basis CPU falls from 4.8112 to 1.7320 seconds. This is a supplied-basis implementation improvement, with no new uniform relation-supply or factoring exponent claim.

The [compact oriented selector](research/compact-oriented-power-search.md) combines quotient tree updates, primitive integer rank tests, and independently replayed unit traces. Its complete 80-bit result has the same 10,240 queried norms and factor as the retained baseline, while the incremental unit certificate shrinks from 3,113,278 to 21,639 JSON bytes. No relation is divided or omitted.

The [larger F16 profile study](research/tensor-f16-capped-search.md) retains complete capped 108-label cubic and three-pair profiles, including every unresolved high-order multiplicity. Their numerical upper pressure fits are 0.710813125460 and 0.712539258890. The [literal low-cell study](research/tensor-literal-low-cell-retractions.md) proves the field-independent floor rho >= 27/38 for the specified 55-label height-15 retraction family, and completes all F4 splits of the earlier height-8 and height-9 cells at 23 labels. These results do not improve the compact pressure frontier or the factoring exponent.

The [prepared selector](tools/prepared_oriented_power_search.py) combines 512-key remainder batches, small-magnitude primitive integer pivots, and compact unit traces. Complete [64/80-bit matched controls](evidence/prepared-oriented-power-search.json) preserve every mathematical search field. The [128-bit input and declared policy](evidence/prepared-oriented-fresh-128-policy.json) fix the 262,144-query cap, 16,384-query checkpoints, 1,024-column extraction cap, 2 GiB address-space cap, and 24 source hashes before the first query. A successful search requires separate saved-artifact replay; this policy is not a result claim.

The [supplied-basis charging theorem](research/oriented-basis-linear-charges.md) bounds orientation refinement by at most the total input logarithmic size in gcds. Combined with Bernstein's classical constructive norm-basis algorithms, it gives S(n + log(S+2))^O(1) bit processing for supplied oriented ideal keys. A separate potential argument proves at most twice the summed input floor-logarithms minus the nonunit-key count in total sparse usage writes for the existing incremental implementation. The executable wrapper and full ideal controls are retained; neither result supplies the missing uniform relation-search bound.

The [complete 112-bit oriented prefix](research/oriented-power-112-prefix.md) factors the certified input after 73,728 queried norms versus 96,752 for the original exact-key selector. Its saved 12.5 MB compressed certificate retains every witness, all 72,687 literal ideals, the complete graph/unit maps, and both independently replayed searches. The oriented method uses fewer queries but more time (259.123 versus 35.585 CPU seconds); this is a finite control, with no uniform exponent improvement.

The [classical quarter-bit terminal](research/circuit-box-ratio-small-roots.md) now has an exact implementation: explicit integer LLL parameters, independently reconstructed short-polynomial certificates, and Sturm isolation avoid scanning exponentially long root intervals. The complete reference algorithm is 2^(n/4) poly(n), handles all factor-size ranges, and is weaker than the established repository theorem. The 255-bit controls are local known-factorization tests, not complete factoring benchmarks.

A [two-pass orientation construction](research/circuit-box-orientation-batch.md) reduces both norm separation and root-sign separation to classical coprime factorization on O(S) integer data. Sparse incidence bookkeeping preserves every literal ideal word, including prime powers and repeated norms. It strengthens the supplied-data processing argument without a relation-generation or score claim.

The [certified 192-bit run](research/batch-smooth-fresh-192.md) returns a proper factor after 48,097,746 generated norms, with all 1,481,611 generated-ahead norms charged. The complete 22.2 MB certificate independently replays 98,968 partial ideal identities, 95,676 prime-residual certificates, and 3,292 integer DAG eliminations. Three checksummed archive segments reconstruct the exact saved certificate. Fixed-base rank is 3,292 of 3,299. This is a finite control with no matched 192-bit speedup or new score claim; the declared 208-bit follow-on started after saved-artifact replay passed.

The [classical coprime engine](research/classical-oriented-basis-engine.md) now constructs all supplied norm data and exact oriented ideal words without a factoring oracle. Complete 64-, 80-, and 96-bit fixtures match their saved bases. This makes the classical supplied-data bound executable; the faster prepared implementation remains the practical search path, and the score exponent remains 1/2.

The [literal full-tail reduction](research/tensor-literal-full-tail-core.md) handles fixed low-mode maps, including maps mixing heights, through an explicit 8k-mode active core and charged spectator factors. Independent literal checks and 78 complete F4 profiles are retained. The tested rank-one families do not improve the current tensor frontier or the factoring score.

[Sparse intermediate graph projections](research/sparse-projected-signed-graph.md) preserve every integer presentation map while avoiding a redundant dense allocation. Independent component and full 64/80 integration replays pass. Both 2-GiB 128-bit failures are retained with their available costs; the second reached rank 418 at 245760 queries. A separate unchanged-source 3-GiB retry is declared, with no success claim.

[Small internal integer pivots](research/small-pivot-integer-row-selection.md) retain exactly the same original relation subsequence. Three paired measurements reduce rank-selection CPU medians from 0.137 to 0.068 seconds on the 96-bit core and from 16.918 to 7.002 seconds on the 112-bit core; these are component measurements.

The [two-pass oriented constructor](research/circuit-box-oriented-batch-controls.md) now uses two explicit classical coprime-basis calls. Unitary preprocessing removes irrelevant second-pass atoms while preserving literal ideal words; the raw-difference option is also tested. Independent replay checks 2932 ideal identities and 3357 prime expansions, including 24 variant comparisons. This improves supplied-data processing without a new relation-supply bound.

[Completed Laurent-frame and F16 searches](research/tensor-weighted-laurent-frame-search.md) retain the full 256-case binary frame census, 152 F16 representatives covering 544 normalized pairs, and finished 115-, 91-, and 55-letter profiles. The Laurent constructor passed independent raw-tensor, endpoint, reciprocal, and Smith checks. These searches did not improve the numerical frontier; capped data and numerical fits are distinguished explicitly.

The [full-rank affine full-tail census](research/tensor-full-tail-affine-search.md) computes both complete profiles for all 256 binary low-mode map pairs at M=16, retaining 131072 active orders and every spectator factor. The constructor passes 524288 independent raw coefficient comparisons. The Gaussian control remains best in this census; no new numerical or factoring frontier is claimed.

[Equivalent-norm stencils](research/equivalent-norm-stencils.md) add exact complementary norms to supplied form powers and replay all class, ideal, graph, and unit identities. In 24 cases with 57600 charged norm tests, the leading-only method retains at least the refined rank of both alternatives at every equal budget. This negative result is retained; the active policy and score are unchanged.

[Compact reconstructible power witnesses](research/compact-reconstructible-power-witnesses.md) preserve the complete query logic and mathematical result while storing only base words, residual indices, and signs. Full 64/80 replays and 25 additional policies pass. The 80-bit witness section falls from 3390155 to 264539 canonical JSON bytes; this is a representation measurement, not a peak-memory claim.

The [declared compact-witness 128-bit policy](evidence/compact-witness-oriented-fresh-128-policy.json) repeats the same certified input and query schedule with a 2-GiB address-space bound. It tests the representation change against the separately retained projected search. This is a charged new attempt; a saved result requires independent replay before certification.

[Compact oriented construction and verification](research/compact-oriented-basis-input-bound.md) now have a deterministic B^(1+o(1)) bit bound in total supplied input size, counting the discriminant once. Batched discriminant validation, two coprime passes, and a balanced arithmetic verifier avoid per-key large-discriminant work and reduced-form output. Controls include 3146 independent HNF identities and constructed large-discriminant inputs. Relation supply and the factoring score remain unchanged.

The [auxiliary-mode Clifford theorem](research/tensor-affine-clifford-profiles.md) recovers complete Smith profiles for Gaussian circuits interleaved with arbitrary invertible one-mode maps over characteristic-two fields. One auxiliary mode gives a matrix of order 2n+2, with scalar valuations and both cut clearings retained. The 18-row implementation matches all 131072 prior binary orders and independent original F4 kernels. A completed 1792-case F4 degree census finds no pressure improvement.
