# Research progress: 30 October 2026

This continuation began at 05:11 UTC. Its goal is to improve the
deterministic factoring bound beyond the existing exponent one half,
without adding Lean formalizations. This note records the mathematical
results and the remaining gap. The
[verification ledger](../evidence/pr23-beyond-half-verification.json) records
the source pins, tested scope, and checkpoint evidence.

## Strongest unconditional result

For an `n`-bit integer `N>=2`, the strongest completed time and workspace
bound is now

```math
T(n)=\exp\!\left(O\!\left(\frac{n}{\sqrt{\log(n+2)}}\right)\right),
\qquad
S(n)=O\!\left(\frac{n}{\sqrt{\log(n+2)}}\right).
```

The [sublinear-space theorem](sublinear-space-half-rate.md) computes
the least prime factor, complete prime factorization, sorted divisors,
and exact divisor-prefix counts. The input and prefix data are read-only;
output is streamed to a write-only tape. Thus `S(n)` measures working
bits. The result does not compress the input itself into fewer than
`n` bits.

The time exponent remains one half. The advance in this continuation
is the simultaneous sublinear workspace bound and the stronger exact
output interface. The earlier
[composition theorem](half-rate-factoring.md) remains the source of
the half-rate time improvement over the original reference bound.

Three storage changes make the new result possible. The
[local factor algorithms](near-linear-factor-space.md) use compressed
partition descriptions and scalar coefficient evaluations in `O(lambda)`
bits. The [top-table queries](streamed-top-space.md) reconstruct input,
middle, and inverse arithmetic messages on demand, using polynomial
time and `O(n/P+Af)` storage. A designated pass accumulates input
weights once, even when arithmetic values are later recomputed.
Finally, the recursive evaluator stores the actual short local-setting
counters, whose lengths sum to `O(n/sqrt(log n))` along a live path.
It does not retain all tables or all settings.

Recovery also avoids storing an `n`-bit candidate. A hash-bucket
descriptor generates either factor's bits on request. One square-zero
tag at a time keeps the coefficient ring at sixteen elements.
Column-by-column multiplication verifies the two virtual factors with
logarithmic arithmetic storage. Successor scans remove duplicate
representations. Divisor-pair product tests identify primes, and
successive prime-power searches recover their multiplicities. All
queries refer to the original read-only `N`; no cofactor register or
factor list is needed.

The [direct exact-counting refinement](direct-sublinear-divisor-counts.md)
also returns an exact prefix count from one weighted evaluator query,
with the same time and sublinear workspace bounds. It combines the
[integral Schur factors](integral-counting-transfer.md) with logarithmic
type blocks and a uniform divisor bound `tau(N)<=2^q`, where
`q=O(n/log n)`. Evaluation modulo `2^(q+1)` returns the integer
count itself. Large residues occupy only a constant number of registers
per recursion frame; their total storage is `O(q loglog n)`, smaller
than the half-rate workspace. Read-only interval endpoints use at most
`2n+1` prefix queries. Counting needs no preliminary factorization or
witness-isolation transform.

## Better local queries and a broader transfer theorem

The [blocked half-rate family](blocked-half-family.md) reduces a local
piece or factor query from `2^O(lambda)` time to

```math
\lambda\,2^{O(\sqrt{\log(\lambda+2)}\log\log(\lambda+4))}
=\lambda^{1+o(1)},
```

with `O(lambda)` workspace. It uses short regular-partition blocks
and retains their full tuples of piece indices. The resulting piece
entropy is `O(lambda/sqrt(log lambda))`. The endpoint of the relaxed
transfer charges that entropy without an extra depth factor, so the
construction retains the exact half exponent and the global sublinear
workspace bound in the same algorithm. It supplies no contraction
below `1/sqrt(2)`.

The [simultaneous time and workspace transfer](simultaneous-time-space-transfer.md)
applies to any supplied fixed-characteristic family with contraction
`rho`, entropy exponent `delta`, and linear-space local queries. With

```math
\beta=\min\{\log_2(1/\rho),\delta\},
```

it gives time `exp(O(n/(log n)^beta))` and workspace
`O(n/(log n)^beta)`. This includes the conditional endpoint
`beta=1`.

Its type-grouped construction can also remove an original local-space
assumption. When a time-only family has polynomially many pieces and
a power-saving width error, logarithmic blocks make each original
oracle small enough to run within the new workspace. Grouping block
words by their multiplicities gives polylogarithmic piece entropy.
The multinomial channel cost is included in the width error. Exact
finite-field coefficient extraction supplies quasipolynomial scalar
queries. These are constructive interface improvements; a better
contracting family is still required to improve the unconditional
time exponent.

The [rational and fixed-number-field transfer](rational-simultaneous-transfer.md)
also removes the fixed-characteristic hypothesis. It collects the
denominators of the finitely many logarithmic-size block factors into
`D`, then reduces their exact identities modulo the odd integer `2D+1`.
All recorded denominators are units. Type sums and the evaluator need
no field division, while binary character recovery only needs two to
be invertible. This gives polynomial coefficient-ring preprocessing
without prime selection or auxiliary factoring. A fixed-degree monic
quotient gives the same result for a fixed number field. The theorem
is conditional on the supplied crossing profile and its uniform
rational-coordinate algorithms.
The pinned Appendix B integer family satisfies those hypotheses, so
the same theorem also gives a second unconditional proof of the
half-rate workspace bound directly from its original time-only
local programs. Their exponential pads are constructed only for
logarithmic-size blocks, and their finite-tape bit cost is included.
Using an odd modulus larger than the divisor bound also gives a
second derivation of the unconditional one-query counting theorem
from that same pinned family.

For direct exact counting, a large coefficient modulus would add a
`loglog n` scalar-storage factor at the conditional endpoint
`alpha=1`. Sequential generalized Chinese remainder reconstruction
removes that factor. Evaluate modulo small integers coprime to the
recorded denominators, retaining only the current residue and its
modulus. An elementary least-common-multiple bound guarantees enough
combined modulus after polynomially many calls. This gives the matched
time and workspace bounds for direct counts at every supplied exponent
`0<alpha<=1`, over rational factors or a fixed number field. It uses
no prime search and no factorization. The unconditional one-call theorem
above remains at `alpha=1/2`.

## Arithmetic networks and faster base evaluators

The [logarithmic-wire network](logarithmic-wire-network.md) gives a
second multiplication reduction with `O(log n)`-bit messages and
`O(n)` total information volume. Bivariate interpolation makes its
top queries polynomial time. Its Fourier form can contract whole
arithmetic blocks and carry segments. The
[mixed linear-and-carry theorem](linear-carry-blocks.md) extends this
to arbitrary allowed subsets of the linear and carry gates. A block
containing the middle multiplication gates together with those other
gates still needs a further contraction argument.

The [faster clique baselines](fast-matrix-clique.md) improve the base
cost inside the recursion. Balanced binary conditioning changes the
leading exhaustive coefficient from `1/2` to `1/3`, while retaining
linear assignment space. A separate three-block matrix product lowers
the arithmetic-operation and table-query coefficient further; an
explicit rectangular split gives approximately `0.310064`. That
version uses exponential matrix storage. Neither change alters the
current logarithmic saving exponent, and the dense-matrix constant
is not asserted as a finite-tape bit-time constant.

## New constructors, and what their certificates establish

The search produced exact
[intersection-matching crossings](intersection-matching-crossing.md),
[balanced-partition crossings](balanced-partition-crossings.md),
and [homogeneous permutation-label families](permutation-homogeneous-labels.md).
Their unit pairings, entry algorithms, and finite profiles are proved
and checked. These constructions enlarge the available search space;
none currently has a proved uniform contraction below `1/sqrt(2)`.

The [six-label split seed](split-permutation-seed.md), for example,
has an exact finite worst rank product of 150. Its attractive
single-piece coefficient does not pay for the number of pieces under
amplification. The
[finite rank-sum criterion](finite-seed-amplification.md) charges that
cost through `sum_i sqrt(r_i c_i)`. A small finite maximum and a
uniform asymptotic improvement are distinct claims.

The following obstructions specify where the present searches stop.
Each has narrower scope than an impossibility theorem for factoring.

| Construction or method | Proved restriction | Routes left outside its scope |
| --- | --- | --- |
| [Full regular homogeneous pieces](regular-fusion-rank.md) | Their actual central ranks attain the half-rate exponent on the full domain at powers of two. | Further label restrictions and recombinations of coefficient pieces. |
| [Coordinate matching support](coordinate-support-rank.md) | A fixed orientation mask within each piece forces `rho>=3/4` when piece entropy is negligible. | Keeping sums across masks together as a new piece. |
| [Permutation-invariant sparse defects](fixed-defect-invariant-barrier.md) | The specified sparse-defect families have pieces with nearly full ranks in both cuts. | Restrictions breaking that symmetry, or other defect regimes. |
| [Full permutation coefficient spans](permutation-span-cycle-obstruction.md) | In characteristic two, connected-cycle blocks exclude the finite half-rate rank-sum target in those spans. The [sharp rank budget](distinguished-piece-rank-budget.md) strengthens the four-point bound. | Dual linear label restrictions, nonidentity kernels, other characteristics, and new collective recombinations of fixed-seed tensor powers. |
| [Conserved convolution kernels](conserved-convolution-crossings.md) | Full group domains have a two-cut block obstruction; exact Smith-envelope certificates also cover the listed small cyclic groups. | Other label embeddings and coefficient-rank arguments beyond that envelope. |
| [Ideal Hankel normalization](ideal-hankel-endpoint.md) | Its regular-unipotent leading coefficient cannot have opposite-cut rank one for `d>=3`. | Parameter-dependent gauges, higher endpoint rank, and other Smith profiles. |

The current [all-field small-seed certificate](minimum-cut-rank-budget.md)
excludes the strict finite rank-sum sufficient criterion through twenty-one
labels. A seed passing that particular criterion therefore needs at
least twenty-two labels. A separate
[two-piece theorem](two-piece-rank-sum-obstruction.md) requires at least
three nonzero pieces at every dimension. These statements do not
exclude a different amplification method on a smaller seed. Their
proofs use swap/product and minimum-cut inequalities with exact
rational rank-budget certificates, rather than a search through all
tensors.

Combining the small-label certificate with surviving cycle blocks also
excludes [every coordinate deletion of the 24-label split-permutation
span](permutation-coordinate-restrictions.md) in characteristic two. General dual linear
embeddings and nonidentity mode kernels remain outside that result.

## Verification level and remaining objective

The results are published on
[`research/beyond-half-20261030`](https://github.com/1337isnot1337/integer-factoring/tree/research/beyond-half-20261030)
and collected in [draft PR #23](https://github.com/1337isnot1337/integer-factoring/pull/23).
Run `make verify` to replay the repository's exact checks. The
[verification ledger](../evidence/pr23-beyond-half-verification.json)
identifies the checked files and records the distinction between the
finite evidence and the all-input proofs.

The substantive claims above have ordinary mathematical proofs, source
audits, and exact finite checks. The checks test the stated formulas,
rank calculations, codecs, recovery procedures, and parameter bounds.
They do not run the full enormous asymptotic evaluator on large inputs,
and finite ranks alone do not establish an infinite-family bound.
There are no new Lean formalizations in this continuation.

Internal review and passing repository checks do not constitute
external peer review, leaderboard acceptance, or a practical factoring
benchmark.

The remaining time objective is precise: obtain an unconditional
bound with `alpha>1/2` in `exp(O(n/(log n)^alpha))`. The current
proofs do not establish that. Plausible routes still outside the
obstructions include different label restrictions, nonidentity
annihilation kernels, collective recombinations across tensor powers,
and a more global evaluator that exploits multiplication-network
structure. A successful construction must include its piece entropy,
uniform generation, scalar-query cost, recursive composition, and
deterministic factor recovery. A better finite rank or an uncharged
piece count is not yet the required result.
