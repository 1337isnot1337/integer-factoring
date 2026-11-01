# Reusable semiprime work in RiemannGaussian

## Source and assessment

This is a source audit of David Sanftenberg's
[`dbsanfte/RiemannGaussian` at `f5fd3ddb407972f9ca65298fbe9be3225fe06635`](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635),
checked on 1 November 2026. The current pinned README differs materially from
the older cached GitHub rendering. No upstream code was modified, no Lean
rebuild was run, and no third-party theorem is claimed independently verified
by this audit. Substantial OpenAI assistant support was used for the review.

The work is relevant to the arithmetic constructor, useful-scalar decoder,
and tensor-contraction backend of this repository. It does not currently
supply a stronger complete factoring exponent than our established 1/2.
The ideas can be studied and used in ordinary proofs without importing its
large Lean development.

## The current target and the runtime boundary

The [pinned README](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/README.md)
targets every binary semiprime in time

```math
2^{C\ell/\sqrt{\log(\ell+2)}}(\ell+1)^d.
```

Here ell is the input bit length. This is the same logarithmic saving
exponent 1/2 as our complete result, rather than a proposed exponent above
1/2. The old sixth-rate milestones are not the latest completion target.

The [current algorithm entry](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1HybridAlgorithm.lean)
states correct factor-byte output for all semiprimes, including squares, and
an eventual bound on the changed program's retained native report. It
explicitly leaves the complete primitive-machine bit-time lowering separate.
The [original call-rate module](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1OriginalCallRate.lean)
defines its input scale as `ell/sqrt(log(ell+2))`.

The missing primitive theorem is not, by itself, evidence that the ordinary
algorithm is wrong. Conversely, a bound on a native report is not already a
bit-operation bound. Reuse as a runtime result requires auditing the complete
cost of source construction, full coefficient words, control, addressing,
memory operations, and report erasure. The source is careful about this
distinction; it should be preserved here even while working without Lean.

## First reuse priority: deterministic boundary contraction

The [installed boundary-hybrid packet](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/docs/research-evidence/semiprime-boundary-hybrid-20261029/README.md)
and its [ordinary research note](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/docs/research-evidence/semiprime-boundary-hybrid-20261029/source-note/RESEARCH.txt)
give a deterministic exact contraction for a loopless multigraph with M
binary edges and maximum degree Delta. The acquired frontier width obeys

```math
W\le M/2+3(M^2\Delta)^{1/3}.
```

With explicitly stored local tables, the time is a polynomial in the graph,
colour count, and complete coefficient-word size, times `2^W`. The
construction pays for the partition and uses sorted records and forward
joins. It includes parallel edges, isolated vertices, full jets, and zero
precision; it does not assume free lookup or multiply two exponential
catalogue sizes. Coefficient arithmetic must be charged at the size of the
whole coefficient vector.

For Delta=o(M), this is a half-leading-exponent leaf algorithm. It could
replace a brute-force leaf or support a sharper explicit cost invoice in
our evaluator. Its effect on the entire recursion must be recomputed: the
leaf saving is not automatically a smaller global volume-contraction ratio
or a logarithmic saving exponent above 1/2. The upstream packet explicitly
does not claim a stronger global complexity class from this improvement.

This is the strongest immediate algorithmic reuse candidate from this audit.

## Second reuse priority: factor moments and honest count recovery

The [packed-network theorem](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PackedNetworkOriginal.lean)
identifies the full acquired multiplication-network contraction with a
useful scalar. For N=pq and B=2^ell, its channels encode

```text
distinct primes:  count = 2, factor sum = p+q;
prime square:     count = 1, factor sum = p.
```

The [arithmetic readout](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PackedMomentSemiprime.lean)
uses a modulus B squared and proves the channel bounds before unpacking.
Unknown primes occur in correctness hypotheses, rather than as program
inputs. Some underlying network identities already apply to every positive N.

This is closely related in purpose to our
[single useful Bezout scalar](leading-bezout-factoring.md) and
[divisor-coefficient decoder](weight-one-divisor-coefficients.md): a precisely
specified arithmetic quantity reveals factors once it has been computed.
The essential remaining question is the complete acquisition/evaluation
cost; a scalar identity alone is not a faster factoring algorithm.

The coefficient ring matters too. The packed readout uses arithmetic
modulo `2^(2 ell)`. Our characteristic-two/F4 crossing interfaces do not
directly provide arithmetic in that ring. Reuse in those interfaces needs
an explicit encoding or recovery argument with its full precision cost.
The generic boundary-contraction algorithm itself works over a commutative
semiring and does not have this particular compatibility restriction.

There is a concrete general-input scope issue. The
[public correctness join](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PublicCorrectness.lean)
uses a mod-4 count readout with semiprime hypotheses. For N=12 the four proper
ordered factor pairs are `(2,6),(3,4),(4,3),(6,2)`, so the ordinary witness
count is zero modulo four. This is a counterexample to a generic
nonzero-count test, not an execution claim about the whole upstream program.
The semiprime completeness theorem cannot simply have its prime hypotheses
deleted.

A possible ordinary-proof extension is to retain the count modulo 2^ell.
An ordinary proper factor pair is determined by its first factor, which is
strictly below N; hence the number of pairs is at most N<2^ell. At that
precision the root count cannot wrap to zero on a composite. Extending the
actual prefix counter, its source compiler, and its full cost bounds to
this precision remains a proof obligation. The
[all-input prefix adapter](sanftenberg-all-input-prefix-adapter.md) adds an
ordinary minimum-factor argument: order the operand-0 wires by descending
numeric significance so exact zero-first prefix search selects the least
divisor. It also computes that a direct whole-root application of the boundary
contractor gets only a linear width bound. A recursive exact-count or
Boolean-existence implementation remains open; this is not an established
all-input runtime reduction.

## Consequence for this campaign

Reuse the deterministic partition/contraction method and the exact
factor-moment interfaces with their source attribution. Compare the
constructor and acquisition obligations against our arithmetic framework.
Keep the semiprime promise, native-report bounds, ordinary cost arguments,
and primitive-machine theorems distinct. None of these source statements
currently warrants increasing the complete score recorded in this repository.
