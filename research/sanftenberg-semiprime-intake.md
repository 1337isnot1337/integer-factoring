# Sanftenberg semiprime source intake

Inspected source: David Sanftenberg's
[`dbsanfte/RiemannGaussian` at f5fd3ddb407972f9ca65298fbe9be3225fe06635](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635),
after the owner supplied his invitation to reuse its semiprime work.
The relevant source headers name Sanftenberg and Apache-2.0. This intake
reuses the mathematical interfaces and independently replays the small
Python prototype. It does not import Lean sources or report a Lean rebuild.

## 1. The endpoint has the same logarithmic saving target

The upstream
[current algorithm](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1HybridAlgorithm.lean)
states semiprime factor-output correctness, including prime squares,
and a bound for the program's own retained report. Its
[boundary-hybrid status](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/docs/research-evidence/semiprime-boundary-hybrid-20261029/README.md)
explicitly leaves the complete primitive-machine bit-time theorem open.
Its target is

```math
2^{C\ell/\sqrt{\log(\ell+2)}}(\ell+1)^d,
```

with ell the input bit length. This has logarithmic saving exponent1/2,
the same exponent as our [complete machine theorem](half-rate-factoring.md).
Halving a leaf's leading constant or tightening a recursion call-count
constant does not increase that saving exponent. No stronger complete
bound is supplied by this intake.

## 2. A genuinely reusable boundary contraction lemma

For a loopless multigraph with M binary edges and maximum degree Delta,
explicit scalar vertex tables over a commutative semiring admit a
deterministically acquired colour partition with contraction width

```math
W\le M/2+3(M^2\Delta)^{1/3}.                         \tag{1}
```

The [ordinary proof and prototype](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635/docs/research-evidence/semiprime-boundary-hybrid-20261029/source-note)
construct the partition by conditional expectations. Boundary tensors
are merged over shared edge variables, with one exponential record
factor, giving poly(|V|,M,b)2^W work when each complete coefficient
value occupies b bits and its arithmetic costs poly(b). Parallel edges,
isolated scalar vertices and nonunits are retained. The current
[closed width theorem](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1BoundaryPartitionWidthRate.lean)
names the acquired partition, rather than assuming a good one as advice.

This is useful for a future evaluator whose explicit graph has few edges
and controlled full coefficient storage. It also gives a concrete model
for charging intermediate records and complete mixed jets.

For our current n-by-n [braid circuit](quantum-block-braid-crossings.md),
the ordinary gate graph has Theta(n^2) edges and bounded local degree.
Thus (1) alone guarantees only 2^O(n^2) contraction. The retained
[path quotient](quantum-braid-path-quotient.md) already permits
2^O(n) entry construction. A useful port would need a smaller graph or
a stronger width estimate. Equation (1) is not the missing favorable
two-cut rank-sum bound for the crossing.

Ordered matrix auxiliaries also require care: the theorem's scalar
semiring is commutative. An auxiliary can be represented by explicit
state wires and scalar entries, but those wires and their storage must
be included in M and b. Treating a noncommutative matrix product as a
free scalar semiring operation would not justify this contraction.

## 3. A useful single-scalar arithmetic interface

The upstream
[packed moment theorem](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PackedMomentSemiprime.lean)
encodes the exact proper-factor count and first factor moment in one
weighted tensor contraction. For N=pq with distinct primes, the channels
are2 and p+q; for a prime square they are1 and p. The
[original-network identity](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeD1PackedNetworkOriginal.lean)
connects that scalar to the actually constructed network.

For distinct primes, a supplied sum s=p+q gives the factors by
(s plus or minus sqrt(s^2-4N))/2 with polynomial bit-cost postprocessing.
This is closely related to our [divisor-moment](divisor-moment-recovery.md)
and [weight-one coefficient](weight-one-divisor-coefficients.md) lanes.
It suggests another precise evaluator target. The efficient computation
of the scalar remains the resource that must be proved; an exact encoding
of factor information does not make that computation cheap.

The [Teichmuller carry backend](teichmuller-cyclic-carries.md) now preserves
2-adic coefficients on an explicit synchronized-polynomial input class,
including a precision-linear residue-root count compiler. It may support
these exact scalar targets if a small source representation is supplied.
This does not provide that missing representation or the full network cost.

## 4. What was independently checked

The two upstream Python files were downloaded at the pinned commit;
their Git blob hashes match exactly. Running `python test_boundary_hybrid.py`
passed48 semantic comparisons,16 exhaustive conditional-expectation checks,
35 arbitrary colourings, complete mixed-jet controls, isolated scalar
controls and malformed-incidence rejection. The larger36-edge case tests
partition acquisition only. The
[local replay receipt](../evidence/sanftenberg-boundary-replay.json) records
source hashes and the exact scope. These finite checks support prototype
correctness; they do not independently verify the Lean build, the general
width theorem, or the outstanding whole-machine runtime claim.

Priority: retain the boundary contractor and packed scalar as attributed
candidate components for new evaluator constructions. The current braid
lane still needs a complete favorable rank sum, and the arithmetic lanes
still need their fast original-index coefficient engines. The external
work is relevant to those designs without discharging either obligation.
