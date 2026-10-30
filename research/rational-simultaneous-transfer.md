# Simultaneous time and workspace from rational crossing families

## Conditional theorem

Fix `1/2<=rho<1` and `0<=theta<1`. Suppose a rational crossing
family is available at every positive binary label length `t`, with:

- A nonzero rational scale, and exact displayed factorizations in
  both cuts.
- At most `C(t+1)^a` pieces, separate binary widths at most `2t`,
  and total width at most `4rho t+C t^theta`.
- Uniform algorithms for piece enumeration, metadata, scales, and
  scalar factor entries in `2^O(t)` **bit time**. Rational output is
  given by binary numerators and nonzero binary denominators.

Then, for `alpha=log_2(1/rho)`, deterministic factoring has the
simultaneous bounds

```math
T(n)=\exp\!\left(O\!\left(\frac n{(\log(n+2))^\alpha}\right)\right),
\qquad
S(n)=O\!\left(\frac n{(\log(n+2))^\alpha}\right).        \tag{1}
```

The conclusions include the least prime factor, complete factorization,
sorted divisors, and exact divisor-prefix counts. The input is read-only
and the output is write-only. Constants depend on the supplied family,
not on `N`.

No fixed characteristic or original local-space bound is assumed.
The algorithm constructs one odd coefficient modulus for the input
length and reduces every rational block factor into that ring. It then applies
the type-grouped construction and the packed evaluator from
[simultaneous-time-space-transfer.md](simultaneous-time-space-transfer.md).
The modulus is `2D+1`, where `D` is a product of the relevant
denominators. It need not be prime. Coefficient-ring preprocessing
takes polynomial time and requires no integer factorization.

This is a conditional transfer theorem. It does not construct a
crossing family with `rho<1/sqrt(2)` or prove a new unconditional time
exponent. It strengthens the earlier
[rational time transfer](rational-height-transfer.md) under the stated
polynomial piece count and power-saving width-error hypotheses.

## 1. Normalize and enumerate the small rational factors

Normalize the scale to one by dividing one factor in each displayed
factorization by the nonzero scale. The argument in the rational time
transfer shows that this preserves `2^O(t)` bit time, ranks, and widths.
We use these normalized tables below.

For one piece and one cut, a scalar factor has two `t`-bit external
labels and at most `2t` channel bits. Thus there are at most `2^(4t)`
entries in each of the four factor tables. Polynomial piece count,
the time bound, and the corresponding output-length bound imply the
following uniform estimate. There is a fixed positive integer `C_*`
such that, for every `m>=1`, all normalized factor entries at every
size `1<=t<=m` can be enumerated, and their positive denominators
multiplied into an integer `D_m`, with

```math
\operatorname{bitlength}(D_m)\le2^{C_*m}.               \tag{2}
```

The enumeration, product construction, and all temporary rational
values use at most `2^(C_*m)` time and bits after increasing `C_*`.
Indeed, there are `2^O(m)` entries, each output has `2^O(m)` bits,
and even schoolbook multiplication of their denominator product is
still `2^O(m)` time. Denominators of zero entries may be taken to be
one. Repeated denominators need not be removed.

The enumeration includes both sides of both factorizations. Complete
piece entries can be obtained from their factor sum, so they need not
be added separately to the denominator list.

## 2. A global cutoff covers every block used by the evaluator

Choose a fixed integer `B_0>=96 C_*`. This choice uses a proved
runtime and height constant for the supplied family. It does not
require an algorithm that infers a time bound from arbitrary source
code. At a label length `lambda`, use the type-grouped blocking length

```math
r(\lambda)=\max\left(1,
 \min\left(\lambda,
 \left\lfloor\frac{\lceil\log_2(\lambda+2)\rceil}{B_0}
 \right\rfloor\right)\right).                          \tag{3}
```

There are full blocks of size `r(lambda)` and at most one shorter
nonempty tail. This differs from the previous type-grouped construction
only in the fixed constant governing its block length.

We need one denominator product and one ring covering every recursion
path. A crude envelope suffices. The separate `2lambda` width caps,
the source grouping bound with effective contraction at most one,
and rounding give

```math
V_{j-1}\le(1+3\zeta_{j-1})V_j+K_{j-1}^2\le4V_j.        \tag{4}
```

Here `zeta<=1/16` and `K_(j-1)^2<=V_j`. The top volume is at most
`1296n`. The schedule satisfies `K_j>=2^(18·2^j)` and its top
`P=K_ell` is at most `A`. The original source parameters give
`log_2 A<=n` for `n>=2`. Consequently `4^ell<=n^2/324` and
every visited volume is at most `1296n·4^ell<=4n^3`.
In particular the explicit integer

```math
\Lambda(n)=2^{16}(n+2)^3                               \tag{5}
```

bounds every visited label length, even before using the sharper
geometric contraction estimate. Inputs too small to have a scheduled
top level use the finite exceptional procedure.

Set

```math
m(n)=\max\left(1,
 \left\lceil
 \frac{\lceil\log_2(\Lambda(n)+2)\rceil}{B_0}
 \right\rceil\right).
```

Every full block and tail in (3) has length at most `m(n)`. From
(2), (5), and the choice of `B_0`, the integer `D=D_(m(n))` satisfies

```math
b(n):=\operatorname{bitlength}(D)
 =O((n+2)^{1/32}).                                     \tag{6}
```

All constants in this estimate are fixed by the supplied family.
One exact integer form of the bound is
`(2^(C_*m(n)))^32 <= 2^(32C_*+6)(n+2)`.
The deliberately loose polynomial envelope (5) avoids any dependence
on the piece choices subsequently taken by the recursion.

## 3. One odd composite modulus makes every denominator invertible

Form the odd integer and coefficient ring

```math
J=2D+1,\qquad R=\mathbb Z/J\mathbb Z.
```

For every denominator `q` in the enumeration, `q` divides `D`, and
therefore

```math
\gcd(q,J)=1,\qquad \operatorname{bitlength}(J)=O(b(n)). \tag{7}
```

Thus each denominator is a unit in `R`, even when `J` is composite.
The extended Euclidean algorithm computes its inverse in polynomial
bit time and `O(b(n))` storage. Forming `J` is linear-time arithmetic
in its bit length. Together with Section 1, all coefficient-ring
preprocessing is polynomial in `n`; no primality test or factoring
subroutine occurs.

Reduce a normalized factor entry `a/q` to `a q^(-1)` in `R`.
Every such denominator divides the product `D`, so its image is a unit.
Reduction is a ring homomorphism from `Z[1/D]` into `R`.
Consequently the exact crossing sum and both displayed factorizations
survive. The old binary channel widths remain valid upper bounds.
Some numerators, entries, or whole pieces may become zero; this does
not damage either identity and does not require changing the modulus.
We retain the displayed factorization dimensions. No matrix-rank
calculation or Gaussian elimination over `R` is needed.

The same `R` is used at every level and for every piece assignment.
Choosing a new modulus independently for individual blocks would not
justify their tensor products; the single product in Section 2 avoids
that problem.

## 4. Type grouping retains the contraction and has linear local space

Apply the type-grouped construction of the simultaneous transfer to
the reduced blocks. Its identities are finite sums and products, so
they hold in the chosen ring. Since `r(lambda)=Theta(log lambda)`,
the total width is

```math
4\rho\lambda+
 O\!\left(\lambda(\log\lambda)^{\theta-1}
       +\lambda\frac{\log\log\lambda}{\log\lambda}\right).
                                                               \tag{8}
```

Its piece-index length is `polylog(lambda)`. In particular its piece
entropy satisfies the required bound for every fixed `delta>0`.
The multinomial word count is retained in both channel dimensions,
and separate widths are capped by the complete-piece identity
factorization when necessary.

Use the explicit word-enumeration algorithm for complete-piece and
capped-factor entries. It requires `2^O(lambda)` ring operations and
`O(lambda)` bits for words, counts, channels, and integer counters.
Only one original rational block oracle is active at a time. By
Section 1 its time and workspace, including its raw output, are at
most `2^(C_*m(n))=O(n^(1/32))`.
Reduction modulo `J` and a ring operation use `O(b(n))` bits.

Every level that actually performs a grid step has

```math
\lambda\ge n^{1/4}/8
```

under the established `n^(3/4)` volume cutoff. Thus all block-oracle
and coefficient-ring scratch in (6) fits within `O(lambda)` for sufficiently
large inputs. Only a constant number of ring elements is live in
one local query. The type tuple, factorial counters, and active block
workspaces therefore fit in `O(lambda)` working bits. Their bit time
is `2^O(lambda)` after charging the polynomial cost of ring arithmetic.

This assertion is needed only at the label sizes actually visited by
the evaluator. It does not claim that a ring element with
`O(n^(1/32))` bits fits in `O(lambda)` for every unrelated small
`lambda`. No local crossing step is executed at an early leaf.

## 5. The larger scalar ring does not multiply every index length

The packed global proof now applies to (8), taking any fixed
`delta>=alpha`. Piece ordinals and binary edge-label counters retain
their previous lengths. A ring value is stored only in a scalar
accumulator or return register; it is not attached to every label bit.
Products and table queries are evaluated sequentially, with no cache
of all scalar values.

There are `O(loglog n)` active grid/query frames. Even charging a
constant number of `O(b(n))`-bit ring values per frame adds only

```math
O(b(n)\log\log(n+2))
=o\!\left(\frac n{(\log(n+2))^\alpha}\right).           \tag{9}
```

The top table reconstructs its arithmetic setting as before. Its input
weights are now in `R[epsilon]/(epsilon^2)`, whose elements occupy
`O(b(n))` bits. That changes neither its uniqueness argument nor the
`O(n/P+Af)` bound for its original arithmetic messages.

All ring operations have polynomial bit cost in `b(n)`. This adds
only polynomial factors to the established evaluator time. The
polynomial coefficient-ring preprocessing is also smaller than the target in
(1). Equations (8)--(9) therefore give the claimed weighted evaluator
time and workspace.

## 6. Recovery over the odd coefficient ring

Use binary evaluation hashes with a control field `F_(2^k)` and
`k=O(n/log(n+2))`, as in the paired-descriptor recovery theorem.
The hash field and the coefficient ring have distinct purposes.
Character phases `+1,-1` belong to `R`; the normalization
`(2^k)^(-1)` exists because `J` is odd. It is computed modulo `J`
without writing the enormous integer `2^k` into a ring register.
For a coordinatewise binary character, each one-coordinate sum is
either `1+1=2` or `1-1=0`. Their product is exactly `2^k` for a
matching bucket and zero otherwise, over this ring as well. Multiplying
by the unit `(2^k)^(-1)` proves the same exact projection identity.
No target-field orthogonality theorem is being assumed.

Hash the first factor's bits and place one square-zero tag on a bit
of either factor. Each actual divisor has a singleton bucket, giving
both virtual factors. Scalar and bit checks, streamed multiplication,
successor scans, prime filtering, and power scans are unchanged.
They use a constant number of short descriptors and accumulators.

Recovery adds `O(n/log(n+2))` to logarithmic time and workspace, plus
the smaller scalar storage in (9). Since `0<alpha<=1`, this is
absorbed by (1), including `alpha=1`. This proves every asserted output
variant.

## 7. A fixed-number-field corollary

The same theorem holds for factors in any **fixed** number field,
provided the supplied algorithms output their rational coordinates in
`2^O(t)` bit time. This is a fixed-degree extension of the interface;
it does not cover algebraic degrees growing with the label length.

Choose once and for all an integral primitive element `theta` for the
number field. Let its monic integer minimal polynomial be `f(X)`, of
fixed degree `e`, and express each value in the basis
`1,theta,...,theta^(e-1)`. Normalize the nonzero scale in the number
field **before** denominator collection. With fixed `e`, exact field
arithmetic and inversion have polynomial bit cost in the rational
coordinate heights, so normalization preserves `2^O(t)` time after
increasing its fixed constant.

Collect the positive denominators of every rational coordinate of
every normalized factor into `D`. The estimates in Sections 1--2
are unchanged after absorbing the fixed factor `e` into `C_*`.
Use the explicit finite coefficient ring

```math
R_K=(\mathbb Z/(2D+1)\mathbb Z)[X]/(f(X)).              \tag{10}
```

The coordinate substitution `theta -> X` defines a homomorphism from
`Z[theta,1/D]` into `R_K`. Every normalized crossing and factorization
identity therefore survives. Because `f` is monic, every element of
`R_K` has a unique remainder of degree less than `e`. Thus `R_K`
is a free module of rank `e` over `Z/(2D+1)`, with an explicit
coefficient-array representation. Neither primality of `2D+1` nor
irreducibility of `f` after reduction is required.

Ring operations use fixed-degree polynomial multiplication and monic
remainder reduction. They take polynomial bit time and `O(e b(n))`
bits. The unit inverses needed for rational coordinates are still
scalar inverses modulo `2D+1`. No new inverse of an arbitrary element
of the quotient ring is requested.

Finally, two remains a unit, so binary character projection works in
`R_K[epsilon]/(epsilon^2)` exactly as in Section 6. The elements zero
and one are distinct, and decoded coordinates are checked against
those two particular elements. Integer pair verification is unchanged.
Since `e` is fixed, scalar storage grows by only a constant factor.
This proves (1) for the stated number-field interface as well.

## 8. An unconditional corollary from the pinned Appendix B interface

The existing integer half-rate family already satisfies the rational
theorem's hypotheses. Thus the theorem also gives a second proof of
the unconditional simultaneous bound

```math
T(n)=\exp\!\left(O\!\left(\frac n{\sqrt{\log(n+2)}}\right)\right),
\qquad
S(n)=O\!\left(\frac n{\sqrt{\log(n+2)}}\right),        \tag{11}
```

with all the output variants stated above. This route starts from the
original time-only local interface. It does not require the separate
compressed-partition local-space algorithms.

All source references in this section use the immutable
`anthropics/formal-math` commit
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.
The definition
[`AppRunningTime.IsPieces`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean)
includes a positive integer scale `d_t`, exact crossing and factorization
identities, at most `(2t)^c` pieces, and total width
`4rho t+c log_2(2t)`. Its `size_d` and `size_f` clauses bound the
scale and factor values by `2^((2t)^c)` in absolute value. The proved
[`GeneralPieces.mayAssume`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/DefGenPieces/MayAssume.lean)
replaces any longer factorization by one with at most `2t` channel
bits. It preserves the scale, pieces, identities, and `IsPieces`
certificates.

The running-time conversion deserves an explicit check. The source
supplies `PolyT` programs for the two factors and the complete-piece
value. Their correctness is guaranteed when their padding input is
a list of at least `2^(ct)` unit symbols. Fix the program and its
constant `c` once. On a size-`t` query, construct a pad of exactly
`2^(ct)` units. After trimming, the external and channel labels use
`O(t)` bits, while a piece ordinal below `(2t)^c` uses `O(log(t+2))`
bits. The complete encoded query therefore has length

```math
L=O(2^{ct}+t+\log(t+2))=2^{O(t)}.                     \tag{12}
```

Constructing that input is included in the cost. By the source's
[`PolyT` definition](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/PolyT.lean),
one fixed finite-symbol program runs in at most `a(L+1)^a` steps for
an absolute integer `a`. Consequently every factor query, including
its explicit pad and binary integer output, takes `2^O(t)` bit time.
The fixed-stack to finite-tape simulation changes this by only a
fixed polynomial factor. The count, scale, and width programs are
already polynomial time on their unary-size inputs. Thus iterating
over the piece ordinals and factor labels also gives the uniform
`2^O(t)` enumeration required in Section 1. This argument materializes
only the logarithmic-size **local-block** pads; it does not
materialize the original global clique evaluator's pad or tables.
The encoding details are recorded in
[half-evaluator-audit.md](half-evaluator-audit.md).

The required family is supplied by the proved theorem
[`AppRunningTime.thmHalfPieces`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Main.lean).
The assembled file passes its proved `AppHalfPieces.pieceFacts` to
`GeneralPieces.thmHalfPieces`; the ring identities are not left as
assumptions. The implementation theorem in
[`GeneralPieces/Half/Main.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Half/Main.lean)
obtains the fixed parameters with `exists_par_half2` and applies
`of_par_realizes` to the named `halfF` programs. It establishes
`IsPieces (1/sqrt(2)) c ...` at every positive label length.

Apply `mayAssume`, divide one factor in each cut by the positive scale,
and take `theta=1/2`, since `log_2(2t)=O(sqrt(t))`. Every hypothesis
of the rational theorem now holds with `rho=1/sqrt(2)`, so
`alpha=log_2(1/rho)=1/2`. Equation (11) follows.

In this particular application denominator collection can be simpler
than Section 1: each normalized denominator divides the scale `d_t`.
It is enough to choose `D=product_(1<=t<=m(n)) d_t`. The source's
size bound gives
`bitlength(D)=O(sum_(t<=m(n)) t^c)=O(m(n)^(c+1))`, which is
polylogarithmic in `n`. This smaller choice is optional; the general
bound (6) already suffices.

The new rational-ring composition and workspace proof are ordinary
mathematical arguments. This source audit adds no Lean formalization
and does not report a fresh build of the pinned project.

## Verification and scope

The companion checker tests denominator enumeration, safe odd moduli,
rational reduction in both cuts, vanishing numerators, composite-ring
character projections, type-grouped crossing identities, and the exact integer cutoff
inequalities. The uniform bounds follow from the proof above, not from
those finite examples.
It also checks reduction of rational-coordinate arithmetic in
`Q(sqrt(2))` into monic quotient rings, including cases where the
polynomial splits and the modulus is composite.

```sh
python3 tools/check_rational_simultaneous_transfer.py \
  --output evidence/rational-simultaneous-transfer.json
```

The result neither assumes a suitable fixed prime for an infinite
rational family nor searches for one. It constructs a sufficient odd
ring for the finitely many logarithmic blocks used at each input length.
