# Weighted divisor-prefix tables over a fixed characteristic

## Result and scope

The multiplication network used in the [half-rate composition](half-rate-factoring.md)
can preserve the identity of each divisor, instead of only its count. This gives
a useful input interface for a clique evaluator in a fixed characteristic.
The construction below is unconditional. A faster modular evaluator is a
separate obligation; this note establishes no factoring exponent above `1/2`.

Fix a prime `p` and an integer `D>=1`. Use the commutative algebra

```math
R_D=\mathbb F_p[Z]/(Z^D-1).
```

Its elements are represented by their `D` coefficients. Let `N>=2` have
binary length `n`, let `P>=2` satisfy `P^40<=n`, and fix a valid divisor-prefix
query `0<=k<=n`, `0<=x<2^(n-k)`. There is one deterministic algorithm that
writes a `P`-vertex clique over `R_D` with common label width `lambda` and
value

```math
\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}Z^{d\bmod D},
\qquad P^2\lambda\le4717n.
\tag{1}
```

It uses the same grouping, edge labels, codec, and setting enumeration as
the existing [top constructor](half-top-constructor.md). Its complete bit
cost and output length are bounded by

```math
2^{O(n/P)}\operatorname{poly}(n,D,\log p).
\tag{2}
```

All polynomial degrees and exponential constants are absolute. In the
fixed-characteristic application, `p` is fixed once for the algorithm.
No finite-field extension, irreducible polynomial, or primitive root is
needed to write these tables. The algebra may have zero divisors; its
coefficient representation and multiplication still work.

The source reduction is Isaac Rudich and Louis-Martin Rousseau's
[factoring repository at c91cc3f](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
This is an ordinary extension of its explicit setting correspondence,
not a new Lean formalization.

## 1. Attach the exponent locally

Retain the source parameters

```math
S=2^s,\qquad A=2S,\qquad h=S,\qquad w=2S^2,
\qquad n\le4S^4=w^2\le16n.
```

Here `w` is the number of bits per input digit. The source names it
`pL s`; it is unrelated to the clique label width. A `faP a` input gate
holds the digits `u_i` at positions `a+A i`, for `0<=i<h`. Define

```math
E_a(u)=\sum_{i=0}^{h-1}u_i\,2^{w(a+A i)}.
\tag{3}
```

Replace only that gate's scalar weight by

```math
\widetilde W(\mathrm{faP}\ a,u)
=W(\mathrm{faP}\ a,u)\,Z^{E_a(u)\bmod D}.
\tag{4}
```

Every other gate retains its original `0/1` weight. All settings, message
functions, message equalities, and edge widths remain the source definitions.
The integer `a+A i` visits each digit position `0,...,w-1` exactly once as
`a` and `i` range over their domains, since `Ah=w`. Consequently

```math
\prod_{a=0}^{A-1}Z^{E_a(u^{(a)})}
=Z^{\sum_{a,i}u_i^{(a)}2^{w(a+A i)}}=Z^d,
\tag{5}
```

where `d` is the first input integer reconstructed from these digits.
The equality holds in `R_D`, with all exponents reduced modulo `D`.

The formula does not require writing an exponentially large monomial.
For a local setting, compute

```math
E_a(u)\bmod D
=\sum_i u_i\,\bigl(2^{w(a+A i)}\bmod D\bigr)\bmod D
```

using binary modular exponentiation. The exponent position is at most
`w^2<=16n`; each digit has at most `w` bits. This computation has cost
polynomial in `n` and `log(D+1)`.

The source digit indexing is `digN_slot`, `digN_dig`, `numOf_dig`, and
`dig_numOf` in
[MulCount.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean).
The gate weights and message definitions are in
[MulNet.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean).

## 2. The value is the divisor polynomial, with coefficient one

Expand the original edge-equality contraction as a sum over all original
gate settings. A setting contributes zero unless it passes the old
`0/1` vertex tests and every message equality. The new monomial factor
cannot turn a previously rejected setting into a contributing one.

For an accepted setting, the source's `good_pinned` theorem determines
all intermediate and carry gates uniquely from the two input digit arrays.
Conversely, `asm_agree` and `good_asm_iff` show that the forced setting
accepts exactly when the reconstructed integers multiply to `N` and the
first one has the requested prefix. The digit maps are inverse bijections.
For every divisor `d` in the query, the second integer is uniquely `N/d`.
There is therefore exactly one accepted setting for each such divisor.

Its new weight is `Z^d` by (5). Summing proves (1). This argument uses
the actual setting bijection, not just the equality of two unweighted
cardinalities. Before exponents are reduced modulo `D`, each divisor
monomial has coefficient exactly one, which is nonzero in every
characteristic. Distinct exponents can collide after reduction; the
collision is retained in the coefficient vector.

To ask whether a prefix contains a divisor greater than one, subtract

```math
[\lfloor1/2^k\rfloor=x]Z
\tag{6}
```

from the contraction, reducing `Z` modulo `Z^D-1`. This removes the
unique divisor-one term. In characteristic two the subtraction is
addition, but the result is still the desired polynomial because its
coefficient was exactly one.

A single ordinary modular count is insufficient for this query. For
example, the nontrivial divisor polynomial of `9` is `Z^3+Z^9` over
`F_2`. It vanishes after setting `Z=1`, and also modulo `Z^3-1`, despite
having nonempty support. Its remainder modulo `Z^5-1` is nonzero.
The [fixed-characteristic transfer](fixed-characteristic-transfer.md)
uses several cyclic moduli to certify that a nonzero sparse polynomial
cannot vanish in all of them.

## 3. Grouped tables and their complete cost

Retain the source grouping chosen in the existing top-constructor proof.
For one prescribed group boundary, enumerate the same local original
settings as before. Reject a setting using the same internal equalities,
old `0/1` vertex predicates, external messages, and zero padding.
For an accepted local setting, add the monomial whose exponent is the
sum of (3) over the `faP` gates in the group.

Thus every table entry is a vector of `D` residues modulo `p`.
It can be initialized to zero and updated by adding one to the indicated
coefficient. The settings need not be retained between updates.
The finite-alphabet tape cost of accessing a coefficient is bounded by
a polynomial in `D` and `log p`; no unit-cost random access is assumed.

The existing proof gives at most `2^H` local settings per group, with

```math
PH\le6768n,\qquad P^2\lambda\le4717n.
```

It also pays for the deterministic grouping search, source arithmetic,
edge ordering, and all `P*2^((P-1)lambda)` table locations. Multiplying
these costs by the coefficient-vector overhead proves (2). In particular,
no enumeration of `p^D` possible algebra elements occurs.

For a group with no original gates, the all-zero boundary entry is the
multiplicative identity of `R_D`, and every other entry is zero. Appended
label bits are still forced to zero. The original off-diagonal codec
merely reorders table locations; it acts separately on the coefficient
vector at each location. Its value-preservation argument therefore holds
over this commutative algebra without modification.

This is a direct construction of tables over `R_D`. It is not a claim
that the existing integer-valued `PInst` type or its compiled program
already accepts coefficient vectors.

## 4. Required evaluator interface

The transfer needs an arithmetic evaluator whose operation count is
independent of the algebra dimension, apart from an explicitly charged
polynomial representation cost. For example, a sufficient form is

```math
2^{A(P,\lambda)}\operatorname{poly}(n,P,\lambda,D)
```

bit operations over the fixed characteristic. Such an evaluator may
perform addition and multiplication in `R_D` and multiply by inverses
of nonzero constants from `F_p`. It must not assume arbitrary nonzero
elements of `R_D` are invertible. A base-field identity involving only
these operations remains valid in the algebra.

Addition of coefficient vectors costs `O(D)` operations in `F_p`.
Schoolbook cyclic convolution costs `O(D^2)` operations. These elementary
algorithms suffice for a polynomial simulation; no multiplication
speedup or irreducible-polynomial construction is needed.

There is a material cost distinction from the old integer evaluator.
If `D=2^{O(n/log n)}`, putting the coefficient-vector length into its
integer bit-bound parameter `B` would create a term of size
`P^2 log B=O(P^2 n/log n)` in its sufficient pad. That term need not fit
the desired factoring time. A modular transfer must establish its
arithmetic-operation bound directly and then charge the polynomial in
`D`; it cannot obtain that bound by this substitution into `PartB`.

If a modular crossing family is defined only over fields, an additional
proof is required before using it over `R_D`. For identities and
factor entries built solely from base-field constants, sums, products,
and base-field constant inverses, the extension is immediate. The
uniform recursion and all of its costs still require their own proof.

## 5. Finite verification

Run

```sh
python3 tools/check_weighted_prefix.py
```

The checker independently compares the source's three-stage Fourier
multiplication formulas with ordinary digit convolution on all 256
small input pairs. It also checks 64 pairs at the next scale, whose
source modulus is the composite integer `2^32+1`, and 40 distributed
digit-weight cases at that scale. It checks carry acceptance and weighted prefix sums
for characteristics `2,3,5` and cyclic dimensions `1,2,3,5,7`.

It also compares direct weighted setting enumeration with contractions
of grouped cyclic-polynomial tables over all 256 binary endpoint-message
maps of a three-gate path, with both group orientations. Additional
checks cover a nontrivial internal equality, an empty group, zero
padding, and the cancellation example above.

These are finite checks of the weighted constructor and codec. They
do not implement a full modular clique evaluator, prove a new exponent,
or replace the all-input setting and cost arguments in this note.
