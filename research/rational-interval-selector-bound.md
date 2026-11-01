# Ordinary interval selectors require a large explicit rational representation

## Scoped result

The new quadratic poles give inexpensive selectors for a specified low-bit
residue class. A single explicitly written rational weight cannot, with
polynomial degree and coefficient bit lengths, select an ordinary upper
half of an exponentially long odd interval.

More precisely, let `K>=1`, and let

\[
 R(a)=A(a)/P(a),\qquad A,P\in\mathbb Z[X],
\]

have no denominator zero on the odd nodes `1<=a<2^(K+1)`. Suppose `R`
is 2-adically integral there and

\[
 R(a)\equiv1_{a\ge2^K}\pmod2.
 \tag{1}
\]

Define

\[
 D_* =\max(1,\deg P),\qquad
 H=\mathrm{bits}\left(\sum_i |P_i|(2^{K+1}-1)^i\right).
 \tag{2}
\]

Then

\[
 \boxed{2^K\le4D_*H.}
 \tag{3}
\]

In particular, the denominator degree and its coefficient bit lengths
cannot both be polynomial in `K`. The bound does not restrict the degree
of the integer numerator. Its coefficients cannot supply missing dyadic
denominators without those denominators being included in `P`.

This is a representation bound for a particular weight. It is not a
factoring lower bound, a lower bound for all uses of the inverse-block
oracle, or a reason to abandon other aggregation methods.

## 1. A unit leaf supplies a low-degree parity polynomial

Use the valuation partition from
[general rational weights](pole-cauchy-rational-weights.md). On each
regular leaf,

\[
 a=\alpha+2^h Z,\qquad
 P(a)=2^s U(Z),\qquad U_0\text{ odd},\quad2^i\mid U_i.
 \tag{4}
\]

The quotient `A(alpha+2^h Z)/U(Z)` has integral series coefficients
whose degree-`i` coefficient is divisible by `2^i`. Retain degrees
through `s` and divide by `2^s`. The resulting rational polynomial
`T(Z)` has degree at most `s` and agrees with `R(a)` modulo two for
every integer `Z` in the leaf's domain.

If `h<=K`, this domain is the complete consecutive range

\[
 0\le Z<2M,\qquad M=2^{K-h}.
\]

The desired upper-half indicator is zero on `Z=0,...,M-1` and one on
`Z=M,...,2M-1`.

If `s<M`, the first `s+1` values `T(0),...,T(s)` are even 2-adic
integers. Every Newton coefficient `Delta^j T(0)` is then even. The
finite identity

\[
 T(Z)=\sum_{j=0}^{s}\Delta^jT(0){Z\choose j}
\]

shows that `T(M)` is even as well, contradicting (1). Thus every leaf
with `h<=K` must satisfy

\[
 s\ge2^{K-h}.
 \tag{5}
\]

This argument uses 2-adic integrality, not ordinary integer-valuedness
outside the finite range.

## 2. The degree budget forces a shallow leaf

Let `h0` be the smallest residue depth of any nonempty unit leaf.
Before the first leaf, every possible odd residue branch is still
active. All such residue classes are nonempty on the full interval as
long as their depth is at most `K+1`.

If `2<=h0<=K`, there are `2^(h0-2)` active primitive nodes at depth
`h0-1`. Each uses at least one unit of the partition's mod-two degree
budget. Hence

\[
 D_*\ge2^{h_0-2}.
\]

The shallow leaf has valuation `s<H` by the explicit value-height bound.
Combining this with (5) yields

\[
 D_*H>2^{h_0-2}2^{K-h_0}=2^{K-2},
\]

which implies (3). If `h0=1`, equation (5) already gives
`H>2^(K-1)`, a stronger conclusion.

If `h0>K`, every odd residue class at depth `K` is still active.
There are `2^(K-1)` such classes, so `D_*>=2^(K-1)`. Since `H>=1`,
(3) follows in this case too.

The argument uses the active degree budget, rather than multiplying
the total number of leaves by their maximum valuation. That distinction
gives the linear product `D_*H` in (3).

## 3. Quadratic coefficients and the residue-selector distinction

For a rational function with integer coefficient pairs in the basis
`1,w`, multiply numerator and denominator by the coefficientwise
quadratic conjugate of the denominator. Its real coordinate has an
integer-polynomial numerator and the integer norm polynomial as
denominator. The norm degree is at most twice the original degree,
and its coefficient and value-height bounds have polynomial growth.
If the quadratic function satisfies (1), its real coordinate does too,
so the same obstruction applies after this explicit conversion.

There is no conflict with the cubic low-bit selector
`(2^h/(a-alpha+2^h w))^3`. Its selected nodes are a single residue class,
and its parity is constant on that entire class. The upper-half target
in (1) has opposite values at `a` and `a+2^K` for every lower-half node.
Every residue class of depth at most `K` therefore contains both values,
which is the property used in (5).

The inverse-block engine computes residue restrictions efficiently.
It does not by itself convert ordinary interval restrictions into a
single small rational weight, nor does it aggregate the product
congruence graph without additional work.

## 4. Controls and a denominator-only certificate

A leaf with `h<=K` and `s<2^(K-h)` is already a certificate that **no
integer-polynomial numerator** can turn that denominator into the
upper-half parity selector. The certificate consists of the exact
leaf identity (4) and the failed necessary inequality (5).

The [checker](../tools/check_rational_interval_selector_bound.py)
enumerates small integer-valued polynomial parity patterns in the
Newton basis and verifies that a consecutive half-step requires degree
at least its first zero-run length. It checks the active-depth budget
on randomly generated explicit denominators and records concrete
denominator-only obstruction witnesses.

It also checks the valid family

\[
 R_K(a)=\frac{\prod_{j=0}^{d-1}(a-(2j+1))}{2^{2d-1}},
 \qquad d=2^{K-1}.
\]

For odd `a=1+2X`, this is an odd integer multiple of `binom(X,d)`.
Its parity is the required upper-half indicator on the stated finite
domain. The exponentially large denominator exponent is explicitly
charged in (2), as required.

The [evidence](../evidence/rational-interval-selector-bound.json)
records the finite checks. No statement about the optimal rational
degree-height tradeoff or the complexity of arbitrary factoring
algorithms is inferred from these controls.

```sh
python3 tools/check_rational_interval_selector_bound.py \
  --output evidence/rational-interval-selector-bound.json
```
