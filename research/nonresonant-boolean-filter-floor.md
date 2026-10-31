# Nonresonant two-power filters are excluded for every Boolean pattern

## Result and scope

The complete homogeneous obstruction extends to many nonhomogeneous
Boolean filters. Let `n` be a power of two, `P` positive and odd,
`m=nP`, and select `r<=n` distinct physical identity modes. Let `g`
be any nonconstant Boolean polynomial of degree `d>=1`, in its unique
multilinear form over `F_2`. Consider

\[
f_k(x)=1+(x^{md}+x^L)g(k),\qquad L>md.
\tag{1}
\]

The complete critical pressure is strictly positive at every real
aggregate charge and every positive scale if either

\[
L\not\equiv0\pmod m\qquad\text{or}\qquad L>2md.
\tag{2}
\]

No homogeneity, read-once, or disjoint-monomial assumption on `g` is
needed. The only values of `L>md` not covered by this extension are

\[
L=m(d+s),\qquad 1\le s\le d.
\tag{3}
\]

Equation (3) identifies the possible integer-power cancellations left
open by the argument. It does not assert that those cancellations give
a negative full pressure or an improved factoring exponent.

## 1. Exact beta orders before the pressure comparison

Decompose the algebraic normal form by degree:

\[
g=\sum_{j=0}^d g_j,
\]

where `g_j` is homogeneous of degree `j`. For an occupation state `k`
of weight `w`, the complete transform gives

\[
\beta_k(x)=x^{mw}+(x^{md}+x^L)
 \sum_{j=0}^{\min(d,w)}x^{m(w-j)}g_j(k).
\tag{4}
\]

This follows by applying the exact monomial identity to every ANF
monomial. No term with `j>w` contributes.

Put `delta=L-md>0`. Every exponent in (4) is at least `mw`, so the
beta order `B_k` always satisfies `B_k>=mw`. If `g_d(k)=0`, the
coefficient at `x^(mw)` stays one and `B_k=mw`.

If `g_d(k)=1`, that coefficient cancels. The second branch in (4) has
a nonzero term at exponent `mw+delta`, and all its other exponents
are larger. The remaining first-branch exponents lie in

\[
\{mw+m,mw+2m,\ldots,mw+md\}.
\tag{5}
\]

When `delta` is not a multiple of `m`, the term at `mw+delta` cannot
cancel with (5). When `delta>md`, it lies above the entire range (5),
so either an earlier term survives or that term is the leading one.
Under either hypothesis in (2), therefore,

\[
\boxed{mw\le B_k\le mw+L-md\quad\text{for every }k.}
\tag{6}
\]

For the narrower range `0<delta<m`, the exact formula is

\[
B_k=mw+\delta\,g_d(k).
\tag{7}
\]

The zero improvement in (7) is controlled by the top homogeneous
component, whereas the infinity cost is supported on the whole
Boolean function `g`. These supports can differ; the proof below does
not identify their probabilities.

## 2. The absolute support loss cap still applies

Write `t=P log s`, and let `p` be the physical support probability
`Pr(g=1)`. The full degree of (1) is exactly `L` on that support and
zero outside it. Thus its infinity factor is

\[
1-p+p e^{2tL/P}.
\]

The order bound (6) makes every term of the zero expectation at least

\[
a=e^{-2t(L-md)/P}.
\]

Consequently the exact joint expectation identity gives

\[
\mathcal R\ge a(1-p)+p e^{2ntd}>p e^{2ntd}.
\tag{8}
\]

The strict inequality uses `0<p<1`, which holds for a nonconstant
Boolean function under a product measure with strictly positive digit
probabilities.

The [biased support lemma](homogeneous-boolean-filter-bounds.md#3-a-biased-boolean-support-lemma)
already applies to arbitrary Boolean polynomials of degree at most `d`.
Distinct physical mode energies therefore give exactly the same two
consequences as before:

\[
\log\mathcal R>-n\log(1+e^{-2nt}),
\tag{9}
\]

and `mathcal R>1` whenever `t>=t_n^*`, where

\[
e^{2nt_n^*}(1-e^{-t_n^*})=1.
\tag{10}
\]

The [complete Gaussian support-floor theorem](global-homogeneous-boolean-filter-floor.md)
proves

\[
G_n(J,t)>n\log(1+e^{-2nt})\qquad(0<t\le t_n^*)
\]

for every real charge and every admissible power-of-two length. This
absorbs (9) on the small-scale branch. Above the threshold, (10) and
the all-charge Gaussian floor suffice. The complete filtered pressure
is therefore strictly positive at every charge and scale, proving the
stated theorem.

## 3. A resonant order-cap counterexample

The upper bound (6) cannot be extended to every integer `L`. Take
`n=2`, `P=1`, `m=2`,

\[
g(k_0,k_1)=k_0+k_0k_1=k_0(1+k_1),\qquad d=2,\qquad L=6.
\tag{11}
\]

The only nonconstant `f` is at state `(1,0)`, where it equals
`1+x^4+x^6`. Direct transformation gives the beta orders, in bit-mask
order, as

\[
(B_{00},B_{10},B_{01},B_{11})=(0,2,2,8).
\tag{12}
\]

At state `(1,1)`, equation (6) would give the incorrect upper bound
`mw+L-md=6`; the actual beta polynomial is `x^8`. The first- and
second-branch terms at exponent six cancel. This example is retained
as an exact negative control. It demonstrates failure of this order
cap, without asserting a favorable full pressure.

## 4. Exact controls

The [checker](../tools/check_nonresonant_boolean_filter_floor.py)
expands literal full-state coefficient polynomials, independently of
the ANF order argument. It checks all nonconstant Boolean polynomials
through three variables and a stated deterministic sample in four
variables, at two admissible ramifications. It verifies every claimed
order interval, the exact near-power formula (7), full degrees,
endpoints, rational support-loss comparisons, and the resonant negative
control. Its [record](../evidence/nonresonant-boolean-filter-floor.json)
contains the exact finite counts and scope.

```sh
python3 tools/check_nonresonant_boolean_filter_floor.py --output /tmp/nonresonant-filter-replay.json
```

The uniform theorem follows from (4)--(10) and the established Gaussian
floor. The finite controls check its transcription and preserve the
exception that guides the subsequent resonant code-coset search.
