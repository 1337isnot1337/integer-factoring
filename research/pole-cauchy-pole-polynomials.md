# Batched proper rational weights and normalized inverse blocks

Date: 2026-11-01. Status: proved constructive primitive with exact checks.

## Result and scope

The inverse of the shifted Cauchy matrix admits a polynomial bit-time
contraction against two explicitly given proper rational functions with one
pole each. The coefficients may lie in the quadratic integer ring
`Z[w]`, where `w^2+w+1=0`. The algorithm assembles only two univariate
product jets. It uses the actual divisibility of every input coefficient
before choosing working precision.

This gives a substantially cheaper implementation of the full-precision
residue blocks in [the residue-block construction](pole-cauchy-residue-blocks.md).
The numerical degrees, coefficient bit lengths, requested precision, and
encoded domain sizes are all charged. The exponentially many matrix entries
are not enumerated. A block sum remains a block sum: this result does not
supply the compressed Smith-unit count needed for a factoring improvement.

Implementation: [pole_cauchy_pole_polynomials.py](../tools/pole_cauchy_pole_polynomials.py).
Checker: [check_pole_cauchy_pole_polynomials.py](../tools/check_pole_cauchy_pole_polynomials.py).
Evidence: [pole-cauchy-pole-polynomials.json](../evidence/pole-cauchy-pole-polynomials.json).

## The generating function

Use the domains and matrix convention of
[the rational-pole formula](pole-cauchy-rational-poles.md). Put `z=N+2^k w`.
For admissible poles `c,d`, its double resolvent contraction is

\[
 F(c,d)=\frac{z}{2^k w(cd-z)}(1-Q(c,d)),
 \qquad
 Q(c,d)=
 \frac{\prod_a(z-ad)\prod_b(z-bc)}
 {z^m\prod_a(a-c)\prod_b(b-d)}.
\]

The poles must avoid their respective finite node domains, and `cd-z` must
be nonzero. The existing affine valuation routines check these conditions
and partition each product into its valuation strata.

Let `q=v_2(Q)` and `v_0=v_2(cd-z)`. Choose an integer `S` larger than
`v_0` and every individual affine-factor valuation, and put `delta=2^S`.
After extracting valuations, the shifted products give integral unit jets

\[
 Q(c+\delta X,d+\delta Y)=2^q U(X)V(Y).
\]

The constant coefficient of

\[
 E(X,Y)=\frac{(c+\delta X)(d+\delta Y)-z}{2^{v_0}}
\]

is a unit. Write `e_ij=[X^iY^j]E^{-1}`. Since `E` has only constant,
linear, and bilinear terms, the rectangular table of these coefficients
is obtained with one three-predecessor recurrence per cell.

## Coefficient normalization and the correlation identity

The requested weights are

\[
 f(a)=\sum_{i=0}^{R-1}\frac{\alpha_i}{(a-c)^{i+1}},
 \qquad
 g(b)=\sum_{j=0}^{T-1}\frac{\beta_j}{(b-d)^{j+1}}.
\]

For a nonzero quadratic integer, its valuation is the minimum valuation
of its two integer coordinates. Define

\[
 A=\min_{\alpha_i\ne0}(v_2(\alpha_i)-iS),\qquad
 B=\min_{\beta_j\ne0}(v_2(\beta_j)-jS),
\]

and integral weights `W_i=alpha_i 2^{-iS-A}` and
`H_j=beta_j 2^{-jS-B}`. The desired contraction is `2^{A+B}` times the
`W,H` coefficient functional of `F(c+delta X,d+delta Y)`.
There are no factorials in this statement: the coefficient of `X^i` in
`1/(a-c-delta X)` is exactly `delta^i/(a-c)^{i+1}`.

Set

\[
 \widehat W_u=\sum_{i\ge u}W_i U_{i-u},\qquad
 \widehat H_v=\sum_{j\ge v}H_j V_{j-v}.
\]

Reordering finite sums gives the exact identity

\[
 \sum_{i,j}W_iH_j[X^iY^j]\bigl(E^{-1}UV\bigr)
 =\sum_{u,v}e_{uv}\widehat W_u\widehat H_v.
\]

Thus the product term uses two one-dimensional correlations and one
rectangular contraction. The constant term uses a second rectangular
contraction. After constructing the product jets, this costs
`O(R^2+T^2+RT)` ring operations, without constructing a separate pair
of product jets for each nonzero pair of input coefficients.

## Precision and bit cost

Put `e=min(0,q)` and

\[
 \sigma=A+B+e-k-v_0,\qquad P=p-\sigma.
\]

The answer is `2^sigma` times the integral expression

\[
 \frac zw\left(
 2^{-e}\sum_{i,j}e_{ij}W_iH_j
 -2^{q-e}\sum_{u,v}e_{uv}\widehat W_u\widehat H_v
 \right).
\]

If `P<=0`, the whole answer vanishes at absolute precision `p`.
Otherwise calculate the denominator table and coefficient weights modulo
`2^P`. Only

\[
 P_Q=\max(0,P-(q-e))
\]

bits of the product jets are needed. If `P_Q=0`, omit both jets and the
entire product contribution. Lifting their coefficients arbitrarily from
precision `P_Q` to precision `P` is harmless because the product term is
multiplied by `2^{q-e}`. In particular, a huge positive `q` never causes
the allocation of a huge integer `2^q`.

The earlier inverse-conditioning theorem gives

\[
 q\ge\min(0,k+v_0-L-u-v),
\]

where `L` bounds the loss in inverse precision and `u,v` bound the pole
denominator valuations on the two domains. Individual factor valuations,
`L,u,v`, and the encoding length of `q` are polynomially bounded in the
encoded input. Together with the explicit degree and coefficient charges,
this bounds `S,A,B,P` and all arithmetic operand lengths polynomially.
The normalized interval-jet product algorithm supplies the two jets in
polynomial bit time. This proves the claimed complexity, including the
zero and skipped-product branches.

## Full-precision residue blocks

Let `P_p(x)=sum_j gamma_j x^j` be the integral Hermite selector of degree
`2p-1`, with `j` ranging from `p` to `2p-1`. Use the quadratic poles
`c=r-2^h w` and `d=s-2^g w`. The coefficient arrays are

\[
 \alpha_{3j-1}=2^{L+3hj}\gamma_j,\qquad
 \beta_{3j-1}=2^{3gj}\gamma_j,
\]

and zero elsewhere. Each array has length `6p-3` and `p` nonzero entries.
The factor `2^L` normalizing the inverse is included before choosing
precision. The result is the sum of the entries of `2^L M^{-1}` in the
specified row and column residue classes, modulo `2^p`.

The original construction evaluated `p^2` repeated-pole contractions.
This implementation replaces those calls by one batched contraction.
It preserves the exact selector and therefore the proof of the block
oracle, including off-class error bounds.

## Exact checks and scale controls

The saved checker run includes 72 dense rational contractions with
quadratic coefficients, 24 dense inverse-block comparisons, and 24
comparisons with the original termwise block algorithm. It checks unequal
degrees, coefficient divisibility, negative output scales, skipped product
jets, the zero-answer branch, five invalid contracts, and the `N=15`
singleton controls that distinguish an inverse entry from a factor count.

Six large controls use `N=2^401+45`, two odd arithmetic progressions of
length `2^200+13`, and `k=414`. The row residues are 5 or 7 modulo 128,
and the column residue is 11 modulo 512. No matrix nodes are enumerated.
For row residue 5, precision 1 needs only two working bits and takes about
0.03 seconds on the test machine. Precision 5 needs 100 working bits and
takes about 1.46 seconds. For row residue 7, precision 5 returns
`[16,16]` modulo 32 with 103 working bits. Precision 9 returns `[80,80]`
modulo 512 with 203 working bits, in about 11.89 seconds. The residues
agree on reduction to the lower checked precisions.

These timings describe this implementation and machine. The exact
identities, recorded precisions, and nonzero large-domain residues are
the reproducible claims. The complete deterministic factoring exponent
is unchanged by this primitive alone.
