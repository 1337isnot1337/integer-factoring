# Resonant Cauchy Smith profiles stabilize at binary input precision

## Result and scope

The [Cauchy content-slope interface](pole-cauchy-content.md) recovers a
rectangle's factor count from two padded determinant valuations. Here the
entire Smith profile stabilizes at precision controlled by the input bit
length, independently of the number of matrix entries. The statement also
works directly for rectangular matrices, without padding.

Let `N` be odd, and let `A,B` be finite sets of distinct positive integers,
with sizes `r,s`. Write `q=min(r,s)`, let

\[
 D=\#\{(a,b)\in A\times B:ab=N\},
\]

and choose `h>=1` such that `|N-ab|<2^h` for every pair. Work over the
unramified quadratic extension of the two-adic field, with valuation ring
`O`, valuation `v(2)=1`, and `w^2+w+1=0`. For every integer `k>=h`, put

\[
 z_k=N+2^kw,\qquad
 M_k(a,b)=\frac{2^kw}{z_k-ab}.
 \tag{1}
\]

Remove from `A` every node that occurs in a factor pair, and remove its
partner from `B`. Call the remaining sets `A_0,B_0`. They have sizes
`r-D,s-D`, and the rational Cauchy matrix

\[
 C_0(a,b)=\frac1{N-ab},\qquad a\in A_0,\ b\in B_0
 \tag{2}
\]

has no zero denominator. Let its `q-D` Smith orders be
`sigma_1<=...<=sigma_(q-D)`. They can be negative since (2) need not be
integral. Then the complete Smith orders of (1) are exactly

\[
 \boxed{\quad
  \underbrace{0,\ldots,0}_{D\text{ terms}},\quad
  k+\sigma_1,\ldots,k+\sigma_{q-D}.
 \quad}
 \tag{3}
\]

Every order after the first `D` is positive. Thus the profile changes by
zero on the factor part and by one everywhere else when `k` increases by
one. There is no matrix-size-dependent stabilization threshold.

This is an ordinary theorem about an implicit matrix. Its proof uses the
unknown factor pairs to identify a unit pivot block; it does not find that
block. No fast global Smith, content, or rank evaluator is supplied, and
the complete factoring score remains `1/2`.

## 1. True factor pivots have unit normalization

Every true factor pair occupies its own row and column. Order those rows
and columns first, with their partners aligned. For a nonzero
`delta=N-ab`, the assumption gives `v(delta)<=h-1<k`. Consequently

\[
 v(z_k-ab)=v(N-ab),\qquad
 v(M_k(a,b))=k-v(N-ab)\ge1.
 \tag{4}
\]

At a true pair the entry is exactly one. The leading `D` by `D` block
therefore reduces to the identity modulo two. It is invertible over `O`.
Integral block row and column operations split off this block and leave
its Schur complement. The leading block contributes `D` zero Smith
orders.

Let `U,V` be the matched node sets. The exact Cauchy Schur formula from
the preceding note gives, for unmatched `a,b`,

\[
 S_k(a,b)=
 \frac{z_k^D\prod_{u\in U}(a-u)\prod_{v\in V}(b-v)}
      {\prod_{v\in V}(z_k-av)\prod_{u\in U}(z_k-ub)}
 M_k(a,b).
 \tag{5}
\]

Each `v` is the partner `N/u` of one `u`, and both are odd. Since `a`
is unmatched, `a-u` is nonzero. Equation (4) and `N=uv` give

\[
 v(z_k-av)=v(N-av)=v\bigl(v(u-a)\bigr)=v(a-u).
 \tag{6}
\]

The analogous column identity is
`v(z_k-ub)=v(b-v)`. Also `z_k` is a unit. The factor multiplying (1) in
(5) therefore separates into one unit for each row and one unit for each
column. Multiplication by these diagonal unit matrices preserves all
Smith orders. Thus the residual Smith profile is exactly that of the
restriction of `M_k` to `A_0` by `B_0`.

This conclusion depends on eliminating a block formed from actual
factor pairs. The earlier example in which an arbitrary leading block
produces a negative-valuation Schur complement is consistent with (5):
that leading block is not a unit block.

## 2. Every residual minor has an exact affine valuation

Take equally sized subsets `R` of `A_0` and `S` of `B_0`, of size `j`.
The Cauchy determinant formula gives

\[
 \det M_k[R,S]
 =
 \frac{(2^kw)^jz_k^{j(j-1)/2}V(R)V(S)}
      {\prod_{a\in R,b\in S}(z_k-ab)}.
 \tag{7}
\]

There are no true pairs in this rectangle. All denominator valuations
in (7) therefore equal those at `z=N`, and both `z_k` and `N` are units.
It follows that

\[
 v\bigl(\det M_k[R,S]\bigr)
   =jk+v\bigl(\det C_0[R,S]\bigr).
 \tag{8}
\]

This identity holds for every minor, not merely a selected pivot chain.
Minimizing over all size-`j` minors gives the valuation of the `j`th
determinantal ideal. Successive differences of these valuations prove
that all residual Smith orders are `k+sigma_i`. Together with the unit
block, this proves (3).

The rational matrix (2) has full possible rank: each square minor is
nonzero by the same Cauchy formula, since nodes are distinct and `N` and
all denominators are nonzero. Every entry of (2) has valuation at least
`1-h`, so its smallest Smith order is at least `1-h`. The residual orders
in (3) are therefore at least one even at `k=h`.

## 3. A rectangular content slope

For an `r` by `s` matrix of full possible rank, define its maximal-minor
content valuation to be the minimum valuation among its `q` by `q`
minors. Equivalently, it is the valuation of the ideal generated by those
minors, or the sum of its `q` Smith orders. Denote this value for `M_k` by
`E_k`. Equation (3) gives

\[
 E_k=(q-D)k+\sum_{i=1}^{q-D}\sigma_i,
 \qquad
 \boxed{D=q-(E_{k+1}-E_k)}.
 \tag{9}
\]

For square matrices this is the existing determinant slope. Equation
(9) is a rectangle-native alternative to padding: only the smaller side
sets the number of Smith orders. It does not imply that minimizing over
the many maximal minors is inexpensive. If one side is empty, take
`q=D=E_k=0`, so (9) remains valid.

For endpoints and `N` of `O(n)` bits, `h=O(n)` can be read from the two
extreme product corners. This proves the small stabilization precision.
The matrix dimensions can still be exponential in `n`. An evaluator
whose cost is polynomial in the numerical dimensions would not beat the
existing factoring bound.

## 4. Exact controls

Run:

```
python3 tools/check_pole_smith_stabilization.py \
  --output /tmp/pole-smith-stabilization.json
cmp /tmp/pole-smith-stabilization.json evidence/pole-smith-stabilization.json
```

The checker performs valuation-minimizing Gaussian elimination over
`Q(w)`, independently of the Cauchy minor formula. It compares every
Smith order at three displacement precisions with (3). It includes all
selected interval starts and rectangular dimensions up to four, both
orientations of an unequal-size control, and the earlier nonintegral
arbitrary-prefix example. Exact case counts and retained witnesses are
in the evidence receipt.

These are controls for the all-size proof. Neither these finite profiles
nor the existential choice of a resonant unit block establish an
efficient factor-pivot search.
