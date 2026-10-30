# A two-product approximation leaves a half-dimensional complement

## Theorem and scope

Let `V=F^d` over any field, and let `E_d` be the swap operator on
`V tensor V`. For arbitrary endomorphisms `P,Q,R,S` of `V`,

```math
\boxed{\quad
\operatorname{rank}\bigl(E_d-P\otimes Q-R\otimes S\bigr)
\ge \frac{d(d-1)}2.
\quad}                                                     \tag{1}
```

The bound is sharp for the class of sums of at most two product
operators: take `P=Q=I` and `R=0`. It holds in characteristic two,
for singular matrices, and when every member of the pencil spanned
by `P,R` is singular. No positivity or invertibility hypothesis is
needed.

For a crossing piece `A` with column-cut rank at most two, its row
matrix is a sum of at most two such products. Therefore

```math
c(A)\le2\quad\Longrightarrow\quad
r(E_d-A)\ge d(d-1)/2.                                  \tag{2}
```

This strengthens the general quotient lower bound
`max_s s(d-2s)`, which is only about `d^2/8` when the column-cut
rank is two. The result is a structural restriction on crossing
decompositions. It does not supply an improved factoring exponent.

## 1. Two maps admit common input and output flags

**Lemma.** Over an algebraically closed field, any two linear maps
between two vector spaces of the same finite dimension admit input
and output bases in which both matrices are upper triangular.

The input and output changes of basis are independent. The statement
asserts simultaneous equivalence, not simultaneous similarity.

Proceed by induction on the common dimension `d`. Dimension zero
is immediate. The homogeneous polynomial

```math
\det(\alpha P+\beta R)
```

has a zero at some projective point `[alpha:beta]`: if it vanishes
identically, any point works; otherwise a nonconstant homogeneous
polynomial of positive degree in two variables has a projective zero
over an algebraically closed field. Choose a nonzero vector `v` in
the kernel of this singular combination. Then `Pv` and `Rv` are
linearly dependent.

Choose a one-dimensional output subspace containing both images.
If both images vanish, choose any output line. The two maps descend
from the quotient by the input line `Fv` to the quotient by that
output line. These quotients both have dimension `d-1`, so the
induction hypothesis supplies their input and output bases. Lift
those bases and place the chosen lines first. Both original matrices
have zeroes below their first diagonal entry, and both lower-right
blocks are upper triangular. This proves the lemma.

This induction covers a wholly singular pencil. Neither map is
inverted, and no rectangular pencil block is discarded. At every
inductive step the two quotient spaces have equal dimension.

## 2. The basis change preserves the swap

Write the simultaneous equivalence from the lemma as

```math
P'=UPV^{-1},\qquad R'=URV^{-1},
```

where `P',R'` are upper triangular and `U,V` are invertible. The
swap identity `E_d(A tensor B)=(B tensor A)E_d` gives

```math
(U\otimes V)E_d(V^{-1}\otimes U^{-1})=E_d.             \tag{3}
```

Consequently the rank in (1) is the rank of

```math
E_d-P'\otimes Q'-R'\otimes S',\qquad
Q'=VQU^{-1},\quad S'=VSU^{-1}.                        \tag{4}
```

Thus the independent left and right basis changes needed for the
first-factor matrices are legitimate. They induce corresponding
changes in the second factors, on which no triangularity condition
is imposed.

For a field that is not algebraically closed, extend scalars first.
Matrix rank is unchanged by field extension, because it is the
largest size of a nonzero minor. The resulting lower bound therefore
also holds over the original field.

## 3. A kernel bound for triangular first factors

A slightly more general statement supplies the rank estimate. If
`T_1,...,T_m` are upper triangular and `U_1,...,U_m` are arbitrary
`d`-by-`d` matrices, then

```math
\dim\ker\left(E_d-\sum_{j=1}^m T_j\otimes U_j\right)
\le \frac{d(d+1)}2.                                   \tag{5}
```

Identify a tensor with its coefficient matrix `X=(x_ab)`. The
kernel equation, in coordinates, is

```math
x_{ba}=\sum_{j=1}^m\sum_{u\ge a}\sum_{v=1}^d
       (T_j)_{au}(U_j)_{bv}x_{uv}.                    \tag{6}
```

Project the kernel linearly onto the last row of `X`. Its image has
dimension at most `d`. If that last row is zero, use `a=d` in (6):
the right side is zero, so the last column is zero as well.

The upper-left `(d-1)`-by-`(d-1)` submatrix then satisfies the
same kernel equation for the principal submatrices of every `T_j`
and `U_j`. The first factors remain upper triangular. Restriction
is injective on the kernel of the last-row projection, since the
last row and column have vanished. Further equations may impose
additional conditions; they can only reduce its dimension. Induction
therefore bounds that kernel by `(d-1)d/2`. Rank-nullity for the
projection yields

```math
d+(d-1)d/2=d(d+1)/2,
```

proving (5). Applying it to (4) gives (1).

For sharpness, `E_d-I` is zero on the diagonal basis tensors. Each
unordered pair `{i,j}` contributes the rank-one matrix
`[[-1,1],[1,-1]]` on the span of `e_i tensor e_j` and
`e_j tensor e_i`. There are `d(d-1)/2` such pairs. This argument
also works in characteristic two.

## 4. Input to the finite rank-budget certificate

In the extremal-piece classification used by
[the residual-budget theorem](weighted-cut-residual-budget.md), let
`N=d^2`, let the chosen piece have ranks `(a,b)`, and let

```math
a_0=(k-1)^2b+1,\qquad a_1=\min(k^2b,N).
```

When `b=2`, equation (2) permits the strengthened residual row
budget

```math
R_0\longleftarrow
\max\left(N-a_1,\max_s s(d-2s),\frac{d(d-1)}2\right). \tag{7}
```

It also strengthens the ordinary total-rank surplus to at least

```math
e_2(d,k)=\max\left(0,a_0+\frac{d(d-1)}2-N\right).      \tag{8}
```

These inequalities concern the sum of the *other* pieces, since
rank subadditivity bounds their row-rank sum below by `r(E_d-A)`.
The existing chord and positive-dual arguments can use (7) and (8)
without further changes. The stronger rank-one refinement remains
available separately when `b=1`.

## 5. Three arbitrary products require a different bound

The same numerical bound fails for three products already at `d=3`.
Let `Omega=sum_i e_i tensor e_i` and set

```math
B=E_3-\Omega\Omega^{\mathsf T}.
```

Its complement has rank one. The opposite cut of
`Omega Omega^T` is the identity matrix, so the opposite cut of `B`
is `E_3-I`, of rank three. Explicitly, with
`A_ij=E_ij-E_ji`,

```math
B=-\sum_{1\le i<j\le3} A_{ij}\otimes A_{ij}.          \tag{9}
```

All identities hold in characteristic two as well. Hence an
unrestricted three-product extension of (1) would assert `1>=3`.
The triangular-first-factor lemma still applies to any number of
terms when a common input/output flag is actually available; such
flags cannot be assumed for arbitrary triples.

## 6. A common block flag limits multiscale constructions

The kernel argument has a quantitative block version. Suppose all first
factors `P_j` admit common input and output flags with matching block
dimensions `b_1,...,b_l`, where `sum b_i=d`. After independent changes
of the two bases, every `P_j` is upper block triangular. Then

```math
\operatorname{rank}\left(E_d-\sum_jP_j\otimes Q_j\right)
\ge\frac{d^2-\sum_i b_i^2}{2}
=\sum_{i<j}b_i b_j.                                   \tag{10}
```

To prove this, project the kernel onto its last `b_l` rows, an image
of dimension at most `db_l`. If those rows vanish, the last `b_l`
columns vanish by the same coordinate equation (6). Restriction
injects into the smaller problem. Induction gives kernel dimension at
most `(d^2+sum b_i^2)/2`, proving (10). The swap-preserving basis
changes are exactly (3).

This block bound is sharp for every partition and every field. Let
`B` act as swap when both input indices belong to the same block,
and as identity when they belong to distinct blocks. Its first-factor
matrix space is block diagonal. The complement `E_d-B` has one
rank-one two-dimensional block for each unordered pair of indices in
distinct blocks, giving equality in (10).

Apply this to the abstract dyadic profile from
[the scalar-relaxation limit](weighted-cut-relaxation-limit.md).
Let `d=2^m`, and call the pairs `(d2^j,d/2^j)` the forward pieces;
the reversed pairs are the backward pieces. If an actual crossing
decomposition had these cut ranks, and `B` were the sum of its
forward pieces, then rank subadditivity for the backward pieces would
give

```math
r(E_d-B)\le\sum_{j=0}^m d/2^j=2d-1.                  \tag{11}
```

Thus all product first factors of the forward pieces cannot share
a complete input/output flag when `d>=8`: (5) would contradict
`d(d-1)/2>2d-1`. In particular, a construction from jointly
triangularizable first factors, including commuting projections,
cannot realize this profile.

More generally a common block flag would have to satisfy

```math
\sum_i b_i^2\ge d^2-4d+2.                            \tag{12}
```

If its largest block has size `b`, then `sum b_i^2<=bd`, so
`b>=d-3`. For `d>=9`, write `d-b=s<=3` and use
`sum b_i^2<=b^2+s^2`. Equation (12) implies `bs<=2d-1`;
the case `s=3` would give `3d-9<=2d-1`, impossible for `d>8`.
Consequently `b>=d-2` for dyadic `d>=16`.

This is a restriction on compatible first-factor subspaces that the
scalar profile alone does not capture. Any realization of that
profile would require first-factor structure that does not split
under a common flag into substantially smaller blocks.

Small-rank perturbations of a flag construction have a quantitative
cost as well. If `B=T+C`, where `T` has the block-flag structure in
(10), rank subadditivity applied to `E_d-T=(E_d-B)+C` gives

```math
\operatorname{rank} C\ge
\frac{d^2-\sum_i b_i^2}{2}-\operatorname{rank}(E_d-B). \tag{13}
```

For a complete flag and the dyadic forward sum, this is at least
`d(d-1)/2-(2d-1)`. Thus a perturbation with row rank of lower order
than `d^2` cannot repair that common-flag construction.

## Verification and provenance

The companion exact checker tests small exhaustive examples, random
examples over several prime fields, swap-preserving basis changes,
and equality for `E_d-I`. It also checks the three-product
counterexample (9), block-flag inequalities and their sharpness,
and the dyadic-profile obstruction. These finite checks support the algebraic
interfaces; Sections 1–3 prove the theorem for all dimensions and
all fields. No numerical optimizer or Lean proof is required.

```sh
python3 tools/check_two_product_complement.py \
  --output evidence/two-product-complement.json
```

The proof was developed jointly in this research session: the
triangular-first-factor kernel bound and the algebraic two-map flag
induction combine to remove the singular-pencil restriction. The
crossing cut convention and quotient comparison are those of
[weighted-cut-rank-budget.md](weighted-cut-rank-budget.md). The work
used substantial OpenAI assistant support.
