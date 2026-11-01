# A rank quotient for arbitrary fixed-target Feistel relabelings

## Statement and scope

Work over a field of characteristic two. Split the `lambda` label bits into
target coordinates `A` and control coordinates `B`. The split is fixed, but
the target set can have any size. For any Boolean function
`q : F_2^B -> F_2^A`, let

```math
\pi_q(a,b)=(a+q(b),b).
```

Let `P(z)=product_i(C_i+z(E_i+C_i))` be the elementary crossing pencil,
where the factors refer to the bit positions. Apply any of the common
four-leg relabelings `pi_q` to this pencil and take the linear span of all
their coefficients. There is a linear quotient from this entire span onto
the coefficient span of the elementary pencil which preserves the crossing
and does not increase either cut rank. Arbitrarily many functions `q` may
be used, of unrestricted algebraic degree.

Consequently this class inherits the elementary entropy lower bound
`rho >= H_2(1/4) = 0.811278124459...`. This is the coefficient quotient and
entropy theorem proved in [the frame-filtration theorem](frame-filtration-barrier.md), applied after the
different filtration established below. The present theorem does not cover
compositions of nonlinear Feistel layers with different target sets.

## 1. A common filtration

For a target bit, identify its two-dimensional label space with the group
algebra `F[F_2]`. Write `a=1+g` for its augmentation generator. Its basis is
`(1,a)=(e_0,e_0+e_1)`, and `a^2=0`. For all target bits use the monomial
basis of these generators; leave the control labels in their original
coordinate basis. Filter the label space by target augmentation degree.

For a fixed control word `b`, translation by `q(b)` acts as multiplication
by `g_(q(b))=product_j(1+a_j)^(q_j(b))`. Thus the permutation operator
`Q_q` is the identity plus terms that strictly increase target degree.
The same assertion holds for its inverse, which equals itself.

On either tensor cut, the common four-leg relabeling acts by

```math
T \longmapsto (Q_q\otimes Q_q) T (Q_q\otimes Q_q)^{-1}.
```

Use the sum of the two target degrees on each cut space. Conjugation by
this operator acts trivially on the associated graded of every
degree-preserving linear operator. The untouched control coordinates
simply have degree zero.

## 2. The row-cut symbol

In the row cut, the local elementary pencil is

```math
I+zN,\qquad N=E+I,
```

where `E` swaps the two local label factors. In the augmentation basis,
both `E` and `I` are homogeneous. Hence every coefficient of the original
row pencil is degree preserving, and its associated graded operator is
itself. Every Feistel conjugate has precisely the same row symbol.

For any linear combination of coefficient tensors, let `a_i` be the sum
of the coefficients placed on degree `i`, across all the relabelings.
The associated graded row operator is exactly `sum_i a_i B_i`, where
`B_i=[z^i]P(z)` in the row cut. Its rank is at most the original rank.

## 3. The column-cut symbol

In the column cut, put

```math
R=\operatorname{vec}(I_2)\operatorname{vec}(I_2)^T.
```

The elementary local pencil is `R+z(E+R)`. Use the two-generator
augmentation basis `(1,a,b,ab)`. Since
`vec(I_2)=a+b+ab`, direct evaluation gives

```math
R(1)=a+b+ab,\quad R(a)=a+b+ab,\quad
R(b)=a+b+ab,\quad R(ab)=0.
```

Thus `R` preserves the degree filtration and its degree-zero symbol
sends `a,b` to `a+b` and kills `1,ab`. This symbol is exactly `N=E+I`.
The local column symbol at a target coordinate is consequently

```math
N+zI.
```

It is strictly equivalent, with matrices independent of `z`, to the
original local column pencil. Indeed `U=E+R` satisfies
`U^2=I` and `UR=R`, so

```math
U(R+z(E+R))=R+zI.
```

Both `R` and `N` are nonzero rank-one square-zero operators on dimension
four. Their Jordan forms are therefore identical, over the original
field. A constant similarity takes `R+zI` to `N+zI`.
Tensor these constant equivalences over the target coordinates and use
the identity on the controls. This shows that every coefficient
combination of the column symbol has exactly the corresponding original
elementary column rank.

The same coefficient sums `a_i` as in the row cut occur. The column rank
of `sum_i a_i B_i` is therefore at most the column rank of the original
mixed Feistel tensor.

## 4. Well-defined quotient and exact crossing

If a linear combination of the relabeled coefficients is the zero tensor,
its row symbol is zero. The elementary coefficients `B_i` are linearly
independent: the tensors in the independent local alphabet `{I,N}` are
linearly independent, and different `B_i` have disjoint word supports.
Therefore all `a_i` vanish. This proves that the assignment

```math
\Pi\left(\sum_{q,i} b_{q,i} B_i^{\pi_q}\right)
  =\sum_i\left(\sum_q b_{q,i}\right)B_i
```

is well defined. It is onto because the identity relabeling is available;
if it is not among the selected pencils, any one selected pencil maps
onto the same elementary span. Every pencil evaluates to `E_d` at `z=1`,
and its quotient is the original pencil, so `Pi(E_d)=E_d`.

For every decomposition inside this span, applying `Pi` piece by piece
produces an elementary decomposition with no larger rank in either cut.
The sharp elementary entropy barrier therefore applies with no bound on
the number of original pieces or relabelings.

At `lambda=7`, a piece whose quotient has nonzero coefficient of degree
three must have row rank at least `1820` and column rank at least `576`.
The resulting product `1048320` exceeds `956^2=913936`, while
`128^(sqrt(2))<956`. Such a piece exists in every crossing decomposition.
Thus every fixed-target Feistel union on 128 labels fails the strict
finite seed criterion, even before summing the costs of other pieces.

## Scope limits

The fixed target set is essential to this proof. Two rounds with
different target sets need not act trivially on any common augmentation
filtration of the form used here. No bound for that larger family is
asserted. The theorem concerns coefficient spans of the elementary
pencil, not arbitrary tensors invariant under a controlled translation.

There is an explicit failure of this rank quotient outside the stated
scope. On four bits, first replace bit zero by `x_0+x_1 x_2`, then replace
bit one by its current value plus the product of the current bit zero
and bit three. This is a two-round permutation `pi`. The tensor
`B_0+B_1+B_1^pi` has cut ranks `(254,18)`, while its proposed coefficient
quotient would be `B_0`, with ranks `(256,1)`. The row rank would increase.
Thus different-target compositions require a different argument, even
for this particular coefficient assignment.

## Replay

```sh
python tools/check_feistel_filtration.py \
  --output evidence/fixed-target-feistel-barrier.json
```

The checker verifies the local column symbol, 120 controlled translations
with arbitrary truth tables, 2,560 dense cut-rank comparisons across
two- and three-pencil mixtures, the exact 128-label bound, and the
different-target rank-quotient counterexample. These finite checks support
the filtration proof; the all-size conclusion follows from that proof.
