# The full sixteen-label CNOT span misses the strict seed target

## Result

The existing [geometric cut-budget theorem in PR #27](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md) already excludes exact finite seeds through 95 labels. This note gives a shorter independent argument for a smaller range; it does not extend that cutoff. A native finite-size exclusion does not by itself exclude asymptotic constructions obtained from tensor powers of a small polynomial source.

Let `V=span{a,b,c,d}` be the four-dimensional local coefficient space
of the [three-CNOT crossing](cnot-crossing-modules.md). Over every field
of characteristic two, every decomposition

```math
E_{16}=\sum_i T_i,\qquad T_i\in V^{\otimes4},
```

satisfies

```math
\boxed{\quad\sum_i\sqrt{r_i c_i}\ge\frac{688}{13}>16^{\sqrt2}.\quad}
\tag{1}
```

This covers arbitrary coefficient recombinations in the full tensor
space and arbitrary characteristic-two extensions. It strengthens
the earlier eight-label result in the
[full-span analysis](cnot-full-tensor-span.md). A constant conjugation
that moves the coefficient space outside `V^tensor(4)`, or a larger
tensor power, is outside this theorem.

The argument combines the existing all-field rank-one classification
with the new [mixed-cut rank budget](mixed-cut-rank-budget.md).
No classification of row-rank-two tensors is required.
This is a restriction on the finite-seed construction, not a factoring
lower bound. The achieved factoring exponent is unchanged.

## 1. Two exceptional rank-one classes

Write `n=256`. The rank-one classification supplies two classes:

- `mathcal U` is the unique column-rank-one line with row rank 256.
  It is spanned by `U=C0^tensor(4)`, whose row matrix is the identity.
- `mathcal V` consists of all row-rank-one tensors with column rank 256.

Every nonzero member of either class costs 16. Every other nonzero
tensor in the space has

```math
1/128\le r/c\le128.
\tag{2}
```

Indeed, if the smaller rank is at least two, the larger rank is at
most 256. If one rank equals one, the other is a power of four;
outside the exceptional classes it is at most 64.

We may assume that at most one piece belongs to each exceptional
class. Combine collinear `mathcal U` pieces. If two `mathcal V`
pieces occur, replace them by their sum. The sum stays in the same
coefficient space, has row rank at most two and column rank at most
256, and therefore costs at most `sqrt(512)<32`, the previous cost.
Delete a zero sum and repeat. This process cannot create a new
`mathcal U` piece, whose row rank would be 256.
It gives a decomposition of no greater cost with the claimed property.

For every remaining, nonexceptional piece, (2) implies

```math
13\sqrt{rc}\ge r+c+11\min(r,c).
\tag{3}
```

To check this, set `t=sqrt(max(r,c)/min(r,c))`. Then `1<=t<12`,
and (3) is `(t-1)(t-12)<=0` after division by the smaller rank.

## 2. Rank left after an exceptional piece

Let `P_q` be the flip on `F^q tensor F^q`. In characteristic two,

```math
\operatorname{rank}(P_q+\alpha I)\ge q(q-1)/2
\tag{4}
```

for every scalar `alpha`. For `alpha!=1` the matrix is invertible.
For `alpha=1`, every unordered pair of distinct basis indices gives
one rank-one two-by-two block, while diagonal basis indices give zero.

We also use the following elementary lemma, valid over every field:

```math
\operatorname{rank}(P_q+\mu M\otimes M)\ge q(q-1)/2
\tag{5}
```

for every square matrix `M` and scalar `mu`. Extend the field to an
algebraic closure and triangularize `M`. Conjugation by `S tensor S`
preserves the flip, so assume `M` is upper triangular.
Choose rows `(i,j)` with `i<j` and columns `(l,k)` with `k<l`,
ordered by the increasing pairs `(i,j)` and `(k,l)`. The flip restricts
to the identity. A nonzero entry `M_(i,l) M_(j,k)` requires `j<=k`,
and hence `k>i`. The added matrix is strictly upper triangular on
this selected minor. Its determinant is one. Field extension preserves
rank, which proves (5) over the original field.

The local row-rank-one conic has original row matrix

```math
xy\,b+x^2c+y^2d=(xX+yY)(xX+yY)^{\mathsf T},
```

where `X=(0,0,1,1)^T` and `Y=(0,1,0,1)^T`. Thus a tensor in
`mathcal V` has row matrix `mu zz^T`, with `z=vec(M)` and `M`+a tensor product of the two-by-two matrices

```math
\begin{pmatrix}0&y\\x&x+y\end{pmatrix}.
```

In the column convention `(x,u)|(y,v)`, its matrix is
`mu M tensor M`, whereas the crossing's column matrix is `P_q`.
Consequently (5) applies to removing any `mathcal V` piece.
The alternate convention `(x,u)|(v,y)` only permutes columns and
preserves the rank statements.

At `q=16`, (4) and (5) give the residual rank lower bound 120.

## 3. The three endpoint cases

The mixed-cut theorem gives

```math
\sum_i\min(r_i,c_i)\ge16
\tag{6}
```

for the whole crossing decomposition. We only need this weaker form
of that theorem. Sum (3) over the nonexceptional pieces.

**No exceptional pieces.** The two cut-rank sums are each at least
256, and the minimum-rank sum is at least 16. Therefore

```math
S\ge(512+11\cdot16)/13=688/13.
```

**Exactly one exceptional piece.** It costs 16 and consumes one
unit of the budget (6). If it is in `mathcal U`, the remaining
row rank is at least 120 by (4), and the remaining column rank
is at least 255 by rank subadditivity. For `mathcal V`, use (5)
and exchange the two cuts. In either case,

```math
S\ge16+(120+255+11\cdot15)/13=748/13.
```

**One piece from each class.** The two pieces cost 32 and leave
minimum-rank budget at least 14. Remove the `mathcal U` piece first:
the row rank is at least 120. Removing the row-rank-one
`mathcal V` piece reduces it by at most one, leaving row rank
at least 119. The sum of the remaining column ranks is at least
their minimum-rank sum, hence at least 14. Thus

```math
S\ge32+(119+14+11\cdot14)/13=703/13.
```

The least of these three bounds is `688/13`. The integer comparisons

```math
688^{70}>16^{99}13^{70},\qquad 99^2>2\cdot70^2
```

prove (1) without numerical logarithms or floating-point ranks.

## 4. Reproduction

```sh
python3 tools/check_cnot_sixteen_span.py \
  --output evidence/cnot-sixteen-span.json
```

The checker verifies the flip-plus-Kronecker minor lemma on all binary
matrices of sides two and three, deterministic larger binary controls,
and independent four-element-field controls. It also regenerates the
local conic and its column reshaping, checks all 81 binary row-rank-one
tensor products at four sites, and verifies the exact three-case bounds.
The all-field conclusion follows from the proof and the prior all-field
classification, not from a finite coefficient search.
