# First-residual rigidity in the rational joint row module

## Statement and scope

Consider the fixed nine-letter ambient crossing and its eight-letter
coordinate embedding from the
[complete rational row model](joint-hot-rational-retraction-module.md).
Work in the original variable w, over any characteristic-two field.
Require the output tensor to be polynomial, with its strict endpoint,
and require both naturally regularized cut matrices to retain all raw
coefficients of orders below k.

For `k = 8, 10, 12, 14`, the next residual matrix in the baseline Smith
elimination has rank one at each cut. If both next residuals still have
rank one after a rational left-retraction change, then:

- the infinity residual is exactly unchanged;
- the zero residual retains its fixed row factor;
- the possible pairs of residual changes form a linear space of dimension
  `3, 6, 15, 0`, respectively.

The dimensions concern residual images, not the entire polynomial parameter
space. Changes invisible at this jet can remain. This theorem does not
exclude higher-order profile improvements, changing earlier raw jets while
preserving their Smith orders, moving the embedding, or other tensor families.
It gives no new factoring exponent.

The result is not an extrapolation from binary enumeration. A quadratic
output-space quotient forces all anchor changes to vanish over every
characteristic-two extension field.

## 1. Complete bounded-jet model

Let A be the raw 64-by-64 matrix after the fixed coordinate restriction,
with row `(x,y)` and column `(v,u)`. Put
`h_i = 2(height(x)+height(y))` for a pair label i. The zero regularization
is

\[
M_{ij}(w)=w^{h_i-32}A_{ij}(w).
\]

The infinity regularization first reshuffles to row `(x,u)` and column
`(y,v)`, then reverses each polynomial at degree `h_(v,u)+32`.
A correction D therefore preserves all raw jets below k at both places
exactly when every nonzero entry satisfies

\[
32+k-h_i\le v_w(D_{ij}),\qquad
\deg D_{ij}\le h_j+32-k.                         \tag{1}
\]

The previous normal-form theorem writes every strict polynomial correction
row as `(w+1) H(w) S(w squared)`, with arbitrary polynomial H. Its independent
leading rows imply a predictable degree. Thus (1) bounds the coefficient
degrees and gives a complete finite linear system, not a chosen pole-order
cutoff. Odd powers of w are included.

## 2. First affected residuals

Apply ordinary DVR row and column elimination to the baseline matrices,
retaining coefficients through order k. Every pivot before order k is
unchanged by (1). Let U and V be the constant parts of the remaining row
and column transformations. The correction to the first residual is

\[
U\,[w^k]D_{\rm reg}\,V.                            \tag{2}
\]

Earlier transformations are regular and invertible. Products of two
order-k perturbations divided by an earlier pivot have order greater
than k, so they do not affect (2). The paired residual correction is
therefore linear in every allowed polynomial coefficient.

Both baseline residuals have anchor entry `(0,0)=1`, and all allowed
corrections have zero at that anchor. At the zero cut, the entire anchor
row is fixed. Consequently, requiring its residual to have rank one is
a linear condition: each row must equal its first entry times that fixed
row. Impose these equations before considering the infinity residual.

The dimensions are as follows:

| k | Paired image after zero rank-one condition | Infinity image rank | Infinity anchor variables | Zero-anchor kernel |
|---:|---:|---:|---:|---:|
| 8 | 406 | 403 | 14 | 389 |
| 10 | 257 | 251 | 12 | 239 |
| 12 | 149 | 134 | 12 | 122 |
| 14 | 57 | 57 | 10 | 47 |

All ranks are obtained by exact binary elimination. They remain unchanged
under any characteristic-two field extension.

## 3. Eliminate the nonlinear anchor motion

A rank-one matrix with anchor entry one is determined by its anchor
column u and anchor row v, as `u v transposed`. Split the infinity image
into a basis for its anchor variations and the subspace with zero anchor
row and column. Let x be the vector of free anchor coordinates.

Subtract a linear image representative with those same anchors, and
reduce the remainder modulo the zero-anchor subspace. The rank-one
condition becomes

\[
Lx+Q(x)=0,                                         \tag{3}
\]

where Q is homogeneous quadratic. Squares `x_i squared` are retained
as genuine quadratic monomials; no binary relation `x_i squared=x_i`
is imposed.

Let W be the linear span of all coefficient vectors of Q. Project (3)
onto the quotient by W. The resulting linear map has full column rank:

| k | Quadratic coefficient rank | Total output rank | Linear rank modulo quadratic coefficients |
|---:|---:|---:|---:|
| 8 | 46 | 60 | 14 |
| 10 | 36 | 48 | 12 |
| 12 | 36 | 48 | 12 |
| 14 | 25 | 35 | 10 |

Hence x=0 over every characteristic-two field. The anchor row and column
cannot change. A rank-one infinity residual with those anchors is the
baseline residual itself. The remaining paired image is the kernel of
projection onto the infinity image; subtracting the second column of
the first table from the first gives dimensions `3,6,15,0`. Conversely,
every such image leaves infinity unchanged and retains the fixed zero
row factor, so both residuals have rank one. This proves the statement.

## 4. Reproduction and controls

The [module](../tools/joint_hot_jet_rankone.py) reconstructs the complete
polynomial row space, both leading-residual maps and every quotient above.
The [checker](../tools/check_joint_hot_jet_rankone.py) repeats all exact
rank identities. It also constructs 32 literal binary polynomial tensor
perturbations and profiles both complete 64-by-64 regularized matrices
with the independent native Smith engine. All 64 controls agree with
the predicted first-residual ranks and earlier-order counts.

The [receipt](../evidence/joint-hot-jet-rankone.json) retains the linear,
quadratic and square coefficient vectors, dimensions, finite controls
and source hashes. It does not claim a full-profile optimum, an arbitrary
rational-gauge obstruction, or a faster factoring algorithm.
