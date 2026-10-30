# A rank-sum obstruction for every binary conservation kernel

## Statement and consequence

Let `F` be any field of characteristic two, let `G=F_2^lambda`, and put
`d=2^lambda`. For a scalar matrix `K:G×G→F`, define

```math
T_K(x,y,u,v)=[x+y=u+v]K(x+y,x+u).                       \tag{1}
```

The crossing is `T_I`, where `I(S,w)=[S=w]`. For any finite exact
crossing decomposition within this entire kernel space,

```math
I=\sum_i K_i,
\qquad
\boxed{\ \sum_i\sqrt{r(T_{K_i})c(T_{K_i})}\ge d^{3/2}.\ } \tag{2}
```

There is no incidence-support condition, position symmetry, degree bound,
projector condition, restriction on coefficients, or bound on piece count.
The space has dimension `d^2=4^lambda`.

For every `d>1`, (2) is strictly larger than the sufficient finite-seed
threshold `d^sqrt(2)`. Thus **no characteristic-two conservation-kernel
seed can pass that criterion at any size**. For polynomially many pieces,
the usual leading width interface requires `rho>=3/4`, giving at most
`log_2(4/3)` through this restricted construction.

This is a theorem about a specified tensor space. It is not a lower
bound on arbitrary crossing tensors or factoring algorithms. The existing
half-rate construction has a broader operator structure, and its proved
factoring exponent remains `1/2`.

## 1. Convolution ranks and Boolean degree

For `f:G→F`, write `C_f(a,b)=f(a+b)`. Permuting the rows and columns
of the two tensor cuts gives the exact block identities

```math
r(T_K)=\sum_{S\in G}\operatorname{rank}C_{K(S,\cdot)},
\qquad
c(T_K)=\sum_{w\in G}\operatorname{rank}C_{K(\cdot,w)}.   \tag{3}
```

For the first cut, fix `S=x+y=u+v` and use `x,u` as free coordinates.
For the second, fix `w=x+u=y+v` and use `x,y`. This second equality
uses the exponent-two label group. Formula (3) was also used in the
[incidence-support theorem](incidence-kernel-rank-sum-barrier.md).

Every function on `G` has a unique multilinear polynomial in the binary
coordinates, with coefficients in `F`. For a nonzero function `f`, let
`deg(f)` be its total polynomial degree. Then

```math
\operatorname{rank}C_f\ge2^{\deg(f)}.                  \tag{4}
```

To prove (4), choose a nonzero monomial of largest degree `j`, with
coordinate support `A`. Restrict the function to the coordinate subgroup
`H=F_2^A`, setting all other coordinates to zero. Its coefficient of
the full `A` monomial is nonzero. Summing the restricted function over
`H` gives exactly this coefficient: every smaller monomial appears an
even number of times.

The group algebra `F[H]` is
`F[z_1,...,z_j]/(z_1^2,...,z_j^2)`, under `z_i=1+g_i`.
An element is a unit precisely when its augmentation is nonzero.
Indeed its augmentation ideal is nilpotent, so the inverse is a finite
geometric sum after removing the nonzero constant term. The convolution
matrix of the restricted function is therefore invertible. It is a
`2^j`-by-`2^j` submatrix of `C_f`, proving (4).

The degree filtration on all functions has dimensions

```math
\dim\{f:\deg f\le j\}=\sum_{a=0}^j {\lambda\choose a}.
                                                               \tag{5}
```

Consequently every basis `f_1,...,f_d` of the function space satisfies

```math
\sum_{a=1}^d\deg f_a
 \ge\sum_{j=0}^\lambda j{\lambda\choose j}
 =\frac{\lambda d}{2}.                                \tag{6}
```

For example, sort the basis degrees. At most the dimension in (5) can
be at most `j`, so every sorted degree is at least the corresponding
sorted monomial degree. Sum these componentwise inequalities. This
argument works over arbitrary extensions of `F_2`.

## 2. Matched counts from an exact determinant expansion

We prove the stronger matrix statement: if `sum_i K_i` is invertible,
then the sum of costs in (2), defined by (3), is at least `d^(3/2)`.
Zero matrices can be omitted.

For each piece, put `q_i=rank K_i`. Choose an actual independent
column set `C_i`, of size `d×q_i`, and an actual independent row set
`R_i`, of size `q_i×d`. There is an invertible `q_i×q_i` matrix `M_i`
such that

```math
K_i=C_i M_i R_i.                                      \tag{7}
```

For completeness, let `J_i` and `L_i` be the selected column and row
indices. The intersection matrix `K_i[L_i,J_i]` is invertible: the
selected columns span all columns, while restricting to the selected
independent rows preserves the full rank. Its inverse is a choice of
`M_i`.

Concatenate all columns into `C`, stack all rows into `R`, and put
`M=diag(M_1,...,M_k)`. Then

```math
\sum_i K_i=CMR.                                       \tag{8}
```

Expand the nonzero determinant of (8) twice by Cauchy--Binet. At least
one term is nonzero. That term selects `d` actual columns from `C`
and `d` actual rows from `R`, with both selected matrices invertible.
Its middle factor is a nonzero `d×d` minor of `M`.

Because `M` is block diagonal, this middle minor can be nonzero only
when it selects the **same number `nu_i` of rows and columns in each
block**. This follows directly from the determinant expansion: every
nonzero permutation must match entries within their original block.
Thus

```math
\sum_i\nu_i=d,
```

and piece `i` supplies `nu_i` selected actual row functions and
`nu_i` selected actual column functions. The total row list is a basis
of the function space, as is the total column list. Equation (6) applies
to both lists.

The matched per-piece counts are the essential extra information.
Choosing unrelated bases separately would not justify the next step.

## 3. Sum the rank costs

Write `a_(i,j)` and `b_(i,j)` for the Boolean degrees of the selected
row and column functions from piece `i`, for `1<=j<=nu_i`.
Since these are actual distinct rows and columns of `K_i`, (3)--(4)
give

```math
r_i\ge\sum_{j=1}^{\nu_i}2^{a_{i,j}},\qquad
c_i\ge\sum_{j=1}^{\nu_i}2^{b_{i,j}}.                  \tag{9}
```

For `nu_i>0`, apply the arithmetic--geometric mean inequality to each
sum in (9):

```math
\sqrt{r_i c_i}
 \ge\nu_i\,2^{D_i/(2\nu_i)},\qquad
D_i=\sum_j(a_{i,j}+b_{i,j}).                           \tag{10}
```

The two basis-degree bounds imply `sum_i D_i>=lambda*d`.
Weighted convexity of the exponential, with weights `nu_i/d`, now gives

```math
\begin{aligned}
\sum_i\sqrt{r_i c_i}
&\ge\sum_{i:\nu_i>0}\nu_i2^{D_i/(2\nu_i)}\\
&\ge d\,2^{(\sum_iD_i)/(2d)}
 \ge d\,2^{\lambda/2}=d^{3/2}.
\end{aligned}                                        \tag{11}
```

Pieces with `nu_i=0` have nonnegative cost and may be omitted from
this lower bound. This completes the all-size proof over every field
of characteristic two within the stated kernel model.

## 4. Supplied formal limits

Suppose the individual kernels have entries in `F((t))`, while their
sum is integral and specializes to the identity. Its determinant has
nonzero constant term, so the sum is invertible over `F((t))`.
Apply the stronger matrix statement above to that field. The same
bound holds for the generic cut ranks, even when individual pieces
have Laurent-series poles.

Every piece must remain in the same fixed conservation-kernel space.
This assertion does not cover arbitrary formal crossings outside it.

## 5. A matrix extension to finite abelian p-groups

There is a broader matrix statement, with an important scope distinction.
Let `G=product_j C_(q_j)` be a finite abelian p-group, where each
`q_j=p^(e_j)>=2`, and work over any field of characteristic `p`.
For a matrix `K` indexed by `G`, define

```math
R(K)=\sum_{x\in G}\operatorname{rank}C_{K(x,\cdot)},
\qquad
C(K)=\sum_{y\in G}\operatorname{rank}C_{K(\cdot,y)}.
```

If `sum_i K_i` is invertible and `d=|G|`, then

```math
\sum_i\sqrt{R(K_i)C(K_i)}\ge d^{3/2}.                 \tag{12}
```

To see this, represent a cyclic coordinate by integers `0<=x_j<q_j`.
The functions

```math
b_a(x)=\prod_j {x_j\choose a_j},\qquad 0\le a_j<q_j,
```

form a basis: their evaluation matrix is triangular, with diagonal one.
They are periodic modulo `q_j` in characteristic `p`. Indeed
`(1+z)^(q_j)=1+z^(q_j)`, so adding `q_j` to an integer upper argument
does not change coefficients of powers below `q_j`. Pascal's identity
therefore gives the cyclic forward-difference rule
`Delta_j b_a=b_(a-e_j)`, using zero for a negative index.

Give a basis index the weight

```math
W(a)=\sum_j\frac{a_j\log_2 q_j}{q_j-1},
```

and order indices by this weight, breaking ties lexicographically.
For a nonzero function `f`, choose its largest supported index `a`.
Every iterated difference `Delta^b f` with `0<=b_j<=a_j` has leading
index `a-b`. These indices are distinct, so these differences are
linearly independent. They all lie in the span of translates of `f`.
Consequently

```math
\operatorname{rank}C_f\ge\prod_j(a_j+1)\ge2^{W(a)}.    \tag{13}
```

The last inequality follows from concavity of `log_2(x+1)` on
`0<=x<=q_j-1`, comparing it with the chord joining its endpoints.
As in (6), any basis of the function space has total leading weight
at least the sum of the weights of all `b_a`. That sum is exactly
`d log_2(d)/2`, since each coordinate index has mean `(q_j-1)/2`.
The matched-count determinant argument and weighted convexity now prove
(12) without any other change.

Equation (12) concerns the displayed **matrix convolution costs**.
The identification of those costs with both cuts of the ordinary
crossing endpoint in (1) uses `G=F_2^lambda`. For other p-groups,
including cyclic groups of order greater than two, this section makes
no new claim about the standard factoring crossing tensor.

## 6. Construction families covered and remaining directions

The local four-feature basis in
[PR26, pinned at `ae5115b3`](https://github.com/1337isnot1337/integer-factoring/blob/ae5115b3d6b4175693eb0899288b71d7b9b4629f/research/four-feature-parity-crossings.md)
is exactly the basis of kernels `K(S,w)` on one label bit. Its four
tensors are

```math
[x+y=g]\,[x+u=h]\,[x+y+u+v=0],\qquad g,h\in F_2.
```

Tensoring over positions gives all `4^lambda` independent kernel delta
functions. Thus (2) covers every characteristic-two four-feature
finite-state or matrix-product coefficient tensor, arbitrary position
dependence, and arbitrary coefficient recombinations in that space.
The prior exact constructions and their finite rank calculations remain
valid; their scalar coefficient pieces cannot give the required rank sum.

It also covers all nonlinear selectors that retain (1), unrestricted
ordinary conjugations of projector kernels, and endpoint corrections
using nilpotent kernels. No assumption about ordinary diagonalizability
or the ranks of `K_i^2` is needed.

The [bi-exterior construction](bi-exterior-crossing-algebra.md) contains
local tensors with odd total parity and is a different space. Its
256-label exclusion remains a separate theorem. Likewise, mixtures of
incompatible physical label bases and nonlinear relabelings that break
the common conservation identity are outside this argument.

## Verification and limitations

The companion exact checker tests the convolution-degree lower bound,
matched basis selection, full-rank kernel decompositions, extension-field
controls, and the actual four-tensor cuts. It includes a control showing
why the characteristic assumption is essential to the stronger matrix
statement: a Walsh matrix over `F_3` is invertible, while each row and
column is a character and has convolution rank one. Its cost is `d`.
This is a counterexample for that stronger matrix statement outside
characteristic two; it is not an odd-characteristic crossing improvement.

```sh
python tools/check_full_conservation_rank_sum.py \
  --output evidence/full-conservation-rank-sum.json
python tools/check_p_group_convolution_degree.py \
  --output evidence/p-group-convolution-degree.json
```

The second checker separately verifies the matrix extension on six
abelian p-groups in their matching characteristics, including mixed
cyclic orders. It exhausts 39,915 nonzero functions and uses exact
rational weights, with the common `log_2(p)` factor suppressed.

The proof, rather than the finite checks, establishes (2). This result
adds no Lean formalization and does not improve the factoring exponent.
