# Formal crossing limits are excluded through sixty-four labels

The [factor-kernel refinement](rank-one-formal-cutoff.md) extends this
formal cutoff to seventy-one labels. The present harmonic certificate
is retained as its independently replayed base case.

## Theorem

Let `F` be any field, let `d>=2`, and let a four-leg tensor
`T(t)` have entries in `F[[t]]`, with `T(0)=E_d`. Suppose

```math
T(t)=\sum_i A_i(t),
\qquad A_i(t)\in\operatorname{Mat}_{d^2}(F((t))).     (1)
```

The individual tensors may have poles. Let `r_i,c_i` be their two
generic cut ranks over `F((t))`. For every `2<=d<=64`,

```math
\boxed{\quad\sum_i\sqrt{r_i c_i}>d^{\sqrt2}.\quad}   (2)
```

Thus a supplied one-parameter formal degeneration passing the strict
rank-sum seed criterion needs at least **sixty-five labels**. This
extends the older formal-degeneration exclusion through forty-one
labels. It is separate from the stronger finite exclusion for exact
crossing decompositions.

The theorem uses the moving-subspace result in
[`geometric-subspace-compression.md`](geometric-subspace-compression.md).
It does not use the exact-only half-rank complement or the stable-rank
swap/product refinement. No better factoring exponent or new Lean
formalization is asserted.

## 1. A complement function valid with poles

For `m>=2`, let `F_geo(d,m)` be the maximum `s q` over the integer
certificates

```math
2\le a\le2m-2,\qquad b=2m-a,
\qquad 2s+aq\le2d,\qquad bs+2q\le2d.                (3)
```

The geometric theorem proves that every `m` first-factor matrices
over an algebraic closure admit an `s`-dimensional input subspace
whose joint images lie in a codimension-`q` output subspace.
Its saturated-lattice proof permits both subspaces to move with `t`.

Define

```math
F_{\rm form}(d,0)=d^2,\qquad
F_{\rm form}(d,1)=\lfloor d^2/4\rfloor,\qquad
F_{\rm form}(d,m)=F_{\rm geo}(d,m)\quad(m\ge2).       (4)
```

The one-product value follows by padding with a zero second product:
the two-map geometric certificate permits `s+q<=d`, whose largest
product is `floor(d^2/4)`.

If `T=A+B` and `B` has opposite-cut rank at most `m`, factor it
over `F((t))` as a sum of `m` products and apply the moving-subspace
projection. All poles of `B` are killed before reduction, giving

```math
r(A)\ge F_{\rm form}(d,m).                          (5)
```

For `m=0`, the tensor `A=T` has full generic row rank `d^2` since
its specialization is the invertible crossing. The same argument
also gives full generic rank for the other cut of `T`.

The distinction at small `m` is deliberate. Equation (4) uses
`floor(d^2/4)` for one or two products, rather than the larger
`d(d-1)/2` theorem available for exact complements of `E_d`.

## 2. Weighted budgets for a formal decomposition

Fix an integer `s` with `0<=s<=d`. For each piece choose the smaller
of `r_i` and `s^2c_i`. Let

```math
J=\{i:r_i>s^2c_i\},\qquad m=\sum_{i\in J}c_i.
```

The sum `B=sum_(i in J) A_i` has opposite-cut rank at most `m`.
Applying (5) to its complement and using rank subadditivity gives

```math
\sum_i\min(r_i,s^2c_i)
 =\sum_{i\notin J}r_i+s^2m
 \ge F_{\rm form}(d,m)+s^2m.                        (6)
```

Consequently, with

```math
G_d(s)=\min_{0\le m\le2d-1}
          \bigl(F_{\rm form}(d,m)+s^2m\bigr),        (7)
```

both weighted budgets hold:

```math
\sum_i\min(r_i,s^2c_i)\ge G_d(s),\qquad
\sum_i\min(s^2r_i,c_i)\ge G_d(s).                   (8)
```

The finite range in (7) is complete. The geometric complement
function vanishes for `m>=2d-1`, and `s^2m` is nondecreasing there.
The cutoffs `d+1` or `d` from weaker complement functions must not
be reused. Ordinary rank subadditivity separately gives

```math
\sum_i r_i\ge d^2,\qquad \sum_i c_i\ge d^2.          (9)
```

Every nonzero piece has `1<=r_i,c_i<=d^2`; zero pieces may be
discarded.

## 3. A rational lower interpolant

For `k>=2`, put

```math
w_1=\frac16,\qquad
w_s=\frac{2}{(2s-1)(2s+1)}\quad(s\ge2),
```

and define

```math
Q_k(r,c)=\frac{r+c}{2k-1}
 +\sum_{s=1}^{k-1}w_s
   \bigl(\min(r,s^2c)+\min(s^2r,c)\bigr).            (10)
```

As proved in
[`weighted-cut-rank-budget.md`](weighted-cut-rank-budget.md),
this positive homogeneous function satisfies

```math
Q_k(r,c)\le\sqrt{rc}
\quad\text{for }k^{-2}\le r/c\le k^2.                (11)
```

The proof is a finite interpolation argument: after dividing by
`c`, the expression is affine between consecutive square and
reciprocal-square knots, equals the square root at those knots,
and stays below the concave square-root graph between them.
The checker verifies the knot equalities exactly.

Apply (11) with `k=d`, sum over all pieces, and use (8)--(9):

```math
\sum_i\sqrt{r_i c_i}\ge
H_d:=\frac{2d^2}{2d-1}
       +2\sum_{s=1}^{d-1}w_sG_d(s).                 (12)
```

The right side is a rational number computed from finite integer
maxima and minima.

## 4. Exact finite certificate

The checker evaluates (3)--(12) for every `2<=d<=64`. For each
rational result `H_d`, it verifies the integer inequality

```math
H_d^{408}>d^{577},\qquad577^2-2\cdot408^2=1.          (13)
```

Thus `H_d>d^(577/408)>d^sqrt(2)` throughout the stated range.
At the final dimension,

```math
H_{64}=\frac{2347092541}{6538455}
       =358.967453473335\ldots,
\qquad64^{\sqrt2}=358.363867829704\ldots.             (14)
```

At dimension sixty-five the same scalar certificate is insufficient:

```math
H_{65}=\frac{259385858696}{708864255}
       =365.917531976555\ldots<65^{\sqrt2}.           (15)
```

The latter comparison is also checked exactly, using
`1393^2-2*985^2=-1` and the lower exponent `1393/985`.
It is a limit of this certificate, not a realizable sixty-five-label
border seed.

## 5. The class-conditioned refinement remains valid

For future refinements, choose an actual piece of maximum rank ratio,
exchange the cuts globally if necessary, and place it in the class

```math
a_0=(k-1)^2b+1\le a\le a_1=\min(k^2b,d^2).
```

The remaining ordinary rank totals satisfy

```math
R\ge\max(d^2-a_1,F_{\rm form}(d,b)),\qquad
C\ge d^2-b.                                        (16)
```

The total row-rank surplus is at least
`max(0,a_0+F_form(d,b)-d^2)`. No exact swap/product or Fitting
refinement is used when `b=1`.

For every `s<k`, the selected piece belongs strictly to the
minimizing subset in (6), because `a>s^2b`. Hence its product
count is at least `b`, and one may replace (7) on that side by

```math
G_d(s;b)=\min_{b\le m\le\max(b,2d-1)}
          \bigl(F_{\rm form}(d,m)+s^2m\bigr).        (17)
```

The same residual chord inequalities as for an exact decomposition
apply to the numerical rank ratios. These statements are recorded
to prevent future refinements from importing an exact-only
small-count assumption. They are not needed for the exclusion
through sixty-four labels.

## 6. Replay and scope

```sh
python tools/check_geometric_border_cutoff.py \
  --output evidence/geometric-border-cutoff.json
```

The standard-library replay checks every claimed dimension, exact
interpolation knots, weighted-minimum closure, and moving-input
and moving-output Laurent-polynomial controls. Those controls use
actual poles of orders one and three, and deliberately include
unsaturated generators whose reductions lose rank. The projected
matrices stay integral because the entire polar part is annihilated.

The theorem concerns the supplied formal family (1). It does not
silently identify an arbitrary abstract closure notion with such a
family, and it does not apply the larger exact decomposition cutoff
to formal limits.

The construction, proof development, and replay used substantial OpenAI assistant support.
