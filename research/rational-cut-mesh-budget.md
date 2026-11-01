# Rational cut meshes extend the exact finite-seed cutoffs

## Result

Every exact crossing decomposition satisfies

\[
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}
\]

for **all fields when `2<=d<=96`**, and for **all characteristic-two fields
when `2<=d<=97`**. The preceding PR27 all-field cutoff was 95. This
successor uses rational cut thresholds to prove the additional all-field
case at 96, then combines those thresholds with binary interior Jacobians
for the additional characteristic-two case at 97.

These are exclusions for completed finite crossing seeds. The achieved
factoring exponent remains `1/2`. A polynomial source on a smaller label
space can still produce a larger useful crossing after taking tensor
powers and grouping coefficients; that possibility is not excluded here.
The separate formal-limit cutoff is unchanged.

The proof reuses the exact subset budgets from the pinned
[PR27 geometric cut theorem](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md)
and the binary Jacobian certificate lemma in
[the interior-Jacobian note](interior-jacobian-budget.md).

## 1. A rational scalar interpolant

Fix an integer `k>=2`, a positive integer `H`, and put `h=1/H`. Let

\[
\mathcal S=\{1,1+h,\ldots,k-h\}.
\]

Define positive rational coefficients

\[
a=\frac1{2k-h},\qquad
w_1=\frac{h}{2(2+h)},\qquad
w_s=\frac{2h}{(2s-h)(2s+h)}\quad(s>1).
\]

For ranks whose ratio lies in `[1/k^2,k^2]`,

\[
\sqrt{rc}\ge a(r+c)
 +\sum_{s\in\mathcal S}w_s
   \bigl(\min(r,s^2c)+\min(s^2r,c)\bigr).
\tag{1}
\]

To prove (1), divide by `c` and put `t=r/c`. The right side is a
piecewise-linear function `L(t)` satisfying `L(t)=tL(1/t)`. On
`[1,k^2]`, its knots are the squares of consecutive grid points. The
slope between `s^2` and `(s+h)^2` is `1/(2s+h)`, which is exactly the
slope of the corresponding chord of `sqrt(t)`. The differences of
adjacent slopes are the positive weights `w_s`. Also
`a+sum_s w_s=1/2`, so `L(1)=1`. Therefore `L` is the chord interpolant
of `sqrt(t)` on this interval and lies below that concave function.
Reciprocal symmetry proves the other half of the interval. All
coefficients and knot identities are rational.

The integer-threshold interpolant from PR27 is exactly `H=1`. The new
finite certificates use `H=4`.

## 2. Conditioning with a fractional threshold

Choose the same maximal-ratio piece and integer class `(k,b)` as in
PR27. Its ranks `(a_piece,b)` satisfy

\[
a_0=(k-1)^2b+1\le a_{\rm piece}\le a_1=\min(k^2b,d^2).
\]

Let `F(m)` be a certified complement table, and define

\[
G(s;b)=\min_{m\ge b}\{F(m)+s^2m\}.
\]

The minimum stops at `max(b,2d-1)` for the tables used here. Let `R,C,e`
be the preceding exact residual row budget, residual column budget,
and extra sum-of-cuts budget. These are unchanged from the
interior-Jacobian note. For a rational threshold `1<=s<k`, put

\[
\mu_s=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_s=\frac{k^2-s^2}{k^4-1}.
\]

The selected piece is guaranteed to lie in the first minimizing subset
only when `a_0>s^2b`. Accordingly define

\[
b_s=\begin{cases}b,&a_0>s^2b,\\0,&a_0\le s^2b.\end{cases}
\]

The valid conditioned bounds are

\[
U_s=\max\{G(s;b_s),\min(a_0,s^2b)+\mu_sR+\nu_sC\},
\]

\[
V_s=\max\{G(s;0),b+\nu_sR+\mu_sC\}.
\tag{2}
\]

The chord bounds on the residual pieces hold for every positive real
threshold, so no integrality premise is changed. The distinction
between `b_s=b` and `b_s=0` is necessary for thresholds in the upper
part of a ratio class; the checker applies it exactly.

Combining (1), (2), and the ordinary rank-sum budget gives

\[
S\ge a(2d^2+e)+\sum_{s\in\mathcal S}w_s(U_s+V_s).
\tag{3}
\]

For rational `s`, every term of (3) is rational. The checker computes
`G(j/H;b)` using the integer numerators `H^2 F(m)+j^2m`, avoiding any
floating-point comparison.

## 3. The all-field case at 96 labels

Use only the existing all-field geometric complement table. Among
15,014 ratio classes, 14,471 satisfy the desired common floor from the
nonnegative base term alone. Another 540 satisfy it using the existing
integer threshold formula. The remaining three receive a complete
quarter-grid evaluation of (3). Every class has cost at least `636`.

The smallest value among the three refined classes occurs at
`(k,b)=(52,3)` and is

\[
\frac{4390766431075362644967128691859121}
     {6903470448235203599488214733030}
=636.023064630893\ldots.
\]

The exact integer comparison

\[
636^{408}>96^{577},\qquad577^2-2\cdot408^2=1
\]

proves `636>96^sqrt(2)`. Together with PR27, this gives the all-field
cutoff through 96. No new field-dependent Jacobian is used for this case.

## 4. Binary interior certificates at 97 labels

The following 13 instances of the Jacobian
`(X,Y) -> (D_iX-YA_i)_i` have full row rank over `F_2`:

| Maps `m` | `s` | `q` | Rows | Columns | Old complement bound | New bound |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 28 | 37 | 4,144 | 4,152 | 1,024 | 1,036 |
| 5 | 23 | 33 | 3,795 | 3,814 | 756 | 759 |
| 6 | 21 | 28 | 3,528 | 3,528 | 576 | 588 |
| 7 | 13 | 36 | 3,276 | 3,288 | 462 | 468 |
| 8 | 10 | 39 | 3,120 | 3,132 | 387 | 390 |
| 9 | 9 | 37 | 2,997 | 3,012 | 328 | 333 |
| 10 | 7 | 42 | 2,940 | 2,940 | 287 | 294 |
| 11 | 7 | 37 | 2,849 | 2,850 | 258 | 259 |
| 12 | 6 | 39 | 2,808 | 2,808 | 230 | 234 |
| 14 | 4 | 48 | 2,688 | 2,724 | 188 | 192 |
| 17 | 3 | 51 | 2,601 | 2,628 | 150 | 153 |
| 18 | 3 | 48 | 2,592 | 2,634 | 141 | 144 |
| 22 | 3 | 38 | 2,508 | 2,524 | 112 | 114 |

The evidence stores every matrix entry as a binary mask and identifies
a full-row nonsingular minor. The chart and proper-projection proof
therefore supplies these complement bounds over every field of
characteristic two, after scalar extension. Dimension feasibility alone
is not used in the theorem.

Using this strengthened table, all 15,340 ratio classes at `d=97` have
cost at least `6453/10`. The base term settles 14,825 classes; the integer
mesh settles 512; and the quarter mesh settles the last three.
The smallest value among those three, again at `(k,b)=(52,3)`, is

\[
\frac{21372032485470851400403906335263}
     {33118065183586580154300445830}
=645.328534955082\ldots.
\]

The exact comparison

\[
6453^{408}>97^{577}\,10^{408}
\]

proves `6453/10>97^sqrt(2)`. Thus the characteristic-two finite-seed
cutoff extends through 97.

## 5. Exact replay

```sh
python tools/check_rational_cut_mesh_budget.py \
  --output evidence/rational-cut-mesh-budget.json
```

Keep `check_interior_jacobian_budget.py` beside the new checker. The replay
uses only the Python standard library and takes about ten seconds in the
recorded workspace. It verifies 434 exact reciprocal-knot interpolation
identities, checks all 13 literal Jacobians and direct evaluations of
their coordinate maps, covers every scalar class at both dimensions,
and proves both final comparisons with integers. It also checks that
`H=1` reproduces the prior formula on representative classes.

Exploratory floating-point scans selected these two exact cases. They
are not part of the proof. The claims above are supported by complete
finite rational coverage and the literal binary rank certificates.
