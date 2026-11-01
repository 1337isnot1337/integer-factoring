# Graded adjugates add upper-tail compatibility constraints

The adjugate mechanism extends to independent physical heights. It gives
two new families of necessary inequalities for scalar strict crossings
over **every field**, without a degree bound. They exclude the previously
retained negative five- and six-letter cell-cap profiles. A nearby integer
negative profile satisfies the stronger inequalities, so this does not
close the general five-letter route or improve the factoring score **1/2**.

This extends [bare adjugate compatibility](adjugate-seven-letter-pressure-floor.md),
using its rank-one reshuffle lemma, and applies the existing
[cell-set](cell-minor-caps.md) and [shared-prefix](shared-minor-profile-relaxation.md)
theorems to the graded adjugate. The new ingredients are the graded
cofactor bookkeeping, the resulting inequalities and their applications
to the gauged relaxation.

## 1. Statement

Use the normalization of the cell-cap note: `q>=2`, `n=q^2`, normalized
increasing complete profiles `r,s` with `r_1=s_1=0`, normalized deficit
`d>=0` and nonnegative paired energy lists `X,Y`, both with minimum zero.
Put

```math
C=(n-1)d+q\left(\sum X+\sum Y\right).
```

Let `TopCell_k` be the sum of the `k` largest of the `n` cell energies
`X_i+Y_j`. Let `TopX_l,TopY_l` be the sums of the `l` largest paired
energies, and `Xmax,Ymax` their maxima. Then

```math
(k-1)s_2+k\sum_{i<n}r_i\le kC-TopCell_k
                  \qquad(1\le k\le n),                 \tag{1}
```

and

```math
(l-1)s_2+l\sum_i r_i-\sum_{i=n-l+1}^n r_i
\le lC-\max(TopX_l+lYmax,TopY_l+lXmax)
                  \qquad(1\le l\le q).                 \tag{2}
```

Both also hold with `r,s` exchanged. Averaging the two orientations gives
the same inequalities for `gamma=(r+s)/2`. They therefore strengthen the
convex profile relaxation before any realizing tensor is known.

Integer heights suffice for the proof. Rational heights follow by
ramification; real heights follow from continuity of the determinantal
prefix functions and their differences. Auxiliary matrix-valued entries,
unequal physical alphabet sizes and endpoints obtained only after type
projection are outside the statement.

## 2. The graded cofactor identity

Write the normalized zero tensor as `R(w)`. Its strict endpoint is the
crossing, and its infinity tensor is

```math
S_{xyuv}(w)=w^{d+c_x(x)+c_y(y)+c_u(u)+c_v(v)}R_{xyuv}(w^{-1}).
```

Here `c_x+c_v=X`, `c_y+c_u=Y`; individual `c` values can have either sign.
Both tensors are polynomial. View their row cuts as ordinary operators
on joint labels, with columns `(v,u)`, so the endpoint is the identity.
Put `p_(x,y)=c_x(x)+c_y(y)` and `a_(v,u)=c_v(v)+c_u(u)`. Then

```math
T(w):=S_{\rm row}(w)=w^dD_p(w)R_{\rm row}(w^{-1})D_a(w).
```

Taking cofactors gives the exact entrywise identity

```math
\operatorname{adj}T(w)
 =w^C D_{-a}(w)\operatorname{adj}R_{\rm row}(w^{-1})D_{-p}(w).
                                                               \tag{3}
```

Indeed each cofactor uses `n-1` factors of `w^d` and all row/column
weights except the omitted two. The sum of all weights is
`q*(sum X+sum Y)`. The inverse row is an original column and the
inverse column an original row; this transpose of the grading is essential.

Let `A` be the strict crossing represented by the polynomial operator
`adj T`. Use its **crossed cut at zero** and its row cut at infinity.
Denote the zero orders by `eta`. Giving the infinity legs the negative
reversed original charges, equation (3) gives its complete infinity
profile exactly:

```math
\theta_i=\sum_j r_j-r_{n+1-i}-C.                         \tag{4}
```

The new unnormalized paired energy lists are `-X,-Y`. After subtracting
their minima they become `Xmax-X,Ymax-Y`.

## 3. Rank-one jets survive adjugation

Since the crossed cut of `T` has profile `s`, it has rank one modulo
`w^(s_2)` when `s_2>0`: a unit pivot and the two-minor relations give
an exact rank-one lift. As in the bare adjugate theorem, such an operator
is a transpose channel `M -> B M^T D`. Its adjugate is another
rank-at-most-one transpose channel, as a polynomial identity over every
field, including singular leading matrices. Consequently

```math
\eta_1\ge0,\qquad \eta_i\ge s_2\quad(i\ge2).             \tag{5}
```

For `s_2=0` the inequalities are automatic. No monomial determinant or
invertible leading operator is assumed.

Apply the shared-line theorem to the new crossing with its actual
infinity profile (4) and its unnormalized negative energies. It gives

```math
\sum_{i\le l}\eta_i+l\sum_j r_j-\sum_{i=n-l+1}^n r_i-lC
 \le-\max(TopX_l+lYmax,TopY_l+lXmax).
```

Substitute (5) to obtain (2).

For the cell-set theorem, put `a=eta_1`. The normalized deficit of the
new crossing is

```math
d'=C-\sum r+r_n-Xmax-Ymax-a.
```

Its cell cap at zero is

```math
\sum_{i\le k}(\eta_i-a)
\le kd'+k(Xmax+Ymax)-TopCell_k.
```

The `a` terms cancel. Substituting (5) gives (1). Reversing the parameter
and exchanging the middle legs proves the opposite orientation.

## 4. Prefix bounds cannot repair an excluded target

Suppose the paired energy minima are zero and the absolute profiles are
`mu,nu`, with `mu_i=r_i+mu_1`, `nu_i=s_i+nu_1` and
`d=-mu_1-nu_1`. Write `S_j` for a first-`j` sum and
`W=q*(sum X+sum Y)`. The left side of (2), with its deficit moved to
the left, becomes

```math
(l-1)S_n(\mu)+S_{n-l}(\mu)
 +(l-1)S_2(\nu)+[l(n-3)+2]S_1(\nu),                     \tag{6}
```

with right side `lW-max(TopX_l+lYmax,TopY_l+lXmax)`.
For (1) the corresponding left side is

```math
kS_{n-1}(\mu)+(k-1)S_2(\nu)+[k(n-3)+2]S_1(\nu).         \tag{7}
```

All prefix coefficients are nonnegative. Thus an excluded row-lattice
target cannot be repaired by larger actual Smith prefixes. The checker
verifies these equivalent forms exactly, including a rational scalar
translation of the two cuts.

## 5. Old targets and a surviving integer target

The exact checker excludes each old cell-cap optimizer at:

- `q=5`, `d=1`, energies `(0,0,2,8/3,3)`;
- `q=5`, `d=7/6`, energies `(0,0,19/9,23/9,37/9)`;
- `q=6`, `d=8/7`, energies `(0,0,15/8,27/10,4,4)`.

They remain negative points of the *older* relaxation, but fail the new
necessary compatibility conditions. This does not exclude nearby
profiles at the same physical energies.

In fact, set `d=6`, both paired energy lists `(0,0,12,16,18)`, and both
normalized target profiles equal to

```text
0,6,9,9,14,19,19,20,22,22,22,22,24,24,24,24,30,34,34,
36,36,37,38,40,41.
```

Every old cell/shared cap and every new graded adjugate inequality holds.
At scale `1/6`, 90-bit outward arithmetic certifies critical pressure
below **`-29/10000`** (approximately `-0.0029956091`). This is a necessary-
inequality target, **not a tensor realizing these profiles**.

The previous integer target ended in `38,40,40,42`; it violates (2),
and (6) excludes its entire proposed row-bound sector. The adaptive affine
run targeting it was stopped with an explicit incomplete checkpoint.
A registered successor uses the compatible new tail and rebuilds local
binary frames in the full degree-42 strict-endpoint space. Its
[manifest](../evidence/adaptive-q5-affine-v2-manifest.json) is registration,
not completed coverage or a contraction certificate.

## 6. Exact controls and scope

Run `make verify-graded-adjugate-caps`. The
[checker](../tools/check_graded_adjugate_caps.py) independently profiles
28 literal tensors over `F_10007` at both cut orientations and checks
864 graded inequalities. Twelve binary-alphabet cases construct full
polynomial adjugates by independently interpolated minors, verify 192
graded cofactor identities, and compute both new profiles to check
(4)--(5). An additional 300 fraction-free characteristic-two orders
check the five-letter source under six declared physical gauges.

The [packet](../evidence/graded-adjugate-caps.json.xz) retains the literal
inputs, exact profiles, cofactor controls, exclusions, new target and
source hashes. The theorem is proved by the argument above; the finite
controls test its implementation and normalization rather than replacing
the all-field argument. Substantial OpenAI assistant support was used.
No new Lean formalization or stronger factoring exponent is claimed.
