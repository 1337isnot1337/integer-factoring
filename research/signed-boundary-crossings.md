# Signed boundary intertwiners and isotropic code reduction

## Result and the constructive target

The [code-quotient construction](code-quotient-crossings.md) used uniform
averaging. There is a larger exact interface: **any normalized boundary
functional fixed by the discarded-space transform preserves the crossing
duality**. Its weights may have signs or finite-field cancellations.
This removes the positivity assumption behind the preceding classical
pressure restriction; it does not itself prove a better rank profile.

This note supplies the intertwining proof, complete common-factor
normalization, explicit effective boundary families, and the isotropic
code special case. Two full receipts are retained, with no critical
pressure witness. The factoring exponent remains `1/2`; no Lean or
independent mathematical review is claimed. Work used substantial OpenAI
Codex assistance.

## 1. The exact intertwining identity

Let `W=Even(n)` for odd `n`, with quadratic form `q(x)=wt(x)/2 mod2`.
Split its polar space orthogonally as `W=V perpendicular D`, both
summands nondegenerate. Write `dim V=2lambda`, `dim D=2s`, so
`K=(n-1)/2=lambda+s`. The regular quadratic transforms factor as
`T_W=T_V tensor T_D`; their Gaussian signs and normalizing dimensions
multiply under orthogonal direct sum.

Let `b` be a boundary vector satisfying

```math
T_D b=b,\qquad \sum_{d\in D}b_d=1.
```

The transforms are symmetric, so these are also the required left
functional identities. For the ambient folded polynomials define

```math
\Phi_v(z)=\sum_{d\in D}b_d f_{v+d}(z).
```

The ambient identity from the [folded-cube note](folded-bell-crossings.md)
gives

```math
T_V\Phi(z)
=(I\otimes b^T)(T_V\otimes T_D)f(z)
=z^K\Phi(1/z),\qquad \Phi(1)=\boldsymbol1.
\tag{1}
```

No positivity or characteristic-zero hypothesis is used. A quadratic
isometry and an input translation convert `T_V` to the physical Bell
transform, as in the parent construction, preserving both ranks and
the endpoint.

There is a general effective source of boundaries. Given a vector `a`
whose coordinate sum is nonzero in the coefficient field, set

```math
b=\frac{a+T_Da}{2\sum_d a_d}.
\tag{2}
```

Regularity and involutivity prove both required identities directly.
For instance, in `F_3` a six-bit cubic block
`a_x=(-1)^(x0*x1*x2+x2*x3*x4)` has integer sum 40. Products of such
blocks have nonzero sums modulo three, so (2) constructs explicit
non-Gaussian boundaries of growing dimension.

## 2. Removing every common scalar cost correctly

Let `a0=min_v ord_0 Phi_v` and let the gcd of all components, normalized
to have value one at `z=1`, be `g(z)=z^a0 d(z)`, where `d(0)!=0`.
Put `delta=deg d`. Equation (1) and invertibility of `T_V` imply

```math
\max_v\deg\Phi_v=K-a0,
\qquad z^\delta d(1/z)=d(z).
```

For the second identity, take the componentwise gcd on both sides of
(1). The nonzero roots, with multiplicities, occur in reciprocal pairs.
The possible scalar is one because evaluation at `z=1` is nonzero.
Dividing every component by `g` therefore gives an exact normalized path

```math
P_v(z)=\Phi_v(z)/g(z),\qquad
T_VP(z)=z^hP(1/z),\qquad P(1)=\boldsymbol1,
\quad h=K-2a0-\delta=K-a0-\deg g.
\tag{3}
```

Thus a common zero of order `a0` at zero removes **twice** that order
from the reciprocal degree. A nonzero reciprocal common factor removes
its ordinary degree. Merely subtracting `deg g` from `K` overcharges
the monomial factor. The retained control with `g=z` has `K=6`, `h=4`.

All operations are effective. For `s=O(lambda)`, both source spaces have
`2^O(lambda)` coordinates. Polynomial evaluation, quadratic transforms,
gcd computation and complete factor listing cost `2^O(lambda)` bit
operations; degrees are `O(lambda)`. A fixed finite seed can be tensor
powered with the same bound. No discarded-state oracle is assumed.

## 3. The isotropic-code special case

Let `R<=W` be totally singular: `q(r)=0` for every `r in R`. Its polar
products vanish as well. Average over `R` and restrict to `R^perp`:

```math
P_{[v]}(z)=|R|^{-1}\sum_{r\in R} f_{v+r}(z),
\qquad [v]\in R^\perp/R.
```

The quotient is nondegenerate and has dimension `2(K-dim R)`.
Averaging the ambient transform over `R` kills terms outside `R^perp`
by character orthogonality. Grouping the surviving terms into cosets
of `R` yields its correctly normalized quotient transform. Thus the
same crossing identity descends, with physical label width `K-dim R`.
This is not the LCD quotient of the earlier note: a whole radical is
removed on both sides, and its physical dimension must be recomputed.

Equivalently extend `R` to a hyperbolic summand `D` of dimension
`2dim R`; the normalized indicator of `R` is fixed by `T_D`. It is a
particular boundary in Section 1. All common scalar factors and any
monomial shift still need the normalization (3). A code's distance or
rate alone is not a crossing-width estimate.

## 4. Mixing the source embedding without losing duality

For any nonzero even vector `u`, the binary map

```math
v\longmapsto v+(v\cdot u)u
```

preserves the polar form. Applying such symplectic transvections to
the bases of `V,D` preserves their nondegeneracy and orthogonality.
It can change the Hamming weights in the ambient polynomials and the
individual quadratic signs. The constructor recomputes these signs;
it does not assume that every transvection preserves the quadratic form.
This allows genuinely different embeddings of the same physical and
boundary dimensions while (1) remains exact.

The implementation evaluates quadratic convolution with a Walsh
transform followed by the exact polar-coordinate permutation. Independent
dense convolution controls check that coordinate conversion, involution
and the all-ones endpoint.

## 5. Evidence and unresolved strict inequality

```sh
python3 tools/check_bell_projection.py
```

The [two explicit receipts](../evidence/bell-signed-projection.json) retain
their seed, every transvection, common divisor, full reduced coefficients,
and exact rational pressure. One checks the monomial degree saving; the
other supplies an explicit normalized signed boundary with pressure rho
approximately `0.864988`. Both preserve all 64 physical eigenvalues and
the exact crossing. Neither meets the critical threshold.

The readback reconstructs every coefficient and checks exact duality,
endpoint, common-factor cost, rational pressure, and the sufficient
critical inequality. It also passes 48 independent dense transform controls.
The receipts claim existence of those paths, not optimality over boundary
functionals or embeddings.

The remaining target is a family of intertwiners whose reduced order
distribution supplies `rho<1/sqrt(2)`, or a direct actual-rank argument
beyond the order/degree pressure. The effective algebraic interface is
available; that strict bound remains open.
