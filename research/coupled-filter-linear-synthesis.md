# Linear synthesis with a certified zero-cost upper bound

## Result and scope

The bounded coupled identity cores support an exact linear search over
arbitrary polynomial joint occupation weights. The search can impose a
sufficient upper bound on the complete zero partition while charging
every infinity degree and preserving the actual crossing endpoint.

This gives a way to search beyond monomial filters without repeatedly
computing large Smith matrices. The certificate is sufficient, so an
infeasible system excludes only its specified degree bounds and chosen
energy split. It does not exclude all polynomial filters, other zero
profiles, or off-diagonal occupation kernels.

The finite regions checked here contain no pressure improvement. In the
stated common-degree boxes, the only solution is the Gaussian multiplied
by a common scalar monomial. Larger degree allowances on the tested rare
infinity states also fail to compensate for the specified improvements
on occupied states. These are exact finite results, with the scope
recorded in the evidence.

## 1. A minor bound from entry valuations

Work over a fixed Puiseux extension of `F_p((x))`, where `x=z^2` and
the valuation of `x` is one. In one zero charge block, the active
Gaussian is diagonal up to constant factorial units and zero units.
Let its diagonal orders, sorted increasingly, be

\[
E_1\le\cdots\le E_N.
\tag{1}
\]

The row and column labels carry the same energy list. Choose a rational
number `0<=theta<=1` and an integer `h`. Suppose a candidate matrix
satisfies

\[
\mathrm{ord}_0 T_{ij}
 \ge h+\theta E_i+(1-\theta)E_j
\qquad\text{for every entry.}
\tag{2}
\]

Zero entries satisfy every order lower bound. When the right side is
not an integer, the coefficient equations simply impose its ceiling.

For any `k` selected rows and columns, factor the prescribed powers
from that minor. The remaining matrix is regular. Its determinant has
nonnegative valuation, so the minor has order at least

\[
kh+\theta\sum_{i\in I}E_i+(1-\theta)\sum_{j\in J}E_j
 \ge kh+\sum_{i=1}^k E_i.
\tag{3}
\]

If `lambda_1<=...<=lambda_N` are the candidate Smith orders, the
determinantal-ideal characterization gives

\[
\sum_{i=1}^k\lambda_i
 \ge kh+\sum_{i=1}^kE_i,\qquad 1\le k\le N.
\tag{4}
\]

No equality of total determinant valuations is required. For every real
`t`, (4) implies

\[
\sum_i\bigl(t-(\lambda_i-h)\bigr)_+
 \le\sum_i(t-E_i)_+.
\tag{5}
\]

Indeed, each side is the maximum over `k` of `kt` minus the corresponding
sum of the `k` smallest orders. Integrating (5) against a positive
exponential kernel proves the convex decreasing partition comparison.
Returning to physical `z` orders, which are twice the `x` orders,

\[
\boxed{\quad Z_0(T;s)\le s^{-2h}Z_0(G;s),\qquad s>1.\quad}
\tag{6}
\]

This argument applies separately to every charge block. The aggregate
physical gauge `J` adds a common order shift within each block, so the
comparison holds for every `J` when the entry constraints were formed
at the convenient base gauge `J=1`.

The split `theta=0` is the condition that `x^{-h}TG^{-1}` be regular
after removing diagonal zero units. The balanced split `theta=1/2`
can admit different entry patterns. The implementation supports every
rational split in the stated interval.

## 2. Polynomial occupation weights give linear equations

Use any initial, terminal, or mixed geometry from the
[coupled identity-mode reduction](coupled-identity-filter-reduction.md).
The latter two geometries retain their stated prime-power hypothesis.
For `r` selected modes, write the joint state as
`d in {0,...,p-1}^r`. The active tensor has the exact linear expansion

\[
T(x)=\sum_{\mathbf d}\alpha_{\mathbf d}(x)M_{\mathbf d}(x),
\tag{7}
\]

where the fixed coefficient tensors `M_d` come from

\[
G\left(\prod_j\mathrm{Exp}_p(-Q_j)\right)
 \prod_j\frac{Q_j^{d_j}}{d_j!}.
\tag{8}
\]

They are computed exactly in the bounded active algebra. Every
coefficient at base charge one has an even `z` exponent, so (7) is
stored using integer Laurent powers of `x=z^2`.

Prescribe nonnegative degree caps `D_d` and write every unknown
`alpha_d` as a polynomial of degree at most `D_d`. Condition (2) says
that specified coefficients of (7) vanish. These are homogeneous
linear equations over `F_p`. The endpoint adds the affine conditions

\[
\alpha_{\mathbf d}(1)=1\qquad\text{for every joint state.}
\tag{9}
\]

The implementation first computes a sparse homogeneous nullspace and
then its image under endpoint evaluation. It returns a solution exactly
when the all-ones endpoint vector lies in that image. Successful
polynomials are replayed against the literal coefficient tensor to
verify (2), their degree caps, and (9).

Condition (9) makes the joint multiplier equal one at the crossing
endpoint. Thus the original polynomial tensor construction applies;
there is no separate repair inferred from a favorable profile. Any
ordinary polynomial degree or common denominator clearing remains
charged by that construction.

## 3. The complete cost bound

Let the selected infinity mode orders be `kappa_j`. At infinity the
joint state `d` has exact order

\[
\sum_j\kappa_jd_j-2\deg\alpha_{\mathbf d}.
\tag{10}
\]

Consequently the degree caps and (6) give the following upper bound on
the cost ratio against the Gaussian at the same physical gauges:

\[
\boxed{\quad
\frac{B_T(s)}{B_G(s)}
\le s^{-2h}
 \frac{\displaystyle\sum_{\mathbf d}
 s^{-\sum_j\kappa_jd_j+2D_{\mathbf d}}}
 {\displaystyle\sum_{\mathbf d}s^{-\sum_j\kappa_jd_j}}.
\quad}
\tag{11}
\]

All Gaussian spectators cancel in this ratio. The physical partition
is unchanged. A rigorously strict pressure inequality from (11) would
therefore yield an explicit polynomial crossing with its endpoint and
both cut bounds already justified.

The constant-degree case `D_d=h` makes the right side one. It permits
the Gaussian scalar solution `alpha_d=x^h` for every state. Lowering
some degree caps while allowing others to rise tests whether rare
infinity states can pay for an improvement on more probable states.
The linear equations decide the stated bounds exactly, although the
sufficient zero certificate can miss a useful filter whose actual
Smith profile does not satisfy (2).

## 4. Reproducible implementation and finite results

The synthesis tool is
[`coupled_filter_linear_synthesis.py`](../tools/coupled_filter_linear_synthesis.py).
For example, the following query returns the common scalar solution:

```bash
python tools/coupled_filter_linear_synthesis.py \
  --prime 2 --n 8 --geometry mixed \
  --zero-shift 1 --degrees 1,1,1,1 --balance 1/2
```

The checker is
[`check_coupled_filter_linear_synthesis.py`](../tools/check_coupled_filter_linear_synthesis.py),
and its complete records are saved in
[`coupled-filter-linear-synthesis.json`](../evidence/coupled-filter-linear-synthesis.json).
It enumerates all degree-at-most-two endpoint-one binary polynomials
on four joint states and compares literal tensors with the linear
oracle. It also records exact ranks in the displayed larger binary and
qutrit degree boxes.

The checked common-cap regions use initial, terminal, and mixed binary
cores at `n=8`, initial, terminal, and mixed qutrit cores at `n=3`, and
the mixed qutrit core at `n=9`. The tested shifts are the distinct values
among `0,1,2,floor(n/2),n`, with splits zero and one half. In every such
box the homogeneous nullspace has dimension one and is generated by
the common scalar `x^h`. The endpoint fixes its coefficient to one.

The compensation tests let the less occupied infinity states have
degree caps `n` or `2n`. They ask to reduce the degree of selected
occupied states from `h` to `h-1`, for `h=1,2`, while leaving the fully
occupied state capped by `h`. All stated systems are infeasible. These
finite computations do not assert an all-degree or all-filter theorem.

A separate positive control admits an explicit nonidentity polynomial
correction that vanishes at the endpoint and has high zero order. Its
large infinity degree is recorded. This checks that the feasibility
condition can admit nonidentity filters; it does not make that filter
competitive.

## 5. Completed mixed-direction monomial screens

The separate portable driver
[`search_coupled_identity_filters.py`](../tools/search_coupled_identity_filters.py)
also completed 729 binary exponent tables at interval eight and 107
qutrit tables at interval nine, using the two extreme identity modes.
The [binary record](../evidence/coupled-mixed-binary-n8-screen.json) and
[qutrit record](../evidence/coupled-mixed-qutrit-n9-screen.json) include
all scale and charge bounds, same-gauge comparisons, and boundary flags.
Neither screen found a same-gauge decrease beyond floating tolerance,
and neither found negative critical pressure. A very small qutrit
pressure at a scale boundary matches its Gaussian background; it is
not an improved optimized rate.
