# Fixed physical heights exclude both supporting axes at every scale

## Statement and scope

Let `T(z)` be a finite Laurent-polynomial tensor over any field, with
physical alphabets `X,Y,Y,X` of sizes `A,B,B,A`, `A,B>=2`, and literal
crossing endpoint `T(1)=[x=v][y=u]`. Fix rational physical heights on all
four legs at each of zero and infinity. Put `n=AB`. The two complete
complementary-cut profiles at unit scale are `r,s`, and the restored paired
energy vectors are `G_X,G_Y`. Choose `r` as the strong channel.
The scalar `t` multiplies these unit-scale orders and energy vectors; the
fixed gauge refers to their monomial exponents, rather than oracle weights
held numerically constant while `t` changes.

There is an `epsilon>0`, depending on this fixed tensor and these fixed
heights, such that for **every scale `t>0`** and every unit positive
supporting pair `a=sqrt(1-b^2)`, `0<b<epsilon`, the critical margin

```math
 M=a\log\sum_i e^{-t r_i/a}+b\log\sum_i e^{-t s_i/b}
       -\log\sum_xe^{-tG_X(x)}-\log\sum_ye^{-tG_Y(y)}
```

is strictly positive. The other supporting axis has the same conclusion
after exchanging the cuts and reversing the Laurent parameter.

Thus a fixed strict polynomial crossing and fixed physical heights cannot
improve the supporting threshold by approaching an axis, **even while
varying the intrinsic scale arbitrarily**. This strengthens the earlier
[fixed-ratio pencil endpoint argument](pencil-uniform-transport-endpoints.md)
in this direction. It gives no uniform neighborhood when heights or the
tensor vary with the angle, no exclusion of interior angles, and no
factoring lower bound. Rational paths with additional finite poles and
general type-only local endpoints are outside the stated hypotheses.

The proof uses the existing [shared-minor principle](shared-minor-smith-floor.md)
and endpoint minor bounds. No new field-specific rank assumption or Lean
formalization is used.

## Endpoint and shared-minor bounds

Clear height denominators. At the weak place, let `m=s_1` and normalize
the physically gauged tensor by `z^-m`. All its entries are polynomials:
`m` is their least valuation. Denote this tensor by `U(z)`; it still has
crossing endpoint at one. At the strong place its profile with the summed
physical leg heights is

```math
 \lambda_i=r_i+m.
```

Let `e_1<=...<=e_n` be the multiset of `G_X(x)+G_Y(y)`. After matching
the strong columns with their endpoint rows, every principal minor of
`U(1)` equals one. For any `k` chosen matched states, the corresponding
minor polynomial has degree at least zero. Its weighted infinity order
is therefore at most the sum of their `e` values. Taking the `k` smallest
states gives every prefix bound

```math
 S_k(\lambda)\le S_k(e)\qquad(1\le k\le n).                 (1)
```

In particular

```math
 D=\operatorname{mean}(e)-\operatorname{mean}(\lambda)\ge0.
```

Write `delta=s_2-s_1`. The two least values of the additive array
`G_X(x)+G_Y(y)` can be chosen to share one coordinate: they are
`X_1+Y_1` and `min(X_2+Y_1,X_1+Y_2)`, including ties. Their strong
principal two-minor is also a weak two-minor, up to transpose. The weak
second determinantal valuation of `U` is `delta`, so this polynomial has
valuation at least `delta` and hence degree at least `delta`. It is nonzero
because its endpoint minor is one. Therefore

```math
 S_2(\lambda)\le S_2(e)-\delta.                            (2)
```

The entry/degree argument applies to Laurent tensors after the stated
normalization. It also applies to a triangular endpoint whenever all the
matched principal minors used in (1) are separately proved nonzero.
The mere existence of a typed crossing endpoint does not supply that fact.

For the coordinate Pascal endpoint `[b=y]*det(P[v,x])`, matching the
columns produces the identity on the `B/Y` alphabet tensored with `Q`,
or its transpose, where `Q` is a principal
restriction of the compound Pascal matrix. It is unit triangular on each
occupation block, so every principal minor is one. Thus the proof also
covers the current [coordinate Pascal constructors](rectangular-pascal-coordinate-alphabets.md),
including the constant nonbinary flags that preserve this endpoint.

## A global positive second-moment floor on the saturated locus

Suppose `D=0` and the weak minimum is simple, so `delta>0`. Put
`d_i=lambda_i-e_i`. Then `sum d_i=0`, every prefix sum is nonpositive,
and the second prefix is at most `-delta`. Abel summation gives
`sum e_i*d_i>=0`. Cauchy--Schwarz on the first two and last `n-2`
coordinates gives

```math
 \sum_i d_i^2\ge\frac{\delta^2}{2}+\frac{\delta^2}{n-2}.
```

Since the means agree,

```math
 \boxed{\operatorname{Var}(r)-\operatorname{Var}(G_X)
       -\operatorname{Var}(G_Y)
       =\operatorname{Var}(\lambda)-\operatorname{Var}(e)
       \ge\frac{\delta^2}{2(n-2)}>0.}                     (3)
```

This holds for every fixed rational physical gauge satisfying the
saturation condition, not only the tested four-letter heights. It rules
out the nonpositive-variance hypothesis of the earlier
[sufficient escape test](saturated-endpoint-second-moments.md) throughout
the present strict polynomial class.

## Uniformity over the scale

Let `Z_v(u)=sum_i exp(-u*v_i)`, and let `Q(u)=Z_(s-m)(u)`. The margin is

```math
 M=a\log Z_\lambda(t/a)-\log Z_e(t)+b\log Q(t/b).           (4)
```

The prefix bounds (1) imply `Z_lambda(u)>=Z_e(u)` for every positive
`u`. One elementary proof integrates the hinge sums
`sum_i(c-v_i)^+=max_k(k*c-S_k(v))` against `u^2*exp(-u*c)`.

If the weak minimum has multiplicity `k>=2`, then `Q>=k`. The derivative
of `c*log Z_e(t/c)` with respect to `c` is a Gibbs entropy between zero
and `log n`. Hence

```math
 M\ge b\log k-(1-a)\log n>0
```

for all small enough `b`, uniformly over every `t>0`.

It remains to consider a simple weak minimum. Equation (2) makes the
partition comparison strict, `Z_lambda(u)>Z_e(u)` for every `u>0`.
For `t` in a fixed compact interval bounded away from zero, (4) is
therefore positive for all sufficiently small `b`, uniformly on that
interval. We check both ends of the scale range explicitly.

For small `t`, uniformly with `a` near one, Taylor expansion gives

```math
 a\log Z_\lambda(t/a)-\log Z_e(t)
 =(a-1)\log n+tD
   +\frac{t^2}{2}\left(\frac{\operatorname{Var}(\lambda)}a
                              -\operatorname{Var}(e)\right)+O(t^3).
```

If `D>0`, choose the fixed small-scale interval so that this expression
is at least `-b^2 log n+tD/2`. For `t>=K*b^2`, a fixed sufficiently
large `K` makes it positive. For `t<K*b^2`, `t/b` tends uniformly to zero,
so the positive term `b*log Q(t/b)` dominates `b^2 log n`.

If `D=0`, (3) supplies a fixed positive variance difference `V`.
On a sufficiently small fixed interval the expression is at least
`-b^2 log n+V*t^2/4`. This is positive for `t>=K*b` with sufficiently
large fixed `K`. For `t<K*b`, `log Q(t/b)>=log Q(K)>0`, and again its
linear factor `b` dominates the possible quadratic loss. Thus every
small scale is covered.

For large `t`, let `ell=min(lambda)` and `u=min(e)`. From (1), `ell<=u`.
If `ell<u`, the positive linear term `t*(u-ell)` eventually dominates
all bounded logarithmic terms, uniformly with `a` near one. If `ell=u`,
the minimum multiplicity of `lambda` is at least that of `e`. If it is
larger, the positive limiting entropy difference gives the same uniform
conclusion.

If these multiplicities are equal, (2) forces both to be one: otherwise
the two least values of both lists would equal `u`, contradicting their
strict second-prefix gap. Put `L=lambda_2-ell` and `E=e_2-u`. Then
`0<L<=E-delta<E`. Choose fixed `a_0` with `L/E<a_0<1`. For `a>=a_0`,
factoring the equal minimum values gives the lower bound

```math
 a\log Z_\lambda(t/a)-\log Z_e(t)
 \ge\frac{a_0}{2}e^{-tL/a_0}-(n-1)e^{-tE}>0
```

for every sufficiently large fixed `t`, uniformly over `a>=a_0`.
The weak-channel term in (4) is nonnegative. This covers all large
scales. Combining the small, compact and large ranges proves the theorem.

## Exact finite controls

The [checker](../tools/check_polynomial_crossing_axis.py) reconstructs
strict binary/nonbinary polynomial crossings, including rectangular
examples, reversed paths and constant-determinant Brauer loops with weak
gaps one, two and three. It uses untruncated fraction-free Smith
elimination and exact rational moments to check (1), (2) and (3).
The [receipt](../evidence/polynomial-crossing-axis-controls.json) states
the actual finite coverage. These controls test the algebraic identities;
the uniform-in-scale conclusion follows from the proof above, not a mesh.

The fresh controls check 2,600 literal endpoints, **496 independent
untruncated orders**, 496 endpoint prefixes, sixty shared second-prefix
bounds and eighteen saturated simple variance floors. A separate replay
checks every retained principal profile: **753,750 endpoint prefixes**,
90,584 shared second-prefix bounds and all 939 saturated simple variance
floors. That replay does not independently reconstruct all native orders.

```sh
python3 tools/check_polynomial_crossing_axis.py --recorded-principal
```

The complete factoring score remains `1/2`. Interior supporting angles,
physical heights varying with the angle, and other admissible path
frameworks remain open. No new Lean formalization is included.

The [full verification receipt](../evidence/verify-1a5f1f9.json) and
[compressed log](../evidence/verify-1a5f1f9.stdout.txt.xz) record clean
exact source `1a5f1f9`: full defaults and three named focused gates passed,
352 Python command blocks and 66 plus five unit tests. The polynomial-axis
gate is part of those defaults. Later completed-search archives, the F4
piece-count proof and the subsequent F8 refinement have their separate
focused checks; this whole-suite pin is not extended to later trees.
