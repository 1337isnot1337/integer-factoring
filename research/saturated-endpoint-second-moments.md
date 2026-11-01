# A second-moment test at a saturated supporting endpoint

The [principal-alphabet experiment](pencil-principal-alphabets.md) finds
support `1.000001337545` on the four-letter unperturbed face
`{0,2,3,4}` on both pairs. It remains above one. This face differs from
the full five-letter source: its first-order supporting endpoint gap can
vanish. The following sufficient condition tests a coupled angle/scale
escape. It gives no new factoring bound without an actual nonpositive
second-moment witness.

## Sufficient condition

Fix an admissible crossing with paired alphabet sizes `A,B>=2`, `n=AB`,
complete unit-scale order lists `r,s` of length `n` in both cuts, and restored paired
energies `G_X,G_Y`. All means and variances below are uniform averages
over their respective lists. Let `m=min(s)`, and define

```math
 D=\operatorname{mean}(G_X)+\operatorname{mean}(G_Y)
       -\operatorname{mean}(r)-m,
 \qquad
 V=\operatorname{Var}(r)-\operatorname{Var}(G_X)-\operatorname{Var}(G_Y).
```

**Claim.** If `D=0`, the minimum of `s` is simple, and `V<=0`, then
some positive supporting weights of Euclidean norm one and some positive
scale have strict supporting pressure below one. Rational data and weights
can be chosen after a strict negative pressure has been obtained, by
continuity. The existing [asymmetric transfer](asymmetric-modular-transfer.md)
still supplies the separate uniform algorithmic conclusion.

Here is a direct proof of the strict pressure assertion. Put
`a=sqrt(1-b^2)`, with `b` positive and tending to zero, and `t=b*x`.
The critical margin is

```math
 M=a\log\sum_i e^{-t r_i/a}+b\log\sum_i e^{-t s_i/b}
       -\log\sum_xe^{-tG_X(x)}-\log\sum_ye^{-tG_Y(y)}.
```

For fixed finite lists, Taylor expansion of the strong and physical
partition sums gives

```math
 M=b x D+b\log\!\left(1+\sum_{s_i>m}e^{-x(s_i-m)}\right)
       +\frac{b^2}{2}\bigl(x^2 V-\log n\bigr)
       +O(b^3x^3+b^4(1+x^2)).
```

The simple minimum is essential to the displayed weak-channel expression.
Let `delta` be its positive gap to the next order, choose any fixed
`c>1/delta`, and set `x=c*log(1/b)`. The weak-channel remainder is
`o(b^2)` and the Taylor remainder is also `o(b^2)`. If `V<0`, divide
by `b^2*x^2`; the limit is `V/2<0`. If `V=0`, divide by `b^2`; the
limit is `-log(n)/2<0`. Thus `M<0` for small enough `b`. Since paired
Gibbs pressure decreases strictly with its rate at positive entropy,
the supporting root is below one. This is an existence argument;
a proposed numerical witness still needs exact profiles and pressure.

## Registered four-letter test

On the fixed four-letter face, take the inherited heights at the weak
place and allow all twelve nonreference symbol heights at the strong
place to vary. Determinant valuation fixes the strong mean's dependence
on those heights. The restored physical mean has the same dependence,
so `D` remains zero. The weak profile stays unchanged with a simple
minimum and gap one. This holds for either choice of strong place.

The inherited point has **exact `V=1`**, so it does not satisfy the
sufficient condition. The [moment evaluator](../tools/pencil_endpoint_moments.py)
retains both complete profiles, both physical energy vectors, the mean
gap and an analytic variance derivative. The
[checker](../tools/check_pencil_endpoint_moments.py) independently verifies
192 untruncated orders and six rational moment identities, and controls
24 derivatives. The [receipt](../evidence/pencil-endpoint-moments-controls.json)
records those checks; it is not a global positive variance certificate.

The [manifest](../evidence/pencil-endpoint-moments-manifest.json) fixes
one hundred starts for each strong-place choice before evaluation,
height bound 32, perturbation scales and time/iteration ceilings.
The [driver](../tools/search_pencil_endpoint_moments.py) searches `V`
directly. A value near zero triggers exact replay before any claim.
Positive numerical minima leave other heights and paths open.

```sh
python3 tools/check_pencil_endpoint_moments.py
OPENBLAS_NUM_THREADS=1 python3 tools/search_pencil_endpoint_moments.py \
  --manifest evidence/pencil-endpoint-moments-manifest.json \
  --output /tmp/pencil-endpoint-moments.json
```

No new Lean formalization or factoring improvement is claimed.

The frozen `1ff76b4` [search archive](../evidence/pencil-endpoint-moments-search-complete.json.xz)
now contains all **200** completed optimizer calls, all with local success,
and **46,224** successful complete-profile queries with zero rejections.
It exhausted the registered restart lists in 19.84 seconds, before either
time ceiling. The lowest observed variance difference is **one at both
axes**, attained at the independently checked natural point. No nonpositive
witness was found. These local searches do not prove that one is a global
minimum, so the sufficient condition remains unfulfilled for this face.

The [recorded-profile filter](../tools/filter_principal_endpoint_moments.py)
has a separate [registration](../evidence/principal-endpoint-moment-filter-manifest.json):
check both endpoint orientations of every one of the 45,292 retained
principal restrictions, using exact integer moment sums and rational
normalization. A negative `D` is also a sufficient endpoint escape, even
when the weak minimum is multiple: for large fixed `x`, the coefficient
`x*D+log(sum exp(-x*(s-min(s))))` is negative. The filter records every
such case and every saturated simple case with `V<=0` for independent
profile replay. It does not vary the recorded heights or certify an
unseen profile. No result of this filter is asserted at registration.

The [completed exact filter](../evidence/principal-endpoint-moment-filter.json)
checks all **90,584** registered endpoint states: 89,645 have positive
first-order gap and 939 have zero gap. Every zero-gap weak minimum is
simple; the least variance difference among them is exactly `1/2`.
No candidate satisfies either sufficient criterion. These are exact
arithmetic statements about the retained inherited profiles, not a
free-height exclusion.

The later [polynomial-axis theorem](polynomial-crossing-axis-obstruction.md)
now proves a global positive variance floor
`V >= (s_2-s_1)^2/(2*(AB-2))` whenever `D=0` and the weak minimum is
simple in a strict Laurent-polynomial crossing. It therefore prevents
the nonpositive-variance condition throughout that class and excludes
both axes at every scale for each fixed physical gauge. Rational paths
with additional poles, general type-only endpoints, moving heights and
interior angles retain their separate obligations.
