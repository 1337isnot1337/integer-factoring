# A profile floor for fixed homogeneous cube projections

## Question and result

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338064)
tests a valid way to use the existing Gaussian cube with two parent axes:
keep its physical leg grouping fixed, and substitute

~~~
x=X^a,    y=Y^a,    z=X^b Y^(a-b),
a>=1,    0<=b<=a,    a,b integers.                     (1)
~~~

This preserves the endpoint and makes the kernel homogeneous of degree
6a in X,Y. It avoids the invalid assignment-dependent leg permutations in
the [rejected regrouping](two-axis-cube-regrouping.md). Each parent label
has two original Gaussian modes, so the earlier one-mode two-axis profile
floor cannot simply be cited without checking its hypotheses.

For the whole family (1), the **common-physical-type, common-scale
Smith-profile envelope** obeys

~~~
rho_profile >= 3/4 > 1/sqrt(2).                        (2)
~~~

The balanced projection b=a/2 minimizes this envelope pointwise in the
scale. An exact degree-64 Bernstein certificate proves the remaining
all-scale inequality. No continuous or finite parameter optimization is
used to infer the floor.

This closes the direct use of (1) with that profile estimate. It does not
bound actual coefficient ranks from below. Additional joint-support
information, different physical types, other parameter maps, other cuts
or a different global evaluator remain outside the conclusion. Generic
Smith data are not a free constructor; selected finite Laurent
factorizations below are controls, not a general constructor theorem.
The complete factoring score remains **1/2**.

## 1. Exact one-particle profiles for every projection

Retain the kernel and augmentation basis of the
[homogeneous cube construction](cube-homogeneous-gaussian.md). The X cut
is multiplication by A+B epsilon+C eta+D epsilon eta, with

~~~
A=X^(6a),
B=X^(3a)Y^(3a),
C=X^(3a+3b)Y^(3a-3b),
D=X^(2a+2b)Y^(4a-2b).
~~~

Its augmentation matrix is

~~~
M = [ A 0 0 0
      B A 0 0
      C 0 A 0
      D C B A ].
~~~

Over F_2(Y)[[X]], the first determinantal valuation is

~~~
s_X=min(3a,2a+2b).
~~~

Every two-minor has valuation at least 6a. The possible nonzero forms are
among A^2, AB, AC, AD+BC, B^2, BC and C^2; their valuations have that
bound. The minor on rows (1,3), columns (0,2) is B^2 and has valuation
exactly 6a. Therefore the first two Smith orders sum to 6a.
Characteristic two gives M^2=A^2 I, so inverse Smith pairing makes
opposite orders sum to 12a. The complete profile is consequently

~~~
(s_X, 6a-s_X, 6a+s_X, 12a-s_X).                       (3)
~~~

The Y cut has the same formula with b replaced by a-b:

~~~
s_Y=min(3a,4a-2b).                                    (4)
~~~

In both cases 2a<=s<=3a. One of s_X,s_Y is 3a, and both equal 3a
precisely at the balanced projection. These identities are all-parameter
minor arguments, not extrapolations from the finite panel.

## 2. The balanced profile is pointwise best in this family

For scale z>0 let

~~~
Z_s(z)=(1+z^s)(1+z^(6a-s))(1+z^(6a+s))(1+z^(12a-s)).
~~~

The function f(u)=log(1+z^u) is convex in u, since its second derivative
is `(log z)^2*z^u/(1+z^u)^2>=0`. Apply Jensen separately to the pairs
(s,6a-s) and (6a+s,12a-s), whose means are 3a and 9a. It follows that

~~~
Z_s(z) >= Z_(3a)(z) = (1+z^(3a))^2(1+z^(9a))^2.        (5)
~~~

Thus both cut partition functions are no smaller than the balanced one
at every common scale. No new profile arrangement within (1) improves
this part of the estimate.

## 3. Normalize the actual four-superleg type correctly

Use the original common occupied fraction p in (0,1) on all eight
Gaussian modes, at word length m. Each original physical letter has
entropy `mH(p)` in natural units. A parent superletter contains two such
letters and therefore has entropy `2mH(p)`.

Every surviving entry has 8pm occupied modes. Kernel degree 6a gives
Gaussian degree `3a` per occupied mode, so the X,Y coefficient degrees
sum to `24apm`. The standard profile estimate, using this conservation
law and a common scale z in (0,1), has normalized envelope

~~~
R_(a,b)(p,z)
 = (log Z_(s_X)(z)+log Z_(s_Y)(z)-24ap log z)/(8H(p)).   (6)
~~~

The denominator is four times the superletter entropy per site; this is
the two-cut width normalization, with lambda'=2lambda. The expression
uses only the profiles and total-degree conservation. It does not infer
that every degree split is actually supported or that the rank upper
bound is attained.

By (5), the right side of (6) is at least its balanced value. Put
r=z^(3a), and write Z(r)=(1+r)^2(1+r^3)^2. Then

~~~
R_(a,b)(p,z) >= (log Z(r)-4p log r)/(4H(p)).             (7)
~~~

We next prove this last expression is at least 3/4 for all p,r.

## 4. An exact all-scale partition-function certificate

The required elementary inequality is

~~~
(1+r)^2(1+r^3)^2 >= (1+r^(4/3))^3,     r>0.            (8)
~~~

For 0<r<=1 put r=u^3. Subtract the two sides. The result is u^3 Q(u),
where

~~~
Q(u)=2-3u+u^3-3u^5+2u^6+3u^9
      +2u^12+u^15+2u^18+u^21.
~~~

Write Q in the degree-64 Bernstein basis. If q_j are its monomial
coefficients, its Bernstein coefficients are the exact rationals

~~~
B_k = sum_(j<=k) q_j * binom(k,j)/binom(64,j).
~~~

All 65 coefficients are positive. Their minimum is exactly

~~~
1589593605121 / 28736878651740 > 0.
~~~

The [receipt](../evidence/fixed-cube-projection-floor.json) retains every
coefficient, and the [checker](../tools/check_fixed_cube_projection_floor.py)
independently reconstructs the monomial polynomial from them. Since every
Bernstein basis function is nonnegative on [0,1], this is a finite exact
certificate for the whole real interval, not a sample grid.

For r>1, use Z(r)=r^8 Z(1/r) and
`(1+r^(4/3))^3=r^4(1+r^(-4/3))^3`. The proven inequality at 1/r and
r^8>=r^4 extend (8) to all positive scales.

Finally, entropy duality gives

~~~
log(1+exp(v))-p*v >= H(p).
~~~

Set v=(4/3)log r and multiply by three. Combining with (8) yields
`log Z(r)-4p log r>=3H(p)`. Equations (6)--(7) prove (2). The target
comparison is exact: `(3/4)^2=9/16>1/2`.

## 5. Controls, provenance and stopping condition

The checker retains six fixed projections: both endpoints of the hidden
weight, one balanced case, two mirrored unequal cases and a common
rescaling. It checks 828 determinant minors, 12 exact finite Laurent
Smith transformations, 12 physical Gaussian cut reconstructions and
1,536 endpoint/homogeneity assignments. These validate conventions;
the all-parameter profile and inequality proofs supply the generality.
Run make verify-fixed-cube-projection.

The kernel, Cayley matrix and finite Laurent machinery retain the
existing cube contributors' credit. This projected-family argument and
its exact certificate were developed with OpenAI assistant support.
There is no new Lean or independent mathematical review.

The registered continuation condition would have required this valid
fixed-bundle projection to give a profile below 1/sqrt(2). The floor
excludes that under the stated method, so no finer p, angle, integer a,b
or larger finite projection sweep is justified. Further work must change
the supported rank information, physical restriction, parameter map or
global evaluator rather than optimize inside this closed envelope.
