# Larger-alphabet targets left open by the shared-minor inequalities

The later [graded adjugate constraints](graded-adjugate-caps.md) exclude
the specified five-letter negative optimizer after adding cell-set caps.
They leave a nearby explicit integer negative target. These updates
concern realizability under stronger necessary conditions; the pressure
computations in this older relaxation are retained as historical results.

## Result and scope

The shared minors that prove the universal binary Smith floor also give
necessary inequalities for larger alphabets. Together with every rectangular
principal subalphabet determinant, these inequalities admit negative-pressure
targets at four and five physical symbols. Exact rational interval arithmetic
certifies that fact.

These are points of a relaxation. No polynomial crossing realizing either
profile is supplied. They do not improve the factoring bound. Their purpose
is to distinguish a genuine opening in the present necessary inequalities
from the binary case, which those inequalities already exclude.

## 1. Normalize a general strict crossing

Let T be a Laurent-polynomial strict crossing on q symbols. Give each
physical leg an integer height at each of zero and infinity. Let H_l be
the sum of its two height functions, and let mu_0 and mu_infinity be the
minimum orders of the two weighted cuts. Write

\[
E_X(x)=H_x(x)+H_v(x),\qquad E_Y(y)=H_y(y)+H_u(y).
\]

Choose a minimum-energy symbol on each paired alphabet and relabel it zero.
Put

\[
c_l(a)=H_l(a)-H_l(0),\quad
c_X=c_x+c_v,\quad c_Y=c_y+c_u,\quad
d=E_X(0)+E_Y(0)-\mu_0-\mu_\infty.
\tag{1}
\]

Both paired lists c_X,c_Y are nonnegative and start at zero. Individual
c_l need not be nonnegative. Multiply an entry by its zero-place clearing
monomial:

\[
R_{xyuv}(w)=w^{-\mu_0+\sum_l h_{0,l}(a_l)}T_{xyuv}(w).
\]

The exact entry support bounds imply

\[
R_{xyuv}\in K[w],\qquad
\deg R_{xyuv}\le d+c_x(x)+c_y(y)+c_u(u)+c_v(v).
\tag{2}
\]

The nonzero matched entry at the four baseline symbols proves d>=0.
The normalized infinity matrix has entries
z^(d+sum c_l(a_l)) R_a(z^-1), so it too is polynomial. Its profile s and
the zero profile r are nonnegative, sorted, and have first entry zero.
For unit scale the complete critical pressure is

\[
G=d+\log\sum_i e^{-r_i}+\log\sum_i e^{-s_i}
-\sqrt2\left(\log\sum_x e^{-c_X(x)/\sqrt2}
              +\log\sum_y e^{-c_Y(y)/\sqrt2}\right).
\tag{3}
\]

A general positive scale is absorbed into d,r,s,c_X,c_Y. Rational heights
are handled by ramification. The same determinantal-prefix inequalities
extend to real heights by continuity.

## 2. Shared lines and rectangular subalphabets

Sort each paired energy list independently, and write X_a and Y_b for the
sums of its a or b smallest entries. Put gamma_i=(r_i+s_i)/2 and
S_k=sum_(i<=k) gamma_i.

Fix y=u at their baseline symbol and select k matched X,V symbols. The
resulting determinant is a minor of both flattenings and is nonzero at one.
Its degree bound is kd+X_k. Its infinity version is its reciprocal with
that same total degree shift. Therefore the sum of its two orders is at
most kd+X_k. The first k Smith orders minimize minor valuations, giving
the shared-line bounds

\[
S_k\le\frac{kd+\min(X_k,Y_k)}2\qquad(1\le k\le q).
\tag{4}
\]

Also select a matched X,V subalphabet of size a and a matched Y,U
subalphabet of size b. In either cut this gives an ab-square minor with
permutation endpoint. Each determinant term uses each X,V symbol b times
and each Y,U symbol a times. Equation (2) thus bounds the degree of either
normalized minor by ab d+bX_a+aY_b. Each profile prefix separately is at
most that bound, and hence

\[
S_{ab}\le ab d+bX_a+aY_b\qquad(1\le a,b\le q).
\tag{5}
\]

In particular S_(q^2)<=q^2 d+q(X_q+Y_q). All these minors are retained.
The two rectangular determinants in (5) are generally different
polynomials; (5) does not incorrectly add their orders as if they were
the same shared minor.

## 3. The strongest exponential bound from these prefixes

For each constrained integer k take the smallest right side in (4)--(5),
and additionally impose S_0=S_1=0. Let C(k) be the greatest convex
minorant of these finitely many cap points on [0,q^2]. Its successive
slopes gamma*_i form a sorted nonnegative list.

Every actual prefix function S is convex, since gamma is sorted. At
successive vertices of C, it lies below the cap; convexity then bounds
it by the intervening chord. Thus S_k<=C(k) for every integer k.
This implies

\[
\sum_i e^{-\gamma_i}\ge\sum_i e^{-\gamma_i^*}.
\tag{6}
\]

One direct proof uses
sum_i (u-gamma_i)_+=max_k(ku-S_k). The prefix inequality makes this
quantity at least its counterpart for gamma*. Integrating against e^-u
on [0,infinity) gives (6). No equality of the total prefix sums is needed.
Convexity of log-sum-exp also gives
log Z(r)+log Z(s)>=2log Z((r+s)/2). Therefore (3) is at least

\[
G_{\rm relaxed}=d+2\log\sum_i e^{-\gamma_i^*}
-\sqrt2\left(\log\sum_x e^{-c_X(x)/\sqrt2}
              +\log\sum_y e^{-c_Y(y)/\sqrt2}\right).
\tag{7}
\]

Within this prefix relaxation the bound is attained by r=s=gamma*.
It need not be attained by any crossing tensor.

## 4. Exact negative targets

Take d=1 and identical paired energy lists. The following lists satisfy
every stated prefix constraint:

| q | Paired energies | Average normalized orders, with multiplicities | Certified pressure |
|---|---|---|---|
| 4 | (0,0,2,8/3) | 0; 1; 3/2 twice; 3 twice; 11/3 three times; 13/3 three times; 29/6 four times | G_relaxed < -1/50 |
| 5 | (0,0,2,8/3,3) | 0; 1; 3/2 twice; 7/3; 11/3 four times; 13/3 three times; 14/3 three times; 26/5 five times; 83/15 five times | G_relaxed < -9/100 |

Their approximate pressures are -0.02315207 and -0.09763798.
The certificate uses 90-bit dyadic outward bounds for exponentials,
logarithms, and sqrt(2); its strict margins are rational comparisons.
Multiplying the four-symbol data by six gives entirely integral target
orders and energies, with the compensating scale 1/6. This removes a
possible denominator obstruction but does not supply a realizing tensor.

As a negative control, omit (4) at q=2, d=0, c_X=c_Y=(0,4). The weaker
relaxation falsely permits (0,16/3,16/3,16/3) with negative pressure.
Its second prefix exceeds the correct shared-line cap 2. This explains
why separate cut determinant budgets were insufficient in the binary
analysis and why their apparent gains must not be interpreted as seeds.

## Reproduction and next construction problem

```sh
make verify-shared-minor-targets
```

[`check_shared_minor_profile_targets.py`](../tools/check_shared_minor_profile_targets.py)
constructs all caps and their exact rational convex hull, checks every
prefix, and proves both pressure margins. The record is
[`shared-minor-profile-targets.json`](../evidence/shared-minor-profile-targets.json).

The outstanding constructive problem is to realize suitable simultaneous
profiles, or obtain other profiles with negative full pressure, while
keeping the strict endpoint and all physical costs. The
[complete-entry affine solver](full-alphabet-two-cut-affine-synthesis.md)
addresses finite choices of frames and degree windows. Failure within one
such frame does not decide the general target. Additional necessary
relations between minors could also rule out these relaxed points; none
is assumed here.
