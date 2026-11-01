# A uniform cold-tail bound for moving physical gauges

## Question, result and scope

The [fixed-gauge axis theorem](polynomial-crossing-axis-obstruction.md)
does not give a uniform neighborhood when physical heights change with
the supporting angle. Can a physical distribution become sufficiently
cold to evade that theorem? The mechanism to test was an angle-dependent
height blow-up, rather than another bounded-height fit. A shared two-minor
instead gives the following **uniform necessary condition** for a
critical candidate, across all strict Laurent-polynomial crossings.

Use physical alphabets X,Y,Y,X of sizes A,B,B,A, with A,B>=2. Let r,s
be the complete complementary-place Smith profiles at unit scale, and
G_X,G_Y the restored physical energy lists. Let a,b>0 satisfy
`a^2+b^2=1`, put `c=a+b`, and let t>0 be the common scale. Write M for
the critical margin of the fixed-gauge theorem. Let E be the gap between
the two least values of the additive list `G_X(x)+G_Y(y)`.

For **arbitrary rational physical heights**, at every angle and scale,

```
M >= c*log(1+exp(-t*E/c))
     -log(sum_x exp(-t*(G_X(x)-min G_X)))
     -log(sum_y exp(-t*(G_Y(y)-min G_Y))).                 (1)
```

There is no assumption that the gauges are fixed as the angle or scale
varies. If E>0, put `C=A+B-2`. A simpler consequence is

```
M > 0  whenever  t*E >= c/(c-1)*log(2*C/c).              (2)
```

Thus every nonpositive-margin candidate with a unique physical ground
pair must lie below that explicit threshold. Near either axis, the
threshold is `O(1/min(a,b))` at fixed alphabet sizes. This excludes a
cold blow-up with `tE` growing faster than the reciprocal small axis
weight. It does **not** exclude all moving gauges, the intermediate
scale regime, a degenerate ground pair E=0, or interior candidates below
(2). The complete factoring score is unchanged at `1/2`.

## 1. Two shared-prefix inequalities

Sort r,s and the physical pair energies e, and set

```
D=e_1-r_1-s_1,   u=r_2-r_1,   v=s_2-s_1.
```

The endpoint entry and the shared principal two-minor from the credited
axis theorem imply, over every field,

```
D>=0,        u+v<=E+2D.                                (3)
```

Indeed that theorem gives `r_1+s_1<=e_1` and
`r_1+r_2+2*s_1<=e_1+e_2-(s_2-s_1)`. Substitute the definitions.
The two least additive physical states can be chosen to share one
coordinate; hence their strong principal two-minor is a weak minor too.
This argument includes repeated cut minima. It requires the literal
crossing endpoint, or separately established matched principal minors;
a type-only endpoint does not provide it automatically.

## 2. Convexity removes the tensor-dependent ground slack

Factor the minimum exponential out of each partition sum. Keeping just
the first two terms at each cut gives

```
M >= tD+a*f(tu/a)+b*f(tv/b)-L_X-L_Y,
f(z)=log(1+exp(-z)),
```

where L_X,L_Y are the two ground-subtracted physical logarithms in (1).
Convexity and monotonicity of f, followed by (3), give

```
a*f(tu/a)+b*f(tv/b) >= c*f(t*(u+v)/c)
                              >= c*f(t*(E+2D)/c).
```

For fixed E>=0, the function `tD+c*f(t*(E+2D)/c)` is nondecreasing
on D>=0: its derivative is

```
t-2*t/(1+exp(t*(E+2D)/c)) >= 0.
```

Its value at D=0 is `c*f(tE/c)`, proving (1). No profile tail is
assigned a fictitious order; deleting all later *positive* partition
terms supplies a lower bound.

## 3. A uniform explicit cold region

When E>0, both physical minima are unique. Every nonminimal individual
physical energy has gap at least E, so `log(1+z)<=z` gives

```
L_X+L_Y <= (A+B-2)*exp(-tE).
```

For z>=0, `log(1+exp(-z)) > exp(-z)/2`. Consequently (1) gives

```
M > (c/2)*exp(-tE/c)-C*exp(-tE).
```

Condition (2) makes this positive. The strict inequality also handles
equality at the stated threshold. Since `c-1` is asymptotic to the
smaller weight at an axis, its reciprocal is the stated scale.

The proof is pointwise in all physical heights; no compact set of gauges,
fixed weak gap, simple Smith minimum, fixed scale/angle ratio, field
extension bound, or numerical optimizer is assumed. Real heights follow
by continuity: each weighted determinantal prefix is the minimum of
finitely many affine height functions, and orders are differences of
successive prefixes. The explicit theorem may equivalently be used only
at rational heights, as required for effective crossing certificates.


## Exact arithmetic and independent controls

The [rational checker](../tools/cold_gap_bound.py) verifies a sufficient
cold exclusion using positive Taylor terms for `exp(tE*(c-1)/c)`. It
compares their rational sum against `2*(A+B-2)/c`; there is no
floating-point logarithm, exponential underflow, or sign tolerance in
that certificate. A failed finite Taylor comparison is inconclusive.
The source crossing and unique-ground hypotheses remain explicit.

The [control receipt](../evidence/cold-gap-bound-controls.json) retains
408 fresh independently computed untruncated orders on constant crossings,
nonbinary Brauer polynomials and random strict polynomial tensors, with
square and unequal physical alphabets and gauges at both places. It
checks all 48 shared-gap instances and 768 numerical lower comparisons.
There are 184 applications of 16 distinct rational weight/alphabet cold
certificates, including angles near both axes. A fabricated spectrum
violating the shared-prefix condition has negative cold margin, and an
insufficient Taylor bound remains explicitly inconclusive. The ordinary
proof above supplies uniformity; finite comparisons do not replace it.

```
make verify-cold-gap
```

This reuses the credited endpoint/shared-minor proof and exact polynomial
Smith engine. The new step is the joint two-gap convexity bound, which
removes the tensor-dependent ground slack and quantifies a cold region
uniformly over moving gauges. Directed by Cinna Davis / 1337isnot1337
with substantial OpenAI assistant support. No new Lean or stronger
complete factoring bound is claimed.
